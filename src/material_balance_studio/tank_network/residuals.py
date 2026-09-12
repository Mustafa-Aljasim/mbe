"""Pure candidate evaluation: existing cumulative balance minus net transfer."""
from datetime import datetime,time,timedelta
from math import fsum,isfinite
from material_balance_studio.domain.models import CumulativeVolumes
from material_balance_studio.domain.validation import finite_value
from material_balance_studio.aquifer import NoAquifer
from material_balance_studio.mbe.oil_balance import evaluate_balance
from material_balance_studio.solver.timestep import cumulative_increments
from .history import volumes_at
from .state import TankState,NetworkState,ConnectionState
from .transfer import transfers


def initial_state(network,settings):
    tanks=[]
    for node in network.tanks:
        tank=node.reservoir
        zero=CumulativeVolumes()
        p=tank.initial_pressure
        balance=evaluate_balance(p,tank,zero,settings.normalization_floor)
        _,observation=volumes_at(node,0.)
        tanks.append(TankState(node.name,p,p,zero,zero,tank.pvt_model.properties_at_pressure(p),balance,
            (tank.aquifer or NoAquifer()).initial_state(p),None,0.,0.,0.,0.,observation,()))
    p={s.name:s.pressure for s in tanks}
    edges=tuple(ConnectionState(e.from_tank,e.to_tank,e.transmissibility,e.enabled,p[e.from_tank]-p[e.to_tank],
        e.effective_transmissibility*(p[e.from_tank]-p[e.to_tank]),0.,0.) for e in network.connections)
    for edge in edges:
        finite_value("Initial connection rate",edge.rate)
    return NetworkState(datetime.combine(network.initial_date,time()),0.,0.,tuple(tanks),edges,0.,0.,0.,0.)


def evaluate_candidate(network,previous,pressures,elapsed_seconds,settings):
    if len(pressures)!=len(network.tanks) or not all(isfinite(p) for p in pressures):
        raise ValueError("Candidate pressure vector must be finite with one value per tank.")
    dt=elapsed_seconds-previous.elapsed_seconds
    if dt<=0:
        raise ValueError("Network candidate time must advance.")
    named={t.name:float(p) for t,p in zip(network.tanks,pressures)}
    edges,contributions,increments=transfers(network,previous,named,dt)
    states=[]
    for node,old,pressure in zip(network.tanks,previous.tanks,pressures):
        if node.name!=old.name:
            raise ValueError("Previous state tank order differs from network.")
        lo,hi=node.pressure_bounds
        if not lo<=pressure<=hi:
            raise ValueError(f"{node.name}: candidate pressure outside bounds.")
        tank=node.reservoir
        cumulative,observation=volumes_at(node,elapsed_seconds)
        aq=(tank.aquifer or NoAquifer()).compute_step(old.aquifer_state,old.pressure,float(pressure),dt)
        balance=evaluate_balance(float(pressure),tank,cumulative,settings.normalization_floor,aquifer_influx=aq.cumulative_influx)
        support=fsum(v for _,v in contributions[node.name])
        residual=balance.residual-support
        relative=abs(residual)/max(abs(balance.withdrawal.net),settings.normalization_floor)
        warnings=list(aq.warnings)
        pvt=tank.pvt_model.properties_at_pressure(float(pressure))
        if cumulative.winj and pvt.bwinj is None:
            warnings.append("Bwinj omitted: direct water injection uses this tank's resident Bw.")
        if cumulative.ginj and pvt.bginj is None:
            warnings.append("Bginj omitted: direct gas injection uses this tank's resident Bg.")
        if pressure>tank.initial_pressure:
            warnings.append("Pressure above initial pressure; verify PVT coverage and compressibility validity.")
        states.append(TankState(node.name,float(pressure),old.pressure,cumulative,cumulative_increments(old.cumulative,cumulative),pvt,
            balance,aq.updated_state,aq,support,increments[node.name],residual,relative,observation,contributions[node.name],tuple(warnings)))
    total=fsum(s.residual for s in states)
    norm=max(fsum(abs(s.components.withdrawal.net) for s in states),settings.normalization_floor)
    return NetworkState(datetime.combine(network.initial_date,time())+timedelta(seconds=elapsed_seconds),elapsed_seconds,dt,tuple(states),edges,
        total,abs(total)/norm,fsum(s.intertank_support for s in states),fsum(s.incremental_support for s in states))


def accepted(state,settings):
    if not all(abs(s.residual)<=settings.absolute_tolerance and s.relative_residual<=settings.relative_tolerance for s in state.tanks):
        return False
    for error,magnitude in ((state.transfer_error,fsum(abs(s.intertank_support) for s in state.tanks)),
                            (state.incremental_transfer_error,fsum(abs(s.incremental_support) for s in state.tanks))):
        if abs(error)>settings.transfer_absolute_tolerance or abs(error)/max(magnitude,settings.normalization_floor)>settings.transfer_relative_tolerance:
            return False
    return abs(state.signed_network_residual)<=len(state.tanks)*settings.absolute_tolerance and state.network_relative_residual<=settings.relative_tolerance
