"""Transparent equation terms; no changes to the Phase 4A drive indices."""
from math import fsum


def tank_terms(state):
    b=state.components
    return dict(pressure=state.pressure,Fprod=b.withdrawal.production,
        water_injection=b.withdrawal.water_injection,gas_injection=b.withdrawal.gas_injection,
        Fnet=b.withdrawal.net,N_Eo=b.oil_expansion_support,N_mEg=b.gas_cap_expansion_support,
        N_Efw=b.rock_water_expansion_support,We=b.aquifer_support,X=state.intertank_support,
        residual=state.residual,relative_residual=state.relative_residual,
        transfer_support_ratio=state.intertank_support/b.withdrawal.production if b.withdrawal.production>1e-12 else None)


def maxima(result):
    states=result.internal_states or (result.initial_state,)
    return dict(maximum_tank_absolute_residual=max(abs(t.residual) for s in states for t in s.tanks),
                maximum_tank_relative_residual=max(t.relative_residual for s in states for t in s.tanks),
                maximum_network_relative_residual=max(s.network_relative_residual for s in states),
                maximum_transfer_conservation_error=max(max(abs(s.transfer_error),abs(s.incremental_transfer_error)) for s in states),
                accepted_substeps=len(result.internal_states),reported_events=len(result.states),
                function_evaluations=sum(s.diagnostics.function_evaluations for s in result.internal_states))
