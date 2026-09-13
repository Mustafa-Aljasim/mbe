"""Explicit engineering-file headers around unchanged canonical import contracts."""
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import re
import pandas as pd
from material_balance_studio.units.conversions import UnitSystem,to_si
from material_balance_studio.units.display import to_display,display_unit
from material_balance_studio.pvt.laboratory import lab_from_frame,lab_to_frame

QUANTITIES={'pressure':'pressure','observed_pressure':'pressure','pressure_sigma':'pressure',
    'np':'oil_volume','gp':'gas_volume','wp':'water_volume','winj':'water_volume','ginj':'gas_volume',
    'bo':'bo','rs':'rs','bg':'bg','bw':'bw','bwinj':'bwinj','bginj':'bginj',
    'oil_viscosity':'viscosity','gas_viscosity':'viscosity','z':'dimensionless'}


def default_upload_units(upload,key,project_units):
    import streamlit as st
    from hashlib import sha256
    signature=None if upload is None else sha256(upload.name.encode()+upload.getvalue()).hexdigest()
    if signature is not None and signature!=st.session_state.get(key+'_source_signature'):
        st.session_state[key]=UnitSystem(project_units).value
    st.session_state[key+'_source_signature']=signature
    st.session_state.setdefault(key,UnitSystem(project_units).value)


def file_label(name,units):
    if name not in QUANTITIES:
        return name
    unit=display_unit(QUANTITIES[name],units)
    if UnitSystem(units)==UnitSystem.SI and name in ('np','gp','wp','winj','ginj'):
        unit='m³'
    return f'{name} [{unit}]'


def engineering_file_frame(canonical,units):
    frame=canonical.copy()
    for name in frame:
        if name in QUANTITIES:
            frame[name]=frame[name].map(lambda v: v if pd.isna(v) else to_display(float(v),QUANTITIES[name],units))
    return frame.rename(columns={n:file_label(n,units) for n in frame})


def normalize_file_frame(frame,units):
    """Tagged SI pressure may be MPa or Pa; plain SI headers retain legacy Pa."""
    units=UnitSystem(units)
    result=frame.copy()
    names=[]
    for column in frame.columns:
        match=re.fullmatch(r'([a-z_]+) \[(.+)\]',str(column))
        name=match.group(1) if match else column
        names.append(name)
        if match:
            label=match.group(2)
            expected=file_label(name,units)
            if name in ('pressure','observed_pressure','pressure_sigma') and units==UnitSystem.SI and label in ('MPa','Pa'):
                if label=='MPa':
                    result[column]=pd.to_numeric(result[column],errors='raise')*1e6
            elif str(column)!=expected:
                raise ValueError(f'{column}: header units conflict with selected {units.value} source units or unsupported header. Expected {expected}.')
    if len(set(names))!=len(names):
        raise ValueError('Duplicate input columns after unit-header normalization.')
    result.columns=names
    return result


def template_frame(kind,units):
    root=Path(__file__).resolve().parents[3]/'examples'
    if kind=='Laboratory':
        raw=pd.DataFrame({'pressure':[10e6,18e6,25e6],'bo':[None,None,None],'rs':[None,None,None],
            'oil_viscosity':[None]*3,'z':[None]*3,'bg':[None]*3,'bw':[None]*3,'gas_viscosity':[None]*3})
    else:
        raw=pd.read_csv(root/(kind.lower()+'.csv'))
        if kind=='History':
            for n in ('winj','ginj'):
                if n not in raw: raw[n]=0.
            raw=raw[['date','np','gp','wp','winj','ginj','observed_pressure']]
    return engineering_file_frame(raw,units)


@dataclass(frozen=True)
class SparseLaboratoryInput:
    lab: object
    gas_viscosity: tuple = ()  # (canonical pressure, Pa s); informational, no validated regression


def sparse_lab_input(frame,units='SI'):
    raw=normalize_file_frame(frame,units)
    if 'gas_viscosity' not in raw:
        return SparseLaboratoryInput(lab_from_frame(raw,units))
    if 'pressure' not in raw:
        raise ValueError('Laboratory pressure is required.')
    extra=[]
    for _,row in raw.iterrows():
        value=row['gas_viscosity']
        if pd.isna(value): continue
        if pd.isna(row['pressure']) or not isfinite(float(value)) or float(value)<=0:
            raise ValueError('Measured gas viscosity needs valid pressure and positive finite viscosity.')
        extra.append((to_si(float(row['pressure']),'pressure',units),to_si(float(value),'viscosity',units)))
    return SparseLaboratoryInput(lab_from_frame(raw.drop(columns='gas_viscosity'),units),tuple(sorted(extra)))


def laboratory_file_frame(lab,gas_viscosity=(),units='SI'):
    raw=lab_to_frame(lab,'SI')
    values=dict(gas_viscosity)
    raw['gas_viscosity']=[values.get(p.pressure) for p in lab.points]
    return engineering_file_frame(raw,units)


def informational_composition(h2s=0.,co2=0.,n2=0.):
    values={'h2s_mole_percent':h2s,'co2_mole_percent':co2,'n2_mole_percent':n2}
    if any(not isfinite(v) or not 0<=v<=100 for v in values.values()) or sum(values.values())>100:
        raise ValueError('Informational gas mole percentages must be between 0 and 100 and total no more than 100%.')
    return values
