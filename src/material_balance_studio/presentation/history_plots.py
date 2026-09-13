"""Paired, sparse-aware history plotting; no interpolation or stored-value changes."""
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.units.display import to_display,display_unit

VARIABLES={'Date':('date',None),'Cumulative Oil Production':('np','oil_volume'),
    'Cumulative Gas Production':('gp','gas_volume'),'Cumulative Water Production':('wp','water_volume'),
    'Cumulative Water Injection':('winj','water_volume'),'Cumulative Gas Injection':('ginj','gas_volume'),
    'Observed Reservoir Pressure':('observed_pressure','pressure'),'Producing ratio Rp = Gp/Np':('rp','rs')}


def history_values(history):
    return pd.DataFrame([dict(date=h.date,**vars(h.cumulative),observed_pressure=h.observed_pressure,
        rp=h.cumulative.gp/h.cumulative.np if h.cumulative.np>0 else None) for h in history])


def axis_values(values,quantity,units):
    if quantity is None: return list(values),'Date'
    shown=[to_display(float(v),quantity,units) for v in values]
    maximum=max((abs(v) for v in shown),default=0.)
    factor,prefix=(1.,'')
    if quantity in ('oil_volume','water_volume','gas_volume'):
        if maximum>=1e6: factor,prefix=1e6,'MM' if str(units) in ('FIELD','UnitSystem.FIELD') else 'million '
        elif maximum>=1e3: factor,prefix=1e3,'M' if str(units) in ('FIELD','UnitSystem.FIELD') else 'thousand '
    return [v/factor for v in shown],prefix+display_unit(quantity,units)


def history_figure(history,x,y,units,plot_type='Auto'):
    if x not in VARIABLES or y not in VARIABLES or plot_type not in ('Auto','Line','Scatter'):
        raise ValueError('Select supported history axes and plot type.')
    frame=history_values(history)
    xc,xq=VARIABLES[x]; yc,yq=VARIABLES[y]
    paired=frame.dropna(subset=list(dict.fromkeys((xc,yc))))
    xs,xunit=axis_values(paired[xc],xq,units)
    ys,yunit=axis_values(paired[yc],yq,units)
    mode='lines+markers' if plot_type=='Line' or (plot_type=='Auto' and xc=='date' and yc!='observed_pressure') else 'markers'
    fig=go.Figure(go.Scatter(x=xs,y=ys,mode=mode,name=y,connectgaps=False,
        hovertemplate=f'{x} [{xunit}]: %{{x}}<br>{y} [{yunit}]: %{{y:.7g}}<extra></extra>'))
    fig.update_layout(title=f'{y} vs {x}',xaxis_title='Date' if xq is None else f'{x} ({xunit})',
                      yaxis_title='Date' if yq is None else f'{y} ({yunit})')
    return fig,paired


def history_plot_workflow(history,units):
    st.caption('Rsb is a scalar fluid input; Rs(P) is dissolved gas in the PVT model; producing ratio Rp = Gp/Np comes only from history.')
    if not st.checkbox('Explore history plots',key='history_plot_enabled') or not history:
        return
    frame=history_values(history)
    options=[name for name,(column,_) in VARIABLES.items() if frame[column].notna().any()]
    x=st.selectbox('History X-axis',options,key='history_plot_x')
    y=st.selectbox('History Y-axis',[n for n in options if n!='Date'],key='history_plot_y')
    kind=st.radio('History plot type',['Auto','Line','Scatter'],horizontal=True,key='history_plot_type')
    fig,paired=history_figure(history,x,y,units,kind)
    if paired.empty: st.info('No valid paired observations for these axes.')
    else: st.plotly_chart(fig,key='history_explorer',width='stretch')
    st.caption('Only rows with both plotted values are used. Missing pressure stays missing. Lines join supplied points visually; no observations are filled or interpolated.')
