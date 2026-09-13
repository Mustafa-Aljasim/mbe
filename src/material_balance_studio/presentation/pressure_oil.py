"""Dedicated pressure/Np diagnostic and unchanged temporal plot in two panels."""
from dataclasses import replace
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from material_balance_studio.aquifer import NoAquifer
from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.display import to_display,display_unit
from .charts import pressure_figure
from .history_plots import axis_values


def without_aquifer_reference(tank,history):
    reference=replace(tank,aquifer=NoAquifer())
    return reference,simulate(reference,tuple(history))


def pressure_oil_figure(active,tank,history,units,reference=None):
    key=tank.aquifer.key if tank.aquifer else 'none'
    if key=='none': reference=active
    elif reference is None: _,reference=without_aquifer_reference(tank,history)
    history=tuple(history)
    all_volumes=[h.cumulative.np for h in history]
    scaled,unit=axis_values(all_volumes,'oil_volume',units)
    factor=next((to_display(v,'oil_volume',units)/s for v,s in zip(all_volumes,scaled) if s),1.)
    def x(v): return to_display(v,'oil_volume',units)/factor
    fig=go.Figure()
    observed=[h for h in history if h.observed_pressure is not None]
    fig.add_scatter(x=[x(h.cumulative.np) for h in observed],y=[to_display(h.observed_pressure,'pressure',units) for h in observed],mode='markers',name='Observed Pressure')
    fig.add_scatter(x=[x(s.cumulative.np) for s in reference.states],y=[to_display(s.pressure,'pressure',units) for s in reference.states],mode='lines+markers',name='Calculated — Without Aquifer')
    if key!='none':
        name=next((n for n,c in MODELS.items() if c.key==key),key)
        fig.add_scatter(x=[x(s.cumulative.np) for s in active.states],y=[to_display(s.pressure,'pressure',units) for s in active.states],mode='lines+markers',name='Calculated — With '+name)
    fig.update_layout(title='Pressure vs Cumulative Oil Production',xaxis_title=f'Cumulative Oil Production, Np ({unit})',yaxis_title=f'Average reservoir pressure ({display_unit("pressure",units)})')
    return fig,reference


def results_pressure_panels(active,tank,history,units):
    temporal=pressure_figure(active,units)
    oil,reference=pressure_oil_figure(active,tank,history,units)
    figure=make_subplots(rows=2,cols=1,vertical_spacing=.2,subplot_titles=('Pressure vs Date','Pressure vs Cumulative Oil Production'))
    for trace in temporal.data: figure.add_trace(trace,row=1,col=1)
    for trace in oil.data: figure.add_trace(trace,row=2,col=1)
    figure.update_xaxes(title_text='Date',row=1,col=1)
    figure.update_xaxes(title_text=oil.layout.xaxis.title.text,row=2,col=1)
    figure.update_yaxes(title_text=temporal.layout.yaxis.title.text,row=1,col=1)
    figure.update_yaxes(title_text=oil.layout.yaxis.title.text,row=2,col=1)
    figure.update_layout(height=850,legend=dict(orientation='h',y=-.2),margin=dict(b=150))
    return figure,reference
