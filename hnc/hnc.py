# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import dash
import matplotlib
from matplotlib import pyplot as plt
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

app=dash.Dash(__name__)
app.css.append_css({
    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
})

df=pd.read_csv('hnc/hnc_WT_meta_data.csv',header=0,index_col=0)

numerical_cols=[]
categorical_cols=[]
for c in df.columns:
    if np.issubdtype(df[c].dtype,np.number):
        numerical_cols.append(c)
    else:
        categorical_cols.append(c)

app.layout = html.Div([

    # header
    html.H1('single-cell transcriptomics of the murine salivary gland'),

    # row 2
    html.Div([

        # first column: selects
        html.Div([
            html.Label('select x coordinate'),
            dcc.Dropdown(
                id='x_select',
                style=dict(width='200px'),
                options=[dict(label=c,value=c) for c in numerical_cols],
                value='tSNE_1'
            ),
            
            html.Label('select y coordinate'),
            dcc.Dropdown(
                id='y_select',
                style=dict(width='200px'),
                options=[dict(label=c,value=c) for c in numerical_cols],
                value='tSNE_2'
            ),
            
            html.Label('select coloring'),
            dcc.Dropdown(
                id='color_select',
                style=dict(width='200px'),
                options=[dict(label=c,value=c) for c in df.columns],
                value='cell.type'
            )], className='three columns', style=dict(height='300px')),

        # second column: plot
        html.Div([
            dcc.Graph(
                id='scatter_plot',
            )], className='nine columns', style=dict(textAlign='top'))
        
    ], className='row')
    ])

@app.callback(
    dash.dependencies.Output(component_id='scatter_plot',component_property='figure'),
    [dash.dependencies.Input(component_id='x_select',component_property='value'),
     dash.dependencies.Input(component_id='y_select',component_property='value'),
     dash.dependencies.Input(component_id='color_select',component_property='value')])

def update_scatterplot(xc, yc, col):

    if col in categorical_cols:
        colvals=df[col].unique()
        nvals=len(colvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20

        traces=[]
        for n,cv in enumerate(colvals):
            tr=go.Scatter(
                x=df[df[col]==cv][xc],
                y=df[df[col]==cv][yc],
                text=df[df[col]==cv].index,
                mode='markers',
                opacity=.7,
                marker=dict(
                    size=5,
                    color=matplotlib.colors.rgb2hex(cm(n)[:3]),
                    line={'width': .1, 'color': 'gray'}
                ),
                name=cv
            )
            traces.append(tr)

    else:
        traces=[go.Scatter(
            x=df[xc],
            y=df[yc],
            text=df.index,
            mode='markers',
            opacity=.7,
            marker=dict(
                size=5,
                color=df[col],
                colorscale='Viridis',
                showscale=True),
            name=col
            )]

    return dict(data=traces,
                layout=go.Layout(
                    xaxis={'title': xc},
                    yaxis={'title': yc},
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                    legend={'x':1.05,'y':1},
                    hovermode='closest'
                ))



if __name__ == '__main__':
    app.run_server(debug=True)
                    
