# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import dash
import scipy.io
import matplotlib
from matplotlib import pyplot as plt
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.tools as tools

############################################################################################################
## get data
############################################################################################################

# load data using an environmental variable that points to a folder where meta and expression data should be
# this can be replaced by some sort of HTTP request directly to SODAR

datadir=os.getenv('DASH_DATADIR')

if True:
    # get about file (markdown) if exists
    about_file=os.path.join(datadir,'about.md')
    if os.path.isfile(about_file):
        print('reading about file from '+about_file)
        about_md=open(about_file).readlines()
    else:
        about_md=["""this is a Dash app for data from {0}""".format(datadir),
                  "",
                  """no further information available at this point"""]

    # get meta data 
    meta_file=os.path.join(datadir,'meta.tsv')
    if os.path.isfile(meta_file):
        print('reading meta data from '+meta_file)
        meta=pd.read_csv(meta_file,header=0,index_col=0,sep='\t')
    else:
        raise Exception('could not read meta data from '+datadir)

    # separate numerical and categorical columns for later
    numerical_meta=[]
    categorical_meta=[]
    for col in meta.columns:
        if np.issubdtype(meta[col].dtype,np.number):
            numerical_meta.append(col)
        else:
            categorical_meta.append(col)

    # get expression data: either as dense table, or as sparse mtx with extra files specifying gene and cell names
    expression_file=os.path.join(datadir,'expression')
    cell_file=os.path.join(datadir,'cells.tsv.gz')
    gene_file=os.path.join(datadir,'genes.tsv.gz')
    if os.path.isfile(expression_file+'.tsv.gz'):
        print('reading dense DGE from '+expression_file+'.tsv.gz')
        DGE=pd.read_csv(expression_file+'.tsv.gz',header=0,index_col=0,sep='\t')
        genes=DGE.index
        cells=DGE.columns
    elif os.path.isfile(expression_file+'.mtx.gz') and \
         os.path.isfile(cell_file) and \
         os.path.isfile(gene_file):
        print('reading sparse DGE from '+expression_file+'.mtx.gz')
        m=scipy.io.mmread(expression_file+'.mtx.gz').todense()
        cells=pd.read_csv(cell_file,header=None,index_col=None).squeeze()
        genes=pd.read_csv(gene_file,header=None,index_col=None).squeeze()
        DGE=pd.DataFrame(m,index=genes,columns=cells).fillna(0)
    else:
        raise Exception('could not read expression data from '+datadir)

    # rescale
    DGE=np.log2(1+DGE)

if False:
    marker_file=os.path.join(datadir,'markers.tsv.gz')
    if os.path.isfile(marker_file):
        print('reading markers from '+markers)
        markers=pd.read_csv(marker_file,header=0,sep='\t')
    else:
        markers=None

##############################################################################################################
## app main layout
##############################################################################################################

app=dash.Dash(datadir)
app.config['suppress_callback_exceptions']=True

# use external css (should be fixed)
#app.css.append_css({
#    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
#})
# Loading screen CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})


# specify layout: header, and then a Div with 3 main tabs
app.layout = html.Div([
    # header
    html.H1('single-cell data from '+datadir),
    html.Hr(),
    html.Div([
        
        dcc.Tabs(id='tabs',
                 value='meta',
                 style=dict(width='40%',height='5%'),
                 children=[
                     dcc.Tab(label='about this dataset', value='about'),
                     dcc.Tab(label='explore meta data', value='meta'),
                     dcc.Tab(label='explore gene expression', value='expression')
                 ]),
        html.Div(id='main-tabs', style=dict(backgroundColor='0x22222'))
    ])],
                      style=dict(marginRight=50,
                                 marginLeft=50,
                                 marginBottom=50,
                                 marginTop=50))

###############################################################################################################
## main tab contents
###############################################################################################################

# callback for main tabs
@app.callback(dash.dependencies.Output('main-tabs','children'),
              [dash.dependencies.Input('tabs','value')])

def render_main_content(tab):

    if tab=='about':

        return html.Div([
            dcc.Markdown('\n'.join(about_md))
        ],
                        style=dict(marginRight=50,
                                   marginLeft=50,
                                   marginBottom=50,
                                   marginTop=50))
    
    elif tab=='meta':

        return html.Div([

            # first column: plot controls

            html.Div([
            # top: determine plot type
                html.Div([
                    html.Label('select plot type'),
                    dcc.RadioItems(
                        id='meta_plot_type',
                        options=[dict(label='scatter',
                                      value='scatter'),
                                 dict(label='violin',
                                      value='violin'),
                                 dict(label='bar',
                                      value='bar')],
                        value='scatter')
                ], style=dict(height='100px',backgroundColor='0xAAAAA')),

                html.Hr(),

                # bottom: plot controls (depend on plot type)
                html.Div(id='meta_plot_controls',
                         style=dict(height='300px',
                                    backgroundColor='0x888888'))
                ], className='three columns', style=dict(height='500px')),
    
            # second column: plot (depends on plot type)
            html.Div(id='meta_plot',
                     className='nine columns', style=dict(textAlign='top',
                                                          marginLeft=25))
        ],className='row',
                        style=dict(marginRight=50,
                                   marginLeft=50,
                                   marginBottom=50,
                                   marginTop=50))

    elif tab=='expression':

        return html.Div([
            
            # first column: plot controls
            
            html.Div([
                # top: determine plot type
                html.Div([
                    html.Label('select plot type'),
                    dcc.RadioItems(
                        id='expression_plot_type',
                        options=[dict(label='scatter',
                                      value='scatter'),
                                 dict(label='violin',
                                      value='violin')],
                        value='scatter')
                ], style=dict(backgroundColor='0xAAAAA')),
                
                html.Label('select gene(s)'),
                dcc.Dropdown(
                    id='expression_select_genes',
                    style=dict(width='200px'),
                    options=[dict(label=g,value=g) for g in genes],
                    value=[],
                    multi=True
                ),
                
                html.Hr(),
                
                # bottom: plot controls (depend on plot type)
                html.Div(id='expression_plot_controls',
                         style=dict(height='300px',
                                    backgroundColor='0x888888'))
            ], className='three columns', style=dict(height='500px')),
            
            # second column: plot (depends on plot type)
            html.Div(id='expression_plot',
                     className='nine columns', style=dict(textAlign='top',
                                                          marginLeft=25))
        ],className='row',
                        style=dict(marginRight=50,
                                   marginLeft=50,
                                   marginBottom=50,
                                   marginTop=50))

###############################################################################################################
## meta plot controls
###############################################################################################################

# callback for meta plot controls
@app.callback(dash.dependencies.Output('meta_plot_controls','children'),
              [dash.dependencies.Input('meta_plot_type','value')])
def render_meta_controls (plot_type):

    if plot_type=='scatter':

        return [html.Label('select x axis and scaling'),
                dcc.Dropdown(
                    id='meta_scatter_select_x',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in numerical_meta],
                    value=meta.columns[0]
                ),
                dcc.RadioItems(
                    id='meta_scatter_scale_x',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label='linear',
                                  value='linear'),
                             dict(label='log',
                                  value='log')],
                    value='linear'),
                html.Label('select y axis and scaling'),
                dcc.Dropdown(
                    id='meta_scatter_select_y',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in numerical_meta],
                    value=meta.columns[1]
                ),
                dcc.RadioItems(
                    id='meta_scatter_scale_y',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label='linear',
                                  value='linear'),
                             dict(label='log',
                                  value='log')],
                    value='linear'),
                
                html.Label('select coloring'),
                dcc.Dropdown(
                    id='meta_scatter_select_color',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in meta.columns],
                    value=meta.columns[2]
                )]

    elif plot_type=='violin':

        return [html.Label('select variable(s) and scaling'),
                dcc.Dropdown(
                    id='meta_violin_select_vars',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in numerical_meta],
                    value=numerical_meta[3:4] if len(numerical_meta) > 2 else None,
                    multi=True
                ),
                dcc.RadioItems(
                    id='meta_violin_select_scale',
                    style=dict(width='200px',
                               marginTop=0,
                               marginBottom=20),
                    options=[dict(label='linear',
                                  value='linear'),
                             dict(label='log',
                                  value='log')],
                    value='linear'),
                html.Label('select grouping'),
                dcc.Dropdown(
                    id='meta_violin_select_group',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=categorical_meta[0]
                ),
                html.Label('select split'),
                dcc.Dropdown(
                    id='meta_violin_select_split',
                    style=dict(width='200px',
                               marginTop=0,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=None
                )]

    elif plot_type=='bar':

        return [html.Label('select grouping'),
                dcc.Dropdown(
                    id='meta_bar_select_group',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=categorical_meta[0]
                ),
                html.Label('select split'),
                dcc.Dropdown(
                    id='meta_bar_select_split',
                    style=dict(width='200px',
                               marginTop=0,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=None
                ),
                html.Label('other properties'),
                dcc.Checklist(
                    id='meta_bar_options',
                    style=dict(marginTop=0,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in ['normalized','stacked']],
                    values=[]
                    )]

###################################################################################################################
## meta plots
###################################################################################################################

# callback for meta plot
@app.callback(dash.dependencies.Output('meta_plot','children'),
              [dash.dependencies.Input('meta_plot_type','value')])

def render_meta_plot (plot_type):
    if plot_type=='scatter':
        return [dcc.Graph(id='meta_scatter_plot')]
    elif plot_type=='violin':
        return [dcc.Graph(id='meta_violin_plot')]
    elif plot_type=='bar':
        return [dcc.Graph(id='meta_bar_plot')]

#####################
# meta scatter plot #
#####################
    
@app.callback(
    dash.dependencies.Output(component_id='meta_scatter_plot',component_property='figure'),
    [dash.dependencies.Input(component_id='meta_scatter_select_x',component_property='value'),
     dash.dependencies.Input(component_id='meta_scatter_scale_x',component_property='value'),
     dash.dependencies.Input(component_id='meta_scatter_select_y',component_property='value'),
     dash.dependencies.Input(component_id='meta_scatter_scale_y',component_property='value'),
     dash.dependencies.Input(component_id='meta_scatter_select_color',component_property='value')])

def get_meta_scatterplot(xc, xscale, yc, yscale, col):

    print('meta scatterplot {0} vs {1}, color {2}'.format(yc,xc,col))
    
    if col in categorical_meta:

        # select color palette 
        colvals=meta[col].unique()
        nvals=len(colvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20

        # plot scatter for each category separately
        traces=[]
        for n,cv in enumerate(colvals):
            tr=go.Scatter(
                x=meta[meta[col]==cv][xc],
                y=meta[meta[col]==cv][yc],
                text=meta[meta[col]==cv].index,
                mode='markers',
                opacity=.7,
                marker=dict(
                    size=5,
                    color=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                    line={'width': .1, 'color': 'gray'}
                ),
                name=cv
            )
            traces.append(tr)

    else:
        
        # for numerical data, plot scatter all at once
        traces=[go.Scatter(
            x=meta[xc],
            y=meta[yc],
            text=meta.index,
            mode='markers',
            opacity=.7,
            marker=dict(
                size=5,
                color=meta[col],
                colorscale='Viridis',
                showscale=True),
            name=col
        )]
        
    return dict(data=traces,
                layout=go.Layout(
                    xaxis=dict(title=xc, type=xscale),
                    yaxis=dict(title=yc, type=yscale),
                    margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                    legend={'x':1.05,'y':1},
                    hovermode='closest'
                ))

#####################
# meta violin plot  #
#####################

@app.callback(
    dash.dependencies.Output(component_id='meta_violin_plot',component_property='figure'),
    [dash.dependencies.Input(component_id='meta_violin_select_vars',component_property='value'),
     dash.dependencies.Input(component_id='meta_violin_select_scale',component_property='value'),
     dash.dependencies.Input(component_id='meta_violin_select_group',component_property='value'),
     dash.dependencies.Input(component_id='meta_violin_select_split',component_property='value')])

def get_meta_violinplot(variables, scale, group, split):

    # select color palette
    if split is None:
        groupvals=meta[group].unique()
        nvals=len(groupvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20
    else:
        splitvals=meta[split].unique()
        nvals=len(splitvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20

    fig=tools.make_subplots(rows=len(variables),
                            cols=1,
                            specs=[[{}] for var in variables],
                            shared_xaxes=True,
                            vertical_spacing=.001)

    for nv, var in enumerate(variables):

        print('meta violinplot {0} against {1}, split {2}'.format(var,group,split))

        traces=[]
        if split is None:
            for n,cv in enumerate(groupvals):
                tr=go.Violin(
                    y=meta[meta[group]==cv][var],
                    name=cv,
                    fillcolor=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                    line=dict(color='gray',width=.5),
                    marker=dict(size=3,color='black'),
                    spanmode='hard',
                    showlegend=nv==0,
                    scalegroup=nv
                )
                fig.append_trace(tr, len(variables)-nv, 1)
        else:
            for n,sv in enumerate(splitvals):
                tr=go.Violin(x=meta[meta[split]==sv][group],
                             y=meta[meta[split]==sv][var],
                             name=sv,
                             offsetgroup=n,
                             scalegroup=nv,
                             showlegend=nv==0,
                             fillcolor=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                             line=dict(color='gray',width=.5),
                             marker=dict(size=3,color='black'),
                             spanmode='hard')
                fig.append_trace(tr, nv+1, 1)

        fig['layout']['yaxis'+str(nv+1)].update(title=var,type=scale)

    fig['layout'].update(xaxis=dict(title=group,tickangle=-45),
                         margin={'l': 50, 'b': 80, 't': 10, 'r': 10},
                         legend={'x':1.05,'y':1},
                         hovermode='closest')

    if split is not None:
        fig['layout'].update(violinmode='group')
    
    return fig

#####################
# meta bar plot  #
#####################

@app.callback(
    dash.dependencies.Output(component_id='meta_bar_plot',component_property='figure'),
    [dash.dependencies.Input(component_id='meta_bar_select_group',component_property='value'),
     dash.dependencies.Input(component_id='meta_bar_select_split',component_property='value'),
     dash.dependencies.Input(component_id='meta_bar_options',component_property='values')])

def get_meta_barplot(group, split, options):

    # select color palette
    if split is None:
        groupvals=meta[group].unique()
        nvals=len(groupvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20
    else:
        splitvals=meta[split].unique()
        nvals=len(splitvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20

    print('meta barplot of {0}, split by {1} (options {2})'.format(group,split, str(options)))

    if split is None:
        tally=meta.groupby(group).size()
        if 'normalized' in options:
            tally=tally.divide(tally.sum())
        traces=[]
        for n,gv in enumerate(tally.index):
            tr=go.Bar(x=[gv],
                      y=[tally[gv]],
                      name=gv,
                      marker=dict(line=dict(color='gray',width=.5),
                                  color=matplotlib.colors.rgb2hex(cm(n%20)[:3])))
            traces.append(tr)
    else:
        tally=meta.groupby([group,split]).size()
        if 'normalized' in options:
            tally=tally.divide(tally.sum(level=0).astype(float),level=0)
        traces=[]
        for n,sv in enumerate(splitvals):
            tr=go.Bar(
                x=tally.xs(sv,level=1).index,
                y=tally.xs(sv,level=1).values,
                name=sv,
                marker=dict(color=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                            line=dict(color='gray',width=.5))
            )
            traces.append(tr)
            
    return dict(data=traces,
                layout=go.Layout(
                    barmode='stack' if 'stacked' in options else 'group',
                    xaxis=dict(title=group,tickangle=-45),
                    yaxis=dict(title='cell frequency' if 'normalized' in options else 'cell number'),
                    margin={'l': 50, 'b': 100, 't': 10, 'r': 10},
                    legend={'x':1.05,'y':1},
                    hovermode='closest'))

###############################################################################################################
## expression plot controls
###############################################################################################################

# callback for expression plot controls
@app.callback(dash.dependencies.Output('expression_plot_controls','children'),
              [dash.dependencies.Input('expression_plot_type','value')])

def render_expression_controls (plot_type):

    if plot_type=='scatter':

        return [html.Label('select x axis and scaling'),
                dcc.Dropdown(
                    id='expression_scatter_select_x',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in numerical_meta],
                    value=meta.columns[0]
                ),
                dcc.RadioItems(
                    id='expression_scatter_scale_x',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label='linear',
                                  value='linear'),
                             dict(label='log',
                                  value='log')],
                    value='linear'),
                html.Label('select y axis and scaling'),
                dcc.Dropdown(
                    id='expression_scatter_select_y',
                    style=dict(width='200px'),
                    options=[dict(label=c,value=c) for c in numerical_meta],
                    value=meta.columns[1]
                ),
                dcc.RadioItems(
                    id='expression_scatter_scale_y',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label='linear',
                                  value='linear'),
                             dict(label='log',
                                  value='log')],
                    value='linear'),
                #html.Button('Render plot', id='expression-scatter-render')
                ]

    elif plot_type=='violin':

        return [html.Label('select grouping'),
                dcc.Dropdown(
                    id='expression_violin_select_group',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=categorical_meta[0]
                ),
                html.Label('select split'),
                dcc.Dropdown(
                    id='expression_violin_select_split',
                    style=dict(width='200px',
                               marginTop=0,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=None),
                #html.Button('Render plot', id='expression-violin-render')
                ]

###################################################################################################################
## expression plots
###################################################################################################################

# callback for expression plot
@app.callback(dash.dependencies.Output('expression_plot','children'),
              [dash.dependencies.Input('expression_plot_type','value')])

def render_expression_plot (plot_type):
    if plot_type=='scatter':
        return [dcc.Graph(id='expression_scatter_plot')]
    elif plot_type=='violin':
        return [dcc.Graph(id='expression_violin_plot')]

###########################
# expression scatter plot #
###########################
    
@app.callback(dash.dependencies.Output(component_id='expression_scatter_plot',component_property='figure'),
    #[dash.dependencies.Input(component_id='expression_scatter_render',component_property='n_clicks')],
    [dash.dependencies.Input(component_id='expression_scatter_select_x',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_scale_x',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_select_y',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_scale_y',component_property='value'),
     dash.dependencies.Input(component_id='expression_select_genes',component_property='value')])

def get_expression_scatterplot(xc, xscale, yc, yscale, genelist):

    ngenes=len(genelist)
    nrow=int(np.floor(np.sqrt(ngenes)))
    ncol=int(np.ceil(ngenes/nrow))
    fig=tools.make_subplots(rows=nrow,
                            cols=ncol,
                            specs=[[{} for c in range(ncol)] for r in range(nrow)],
                            shared_xaxes=True,
                            shared_yaxes=True,
                            vertical_spacing=.01,
                            horizontal_spacing=.01,
                            subplot_titles=genelist)

    # rescale all expression values to the same min and max
    maxval=max(1,DGE.loc[genelist].max().max())

    for ng, gene in enumerate(genelist):

        print('expression scatterplot {0} against {1} for {2}'.format(yc, xc, gene))

        # for numerical data, plot scatter all at once
        trace=go.Scatter(
            x=meta[xc],
            y=meta[yc],
            text=cells,
            mode='markers',
            opacity=.7,
            marker=dict(
                size=5,
                color=maxval*DGE.loc[gene]/max(1,DGE.loc[gene].max()),
                colorscale='Viridis',
                showscale=(ng==0)))

        nr=ng%nrow+1
        nc=ng//nrow+1
        fig.append_trace(trace, nr, nc)

    # fix axes labels and scales
    for k in dir(fig['layout']):
        if k.startswith('yaxis'):
            fig['layout'][k].update(title=yc,type=yscale)
        elif k.startswith('xaxis'):
            fig['layout'][k].update(title=xc,type=xscale)
            
    fig['layout'].update(margin={'l': 40, 'b': 40, 't': 40, 'r': 40},
                         showlegend=False,
                         hovermode='closest')

    return fig

###########################
# expression violin plot  #
###########################

@app.callback(
    dash.dependencies.Output(component_id='expression_violin_plot',component_property='figure'),
#    [dash.dependencies.Input(component_id='expression_violin_render',component_property='n_clicks')],
    [dash.dependencies.Input(component_id='expression_select_genes',component_property='value'),
     dash.dependencies.Input(component_id='expression_violin_select_group',component_property='value'),
     dash.dependencies.Input(component_id='expression_violin_select_split',component_property='value')])

def get_expression_violinplot(genelist, group, split):

    # select color palette
    if split is None:
        groupvals=meta[group].unique()
        nvals=len(groupvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20
    else:
        splitvals=meta[split].unique()
        nvals=len(splitvals)
        if nvals < 10:
            cm=plt.cm.tab10
        else:
            cm=plt.cm.tab20

    fig=tools.make_subplots(rows=len(genelist),
                            cols=1,
                            specs=[[{}] for gene in genelist],
                            shared_xaxes=True,
                            vertical_spacing=.001)

    for ng, gene in enumerate(genelist):

        print('expression violinplot {0} against {1}, split {2}'.format(gene,group,split))

        if split is None:
            for n,cv in enumerate(groupvals):
                tr=go.Violin(
                    y=DGE.loc[gene,meta[group]==cv],
                    name=cv,
                    fillcolor=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                    line=dict(color='gray',width=.5),
                    marker=dict(size=3,color='black'),
                    spanmode='hard',
                    showlegend=ng==0,
                    scalegroup=ng
                )
                fig.append_trace(tr, len(genelist)-ng, 1)
        else:
            for n,sv in enumerate(splitvals):
                tr=go.Violin(x=meta[meta[split]==sv][group],
                             y=DGE.loc[gene,meta[split]==sv],
                             name=sv,
                             offsetgroup=n,
                             scalegroup=nv,
                             showlegend=ng==0,
                             fillcolor=matplotlib.colors.rgb2hex(cm(n%20)[:3]),
                             line=dict(color='gray',width=.5),
                             marker=dict(size=3,color='black'),
                             spanmode='hard')
                fig.append_trace(tr, len(genelist)-ng, 1)

    for ng in range(len(genelist)):
        if len(genelist)-ng==1:
            fig['layout']['yaxis'].update(title=gene)
        else:
            fig['layout']['yaxis'+str(len(genelist)-ng)].update(title=gene)

    fig['layout'].update(xaxis=dict(title=group,tickangle=-45),
                         margin={'l': 50, 'b': 80, 't': 10, 'r': 10},
                         legend={'x':1.05,'y':1},
                         hovermode='closest')

    if split is not None:
        fig['layout'].update(violinmode='group')
    
    return fig
    

if __name__ == '__main__':
    app.run_server(debug=True)
                    
