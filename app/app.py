# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import dash
import scipy.io
import anndata
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.tools as tools

##############################################################################
## define colors
##############################################################################

#Set3=['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69',
#      '#fccde5','#d9d9d9','#bc80bd']

tab10=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2',
       '#7f7f7f', '#bcbd22', '#17becf']
tab20=["#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a","#d62728",
       "#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2","#f7b6d2",
       "#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"]
tab40=["#393b79","#5254a3","#6b6ecf","#9c9ede","#637939","#8ca252","#b5cf6b",
       "#cedb9c","#8c6d31","#bd9e39","#e7ba52","#e7cb94","#843c39","#ad494a",
       "#d6616b","#e7969c","#7b4173","#a55194","#ce6dbd","#de9ed6","#3182bd",
       "#6baed6","#9ecae1","#c6dbef","#e6550d","#fd8d3c","#fdae6b","#fdd0a2",
       "#31a354","#74c476","#a1d99b","#c7e9c0","#756bb1","#9e9ac8","#bcbddc",
       "#dadaeb","#636363","#969696","#bdbdbd","#d9d9d9"]

def get_cm(values):
    if len(values) < 10:
        return tab10
    elif len(values) < 20:
        return tab20
    else:
        return tab40

my_gradients=[[[0,'rgb(180,180,180)'], # grey to blue
               [1,'rgb(0,0,240)']],
              [[0,'rgb(180,180,180)'], # grey to red
               [1,'rgb(240,0,0)']],
              [[0,'rgb(180,180,180)'], # grey to green
               [1,'rgb(0,200,0)']],
              [[0,'rgb(180,180,180)'], # grey to black
               [1,'rgb(0,0,0)']]]


##############################################################################
## get data
##############################################################################

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

ad_file=os.path.join(datadir,'data.h5ad')
if os.path.isfile(ad_file):
    print('reading all data from '+ad_file)
    ad=anndata.read_h5ad(ad_file)
    meta=ad.obs
    DGE=ad.to_df().T

# separate numerical and categorical columns for later
numerical_meta=[]
categorical_meta=[]
for col in meta.columns:
    if pd.api.types.is_numeric_dtype(meta[col]):
        numerical_meta.append(col)
    else:
        categorical_meta.append(col)

genes=DGE.index
cells=DGE.columns

if False:
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

###############################################################################
## app main layout
###############################################################################

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

###############################################################################
## main tab contents
###############################################################################

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

                html.Label('select cell sample'),
                dcc.Dropdown(
                    id='meta_select_cell_sample',
                    style=dict(width='200px'),
                    options=[dict(label=g,value=g) for g in [100,1000,5000] if g < meta.shape[0]]+[dict(label='all',
                                                                                                        value='all')],
                    value=min(1000,meta.shape[0]),
                ),

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
                                      value='violin'),
                                 dict(label='dot',
                                      value='dot')],
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

                html.Label('select cell sample'),
                dcc.Dropdown(
                    id='expression_select_cell_sample',
                    style=dict(width='200px'),
                    options=[dict(label=g,value=g) for g in [100,1000,5000] if g < meta.shape[0]]+[dict(label='all',
                                                                                                        value='all')],
                    value=min(1000,meta.shape[0]),
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

###############################################################################
## meta plot controls
###############################################################################

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

@app.callback(dash.dependencies.Output('meta_select_cell_sample','style'),
              [dash.dependencies.Input('meta_plot_type','value')])

def toggle_meta_sample (plot_type):
    if plot_type=='bar':
        return dict(display='none')
    else:
        return dict(display='block')

###############################################################################
## meta plots
###############################################################################

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
     dash.dependencies.Input(component_id='meta_scatter_select_color',component_property='value'),
     dash.dependencies.Input(component_id='meta_select_cell_sample',component_property='value')])

def get_meta_scatterplot(xc, xscale, yc, yscale, col, sample_size):

    if sample_size!='all':
        meta_here=meta.sample(n=sample_size)
    else:
        meta_here=meta

    if col in categorical_meta:

        # select color palette 
        colvals=meta_here[col].unique()
        cm=get_cm(colvals)

        # plot scatter for each category separately
        traces=[]
        for n,cv in enumerate(colvals):
            tr=go.Scatter(
                x=meta_here[meta_here[col]==cv][xc],
                y=meta_here[meta_here[col]==cv][yc],
                text=meta_here[meta_here[col]==cv].index,
                mode='markers',
                opacity=.7,
                marker=dict(
                    size=5,
                    color=cm[n%40],
                    line={'width': .1, 'color': 'gray'}
                ),
                showlegend=False,
                name=cv
            )
            traces.append(tr)
            # extra trace for legend with bigger marker
            traces.append(go.Scatter(x=[None],
                                     y=[None],
                                     mode='markers',
                                     marker=dict(size=20,
                                                 color=cm[n%40],
                                                 line={'width': .1, 'color':'gray'}),
                                     showlegend=True,
                                     visible='legendonly',
                                     name=cv))

    else:
        
        # for numerical data, plot scatter all at once
        traces=[go.Scatter(
            x=meta_here[xc],
            y=meta_here[yc],
            text=meta_here.index,
            mode='markers',
            opacity=.7,
            marker=dict(
                size=5,
                color=meta_here[col],
                colorscale='Viridis',
                colorbar=dict(title=col,
                              titleside='right'),
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
     dash.dependencies.Input(component_id='meta_violin_select_group',component_property='value'),
     dash.dependencies.Input(component_id='meta_violin_select_split',component_property='value'),
     dash.dependencies.Input(component_id='meta_select_cell_sample',component_property='value')])

def get_meta_violinplot(variables, group, split, sample_size):

    if sample_size!='all':
        meta_here=meta.sample(n=sample_size)
    else:
        meta_here=meta

    # select color palette
    if split is None:
        groupvals=meta_here[group].unique()
        cm=get_cm(groupvals)
    else:
        splitvals=meta_here[split].unique()
        cm=get_cm(splitvals)
            
    nvar=len(variables)

    fig=tools.make_subplots(rows=nvar,
                            cols=1,
                            specs=[[{}] for var in variables],
                            shared_xaxes=True,
                            vertical_spacing=.001)

    sg=0
    for nv, var in enumerate(variables):

        traces=[]
        if split is None:
            for n,cv in enumerate(groupvals):
                y=meta_here[meta_here[group]==cv][var]
                tr=go.Violin(
                    y=y,
                    name=cv,
                    fillcolor=cm[n%40],
                    line=dict(color='gray',width=.5),
                    marker=dict(size=1),
                    spanmode='manual',
                    span=[y.min(),None],
                    showlegend=nv==0,
                    scalegroup=sg
                )
                fig.append_trace(tr, nvar-nv, 1)
                sg+=1
        else:
            for n,sv in enumerate(splitvals):
                y=meta_here[meta_here[split]==sv][var]
                tr=go.Violin(x=meta_here[meta_here[split]==sv][group],
                             y=y,
                             name=sv,
                             offsetgroup=n,
                             scalegroup=sg,
                             showlegend=nv==0,
                             fillcolor=cm[n%40],
                             line=dict(color='gray',width=.5),
                             marker=dict(size=1),
                             spanmode='manual',
                             span=[y.min(),None])
                fig.append_trace(tr, nvar-nv, 1)
                sg+=1

    for nv,var in enumerate(variables):
        if len(variables)-nv==1:
            fig['layout']['yaxis'].update(title=var)
        else:
            fig['layout']['yaxis'+str(nvar-nv)].update(title=var)

    fig['layout'].update(xaxis=dict(title=group,tickangle=-45),
                         margin={'l': 50, 'b': 80, 't': 10, 'r': 10},
                         legend={'x':1.05,'y':1},
                         hovermode='closest')

    if split is not None:
        fig['layout'].update(violinmode='group')
    
    return fig

#####################
# meta bar plot     #
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
        cm=get_cm(groupvals)
    else:
        splitvals=meta[split].unique()
        cm=get_cm(splitvals)

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
                                  color=cm[n%40]))
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
                marker=dict(color=cm[n%40],
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

###############################################################################
## expression plot controls
###############################################################################

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

    elif plot_type=='dot':
            
        return [html.Label('select grouping'),
                dcc.Dropdown(
                    id='expression_dot_select_group',
                    style=dict(width='200px',
                               marginTop=10,
                               marginBottom=20),
                    options=[dict(label=c,value=c) for c in categorical_meta],
                    value=categorical_meta[0]
                ),
                html.Label('select split'),
                dcc.Dropdown(
                    id='expression_dot_select_split',
                    style=dict(width='200px',
                               marginTop=0,
                              marginBottom=20),
                    # get at most 4 splits otherwise it gets to messy
                    options=[dict(label=c,value=c) for c in categorical_meta if len(meta[col].unique()) < 4],
                    value=None),
                #html.Button('Render plot', id='expression-dot-render')
                ]

@app.callback(dash.dependencies.Output('expression_select_cell_sample','style'),
              [dash.dependencies.Input('expression_plot_type','value')])

def toggle_expression_sample (plot_type):
    if plot_type=='bar':
        return dict(display='none')
    else:
        return dict(display='block')

###############################################################################
## expression plots
###############################################################################

# callback for expression plot
@app.callback(dash.dependencies.Output('expression_plot','children'),
              [dash.dependencies.Input('expression_select_genes','value'),
               dash.dependencies.Input('expression_plot_type','value')])

def render_expression_plot (genelist, plot_type):
    if len(genelist)==0:
        return []
    if plot_type=='scatter':
        return [dcc.Graph(id='expression_scatter_plot')]
    elif plot_type=='violin':
        return [dcc.Graph(id='expression_violin_plot')]
    elif plot_type=='dot':
        return [dcc.Graph(id='expression_dot_plot')]

###########################
# expression scatter plot #
###########################
    
@app.callback(dash.dependencies.Output(component_id='expression_scatter_plot',component_property='figure'),
    #[dash.dependencies.Input(component_id='expression_scatter_render',component_property='n_clicks')],
    [dash.dependencies.Input(component_id='expression_scatter_select_x',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_scale_x',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_select_y',component_property='value'),
     dash.dependencies.Input(component_id='expression_scatter_scale_y',component_property='value'),
     dash.dependencies.Input(component_id='expression_select_genes',component_property='value'),
     dash.dependencies.Input(component_id='expression_select_cell_sample',component_property='value')])

def get_expression_scatterplot(xc, xscale, yc, yscale, genelist, sample_size):

    if sample_size!='all':
        DGE_here=DGE.sample(n=sample_size,axis=1)
        meta_here=meta.loc[DGE_here.columns]
    else:
        DGE_here=DGE
        meta_here=meta

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
    maxval=max(1,DGE_here.loc[genelist].max().max())

    for ng, gene in enumerate(genelist):
        
        # for numerical data, plot scatter all at once
        trace=go.Scatter(
            x=meta_here[xc],
            y=meta_here[yc],
            text=cells,
            mode='markers',
            marker=dict(
                size=3,
                color=maxval*DGE_here.loc[gene]/max(1,DGE_here.loc[gene].max()),
                colorscale='Viridis',
                colorbar=dict(title='expression',
                              titleside='right'),
                showscale=(ng==0)))

        nc=ng%ncol+1
        nr=ng//ncol+1
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
     dash.dependencies.Input(component_id='expression_select_cell_sample',component_property='value'),
     dash.dependencies.Input(component_id='expression_violin_select_group',component_property='value'),
     dash.dependencies.Input(component_id='expression_violin_select_split',component_property='value')])

def get_expression_violinplot(genelist, sample_size, group, split):

    if sample_size!='all':
        DGE_here=DGE.sample(n=sample_size,axis=1)
        meta_here=meta.loc[DGE_here.columns]
    else:
        DGE_here=DGE
        meta_here=meta

    # select color palette
    if split is None:
        groupvals=meta_here[group].unique()
        cm=get_cm(groupvals)
    else:
        splitvals=meta_here[split].unique()
        cm=get_cm(splitvals)

    ngenes=len(genelist)
    
    fig=tools.make_subplots(rows=ngenes,
                            cols=1,
                            specs=[[{}] for gene in genelist],
                            shared_xaxes=True,
                            vertical_spacing=.01)
    sg=0
    for ng, gene in enumerate(genelist):

        if split is None:
            for n,cv in enumerate(groupvals):
                y=DGE_here.loc[gene,meta_here[group]==cv]
                tr=go.Violin(
                    y=y,
                    name=cv,
                    fillcolor=cm[n%40],
                    line=dict(color='gray',width=.5),
                    marker=dict(size=1),
                    spanmode='manual',
                    span=[y.min(),None],
                    showlegend=ng==0,
                    scalegroup=sg
                )
                fig.append_trace(tr, ngenes-ng, 1)
                sg+=1
        else:
            for n,sv in enumerate(splitvals):
                y=DGE_here.loc[gene,meta_here[split]==sv]
                tr=go.Violin(x=meta_here[meta_here[split]==sv][group],
                             y=y,
                             name=sv,
                             offsetgroup=n,
                             scalegroup=sg,
                             showlegend=ng==0,
                             fillcolor=cm[n%40],
                             line=dict(color='gray',width=.5),
                             spanmode='manual',
                             marker=dict(size=1),
                             span=[y.min(),None])
                fig.append_trace(tr, ngenes-ng, 1)
                sg+=1

    for ng,gene in enumerate(genelist):
        if ngenes-ng==1:
            fig['layout']['yaxis'].update(title=gene)
        else:
            fig['layout']['yaxis'+str(ngenes-ng)].update(title=gene)

    fig['layout'].update(xaxis=dict(title=group,tickangle=-45),
                         margin={'l': 50, 'b': 80, 't': 10, 'r': 10},
                         legend={'x':1.05,'y':1},
                         hovermode='closest')

    if split is not None:
        fig['layout'].update(violinmode='group')
    
    return fig
    
##############################
# expression dot/bubble plot #
##############################

@app.callback(
    dash.dependencies.Output(component_id='expression_dot_plot',component_property='figure'),
#    [dash.dependencies.Input(component_id='expression_dot_render',component_property='n_clicks')],
    [dash.dependencies.Input(component_id='expression_select_genes',component_property='value'),
     dash.dependencies.Input(component_id='expression_dot_select_group',component_property='value'),
     dash.dependencies.Input(component_id='expression_dot_select_split',component_property='value')])

def get_expression_dotplot(genelist, group, split):
    
    groupvals=meta[group].unique()
    ngroup=len(groupvals)
    ngenes=len(genelist)
    
    if split is None:
        means=np.array([DGE.loc[genelist,meta[group]==cv].mean(axis=1).fillna(0) for cv in groupvals])
        sizes=np.array([(DGE.loc[genelist,meta[group]==cv] > 0).mean(axis=1).fillna(0) for cv in groupvals])
        xv,yv=np.meshgrid(np.arange(ngenes),np.arange(ngroup))        
        traces=[go.Scatter(x=xv.ravel(),
                           y=yv.ravel(),
                           mode='markers',
                           showlegend=False,
                           marker=dict(size=20*sizes.ravel(),
                                       color=means.ravel(),
                                       colorscale=my_gradients[0]))]
        # extra invisible traces for legend
        traces+=[go.Scatter(x=[0],y=[0], mode='markers',
                            marker=dict(size=s,color='black'),
                            name=p,visible='legendonly',
                            legendgroup='size',showlegend=True) 
                 for s,p in zip([0,10,20],['0% cells','50% cells','100% cells'])]
        traces+=[go.Scatter(x=[0],y=[0], mode='markers',
                            marker=dict(size=10,color=my_gradients[0][c][1],symbol='square'),
                            name=p,visible='legendonly',
                            legendgroup='color',showlegend=True) 
                 for c,p in zip([0,1],['low','high'])]
        layout=go.Layout(xaxis=go.layout.XAxis(showgrid=False,
                                                zeroline=False,
                                                showline=False,
                                                tickvals=np.arange(ngenes),
                                                ticktext=genelist,
                                                tickangle=-90),
                          yaxis=go.layout.YAxis(showgrid=False,
                                                zeroline=False,
                                                showline=False,
                                                tickvals=np.arange(ngroup),
                                                ticktext=groupvals),
                          margin={'l': 100, 'b': 100, 't': 40, 'r': 40},
                          showlegend=True,
                          hovermode='closest')
    else:
        splitvals=meta[split].unique()
        nsplit=len(splitvals)
        traces=[]
        for n,sv in enumerate(splitvals):
            means=np.array([DGE.loc[genelist,(meta[split]==sv) & 
                                    (meta[group]==cv)].mean(axis=1).fillna(0) for cv in groupvals])
            sizes=np.array([(DGE.loc[genelist,(meta[split]==sv) & 
                                     (meta[group]==cv)] > 0).mean(axis=1).fillna(0) for cv in groupvals])
            xv,yv=np.meshgrid(np.arange(ngenes),(nsplit+1)*np.arange(ngroup))
            tr=go.Scatter(x=xv.ravel(),
                          y=yv.ravel()+n,
                          mode='markers',
                          showlegend=False,
                          marker=dict(size=10*sizes.ravel(),
                                      color=means.ravel(),
                                      colorscale=my_gradients[n]))
            traces.append(tr)
            
        # extra invisible traces for legend
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=s,color='black'),
                            name=p,visible='legendonly',
                            legendgroup='size',showlegend=True) 
                 for s,p in zip([0,5,10],['0% cells', '50% cells','100% cells'])]
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=10,symbol='square',
                                        color=my_gradients[-1][c][1]),
                            name=p,visible='legendonly',
                            legendgroup='color',showlegend=True) 
                 for c,p in zip([0,1],['low','high'])]
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=10,symbol='square',
                                        color=my_gradients[n][1][1]),
                            name=sv,visible='legendonly',
                            legendgroup='split',showlegend=True) 
                 for n,sv in enumerate(splitvals)]

        layout=go.Layout(xaxis=go.layout.XAxis(showgrid=False,
                                               zeroline=False,
                                               showline=False,
                                               tickvals=np.arange(ngenes),
                                               ticktext=genelist,
                                               tickangle=-90),
                         yaxis=go.layout.YAxis(showgrid=False,
                                               zeroline=False,
                                               showline=False,
                                               tickvals=(nsplit+1)*np.arange(ngroup)+(nsplit-1)/2.,
                                               ticktext=groupvals),
                          margin={'l': 100, 'b': 100, 't': 40, 'r': 40},
                          showlegend=True,
                          hovermode='closest')
        
    return dict(data=traces,layout=layout)


if __name__ == '__main__':
    app.run_server(host='0.0.0.0',debug=True)
                    
