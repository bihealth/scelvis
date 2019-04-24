# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import dash
import urllib.parse
import anndata
import dash_core_components as dcc
import dash_bootstrap_components as dbc
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
               [1,'rgb(31,119,180)']],
              [[0,'rgb(180,180,180)'], # grey to red
               [1,'rgb(214,39,40)']],
              [[0,'rgb(180,180,180)'], # grey to green
               [1,'rgb(44,160,44)']],
              [[0,'rgb(180,180,180)'], # grey to black
               [1,'rgb(0,0,0)']]]


##############################################################################
## get data
##############################################################################

# load data using an environmental variable that points to a folder where meta and expression data should be
# this can be replaced by some sort of HTTP request directly to SODAR

datadir=os.getenv('DASH_DATADIR')
about_file=os.path.join(datadir,'about.md')
ad_file=os.path.join(datadir,'data.h5ad')
if datadir is None or not os.path.isfile(ad_file):
    datadir='datasets/pbmc'
    about_file=os.path.join(datadir,'about.md')
    ad_file=os.path.join(datadir,'data.h5ad')

# get about file (markdown) if exists
if os.path.isfile(about_file):
    print('reading about file from '+about_file)
    about_md=open(about_file).readlines()
else:
    about_md=["""this is a Dash app for data from {0}""".format(datadir),
              "",
              """no further information available at this point"""]

if os.path.isfile(ad_file):
    print('reading all data from '+ad_file)
    ad=anndata.read_h5ad(ad_file)
    coords={}
    for k in ['X_tsne','X_umap']:
        if k in ad.obsm.keys():
            coords[k]=pd.DataFrame(ad.obsm[k],
                                   index=ad.obs.index,
                                   columns=[k[2:].upper()+str(n+1) for n in range(ad.obsm[k].shape[1])])
    meta=pd.concat(coords.values(),axis=1).join(ad.obs)
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
    marker_file=os.path.join(datadir,'markers.tsv.gz')
    if os.path.isfile(marker_file):
        print('reading markers from '+markers)
        markers=pd.read_csv(marker_file,header=0,sep='\t')
    else:
        markers=None

###############################################################################
## meta plot controls
###############################################################################

meta_plot_subcontrols=dict(scatter=[html.Label('select x axis'),
                                    dcc.Dropdown(
                                        id='meta_scatter_select_x',
                                        options=[dict(label=c,value=c) for c in numerical_meta],
                                        value=meta.columns[0]
                                    ),
                                    html.Label('select y axis'),
                                    dcc.Dropdown(
                                        id='meta_scatter_select_y',
                                        options=[dict(label=c,value=c) for c in numerical_meta],
                                        value=meta.columns[1]
                                    ),
                                    
                                    html.Label('select coloring'),
                                    dcc.Dropdown(
                                        id='meta_scatter_select_color',
                                        options=[dict(label=c,value=c) for c in meta.columns],
                                        value=categorical_meta[0]
                                    )],
                           violin=[html.Label('select variable(s) and scaling'),
                                   dcc.Dropdown(
                                       id='meta_violin_select_vars',
                                       options=[dict(label=c,value=c) for c in numerical_meta],
                                       value=None,
                                       multi=True
                                   ),
                                   html.Label('select grouping'),
                                   dcc.Dropdown(
                                       id='meta_violin_select_group',
                                       options=[dict(label=c,value=c) for c in categorical_meta],
                                       value=categorical_meta[0]
                                   ),
                                   html.Label('select split'),
                                   dcc.Dropdown(
                                       id='meta_violin_select_split',
                                       options=[dict(label=c,value=c) for c in categorical_meta],
                                       value=None
                                   )],
                           bar=[html.Label('select grouping'),
                                dcc.Dropdown(
                                    id='meta_bar_select_group',
                                    options=[dict(label=c,value=c) for c in categorical_meta],
                                    value=categorical_meta[0]
                                ),
                                html.Label('select split'),
                                dcc.Dropdown(
                                    id='meta_bar_select_split',
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
                                )])

meta_plot_controls=[dbc.Row([
    
    # first column: plot controls
    dbc.Col([
        
        # top: determine plot type
        dbc.Container([
            html.Label('select plot type'),
            dcc.RadioItems(
                id='meta_plot_type',
                options=[dict(label='scatter plot',value='scatter'),
                         dict(label='violin plot',value='violin'),
                         dict(label='bar plot',value='bar')],
                value='scatter',
                labelStyle={'display': 'block'}),

            html.Hr(),

            dbc.Container(id='meta_plot_subcontrols'),
            
            html.Hr(),
            
            html.Label('select cell sample'),
            dcc.Dropdown(
                id='meta_select_cell_sample',
                options=[dict(label=g,value=g) for g in [100,1000,5000]
                         if g < meta.shape[0]]+[dict(label='all',value='all')],
                value=min(1000,meta.shape[0]),
                disabled=False
            )]),
        ],className='col-3'),
    
    # second column: plot (depends on plot type)
    dbc.Col([dcc.Loading(id='meta_plot',
                         type='circle')],
            className='col-9')])
]

###############################################################################
## expression plot controls
###############################################################################

# plot-type-specific controls
expression_plot_subcontrols=dict(scatter=[html.Label('select x axis'),
                                          dcc.Dropdown(
                                              id='expression_scatter_select_x',
                                              options=[dict(label=c,value=c) for c in numerical_meta],
                                              value=meta.columns[0]
                                          ),
                                          html.Label('select y axis'),
                                          dcc.Dropdown(
                                              id='expression_scatter_select_y',
                                              options=[dict(label=c,value=c) for c in numerical_meta],
                                              value=meta.columns[1]
                                          )],
                                 violin=[html.Label('select grouping'),
                                         dcc.Dropdown(
                                             id='expression_violin_select_group',
                                             options=[dict(label=c,value=c) for c in categorical_meta],
                                             value=categorical_meta[0]
                                         ),
                                         html.Label('select split'),
                                         dcc.Dropdown(
                                             id='expression_violin_select_split',
                                             options=[dict(label=c,value=c) for c in categorical_meta],
                                             value=None)],
                                 dot=[html.Label('select grouping'),
                                      dcc.Dropdown(
                                          id='expression_dot_select_group',
                                          options=[dict(label=c,value=c) for c in categorical_meta],
                                          value=categorical_meta[0]
                                      ),
                                      html.Label('select split'),
                                      dcc.Dropdown(
                                          id='expression_dot_select_split',
                                          options=[dict(label=c,value=c) for c in categorical_meta
                                                   if len(meta[c].unique()) < 4],
                                          value=None)]
                                 )


expression_plot_controls=[dbc.Row([
    
    # first column: plot controls
    
    dbc.Col([
        
        # top: determine plot type
        dbc.Container([

            html.Label('select plot type'),
            
            dcc.RadioItems(
                id='expression_plot_type',
                options=[dict(label='scatter plot',value='scatter'),
                         dict(label='violin plot',value='violin'),
                         dict(label='dot plot',value='dot')],
                value='scatter',
                labelStyle={'display': 'block'}),
            
            html.Hr(),
            
            dbc.Container(id='expression_plot_subcontrols'),
            
            
            html.Hr(),
            
            html.Label('select gene(s)'),
            dcc.Dropdown(
                id='expression_select_genes',
                options=[dict(label=g,value=g) for g in genes],
                value=[],
                multi=True
            ),
            
            html.Label('select cell sample'),
            dcc.Dropdown(
                id='expression_select_cell_sample',
                options=[dict(label=g,value=g) for g in [100,1000,5000]
                         if g < meta.shape[0]]+[dict(label='all',value='all')],
                value=min(1000,meta.shape[0]),
                disabled=False)
        ]),
    ],className='col-3'),
    
    # second column: plot (depends on plot type)
    dbc.Col([dcc.Loading(id='expression_plot',
                         type='circle')],
            className='col-9')])
]

                                      
###############################################################################
## app layout
###############################################################################

app=dash.Dash(datadir,
              external_stylesheets=[dbc.themes.BOOTSTRAP,
                                    dbc.themes.CERULEAN])
app.config['suppress_callback_exceptions']=True

# use external css (should be fixed)
#app.css.append_css({
#    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
#})
# Loading screen CSS
#app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})


# specify layout: header, and then a Div with 3 main tabs
app.layout = html.Div([

    ############################################################
    ## header
    ############################################################
    html.H1('interactive single-cell data analysis and visualization'),
    html.Hr(),

    dcc.Tabs(id='main_tabs',
             value='meta',
             style=dict(width='40%',height='5%'),
             children=[
                 
                 # "about" tab: displays markdown
                 dcc.Tab(label='about this dataset',
                         value='about',
                         children=[
                             dcc.Markdown('\n'.join(about_md))
                         ]),

                 # "cell  annotation" tab: scatter, violin or bar plots
                 dcc.Tab(label='explore cell annotation',
                         value='meta',
                         children=meta_plot_controls
                 ),

                 # "gene expression" tab: scatter, violin or dot plots
                 dcc.Tab(label='explore gene expression',
                         value='expression',
                         children=expression_plot_controls
                 ),
             ])],
                      style=dict(marginRight=50,
                                 marginLeft=50,
                                 marginBottom=50,
                                 marginTop=50))
            

###############################################################################
## meta plots
###############################################################################

# callback for meta plot controls
@app.callback([dash.dependencies.Output('meta_plot_subcontrols','children')],
              [dash.dependencies.Input('meta_plot_type','value')])
def get_meta_plot_subcontrols(plot_type):
    return [meta_plot_subcontrols[plot_type]]

# callback for meta plot (disable sampling for dot plot)
@app.callback([dash.dependencies.Output('meta_plot','children'),
               dash.dependencies.Output('meta_select_cell_sample','disabled')],
              [dash.dependencies.Input('meta_plot_type','value')])

def render_meta_plot (plot_type):
    if plot_type=='scatter':
        return [dcc.Graph(id='meta_scatter_plot'),
                html.A('download data for this plot',
                       id='meta_scatter_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],False
    elif plot_type=='violin':
        return [dcc.Graph(id='meta_violin_plot'),
                html.A('download data for this plot',
                       id='meta_violin_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],False
    elif plot_type=='bar':
        return [dcc.Graph(id='meta_bar_plot'),
                html.A('download data for this plot',
                       id='meta_bar_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],True

@app.callback(dash.dependencies.Output('meta_bar_options','options'),
              [dash.dependencies.Input('meta_bar_select_split','value')])

def toggle_meta_bar_options(split):
    if split is None:
        return [dict(label='normalized',value='normalized')]
    else:
        return [dict(label=c,value=c) for c in ['normalized','stacked']]
        


#####################
# meta scatter plot #
#####################
    
@app.callback(
    [dash.dependencies.Output('meta_scatter_plot','figure'),
     dash.dependencies.Output('meta_scatter_download','href')],
    [dash.dependencies.Input('meta_scatter_select_x','value'),
     dash.dependencies.Input('meta_scatter_select_y','value'),
     dash.dependencies.Input('meta_scatter_select_color','value'),
     dash.dependencies.Input('meta_select_cell_sample','value')])

def get_meta_scatterplot(xc, yc, col, sample_size):

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
            tr=go.Scattergl(
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
        traces=[go.Scattergl(
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

    data=meta[[xc,yc,col]]
    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data.to_csv(index=True, header=True, encoding='utf-8'))

    fig=dict(data=traces,
             layout=go.Layout(
                 xaxis=dict(title=xc),
                 yaxis=dict(title=yc),
                 margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                 legend={'x':1.05,'y':1},
                 hovermode='closest'
             ))
    return fig, csv_string

#####################
# meta violin plot  #
#####################

@app.callback(
    [dash.dependencies.Output('meta_violin_plot','figure'),
     dash.dependencies.Output('meta_violin_download','href')],
    [dash.dependencies.Input('meta_violin_select_vars','value'),
     dash.dependencies.Input('meta_violin_select_group','value'),
     dash.dependencies.Input('meta_violin_select_split','value'),
     dash.dependencies.Input('meta_select_cell_sample','value')])

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
            data=meta[variables+[group]]
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
            data=meta[variables+[group,split]]

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

    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data.to_csv(index=True, header=True, encoding='utf-8'))

    return fig,csv_string

#####################
# meta bar plot     #
#####################

@app.callback(
    [dash.dependencies.Output('meta_bar_plot','figure'),
     dash.dependencies.Output('meta_bar_download','href')],
    [dash.dependencies.Input('meta_bar_select_group','value'),
     dash.dependencies.Input('meta_bar_select_split','value'),
     dash.dependencies.Input('meta_bar_options','values')])

def get_meta_barplot(group, split, options):

    if split is None:
        groupvals=meta[group].unique()
        cm=get_cm(groupvals)
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
        splitvals=meta[split].cat.categories
        cm=get_cm(splitvals)

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
            
    fig=dict(data=traces,
             layout=go.Layout(
                 barmode='stack' if 'stacked' in options else 'group',
                 xaxis=dict(title=group,tickangle=-45),
                 yaxis=dict(title='cell frequency' if 'normalized' in options else 'cell number'),
                 margin={'l': 50, 'b': 100, 't': 10, 'r': 10},
                 legend={'x':1.05,'y':1},
                 hovermode='closest'))

    data=tally
    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data.to_csv(index=True, header=True, encoding='utf-8'))

    return fig,csv_string
    

    

###############################################################################
## expression plots
###############################################################################

# callback for expression plot controls
@app.callback([dash.dependencies.Output('expression_plot_subcontrols','children')],
              [dash.dependencies.Input('expression_plot_type','value')])
def get_expression_plot_subcontrols(plot_type):
    return [expression_plot_subcontrols[plot_type]]


# callback for expression plot (disable sampling for dot plot)
@app.callback([dash.dependencies.Output('expression_plot','children'),
               dash.dependencies.Output('expression_select_cell_sample','disabled')],
              [dash.dependencies.Input('expression_select_genes','value'),
               dash.dependencies.Input('expression_plot_type','value')])

def render_expression_plot (genelist, plot_type):
    if len(genelist)==0:
        return [],False
    if plot_type=='scatter':
        return [dcc.Graph(id='expression_scatter_plot'),
                html.A('download data for this plot',
                       id='expression_scatter_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],False,
    elif plot_type=='violin':
        return [dcc.Graph(id='expression_violin_plot'),
                html.A('download data for this plot',
                       id='expression_violin_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],False
    elif plot_type=='dot':
        return [dcc.Graph(id='expression_dot_plot'),
                html.A('download data for this plot',
                       id='expression_dot_download',
                       download='plot_data.csv',
                       href="",
                       target="_blank")],True

# callback for gene selection using Input
#@app.callback(dash.dependencies.Output('expression_select_genes','value'),
#              [dash.dependencies.Input('expression_input_genes','n_submit'),
#               dash.dependencies.Input('expression_input_genes','value')])
#
#def update_gene_selection(n_submit, text):
#
#    geneset=set()
#    for g in re.split('[,\s]',text):
#        if g in genes:
#            geneset.add(g)
#
#    return list(geneset)

###########################
# expression scatter plot #
###########################
    
@app.callback([dash.dependencies.Output('expression_scatter_plot','figure'),
               dash.dependencies.Output('expression_scatter_download','href')],
    [dash.dependencies.Input('expression_scatter_select_x','value'),
     dash.dependencies.Input('expression_scatter_select_y','value'),
     dash.dependencies.Input('expression_select_genes','value'),
     dash.dependencies.Input('expression_select_cell_sample','value')])

def get_expression_scatterplot(xc, yc, genelist, sample_size):

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
            fig['layout'][k].update(title=yc)
        elif k.startswith('xaxis'):
            fig['layout'][k].update(title=xc)
            
    fig['layout'].update(margin={'l': 40, 'b': 40, 't': 40, 'r': 40},
                         showlegend=False,
                         hovermode='closest')

    data=DGE.loc[genelist].T.to_csv(index=True, header=True, encoding='utf-8')
    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data)

    return fig,csv_string

###########################
# expression violin plot  #
###########################

@app.callback(
    [dash.dependencies.Output('expression_violin_plot','figure'),
     dash.dependencies.Output('expression_violin_download','href')],
    [dash.dependencies.Input('expression_select_genes','value'),
     dash.dependencies.Input('expression_select_cell_sample','value'),
     dash.dependencies.Input('expression_violin_select_group','value'),
     dash.dependencies.Input('expression_violin_select_split','value')])

def get_expression_violinplot(genelist, sample_size, group, split):

    if sample_size!='all':
        DGE_here=DGE.sample(n=sample_size,axis=1)
        meta_here=meta.loc[DGE_here.columns]
    else:
        DGE_here=DGE
        meta_here=meta

    # select color palette
    if split is None:
        groupvals=meta_here[group].cat.categories
        cm=get_cm(groupvals)
    else:
        splitvals=meta_here[split].cat.categories
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
            data=DGE.loc[genelist].T.join(meta[group])
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
            data=DGE.loc[genelist].T.join(meta[[group,split]])

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
    
    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data.to_csv(index=True, header=True, encoding='utf-8'))

    return fig,csv_string
    
##############################
# expression dot/bubble plot #
##############################

@app.callback(
    [dash.dependencies.Output('expression_dot_plot','figure'),
     dash.dependencies.Output('expression_dot_download','href')],
    [dash.dependencies.Input('expression_select_genes','value'),
     dash.dependencies.Input('expression_dot_select_group','value'),
     dash.dependencies.Input('expression_dot_select_split','value')])

def get_expression_dotplot(genelist, group, split):
    
    groupvals=meta[group].cat.categories
    ngroup=len(groupvals)
    ngenes=len(genelist)

    def my_agg (x):
        return pd.DataFrame([x.mean(),
                             (x > 0).mean()],
                            index=['expression','pct_cells'])
                             
    
    if split is None:
        data=(DGE.loc[genelist].T
              .join(meta[group])
              .groupby(group)
              .apply(my_agg).unstack(level=1))
        means=data.xs('expression',axis=1,level=1).values
        sizes=data.xs('pct_cells',axis=1,level=1).values
        xv,yv=np.meshgrid(np.arange(ngenes),np.arange(ngroup))        
        traces=[go.Scatter(x=xv.ravel(),
                           y=yv.ravel(),
                           mode='markers',
                           showlegend=False,
                           opacity=1,
                           marker=dict(size=20*np.sqrt(sizes.ravel()),
                                       color=means.ravel(),
                                       colorscale=my_gradients[0]))]
        # extra invisible traces for legend
        traces+=[go.Scatter(x=[0],y=[0], mode='markers',
                            marker=dict(size=20*np.sqrt(s),color='black'),
                            name=p,visible='legendonly',
                            legendgroup='size',showlegend=True) 
                 for s,p in zip([0,.5,1],['0% cells','50% cells','100% cells'])]
        traces+=[go.Scatter(x=[0],y=[0], mode='markers',
                            marker=dict(size=20,color=my_gradients[0][c][1],symbol='square'),
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
        splitvals=meta[split].cat.categories
        nsplit=len(splitvals)
        traces=[]
        data=(DGE.loc[genelist].T
              .join(meta[[group,split]])
              .groupby([group,split])
              .apply(my_agg).unstack(level=2))
        for n,sv in enumerate(splitvals):
            means=data.xs('expression',axis=1,level=1).xs(sv,axis=0,level=1).values
            sizes=data.xs('pct_cells',axis=1,level=1).xs(sv,axis=0,level=1).values
            xv,yv=np.meshgrid(np.arange(ngenes),(nsplit+1)*np.arange(ngroup))
            tr=go.Scatter(x=xv.ravel(),
                          y=yv.ravel()+n,
                          mode='markers',
                          showlegend=False,
                          opacity=1,
                          marker=dict(size=20*np.sqrt(sizes.ravel()),
                                      color=means.ravel(),
                                      colorscale=my_gradients[n]))
            traces.append(tr)
            
        # extra invisible traces for legend
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=20*np.sqrt(s),color='black'),
                            name=p,visible='legendonly',
                            legendgroup='size',showlegend=True) 
                 for s,p in zip([0,.5,1],['0% cells', '50% cells','100% cells'])]
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=20,symbol='square',
                                        color=my_gradients[-1][c][1]),
                            name=p,visible='legendonly',
                            legendgroup='color',showlegend=True) 
                 for c,p in zip([0,1],['low','high'])]
        traces+=[go.Scatter(x=[None],y=[None], mode='markers',
                            marker=dict(size=20,symbol='square',
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
        
    csv_string='data:text/csv;charset=utf-8,'+\
        urllib.parse.quote(data.to_csv(index=True, header=True, encoding='utf-8'))

    fig=dict(data=traces,layout=layout)
    return fig,csv_string


if __name__ == '__main__':
    app.run_server(host='0.0.0.0',debug=True)
                    
