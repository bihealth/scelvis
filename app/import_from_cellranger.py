import os
import sys
import anndata
import scanpy.api as sc
import numpy as np
import pandas as pd
from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('-i','--indir',dest='indir',
                    help="""input directory (=cellranger output directory)""")
parser.add_argument('-o','--outdir',dest='outdir',
                    help="""output directory""")
parser.add_argument('--use_raw',dest='use_raw',
                    default=False,action='store_true',
                    help="""do not normalize DGE (use raw counts)""")
parser.add_argument('--split_species',dest='split_species',
                    default=False,action='store_true',
                    help="""split species""")
parser.add_argument('--nmarkers',dest='nmarkers',
                    default=10,type=int,
                    help="""save top n markers per cluster in markers.file [10]""")

args=parser.parse_args()

tsne_file=os.path.join(args.indir,'outs','analysis','tsne','2_components','projection.csv')
if not os.path.isfile(tsne_file):
    raise Exception('cannot find tSNE output at '+tsne_file)
else:
    print('reading tSNE output from '+tsne_file,file=sys.stderr)
    tsne=pd.read_csv(tsne_file,header=0,index_col=0)

pca_file=os.path.join(args.indir,'outs','analysis','pca','10_components','projection.csv')
if not os.path.isfile(pca_file):
    raise Exception('cannot find PCA output at '+pca_file)
else:
    print('reading PCA output from '+pca_file,file=sys.stderr)
    pca=pd.read_csv(pca_file,header=0,index_col=0)

clustering_file=os.path.join(args.indir,'outs','analysis','clustering','graphclust','clusters.csv')
if not os.path.isfile(clustering_file):
    raise Exception('cannot find clustering output at '+clustering_file)
else:
    print('reading clustering output from '+clustering_file,file=sys.stderr)
    clustering=pd.read_csv(clustering_file,header=0)
    nc=pd.DataFrame(dict(nclust=clustering['Cluster'].value_counts()))
    clustering=clustering.set_index('Cluster')
    clustering['nclust']=nc
    clustering['Cluster']='Cluster_'+clustering.index.astype(str)
    clustering=clustering.set_index('Barcode')['Cluster']
    clustering=clustering.astype('category')

diffexp_file=os.path.join(args.indir,'outs','analysis','diffexp','graphclust','differential_expression.csv')
if not os.path.isfile(diffexp_file):
    raise Exception('cannot find differential expression output at '+diffexp_file)
else:
    print('reading differential expression output from '+diffexp_file,file=sys.stderr)
    diffexp=pd.read_csv(diffexp_file,header=0,index_col=[0,1])
    diffexp.columns=(diffexp.columns.str.replace('Cluster ','Cluster_')
                     .str.split(n=1,expand=True)
                     .rename(['Cluster','Obs']))
    diffexp=diffexp.stack(level=0).reset_index()
    diffexp.columns=['GeneID','gene','Cluster','p_adj','log2_fc','mean_counts']

expression_file=os.path.join(args.indir,'outs','filtered_feature_bc_matrix.h5')
if not os.path.isfile(expression_file):
    raise Exception('cannot find expression file at '+expression_file)
else:
    print('reading gene expression from '+expression_file,file=sys.stderr)
    ad=sc.read_10x_h5(expression_file)
    ad.var_names_make_unique()

print('combining meta data',file=sys.stderr)
ad.obs['cluster']=clustering
ad.obs['n_counts']=ad.X.sum(1).A1
ad.obs['n_genes']=(ad.X > 0).sum(1).A1

print('adding coordinates', file=sys.stderr)
ad.obsm['X_tsne']=tsne.values
ad.obsm['X_pca']=pca.values

if False:
    # collapse DGE for gene names
    print('collapsing DGE by gene name')
    genes,rev=np.unique(features[1],return_inverse=True)
    mask=scipy.sparse.csr_matrix((np.ones(len(rev),int),
                                  (rev,np.arange(len(rev)))),
                                 shape=(len(genes),len(rev)))
    dge=mask.dot(DGE)

if args.split_species:
    print('determining species mixing',file=sys.stderr)
    species=ad.var_names.str.split('_',n=1).str[0]
    for sp in species.unique():
        ad.obs['nUMI_'+sp]=ad.X[:,species==sp].sum(1).A1
    cols='nUMI_'+species.unique()
    ratio=ad.obs[cols].divide(ad.obs[cols].sum(axis=1),axis=0)
    ad.obs['species']=(ratio > .95).idxmax(axis=1).str.split('_').str[1]
    ad.obs['species'][ratio.max(axis=1) < .95]='other'

if not args.use_raw:
    print('normalizing and filtering DGE')
    sc.pp.filter_cells(ad, min_genes=100)
    sc.pp.filter_genes(ad, min_cells=5)
    sc.pp.normalize_per_cell(ad,counts_per_cell_after=1.e4)
    sc.pp.log1p(ad)

if not os.path.isdir(args.outdir):
    raise Exception('output directory '+args.outdir+' does not exist!')

out_file=os.path.join(args.outdir,'data.h5ad')
print('saving anndata object to'+out_file,file=sys.stderr)
ad.write(out_file)

marker_file=os.path.join(args.outdir,'markers.csv')
print('saving top '+str(args.nmarkers)+' markers per cluster to '+marker_file,file=sys.stderr)
diffexp[(diffexp['p_adj'] < .05) & (diffexp['log2_fc'] > 0)]\
    .drop('GeneID',axis=1)\
    .sort_values(['Cluster','p_adj'])\
    .groupby('Cluster')\
    .head(args.nmarkers)\
    .to_csv(marker_file,header=True,index=True)

if False:
    DGE_file=os.path.join(args.outdir,'expression.mtx.gz')
    print('saving DGE to '+DGE_file,file=sys.stderr)
    scipy.io.mmwrite(gzip.open(DGE_file,'wb'),dge)

    cell_file=os.path.join(args.outdir,'cells.tsv.gz')
    print('saving cells to '+cell_file,file=sys.stderr)
    barcodes.to_csv(cell_file,header=None,index=0)

    gene_file=os.path.join(args.outdir,'genes.tsv.gz')
    print('saving genes to '+gene_file,file=sys.stderr)
    pd.Series(genes).to_csv(gene_file,header=None,index=0)
