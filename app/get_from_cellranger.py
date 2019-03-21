import os
import sys
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import gzip
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
    clustering['Cluster']='Cluster_'+clustering.index.astype(str)+'_(n='+clustering['nclust'].astype(str)+')'
    clustering=clustering.set_index('Barcode')['Cluster']

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
    diffexp.columns=['GeneID','Gene','Cluster','p_adj','log2_fc','mean_counts']

expression_file=os.path.join(args.indir,'outs','filtered_feature_bc_matrix','matrix.mtx.gz')
if not os.path.isfile(expression_file):
    raise Exception('cannot find expression file at '+expression_file)
else:
    print('reading DGE from '+expression_file,file=sys.stderr)
    DGE=scipy.io.mmread(expression_file).tocsr()

barcode_file=os.path.join(args.indir,'outs','filtered_feature_bc_matrix','barcodes.tsv.gz')
if not os.path.isfile(barcode_file):
    raise Exception('cannot find barcode file at '+barcode_file)
else:
    print('reading barcodes from '+barcode_file,file=sys.stderr)
    barcodes=pd.read_csv(barcode_file,header=None,index_col=None).squeeze()

feature_file=os.path.join(args.indir,'outs','filtered_feature_bc_matrix','features.tsv.gz')
if not os.path.isfile(feature_file):
    raise Exception('cannot find feature file at '+feature_file)
else:
    print('reading features from '+feature_file,file=sys.stderr)
    features=pd.read_csv(feature_file,header=None,index_col=None,sep='\t')

print('combining meta data',file=sys.stderr)
meta=(tsne.join(clustering)
      .join(pd.DataFrame(dict(nUMI=DGE.sum(0).A1,
                              nGene=(DGE > 0).sum(0).A1),
                         index=barcodes))
      .join(pca[['PC-1','PC-2','PC-3']]))

# collapse DGE for gene names
print('collapsing DGE by gene name')
genes,rev=np.unique(features[1],return_inverse=True)
mask=scipy.sparse.csr_matrix((np.ones(len(rev),int),
                              (rev,np.arange(len(rev)))),
                             shape=(len(genes),len(rev)))
dge=mask.dot(DGE)

if not args.use_raw:
    print('normalizing DGE')
    dge=dge.dot(scipy.sparse.diags(1.e4/DGE.sum(0).A.ravel()))

if not os.path.isdir(args.outdir):
    raise Exception('output directory '+args.outdir+' does not exist!')

if args.split_species:
    print('determining species mixing',file=sys.stderr)
    species=features[0].str.split('_',n=1).str[0]
    nUMI={}
    for sp in species.unique():
        nUMI['nUMI_'+sp]=DGE[np.where(species==sp)[0]].sum(0).A1
    nUMI=pd.DataFrame.from_dict(nUMI)
    nUMI.index=barcodes
    ratio=nUMI.divide(nUMI.sum(axis=1),axis=0)
    nUMI['species']=(ratio > .95).idxmax(axis=1).str.split('_').str[1]
    nUMI['species'][ratio.max(axis=1) < .95]='unclear'
    meta=meta.join(nUMI)

meta_file=os.path.join(args.outdir,'meta.tsv')
print('saving meta data to '+meta_file,file=sys.stderr)
meta.loc[barcodes].to_csv(meta_file,sep='\t',header=True,index=True)

marker_file=os.path.join(args.outdir,'markers.tsv')
print('saving markers to '+marker_file,file=sys.stderr)
diffexp.drop('GeneID',axis=1).to_csv(marker_file,sep='\t',header=True,index=False)

DGE_file=os.path.join(args.outdir,'expression.mtx.gz')
print('saving DGE to '+DGE_file,file=sys.stderr)
scipy.io.mmwrite(gzip.open(DGE_file,'wb'),dge)

cell_file=os.path.join(args.outdir,'cells.tsv.gz')
print('saving cells to '+cell_file,file=sys.stderr)
barcodes.to_csv(cell_file,header=None,index=0)

gene_file=os.path.join(args.outdir,'genes.tsv.gz')
print('saving genes to '+gene_file,file=sys.stderr)
pd.Series(genes).to_csv(gene_file,header=None,index=0)
