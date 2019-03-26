# single-cell visualization in SODAR using dash + docker

prototype still based on locally provided datasets 

* our own mouse salivary gland tumor data (`export DASH_DATADIR=datasets/hnc`)
* public 10X data from their website
  * 10k mouse heart cells (`export DASH_DATADIR=datasets/heart_10k`)
  * 10k mouse brain cells (`export DASH_DATADIR=datasets/neuron_10k`)
  * 1k human-mouse mix (`export DASH_DATADIR=datasets/hgmm_1k`)
* published PBMCs (stimulated and control) from Kang et al. Nat Biotech 2017 (`export DASH_DATADIR=datasets/pbmc`)

## conda

the dash environment is described in `dash_env.hml` and should install everything that's needed; `requirements.txt` can be used for the `pip` call in docker

## add your own data

the app uses a hdf5 file called `data.h5ad` in the folder specified by `DASH_DATADIR`. this is an `anndata` object (read the docs [here](https://anndata.readthedocs.io/en/latest/index.html)) that stores gene expression (sparse CSR matrix) and meta data with very fast read access. use [`scanpy`](https://scanpy.readthedocs.io/en/latest/index.html) to create such an object with your own data or export from Seurat using `Convert(sobj, to='anndata', filename='data.h5ad')`

also, you should write a small markdown called `about.md` that describes this dataset

## get data from 10X runs

use `python import_from_cellranger -i /path/to/10X/output -o datasets/new_data` to get an `anndata` object for raw 10X data. this does no further processing except log-normalization and uses PCA, tSNE and clustering performed by `cellranger` 

## first version of full app

the app reads the environment variable `DASH_DATADIR` (see above) and then uses data from that directory

```
cd app
export DASH_DATADIR=datasets/hgmm_1k
python app.py
```

## docker

edit the variable `DASH_DATADIR` in the project-specific Dockerfile and then run docker with

```
docker build -t ${project} -f datasets/${project}/Dockerfile .
```

and finally run the docker image 

```
docker run -p 8050:8050 ${project}
```

