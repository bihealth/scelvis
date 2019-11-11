=======================================
SCelVis: Easy Single-Cell Visualization
=======================================

.. image:: https://img.shields.io/conda/dn/bioconda/scelvis.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/scelvis/README.html

.. image:: https://img.shields.io/pypi/pyversions/scelvis.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/scelvis.svg
    :target: https://pypi.python.org/pypi/scelvis

.. image:: https://api.codacy.com/project/badge/Grade/9ee0ec1424c143dfad9977a649f917f7
    :target: https://www.codacy.com/app/bihealth/scelvis?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/scelvis&amp;utm_campaign=Badge_Grade

.. image:: https://api.codacy.com/project/badge/Coverage/9ee0ec1424c143dfad9977a649f917f7
    :target: https://www.codacy.com/app/bihealth/scelvis?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/scelvis&amp;utm_campaign=Badge_Coverage

.. image:: https://travis-ci.org/bihealth/scelvis.svg?branch=master
    :target: https://travis-ci.org/bihealth/scelvis

.. image:: https://zenodo.org/badge/185944510.svg
    :target: https://zenodo.org/badge/latestdoi/185944510

|

.. image:: tutorial/scelvis_movie.gif
    :height: 400px
    :align: center

------------
Installation
------------

The only prerequisite is Python 3, everything else will be installed together with the ``scelvis`` package.

You can install SCelVis and its dependencies using ``pip`` or through ``conda``:

.. code-block:: shell

    $ pip install scelvis
    # OR
    $ conda install scelvis

A Docker container is also available via `Quay.io/Biocontainers <https://quay.io/organization/biocontainers>`_.

.. code-block:: shell

    $ docker run quay.io/biocontainers/scelvis:TAG scelvis --help
    $ docker run -p 8050:8050 -v data:/data quay.io/biocontainers/scelvis:TAG scelvis run --data-source /data

look up the latest ``TAG`` to use at `here <https://quay.io/repository/biocontainers/scelvis?tab=tags>`_, e.g.,

.. code-block:: shell

    $ docker run quay.io/biocontainers/scelvis:0.7.0--py_0 scelvis --help
    $ docker run -p 8050:8050 -v data:/data quay.io/biocontainers/scelvis:0.7.0--py_0 scelvis run --data-source /data


--------
Tutorial
--------

explore 1000 cells from a 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells (10X v3 chemistry) or a published dataset of ~14000 IFN-beta treated and control PBMCs from 8 donors (`GSE96583 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_; see `Kang et al. <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_)

.. code-block:: shell

    $ scelvis run --data-source /path/to/scelvis/examples/hgmm_1k.h5ad
    $ scelvis run --data-source https://files.figshare.com/18037739/pbmc.h5ad


and then point your browser to http://0.0.0.0:8050/ or http://localhost:8050/


-------------------
Preparing Your Data
-------------------

Data sets are provided as HDF5 files (`anndata <https://anndata.readthedocs.io/en/latest/index.html>`__ objects) that store gene expression (sparse CSR matrix) and meta data with very fast read access.  

For the input you can either specify one HDF5 file or a directory containing multiple such files.

You can use `scanpy <http://scanpy.rtfd.io>`__ to create this HDF5 file directly or use the ``scelvis convert`` command for converting your single-cell pipeline output.

HDF5 Input
----------

for HDF5 input, you can do your analysis with `scanpy <http://scanpy.rtfd.io>`__ to create an anndata object ``ad``. SCelVis will use embedding coordinates from ``ad.obsm``, cell annotation from ``ad.obs`` and expression data directly from ``ad.X`` (this should contain normalized and log-transformed expression values for all genes). If present, information about the dataset will be extracted from strings stored in ``ad.uns['about_title']``, ``ad.uns['about_short_title']`` and ``ad.uns['about_readme']`` (assumed to be Markdown). Information about marker genes will be taken either from the ``rank_genes_groups`` slot in ``ad.uns`` or from entries starting with ``marker_`` in ``ad.uns``: entries called ``marker_gene`` (required!), ``marker_cluster``, ``marker_padj``, ``marker_LFC`` will create a table with the columns ``gene``, ``cluster``, ``padj``, and ``LFC``.

If you prepared your data with ``Seurat`` (v2), you can use ``Convert(from = sobj, to = "anndata", filename = "data.h5ad")`` to get an HDF5 file.

Text Input
----------

For "raw" text input, you need to prepare at least three files in the input directory:

- ``expression.tsv.gz``, a tab-separated file with normalized expression values for each gene (rows) in each cell (columns), e.g., like this::

        .       cell_1   cell_2   cell_3  ...
        gene_1  0.13     0.0      1.5     ...
        gene_2  0.0      3.1      0.3     ...
        gene_3  0.0      0.0      0.0     ...

- ``annotation.tsv``, a tab-separated file with annotations for each cell, e.g., like this::

        .         cluster     genotype  ...
        cell_1    cluster_1   WT        ...
        cell_2    cluster_2   KO        ...


- ``coords.tsv``, a tab-separated file with embedding coordinates for each cell, e.g., like this::

        .         tSNE_1   tSNE_2   UMAP_1  UMAP_2  ...
        cell_1    20.53    -10.05   3.9     2.4     ...
        cell_2    -5.34    13.94    -1.3    3.4     ...

- ``markers.tsv``, an optional tab-separated file with marker genes and **it needs to have a column named ``gene``**, e.g., like this::

        gene    cluster     log2FC   adj_pval   ...
        gene_1  cluster_1   3.4      1.5e-6     ...
        gene_2  cluster_1   1.3      0.00004    ...
        gene_3  cluster_2   2.1      5.3e-9     ...

- a markdown file (e.g., ``text_input.md``) with information about this dataset::

        ----
        title: An Optional Long Data Set Title
        short_title: optional short title
        ----

        A verbose description of the data in Markdown format.

.. code-block:: shell

    $ scelvis convert --input-dir text_input --output data/text_input.h5ad --about-md text_input.md

in ``examples/dummy_raw.zip`` and ``examples/dummy_about.md`` we provide raw data for a simulated dummy dataset.

Loom Input
----------

for `loompy <http://loompy.org>`__ or `loomR <https://github.com/mojaveazure/loomR>`__ input, you can convert your data like this:

.. code-block:: shell

    $ scelvis convert --i input.loom -m markers.tsv -a about.md -o loom_input.h5ad 

if you prepared your data with ``Seurat`` (v3), you can use ``as.loom(sobj, filename="output.loom")`` to get a ``.loom`` file and then convert to ``.h5ad`` with the above command.

CellRanger Input
----------------

Alternatively, the output directory of ``CellRanger`` can be used. This is the directory called ``outs`` containing either a file called ``filtered_gene_bc_matrices_h5.h5`` (version 2) or a file called ``filtered_feature_bc_matrix.h5`` (version 3), and a folder ``analysis`` with clustering, embedding and differential expression results. This will not no any further processing except log-normalization. Additionally, a markdown file provides meta information about the dataset (see above)

.. code-block:: shell

    $ mkdir -p data
    $ cat <<EOF > data/cellranger.md
    ----
    title: My Project
    short_title: my_project
    ----

    This is my project data.
    EOF
    $ scelvis convert --input-dir cellranger-out --output data/cellranger_input.h5ad --about-md cellranger.md

In ``examples/hgmm_1k_raw`` we provide ``CellRanger`` output for the 1k 1:1 human mouse mix. Specifically, from the ``outs`` folder we selected

- ``filtered_feature_bc_matrix.h5``
- tSNE and PCA projections from ``analysis/tsne`` and ``analysis/pca``
- clustering from ``analysis/clustering/graphclust`` and
- markers from ``analysis/diffexp/graphclust``

``examples/hgmm_1k_about.md`` contains information about this dataset

---------------------
Visualizing Your Data
---------------------

.. code-block:: shell

    $ tree data
    data
    ├── text_input.h5ad
    └── cellranger_input.h5ad

    $ scelvis run --data-source data/cellranger_input.h5ad
    # OR
    $ scelvis run --data-source data

------------
Data Sources
------------

Data sources can be:

- paths, e.g., ``relative/paths`` or ``/absolute/paths`` or ``file://url/paths``
- SFTP URLs, e.g., ``sftp://user:password@host/path/to/data``
- FTP URLs, e.g., ``ftp://user:password@host/path/to/data`` (sadly encryption is not supported by the underlying library `PyFilesystem2 <https://github.com/PyFilesystem/pyfilesystem2>`__.
- iRODS URLS, e.g., ``irods://user:password@host/zoneName/path/to/data``
    - Enable SSL via ``irods+ssl``
    - Switch to PAM authentication with ``irods+pam`` (you can combine this with ``+ssl`` in any order)
    - Enable ticket access by appending ``?ticket=TICKET``.
- HTTP(S) URLs, e.g., ``https://user:password@host/path/to/data``.
- S3 URLs, e.g., ``s3://bucket/path``, optionally ``s3://key:token@bucket/path``.

Data sources can either point to HDF5 files directly or to directories containing multiple HDF5 files.
The only exception is iRODS with ticket-based access.
Because of technical restrictions, you have to assign a unique ticket for each data set and specify the data sets individually.

---------------------
Environment Variables
---------------------

You can use the following environment variables to configure the server.

- ``SCELVIS_DATA_SOURCES`` -- semicolon-separated list of data sources
- ``SCELVIS_HOST`` -- host specification for web server to listen on
- ``SCELVIS_PORT`` -- port for web server to listen on
- ``SCELVIS_CACHE_DIR`` -- directory to use for the cache (default is to create a temporary directory)
- ``SCELVIS_CACHE_REDIS_URL`` -- enable caching with REDIS and provide connection URL
- ``SCELVIS_CACHE_DEFAULT_TIMEOUT`` -- cache lifetime coverage
- ``SCELVIS_CACHE_PRELOAD_DATA`` -- will preload all data at startup
- ``SCELVIS_UPLOAD_DIR`` -- the directory to store uploaded data sets in (default is to create a temporary directory)
- ``SCELVIS_UPLOAD_DISABLED`` -- set to "0" to disable upload feature
- ``SCELVIS_CONVERSION_DISABLED`` -- set to "0" to disable the conversion feature
- ``SCELVIS_URL_PREFIX`` -- set if you want to run scelvis below a non-root path (e.g., behind a reverse proxy)

---------------
Developer Setup
---------------

The prerequisites are:

- Python 3, either
    - system-wide installation with ``virtualenv``, or
    - installed with `Conda <https://docs.conda.io/en/latest/>`__.
- `Git LFS <https://git-lfs.github.com/>`__ must be installed for obtaining the example data files.

For ``virtualenv``, first create a virtual environment and activate it.

.. code-block:: shell

    $ virtualenv -p venv
    $ source venv/bin/activate

For a Conda-based setup create a new environment and activate it.

.. code-block:: shell

    $ conda create -y -n scelvis 'python>=3.6'
    $ conda activate scelvis

Next, clone the repository and install the software as editable (``-e``).
Also install the development requirements to get helpers such as black.
(Note that you must have `Git LFS <https://git-lfs.github.com/>`__ installed to actually obtain the data files).

.. code-block:: shell

    $ git clone git@github.com:bihealth/scelvis.git
    $ cd scelvis
    $ pip install -e .
    $ pip install -r requirements/develop.txt

Afterwards, you can run the visualization web server as follows:

.. code-block:: shell

    $ scelvis run --data-source path/to/data/dir

To explore the datasets provided in the git repository, use ``git lfs fetch`` to download

Releasing Packages
------------------

For the `PyPi package <https://pypi.org/project/scelvis/>`__:

.. code-block:: shell

    $ python setup.py sdist
    $ twine upload --repository-url https://test.pypi.org/legacy/ dist/scelvis-*.tar.gz
    $ twine upload dist/scelvis-*.tar.gz

For the Bioconda package, see `the great documentation <http://bioconda.github.io/updating.html>`__.
The Docker image will automatically be created as a BioContainer when the Bioconda package is built.
