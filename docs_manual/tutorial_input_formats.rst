.. _tutorial_input_formats:

---------------------
Source Format Details
---------------------

HDF5 Input
==========

For HDF5 input (no conversion necessary), you can do your analysis with `scanpy <http://scanpy.rtfd.io>`__ to create an anndata object ``ad``. SCelVis will use embedding coordinates from ``ad.obsm``, cell annotation from ``ad.obs`` and expression data directly from ``ad.raw.X`` (if present) or ``ad.X`` (this should contain normalized and log-transformed expression values for all genes and should be sparse, otherwise performance will suffer)
If present, information about the dataset will be extracted from strings stored in ``ad.uns['about_title']``, ``ad.uns['about_short_title']`` and ``ad.uns['about_readme']`` (assumed to be Markdown).
Information about marker genes will be taken from entries starting with ``marker_`` in ``ad.uns``: entries called ``marker_gene`` (required!), ``marker_cluster``, ``marker_padj``, ``marker_LFC`` will create a table with the columns ``gene``, ``cluster``, ``padj``, and ``LFC``. ``SCelVis`` will also extract marker information from ``ad.uns['rank_genes_groups']``. However, certain datatypes in ``ad.uns`` together with version mismatches of ``scanpy``, ``h5py`` and ``anndata`` can lead to ``.h5ad`` files that are not readable by ``SCelVis`` (see `issue #832 <https://github.com/theislab/scanpy/issues/832>`__). To be on the safe side, it's recommended to delete unneccessary slots in ``ad.uns`` (e.g., ``del ad.uns['rank_genes_groups']``). Also, ``ad.obsm['X_pca']``, ``ad.varm['PCs']`` and entries in ``ad.obsp`` are likely dispensable.

If you prepared your data with ``Seurat`` (v2), you can use ``Convert(from = sobj, to = "anndata", filename = "data.h5ad")`` to get an HDF5 file.

Alternatively, you can use `sceasy <https://github.com/cellgeni/sceasy>`__ to convert your objects into ``anndata`` HDF5 format.

Text Input
==========

For "raw" text input, you need to prepare a file with expression values, cell meta data and coordinates, and potentially information about this dataset as well as cluster markers.

- normalized expression values for each gene (rows) in each cell (columns) can be given either as tab-separated file (dense) or in matrix-market format:

  - if your file is called ``expression.tsv.gz``, ``SCelVis`` expects a tab-separated file , e.g., like this::
  
          .       cell_1   cell_2   cell_3  ...
          gene_1  0.13     0.0      1.5     ...
          gene_2  0.0      3.1      0.3     ...
          gene_3  0.0      0.0      0.0     ...

  - if your file is called ``expression.mtx``, ``SCelVis`` expects this to be a sparse matrix-market file and additional files called ``barcodes.tsv`` (containing a list of cell names / barcodes, one per line, no header or row names) and ``genes.tsv`` (containing a list of gene names, one per line, no header or row names) to be present.

- annotations for each cell can be provided  in a tab-separated file called ``annotation.tsv``, e.g., like this::

        .         cluster     genotype  ...
        cell_1    cluster_1   WT        ...
        cell_2    cluster_2   KO        ...


- embedding coordinates for each cell can be provided in a tab-separated file called ``coords.tsv``, e.g., like this::

        .         tSNE_1   tSNE_2   UMAP_1  UMAP_2  ...
        cell_1    20.53    -10.05   3.9     2.4     ...
        cell_2    -5.34    13.94    -1.3    3.4     ...

- an optional tab-separated file called ``markers.tsv`` can contain information about marker genes. **it needs to have a column named ``gene``**, e.g., like this::

        gene    cluster     log2FC   adj_pval   ...
        gene_1  cluster_1   3.4      1.5e-6     ...
        gene_2  cluster_1   1.3      0.00004    ...
        gene_3  cluster_2   2.1      5.3e-9     ...

- finally, a markdown file (e.g., ``text_input.md``) can provide information about this dataset::

        ----
        title: An Optional Long Data Set Title
        short_title: optional short title
        ----

        A verbose description of the data in Markdown format.

conversion to ``.h5ad`` is then performed like so:

.. code-block:: shell

    $ scelvis convert --input-dir text_input --output data/text_input.h5ad --about-md text_input.md

in ``examples/dummy_raw.zip`` and ``examples/dummy_about.md`` we provide raw data for a simulated dummy dataset.

if you prepared you data with ``Seurat``, you can export to raw text like this

.. code-block:: r

    writeMM(sobj@assays$RNA@data, file = 'expression.mtx')
    write.table(Cells(sobj), file = 'barcodes.tsv', col.names = FALSE, row.names = FALSE, sep = ',')
    write.table(row.names(sobj@assays$RNA@data), file = 'genes.tsv', col.names = FALSE, row.names = FALSE, sep = ',')
    sobj@meta.data$cluster <- paste0('cluster_', sobj@meta.data$seurat_clusters)
    write.table(sobj@meta.data, file = 'annotation.tsv', sep = '\t')
    write.table(sobj@reductions$umap@cell.embeddings, file = 'coords.tsv', sep = '\t')


Loom Input
==========

for `loompy <http://loompy.org>`__ or `loomR <https://github.com/mojaveazure/loomR>`__ input, you can convert your data like this:

.. code-block:: shell

    $ scelvis convert --i input.loom -m markers.tsv -a about.md -o loom_input.h5ad

if you prepared your data with ``Seurat`` (v3), you can use ``as.loom(sobj, filename = "output.loom")`` to get a ``.loom`` file and then convert to ``.h5ad`` with the above command (this is quite slow, however, and exact format specifications for ``.loom`` and ``.h5ad`` are not always compatible between versions)

CellRanger Input
================

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

In ``examples/hgmm_1k_raw`` we provide ``CellRanger`` output for the 1k 1:1 human mouse mix.
Specifically, from the ``outs`` folder we selected

- ``filtered_feature_bc_matrix.h5``
- tSNE and PCA projections from ``analysis/tsne`` and ``analysis/pca``
- clustering from ``analysis/clustering/graphclust`` and
- markers from ``analysis/diffexp/graphclust``

``examples/hgmm_1k_about.md`` contains information about this dataset.
