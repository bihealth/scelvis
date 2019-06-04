## ScelVis Quickstart

SCelVis is a web app for the visualization and interactive exploration of single-cell transcriptomic data.

**Start by selecting the dataset you want to work with from the `Go To` menu on the top right**:
choose an existing dataset, or upload a new one.

Conversion from raw text or `cellranger` output is also possible:
upload a zip or tarball of the `outs` directory, containing at least the `filtered_feature_bc_matrix.h5` and the `analysis` subfolder.

Datasets can be explored using either a cell-centric (i.e., cell annotation such as clustering and some QC) or a gene-centric view (i.e., gene expression, markers, etc.).

**In the cell-centric view, you can create**

* scatter plots of annotations
* violin plots of cell-based statistics (e.g., number of genes per cell)
* bar plots summarizing cell numbers or proportions

![cell figure](assets/cells.png)

**In the gene-centric view, you can create**

* scatter plots showing expression of (multiple) genes
* violin plots showing gene expression per cluster, possibly sub-grouped ("split") by another factor
* dot plots summarizing expression of many genes per cluster, possibly sub-grouped

Genes can be selected from / entered into the dropdown menu or taken from a list of marker genes if available

![gene figure](assets/genes.png)
