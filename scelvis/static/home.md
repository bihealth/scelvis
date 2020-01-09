### SCelVis

SCelVis is a web app for the visualization and interactive exploration of single-cell transcriptomic data.

#### Quickstart

- Start by selecting the dataset you want to work with from the `Go To` menu on the top right: choose an existing dataset, or upload a new one.
- check out our [Documentation](https://scelvis.readthedocs.io/en/latest/index.html), in particular the [Tutorial](https://scelvis.readthedocs.io/en/latest/tutorial_analysis.html), or the movie on our [github page](https://github.com/bihealth/scelvis#Movie)
- some example data to upload: [PBMCs](https://files.figshare.com/18037739/pbmc.h5ad), [human-mouse mix](https://github.com/bihealth/scelvis/raw/master/examples/hgmm_1k.h5ad).

#### Main user interface

![go to https://github.com/bihealth/scelvis for an animated tutorial](https://raw.githubusercontent.com/bihealth/scelvis/master/scelvis/static/scelvis_screenshot.png)

1. switch between tabs:

   - `About`: information about dataset
   - `Cell Annotation`: browse cell meta data
   - `Gene Expression`: explore gene expression

2. select plot type

3. open control panel to filter cells

4. choose x- and y-axis of scatter plot (usually embeddings like tSNE or UMAP)

5. open control panel to perform differential expression analysis between groups of cells selected on the scatter plot

6. hover over points on the plot to get more information

7. hover here to access plot controls (download png, zoom, select cells, etc.)

8. download csv file with underlying data



