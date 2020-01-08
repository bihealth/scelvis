.. _tutorial_analysis:

=================
Analysis Tutorial
=================

This tutorial describes the basics of performing the analysis of data using SCelVis.
For this tutorial, we will use the public demo instance at https://scelvis-demo.bihealth.org.

.. note::

    Data to be visualized can either be uploaded into a SCelVis server or it can be defined when the SCelVis server starts.
    When using a remote SCelVis server such as the public instance at `scelvis-demo.bihealth.org <https://scelvis-demo.bihealth.org>`_, you will most likely upload your data as shown below.
    However, the server can also be started with the path or URL to the data.
    This way, computational staff can provide predefined datasets to non-computational staff.
    See :ref:`tutorial_cli` for more information.

First, download the file `pbmc.h5ad <https://files.figshare.com/18037739/pbmc.h5ad>`_, which is a published dataset of ~14000 IFN-beta treated and control PBMCs from 8 donors (`GSE96583 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_; see `Kang et al. <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_). For a simpler dataset, you could also use  `hgmm_1k.h5ad <https://github.com/bihealth/scelvis/raw/master/examples/hgmm_1k.h5ad>`_, containing data for 1000 cells from a 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells (10X v3 chemistry).
    
-----------
Upload Data
-----------

Then, use the menu on the top right to access the data upload screen :guilabel:`Go To --> Upload Data`.

.. figure:: figures/scelvis_goto_upload.png
    :width: 80%
    :align: center

    Accessing the upload screen.

On the next screen click :guilabel:`Choose File`, select the *pbmc.h5ad* file, and click :guilabel:`Upload`.
The file upload will take a while and return a link to the data analysis screen shown in the next section.

.. note::

    Alternatively, you can use the example data set available from the top-right menu: :guilabel:`Go To --> Data Sets --> PBMC`.

------------------------
Cell Annotation Analysis
------------------------

In the beginning of each analysis, the cell annotation screen from the figure below is shown.

.. figure:: figures/scelvis_cell_scatter.png
    :width: 80%
    :align: center

    The cell annotation scatter plot.

1. **Analysis selection tab.**
   This control allows to switch between

    About
        A textual description attached to the dataset.
    Cell Annotation
        The cell annotation analysis screen that you are looking at.
    Gene Expression
        The gene expression analysis screen.

2. **Plot type selection.**
   This allows to select the plot type for the cell-based analysis.
   Note that changing the plot type will change the subsequent controls.

3. **Cell filter button.**
   Using this control, you can filter cells based on various properties, see Section :ref:`tutorial_cell_filtering`.

4. **Axis and color control.**
   This allows you to select the dimensions to display along the horizontal and vertical axes as well as the
   colouring. 

5. **Differential expression button.**
   This allows you to run a differential expression between two groups of cells, see Section :ref:`tutorial_cell_differential_expression`.

6. **The cell scatter plot.**
   This is a dynamic scatter plot. Each dot corresponds to one cell, colored by what is selected from the :guilabel:`select coloring` list.

7. **Plot controls.**
   When you move your mouse cursor over the plot then some buttons will appear on the top right.
   These are explained in detail in Section :ref:`tutorial_plot_ui_controls` together with some other tricks.

8. **Download data for this plot.**
   Download a CSV file with the data for reproducing the plot outside of SCelVis.

Scatter Plot
============

Usually, you would choose embedding coordinates (e.g., tSNE_1 and tSNE_2) for :guilabel:`select x axis` and :guilabel:`select y axis` to create a standard tSNE or UMAP plot. :guilabel:`select coloring` allows to color cells by different cell annotation attributes, e.g., *cluster* for the cluster annotation or *n_counts* for the # of UMIs per cell. The available choices depend on how the dataset was preprocessed. Categorical variables will be shown with a discrete color scale, numerical variables with a gradient. 

Alternatively, you could also plot, e.g., # of UMIs vs. # of genes for QC.

Violin and Box Plot
===================

When selecting :guilabel:`violin plot` or :guilabel:`box plot` in the :guilabel:`select plot type` control you can draw violin or box plots. For example, selecting *n_counts* and *n_genes* in the :guilabel:`select variable(s)` and :*orig.ident* in :guilabel:`select coloring` will display the *n_genes* value distribution in the upper panel, and the *n_counts* value distribution in the lower panel for the individual samples of this dataset.

You can use the :guilabel:`select split` list to select whether you want to further split the grouping, e.g., by cluster identity. Hovering the mouse over the violin or box shapes will show you various statistical summaries of the distribution.

Bar Plots
=========

With :guilabel:`bar plot`, you can display summary statistics, such as the number of cells per cluster by selecting *cluster* in :guilabel:`select grouping`. With :guilabel:`select split`, you can further investigate how clusters are populated in the different samples or the different donors. Checking :guilabel:`normalized` will switch from cell numbers to fractions, :guilabel:`stacked` will use stacked bars.


------------------------
Gene Expression Analysis
------------------------

When clicking the :guilabel:`Gene Expression` tab, you can investigate gene expression.
As for the :guilabel:`Cell Annotation` analysis, it starts with the :guilabel:`scatter plot` type in **(1)**. The main difference is that all plot types will use the same list of genes selected in :guilabel:`select gene(s)` in **(4)**.

.. figure:: figures/scelvis_gene_scatter.png
    :width: 80%
    :align: center

    The gene expression scatter plot for a monocyte (CCL2) and a T cell marker (SELL).

1. **Plot type selection.**
   This allows to select the plot type for the gene expression analysis.
   Note that changing the plot type will change the subsequent controls.

2. **Cell filter button.**
   Using this control, you can filter cells based on various properties, see Section :ref:`tutorial_cell_filtering`.

3. **Axis control.**
   This allows you to select the dimensions to display along the horizontal and vertical axes.

4. **Selecting genes.**
   Select one or multiple genes from this list or enter them by hand.

5. **Show tables.**
   Check these boxes to display tables with log2-fold change, p-values and other information for marker genes or differential expression results (if available). 

6. **Table selection.**
   Genes can be selected from the table by checking the boxes to the left and clicking :guilabel:`use selected genes` to add them to the list.


Scatter Plot
============

Scatter plots for one or multiple genes will be shown in a grid, with expression values rescaled to the same range.

Violin or Box Plot
==================

Violin and Box plots will show one gene per row, with one violin or box per category selected in :guilabel:`select grouping`, or multiple violins/boxes if :guilabel:`select split` is used.

Dot Plot
========

Dot plots summarise gene expression by size (fraction of expressing cells) and color (expression value), with genes in columns and groups (use :guilabel:`select grouping`) in rows. Dots can be subdivided using :guilabel:`select split`.

.. _tutorial_plot_ui_controls:

Plot Interface Controls
=======================

.. figure:: figures/scelvis_movie_plot_controls.gif
    :width: 80%
    :align: center

    A short demonstration of the plot controls in SCelVis.

The buttons on the top right of the plot (standard features of the `Plotly <https://plot.ly/>`_ library) are as follows:

Save plot as image
    The plot will be saved in PNG format.

Zoom
    After clicking this button, you can select a rectangular area to zoom into.

Pan
    Drag and drop the drawing area to move around in the plot.

Box Select, Lasso Select
    After clicking this button, you can select a rectangular area on the plot or draw a free from shape to select dots. This will be useful for differential expression analysis.

Zoom In / Zoom Out
    Zoom into or out of plot.

Autoscale / Reset Axes
    This will reset the scaling to the original area.
    You can obtain the same behaviour by double-clicking on a white spot in the plot.

Toggle Spike Lines
    Enable horizontal and vertical lines when hovering over a point.

Show Closest / Compare Data on Hover
    Change the spike lines behaviour.

Note that you can enable/disable individual groups by clicking their label in the legend.
You can disable all but one group by double-clicking the label.

.. _tutorial_cell_filtering:

Cell Filtering
==============

.. figure:: figures/scelvis_movie_filtering.gif
    :width: 80%
    :align: center

    A short demonstration of cell filtering in SCelVis.

The cell filtering works as follows:

1. Click the :guilabel:`filter cells` button to open the control panel.

2. Select a criterion by which cells should be filtered. Depending on how the data were preprocessed, the list will include cluster annotation, # of UMIs or genes per cell, etc. It is also possible to filter cells by expression of genes.

3. For categorical variables (e.g., cluster identity), checkboxes will appear and specific clusters can be (un)checked in order to include or exlude them from the analysis. For numerical variables (e.g., n_counts or gene expression), a range slider will appear: move the big circles inwards to remove cells outside the selected range.

4. Hit :guilabel:`update plot` to apply these filters to the current plot   

5. Filters will be combined with AND logic; active filters are listed above the :guilabel:`update plot` button

6. Click :guilabel:`reset filters` to reset all filters and :guilabel:`update plot` to include all cells in the current plot

7. Note that current filter criteria will be applied to all subsequent plots of the current datasets, both in the :guilabel:`Cell Annotation` and the :guilabel:`Gene Expression` tabs
   
.. _tutorial_cell_differential_expression:

Differential Expression Analysis
================================

.. figure:: figures/scelvis_movie_differential_expression.gif
    :width: 80%
    :align: center

    A short demonstration of differential expression analysis in SCelVis.

The differential expression analysis is available only when a scatter plot is displayed in the :guilabel:`Cell Annotation` tab. It works as follows:

1. Click the :guilabel:`differential expression` button, opening the controls panel.

2. Then, use either the box or lasso select tool of the plot for selecting cells in the scatter plot.
   For example, click the lasso select button in the top of the right of the plot.
   Move your mouse cursor to the position you want to start selecting at.
   Keep the left mouse button pressed and draw a shape around the cells that you are interested in.
   Release the mouse button. then Click :guilabel:`group A`.

3. Repeat step 2 but click :guilabel:`group B`.

4. Click :guilabel:`run` to perform the analysis.

The result could read something like *200 DE genes at 5% FDR* (a maximum of 200 genes will be displayed).
You can click :guilabel:`view groups` to show the groups in the scatter plot, or click :guilabel:`table` to see the resulting DE genes in the :guilabel:`Gene Expression` tab table.
You can also download the :guilabel:`results` or the :guilabel:`parameters` that were used for the DE gene analysis.

Clicking :guilabel:`reset` allows you to start a new DE gene analysis.

-------
The End
-------

This is the end of the data analysis tutorial.
You might want to learn about the conversion of data into the HDF5 format next by reading Section :ref:`tutorial_convert`.
