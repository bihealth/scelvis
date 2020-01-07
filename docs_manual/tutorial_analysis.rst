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

First, obtain the file `hgmm_1k.h5ad <https://github.com/bihealth/scelvis/raw/master/examples/hgmm_1k.h5ad>`_ and download it to your computer.
This file contains data for 1000 cells from a 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells (10X v3 chemistry).
For a larger dataset, you could also use `pbmc.h5ad <https://files.figshare.com/18037739/pbmc.h5ad>`_ which is a published dataset of ~14000 IFN-beta treated and control PBMCs from 8 donors (`GSE96583 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_; see `Kang et al. <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583>`_).

-----------
Upload Data
-----------

Then, use the menu on the top right to access the data upload screen :guilabel:`Go To --> Upload Data`.

.. figure:: figures/scelvis_goto_upload.png
    :width: 80%
    :align: center

    Accessing the upload screen.

On the next screen click the :guilabel:`Drag and Drop or Click and Select files.` and select the *hgmm_1k.h5ad* file for upload.
The file upload will take a while and you will be redirected to the data analysis screen shown in the next section.

.. note::

    Alternatively, you can use the example data set available through the top-right menu as :guilabel:`Go To --> Data Sets --> PBMC`.

------------------------
Cell Annotation Analysis
------------------------

In the beginning of each analysis, the cell annotation screen from the figure below is shown.

Scatter Plot
============

.. note:: @Benedikt: here we would need some more pointers towards the biological interpretation?

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
        The gene expression analsysis screen.

2. **Plot type selection.**
   This allows tos elect the plot type for the cell-based analysis.
   Note that changing the plot type will change the subsequent controls.

3. **Cell filter button.**
   Using this control, you can filter cells based on various properties.
   These includes the cluster, species and other properties.

   NB: The values available here depend on the configuration of the analysis software generating the HDF5 file.

4. **Axis and color control.**
   This allows you to select the dimensions to display along the horizontal and vertical axes as well as the colouring.

5. **Differential expression button.**
   This allows you to run a differential expression between two groups of cells.
   This will be explained in detail below in Section :ref:`tutorial_cell_differential_expression`.

6. **The cell scatter plot.**
   This is a dynamic scatter plot.
   When doing cell annotation analysis, each dot corresponds to one cell.
   The cells are colored by cluster by default but you can change this with the :guilabel:`select coloring` list.

7. **Plot controls.**
   When you move your mouse cursor over the plot then some buttons will appear on the top right.
   These are explained in detail in Section :ref:`tutorial_plot_ui_commands` together with some other tricks.

8. **Download data for this plot.**
   Download a CSV file with the data for reproducing the plot outside of SCelVis.

.. _tutorial_plot_ui_commands:

Plot Interface Commands
=======================

The buttons on the top right of the plot are as follows.
NB: these are standard features of the `Plotly <https://plot.ly/>`_ library.

Save plot as image
    The plot will be saved in PNG format.

Zoom
    After clicking this button, you can select a rectangular area to zoom into.

Pan
    Drag and drop the drawing area to move around in the plot.

Box Select, Lasso Select
    After clicking this button, you can select a rectangular area on the plot or draw a free from shape to select dots.

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

.. _tutorial_cell_differential_expression:

Differential Expression Analysis
================================

The differential expression analysis works as follows:

1. Click the :guilabel:`differential expression` button.
   Four buttons appear: :guilabel:`group A`, :guilabel:`group B`, :guilabel:`reset`, :guilabel:`run`.

2. Click :guilabel:`group A`.
   Then, use either the box or lasso select tool of the plot for selecting some points.
   For example, click the lasso select button in the top of the right of the plot.
   Move your mouse cursor to the position you want to start selecting at.
   Keep the left mouse button pressed and draw a shape around the cells that you are interested in.
   Release the mouse button.

3. Click :guilabel:`group B` and repeat step 2.

4. Click :guilabel:`run` to perform the analysis.

The result could read something like *200 DE genes at 5% FDR*.
You can then click :guilabel:`view groups` to show the groups in the scatter plot, click :guilabel:`table` to see the resulting DE genes in the :guilabel:`Gene Expression` tab table.
You can also download the :guilabel:`results` or the :guilabel:`parameters` that were used for the DE gene analysis.

Clicking :guilabel:`reset` allows you to start a new DE gene analysis.

.. note:: @Benedikt: here we would need a specific example with some possible interpretation.

Violin Plot
===========

When selecting :guilabel:`violin plot` in the :guilabel:`select plot type` control you can draw violin plots.
For example, selecting :guilabel:`n_counts` and :guilabel:`nUMI_hg19` in the :guilabel:`select variable(s) and scaling` shows the following plot.
The top half shows the *nUMI_hg19* value distribution and the lower half shows the *n_counts* value distribution for the individual clusters.

.. figure:: figures/scelvis_cell_violin.png
    :width: 80%
    :align: center

    The cell annotation scatter plot.

Use the :guilabel:`grouping` list to select the dimension to display along the horizontal axis.
You can use the :guilabel:`select split` list to select whether you want to further split the grouping.

Note that you can use a subset of the plot controls as described in Section :ref:`tutorial_plot_ui_commands` and also download the data for the violin plot.
Hovering the mouse over the violin shapes will show you various statistical summaries of the distribution.

Bar and Box Plots
=================

More visualization options are available by selecting the :guilabel:`bar plot` and :guilabel:`box plot` plot types.

------------------------
Gene Expression Analysis
------------------------

.. note:: @Benedikt: what about the show diffexp results?

When clicking the :guilabel:`Gene Expression` tab, you can perform a gene expression analysis.
As for the :guilabel:`Cell Annotation` analysis, it starts with the :guilabel:`scatter plot` type in (1).
However, first select some genes, for example *mm10_Hes1* and *hg19_HES1* in (2).

Scatter Plot
============

.. note:: @Benedikt: here we would need some more pointers towards the biological interpretation?

.. figure:: figures/scelvis_gene_scatter.png
    :width: 80%
    :align: center

    The gene expression scatter plot.

You can also display a table with the genes including their log2-fold change, counts and other information usin gthe :guilabel:`show marker table` button.
Using the check boxes on the left-hand-side (below; 1), you can select one or more genes and use the :guilabel:`use selected genes` (below; 2) button to add them to the :guilabel:`select gene(s)` field:

.. figure:: figures/scelvis_gene_table.png
    :width: 80%
    :align: center

    The gene expression scatter plot.

.. note:: @Benedikt: can you do some selection here and offer some biological interpretation?

-------
The End
-------

This is the end of the data analysis tutorial.
You might want to learn about the conversion of data into the HDF5 format next by reading Section :ref:`tutorial_convert`.
