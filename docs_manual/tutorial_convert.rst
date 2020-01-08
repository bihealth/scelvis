.. _tutorial_convert:

===================
Conversion Tutorial
===================

For visualization, data sets are provided as HDF5 files (`anndata <https://anndata.readthedocs.io/en/latest/index.html>`__ objects) that store gene expression (sparse CSR matrices) and meta data with very fast read access.
You can use `scanpy <http://scanpy.rtfd.io>`__ to create these HDF5 files directly or use the ``scelvis convert`` command for converting your single-cell pipeline output (see Section :ref:`tutorial_cli`).

------------------
Web File Converter
------------------

You can use the :guilabel:`Go To --> Convert Data` menu entry to access the file conversion screen.

.. figure:: figures/scelvis_goto_upload.png
    :width: 80%
    :align: center

    Accessing the conversion screen.

Here, you can enter a title, short title and a description of your dataset, and upload a .zip or .tar.gz file containing the data with :guilabel:`Choose File`. Allowed formats are (see also :ref:`tutorial_input_formats`)

- raw text (use `this file <https://github.com/bihealth/scelvis/raw/master/examples/dummy_raw.zip>`_ as an example)
- CellRanger output (zip `this directory <https://github.com/bihealth/scelvis/raw/master/examples/hgmm_1k.raw>`_ as an example)
- loom   

Hitting :guilabel:`Upload` will convert your data to HDF5 and take you to a screen where you can either directly view the converted dataset or download the resulting HDF5 file.

