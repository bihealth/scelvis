.. _manual-main:

=======
SCelVis
=======

SCelVis is a web-based visualization tool for single-cell sequencing data.

---------------
Getting Started
---------------

- you can find the **publically available web app** at https://scelvis-demo.bihealth.org.
- to get started with your analysis, follow the :ref:`tutorial_analysis`.
- to learn how to get your own data into SCelVis, read about :ref:`tutorial_input_formats`.
- if you are a bioinformatician and interested in running your own SCelVis server, follow the :ref:`tutorial_cli`.


------------------------
Running on your Computer
------------------------

The recommended way to run SCelVis is using Docker:

.. code-block:: shell

    $ docker run quay.io/biocontainers/scelvis:TAG scelvis --help
    $ docker run -p 8050:8050 -v data:/data quay.io/biocontainers/scelvis:TAG scelvis run --data-source /data

In the above, replace ``TAG`` with the latest version that you can find on the corresponding `Quay.io <https://quay.io/repository/biocontainers/scelvis?tag=latest&tab=tags>`_ project.
For example, use ``0.7.3--py_0``.

Alternatively, you can install it using ``pip`` ...

.. code-block:: shell

    $ pip install scelvis
    $ scelvis run --data-source ./data

... or using `Bioconda <http://bioconda.github.io/user/install.html#>`_:


.. code-block:: shell

    $ conda install scelvis
    $ scelvis run --data-source ./data


.. toctree::
    :maxdepth: 1
    :caption: Tutorial
    :name: tutorial
    :hidden:
    :titlesonly:

    tutorial_analysis
    tutorial_convert
    tutorial_input_formats
    tutorial_cli

.. toctree::
    :maxdepth: 1
    :caption: Developer Info
    :name: developer
    :hidden:
    :titlesonly:

    developer_setup
    developer_release

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
