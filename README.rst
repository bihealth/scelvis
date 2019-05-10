====================================
Single-Cell Visualization using Dash
====================================

------------
Installation
------------

The only prerequisite is Python 3, everything else will be installed together with the ``scelvis`` package.

You can install SCelVis and its dependencies using ``pip`` or through ``conda``:

.. code-block:: shell

    $ pip install scelvis
    # OR
    $ conda install scelvis

A Docker container is also available:

.. code-block:: shell

    $ docker run bihealth/scelvis:latest --help
    $ docker run -p 8050:8050 -v data:/data bihealth/scelvis:latest run --data-dir /data

-------------------
Preparing Your Data
-------------------

Each data set consists of an HDF5 file called ``data.h5ad`` and a dataset description file ``about.md``.
The HDF5 file is an `anndata <https://anndata.readthedocs.io/en/latest/index.html>`_ object that stores gene expression (sparse CSR matrix) and meta data with very fast read access.
You can use the ``scelvis convert`` command for converting your single-cell pipeline output into an appropriate HDF5 file.
The ``about.md`` file should look as follows:

::

    ----
    title: An Optional Long Data Set Title
    short_title: optional short title
    ----

    A verbose description of the data in Markdown format.

A directory containing both an ``data.h5ad`` and an ``about.md`` file is a **dataset directory**.
For the input you can either specify one dataset directory or a **data directory** containing multiple dataset directories.

You can convert your single-cell transcriptome analysis pipeline as follows.
This does no further processing except log-normalization and uses PCA, tSNE, and clustering performed by ``cellranger``

.. code-block:: shell

    $ mkdir -p data/project
    $ scelvis convert --input-dir cellranger-out --output-dir data/project
    $ cat <<EOF
    ----
    title: My Project
    ----

    This is my project data.
    EOF

    $ tree data
    data
    ├── other
    │   ├── about.md
    │   └── data.h5ad
    └── project
        ├── about.md
        └── data.h5ad

Note that right now only CellRanger output is supported.

---------------------
Visualizing Your Data
---------------------

.. code-block::

    $ tree data
    data
    ├── other
    │   ├── about.md
    │   └── data.h5ad
    └── project
        ├── about.md
        └── data.h5ad

    $ scelvis run --data-dir data/project
    # OR
    $ scelvis run --data-dir data

---------------
Developer Setup
---------------

The prerequisites are:

- Python 3, either
    - system-wide installation with ``virtualenv``, or
    - installed with `Conda <https://docs.conda.io/en/latest/>`_.

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

.. code-block:: shell

    $ git clone git@github.com:bihealth/scelvis.git
    $ cd scelvis
    $ pip install -e .
    $ pip install -r requirements/develop.txt

Afterwards, you can run the visualization web server as follows:

.. code-block:: shell

    $ scelvis run --data-dir path/to/data/dir
