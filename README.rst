=======================================
SCelVis: Easy Single-Cell Visualization
=======================================


.. image:: https://img.shields.io/conda/dn/bioconda/scelvis.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/scelvis/README.html

.. image:: https://img.shields.io/pypi/pyversions/scelvis.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/scelvis.svg
    :target: https://pypi.python.org/pypi/scelvis


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
    $ docker run -p 8050:8050 -v data:/data bihealth/scelvis:latest run --data-source /data

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

    $ scelvis run --data-source data/project
    # OR
    $ scelvis run --data-source data

------------
Data Sources
------------

Data sources can be:

- paths, e.g., ``relative/paths`` or ``/absolute/paths`` or ``file://url/paths``
- SFTP URLs, e.g., ``sftp://user:password@host/path/to/data``
- FTP URLs, e.g., ``ftp://user:password@host/path/to/data`` (sadly encryption is not supported by the underlying library `PyFilesystem2 <https://github.com/PyFilesystem/pyfilesystem2>`_.
- iRODS URLS, e.g., ``irods://user:password@host/zoneName/path/to/data``
    - Enable SSL via ``irods+ssl``
    - Switch to PAM authentication with ``irods+pam`` (you can combine this with ``+ssl`` in any order)
    - Enable ticket access by appending ``?ticket=TICKET``.

Data sources can either point to directories that contain the ``about.md`` string directly (data sets) contain multiple data set directories.
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
- ``SCELVIS_UPLOAD_DIR`` -- the directory to store uploaded data sets in (default is to create a temporary directory)
- ``SCELVIS_UPLOAD_DISABLED`` -- set to "0" to disable upload feature
- ``SCELVIS_CONVERSION_DISABLED`` -- set to "0" to disable the conversion feature

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

    $ scelvis run --data-source path/to/data/dir

Releasing Packages
==================

For the `PyPi package <https://pypi.org/project/scelvis/>`_:

.. code-block:: shell

    $ python setup.py sdist
    $ twine upload --repository-url https://test.pypi.org/legacy/ dist/scelvis-*.tar.gz
    $ twine upload dist/scelvis-*.tar.gz

For the Bioconda package, see `the great documentation <http://bioconda.github.io/updating.html>`_.
The Docker image will automatically be created as a BioContainer when the Bioconda package is built.
