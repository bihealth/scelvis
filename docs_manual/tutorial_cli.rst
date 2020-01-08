.. _tutorial_cli:

=====================
Command Line Tutorial
=====================

This tutorial explains the installation of SCelVis on your own computer.

------------
Installation
------------

A Docker container is also available via `Quay.io/Biocontainers <https://quay.io/organization/biocontainers>`_.

.. code-block:: shell

    $ docker run quay.io/biocontainers/scelvis:TAG scelvis --help
    $ docker run -p 8050:8050 -v data:/data quay.io/biocontainers/scelvis:TAG scelvis run --data-source /data

look up the latest ``TAG`` to use at `here <https://quay.io/repository/biocontainers/scelvis?tab=tags>`_, e.g.,

.. code-block:: shell

    $ docker run quay.io/biocontainers/scelvis:0.7.3--py_0 scelvis --help
    $ docker run -p 8050:8050 -v data:/data quay.io/biocontainers/scelvis:0.7.0--py_0 scelvis run --data-source /data

For installation, the only prerequisite is Python 3, everything else will be installed together with the ``scelvis`` package.
You can install SCelVis and its dependencies using ``pip`` or through ``conda``:

.. code-block:: shell

    $ pip install scelvis
    # OR
    $ conda install scelvis

The most robust way is to use Docker, though.

For the sake of simplicity, we will give the executable as ``scelvis``.
When using Docker, use the corresponding prefix.

---------------
Running SCelVis
---------------


You can run the SCelVis Web server with ``scelvis run``.

.. code-block:: shell

    $ scelvis run --data-source /path/to/scelvis/examples/hgmm_1k.h5ad
    $ scelvis run --data-source https://files.figshare.com/18037739/pbmc.h5ad


and then point your browser to http://0.0.0.0:8050/ or http://localhost:8050/.

We provide the following two example HDF5 files:

- `hgmm_1k.h5ad <https://github.com/bihealth/scelvis/raw/master/examples/hgmm_1k.h5ad>`_
- `pbmc.h5ad <https://files.figshare.com/18037739/pbmc.h5ad>`_

The first command will make SCelVis directly serve the given ``hgmm_1k.h5ad`` file while the second command will first download the file ``pbmc.h5ad`` from the given URL and then complete web server startup.
The files given as ``--data-source`` on server startup will be displayed in the top right :guilabel:`Go To` menu.

------------
Data Sources
------------

Data sources can be:

paths
    e.g., ``relative/paths`` or ``/absolute/paths`` or ``file://url/paths``

HTTP(S) URLs
    e.g., ``https://user:password@host/path/to/data``.

S3 URLs
    ``s3://bucket/path``, optionally ``s3://key:token@bucket/path``.

SFTP URLs
    e.g., ``sftp://user:password@host/path/to/data``

FTP URLs
    e.g., ``ftp://user:password@host/path/to/data`` (sadly encryption is not supported by the underlying library `PyFilesystem2 <https://github.com/PyFilesystem/pyfilesystem2>`__.

iRODS URLS
    e.g., ``irods://user:password@host/zoneName/path/to/data``

    - Enable SSL via ``irods+ssl``
    - Switch to PAM authentication with ``irods+pam`` (you can combine this with ``+ssl`` in any order)
    - Enable ticket access by appending ``?ticket=TICKET``.

Data sources can either point to HDF5 files directly or to directories containing multiple HDF5 files.

-------------
Configuration
-------------

You can configure the web server by passing command line arguments (run ``scelvis run --help`` for all available options).
Also, you can use the following environment variables:

---------------------
Environment Variables
---------------------

You can use the following environment variables to configure the server.

``SCELVIS_DATA_SOURCES``
    semicolon-separated list of data sources

``SCELVIS_HOST``
    host specification for web server to listen on

``SCELVIS_PORT``
    port for web server to listen on

``SCELVIS_CACHE_DIR``
    directory to use for the cache (default is to create a temporary directory)

``SCELVIS_CACHE_REDIS_URL``
    enable caching with REDIS and provide connection URL

``SCELVIS_CACHE_DEFAULT_TIMEOUT``
    cache lifetime coverage

``SCELVIS_CACHE_PRELOAD_DATA``
    will preload all data at startup

``SCELVIS_UPLOAD_DIR``
    the directory to store uploaded data sets in (default is to create a temporary directory)

``SCELVIS_UPLOAD_DISABLED``
    set to "0" to disable upload feature

``SCELVIS_CONVERSION_DISABLED``
    set to "0" to disable the conversion feature

``SCELVIS_URL_PREFIX``
    set if you want to run scelvis below a non-root path (e.g., behind a reverse proxy)

