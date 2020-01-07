.. _developer_setup:

===============
Developer Setup
===============

The prerequisites are:

- Python 3, either

    - system-wide installation with ``virtualenv``, or
    - installed with `Conda <https://docs.conda.io/en/latest/>`__.

- `Git LFS <https://git-lfs.github.com/>`__ must be installed for obtaining the example data files.

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
(Note that you must have `Git LFS <https://git-lfs.github.com/>`__ installed to actually obtain the data files).

.. code-block:: shell

    $ git clone git@github.com:bihealth/scelvis.git
    $ cd scelvis
    $ pip install -e .
    $ pip install -r requirements/develop.txt

Afterwards, you can run the visualization web server as follows:

.. code-block:: shell

    $ scelvis run --data-source path/to/data/dir

To explore the datasets provided in the git repository, use ``git lfs fetch`` to download
