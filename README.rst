=======================================
SCelVis: Easy Single-Cell Visualization
=======================================

.. image:: https://img.shields.io/conda/dn/bioconda/scelvis.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/scelvis/README.html

.. image:: https://img.shields.io/pypi/pyversions/scelvis.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/scelvis.svg
    :target: https://pypi.python.org/pypi/scelvis

.. image:: https://api.codacy.com/project/badge/Grade/9ee0ec1424c143dfad9977a649f917f7
    :target: https://www.codacy.com/app/bihealth/scelvis?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/scelvis&amp;utm_campaign=Badge_Grade

.. image:: https://api.codacy.com/project/badge/Coverage/9ee0ec1424c143dfad9977a649f917f7
    :target: https://www.codacy.com/app/bihealth/scelvis?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=bihealth/scelvis&amp;utm_campaign=Badge_Coverage

.. image:: https://readthedocs.org/projects/scelvis/badge/?version=latest
    :target: https://scelvis.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://travis-ci.org/bihealth/scelvis.svg?branch=master
    :target: https://travis-ci.org/bihealth/scelvis

.. image:: https://zenodo.org/badge/185944510.svg
    :target: https://zenodo.org/badge/latestdoi/185944510

-------------
Documentation
-------------

The SCelVis documentation at `ReadTheDocs.org <https://scelvis.readthedocs.org>`_ contains comprehensive information including a `Tutorial <https://scelvis.readthedocs.io/en/latest/tutorial_analysis.html>`_.

-----
Movie
-----

.. image:: docs_manual/figures/scelvis_movie.gif
    :height: 400px
    :align: center

------------
Installation
------------

You can install with Pip and Bioconda or run directly via Docker.

.. code-block:: shell

    $ pip install scelvis
    # OR
    $ conda install scelvis
    # OR
    $ docker run ghcr.io/bihealth/scelvis:0.8.8-0 scelvis --help

Look up the latest version (instead of ``0.8.8-0`` to use at `here <ghcr.io/bihealth/scelvis>`_)

---------------------
Building Docker Image
---------------------

1. Push changes to Github (we need a [git tree-ish](https://stackoverflow.com/questions/4044368/what-does-tree-ish-mean-in-git), so a tag works as well as a feature branch).
2. Call ``GIT_TAG=tag-or-branch bash docker/build-docker.sh``.
3. Run with ``docker run ghcr.io/bihealth/scelvis:[version]-0``.
4. Push with ``docker push ghcr.io/bihealth/scelvis:[version]-0``. 
