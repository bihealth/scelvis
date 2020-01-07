.. _developer_release:

==================
Releasing Packages
==================

For the `PyPi package <https://pypi.org/project/scelvis/>`__:

.. code-block:: shell

    $ python setup.py sdist
    $ twine upload --repository-url https://test.pypi.org/legacy/ dist/scelvis-*.tar.gz
    $ twine upload dist/scelvis-*.tar.gz

For the Bioconda package, see `the great documentation <http://bioconda.github.io/updating.html>`__.
The Docker image will automatically be created as a BioContainer when the Bioconda package is built.
