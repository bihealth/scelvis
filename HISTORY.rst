=======
History
=======

------
v0.8.7
------

- fixing problem with ports in HTTP URLs

------
v0.8.6
------

- bug fixes for cell filtering

------
v0.8.5
------

- allow text import with mtx files
- display only most frequent categories
- use `ad.raw` if present
- enable downsampling during conversion
- link to publication

------
v0.8.4
------

- fix bug in custom static folder route
- fix iRODS authentication defaults

------
v0.8.3
------

- allow custom markdowns for home screen
- enable gene-gene plots and change colormap for gene scatter

------
v0.8.2
------

- fix iRODS authentication issues
- update docs for .h5ad input

------
v0.8.1
------

- fixing iframe sources (important if running behind reverse proxy with HTTPS)
- link documentation on home screen

------
v0.8.0
------

- adding sphinx documentation for ReadTheDocs
- upload and convert with IFrames; increase size limit
- add box plots

------
v0.7.3
------

- fix cache timeout error

------
v0.7.2
------

- fix ReferenceError
- updated README and tutorialmovie

------
v0.7.1
------

- improved cache handling
- improved user feedback for filtering & differential expression

------
v0.7.0
------

- added conversion from .loom files
- cell filtering also supports downsampling
- added PBMC dataset hosted on figshare
- added demo movie

------
v0.6.0
------

- cell filtering
- differential expression

------
v0.5.0
------

- upgrades to Dash v1
- fixes to UI, upload and conversion
- avoid creation of dense matrices

------
v0.4.1
------

- Fixing bug with specifying single ``.h5ad`` file as data source.
- Adding ``Dockerfile`` for building Docker images from intermediate versions.

------
v0.4.0
------

- Adding support for HTTP(S) data sources.
- Embedding ``about.md`` information in Anndata file.
- Adding support for passing

------
v0.3.0
------

- Adding example data set.
- Adding nice introduction to start page.
- Adding functionality for creating simple fake data set.
- Making import of ``ruamel_yaml`` more robust still.
- Adding tests.
- Adding Travis CI--based continuous integration tests.

------
v0.2.1
------

- Fixing SFTP support.
- Fixing import of ``ruamel_yaml``.

------
v0.2.0
------

- More refactorization.
- Fixing dependency on ``ruamel-yaml`` to ``ruamel.yaml``.
- Adding conversion feature.
- Adding upload feature.
- Adding support to load from SSHFS, FTP through pyfilesystem (no FTPS support).
- Adding support to load from iRODS, also works via tickets (pass ``?ticket=TICKET`` to the query parameters).

------
v0.1.0
------

Initial release.

- Everything is new!
