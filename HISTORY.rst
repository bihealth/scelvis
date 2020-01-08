=======
History
=======

------
v0.8
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
