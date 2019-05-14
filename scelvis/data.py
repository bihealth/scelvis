"""Code for the definition of the data and meta data and loading them from disk.

The types of this modules are used commonly in the scelvis app but the actual ``load_`` functions
should only be accessed through the ``.store`` module such that they can be cached by the app.
This module is unaware of the Dash app.
"""

import os
import os.path
import contextlib
from urllib.parse import urlunparse as _urlunparse
import shutil
import ssl
import typing
from urllib.parse import parse_qs

import anndata
import attr
import fs.path
import fs.tools
from fs.sshfs import SSHFS as SSHFS
from fs.ftpfs import FTPFS
from fs.osfs import OSFS
from fs.tempfs import TempFS
from irods.session import iRODSSession
from irods.ticket import Ticket
from logzero import logger
import pandas as pd
from ruamel.yaml import YAML

from . import settings
from .exceptions import ScelVisException


def redacted_urlunparse(url, redact_with="***"):
    """``urlunparse()`` but redact password."""
    netloc = []
    if url.username:
        netloc.append(url.username)
    if url.password:
        netloc.append(":")
        netloc.append(redact_with)
    if url.hostname:
        if netloc:
            netloc.append("@")
        netloc.append(url.hostname)
    url = url._replace(netloc="".join(netloc))
    return _urlunparse(url)


#: Schemes supported through PyFilesystem
PYFS_SCHEMES = ("file", "ftp", "sftp")


def make_osfs(url):
    """Construct OSFS from url."""
    if url.scheme != "file":
        raise ValueError("Scheme must be == 'file'")
    return OSFS("/")


def make_ftpfs(url):
    """Construct FTPFS from url."""
    if url.scheme != "ftp":
        raise ValueError("Scheme must be == 'ftp'")
    return FTPFS(host=url.hostname, user=url.username, passwd=url.password, port=(url.port or 21))


def make_ssfs(url):
    """Construct SSHFS from url."""
    if url.scheme != "sftp":
        raise ValueError("Scheme must be == 'sftp'")
    return SSHFS(host=url.hostname, user=url.username, passwd=url.password, port=(url.port or 22))


def make_fs(url):
    """Create PyFilesystem FS for the given url."""
    factories = {"file": make_osfs, "ftp": make_ftpfs, "sftp": make_ssfs}
    if url.scheme not in factories:
        raise ValueError("Invalid scheme '%s'" % url.scheme)
    else:
        return factories[url.scheme](url).opendir(url.path)


def create_irods_session(url):
    """Create an iRODS session."""
    if "ssl" in url.scheme:
        ssl_settings = {
            "ssl_context": ssl.create_default_context(
                purpose=ssl.Purpose.SERVER_AUTH, cafile=None, capath=None, cadata=None
            )
        }
    else:
        ssl_settings = {}

    if 'pam' in url.scheme:
        irods_authentication_scheme = 'pam'
    else:
        irods_authentication_scheme = 'native'

    logger.info("Creating iRODS session")
    session = iRODSSession(
        host=url.hostname,
        port=(url.port or 1247),
        user=url.username,
        password=url.password or '',
        irods_authentication_scheme=irods_authentication_scheme,
        zone=(url.path or '/').split('/')[1],
        **ssl_settings
    )
    query = parse_qs(url.query)
    if 'ticket' in query:
        logger.info("Setting ticket for session")
        ticket = Ticket(session, query['ticket'][0])
        ticket.supply()
    return session


@attr.s(auto_attribs=True)
class MetaData:
    """Class to bundle the data loaded for SCelVis."""

    #: ID (= folder name) of the dataset
    id: str
    #: Title of the data set
    title: str
    #: Short title of the data set
    short_title: str
    #: String with Markdown-formatted README.
    readme: str


@attr.s(auto_attribs=True)
class Data:
    """Class to bundle the data loaded for SCelVis."""

    #: String with Markdown-formatted "about" text.
    metadata: MetaData

    #: The raw ad file content.
    ad: anndata.AnnData
    #: The coordinates.
    coords: typing.Dict
    #: The meta information.
    meta: pd.DataFrame
    #: The DGE data
    DGE: pd.DataFrame
    #: The genes data
    genes: pd.Index
    #: The cells data
    cells: pd.Index
    #: The markers data
    markers: pd.DataFrame
    #: Numerical data
    numerical_meta: object
    #: Categorical data
    categorical_meta: object


@contextlib.contextmanager
def download_file(url, path, *more_components):
    """Download file (if necessary) and yield filename if necessary."""
    path = fs.path.join(path, *more_components)
    full_path = fs.path.join(url.path, path)
    logger.info("Attempting to download file %s", url._replace(path=full_path))
    if url.scheme == "file":
        logger.info("No need for download, just use %s", full_path)
        if os.path.exists(full_path):
            yield full_path
        else:
            yield None
    else:
        basename = fs.path.basename(full_path)
        if url.scheme in PYFS_SCHEMES:
            src_fs = make_fs(url)
            # Download file if it exists.
            if not src_fs.exists(path):
                yield None
            else:
                with TempFS() as tmpfs:
                    logger.info("Downloading file %s from %s" % (path, redacted_urlunparse(url)))
                    if url.scheme == "ftps":
                        # NB: download() too slow with built-in methods of fs.sshfs
                        src_fs._wrap_fs._sftp.get(full_path, tmpfs.getospath(basename))
                    else:
                        with open(tmpfs.getospath(basename), "wb") as tmpf:
                            src_fs.download(path, tmpf)
                    logger.info("Download complete.")
                    yield tmpfs.getospath(basename)
                    logger.info("Releasing %s" % tmpfs)
        elif url.scheme.startswith("irods"):
            with create_irods_session(url) as irods_session:
                logger.info("Downloading file...")
                with TempFS() as tmpfs:
                    logger.info("Downloading file %s from %s" % (path, redacted_urlunparse(url)))
                    path_tmp = tmpfs.getospath(basename)
                    collection = irods_session.collections.get(fs.path.dirname(full_path))
                    name = fs.path.basename(full_path)
                    # import ipdb; ipdb.set_trace()
                    for data_object in collection.data_objects:
                        if data_object.name == name:
                            with data_object.open() as inputf:
                                with open(tmpfs.getospath(basename), "wb") as outputf:
                                    shutil.copyfileobj(inputf, outputf, settings.MAX_UPLOAD_DATA_SIZE)
                                    break
                    else:
                        raise ScelVisException("Could not find %s in %s" % (full_path, redacted_urlunparse(url)))
                    logger.info("Download complete.")
                    yield path_tmp
                    logger.info("Releasing %s" % tmpfs)
        else:
            raise ScelVisException("Invalid URL scheme: %s" % url.scheme)


def load_metadata(data_source, identifier):
    """Load metadata from a data_source URL and identifier."""
    logger.info("Loading metadata for %s from %s", identifier, redacted_urlunparse(data_source))
    with download_file(data_source, identifier, "about.md") as about_path:
        with open(about_path, "rt") as inputf:
            header = []
            lines = [line.rstrip() for line in inputf.readlines()]

            # Load meta data, if any.
            if lines and lines[0] and lines[0].startswith("----"):
                for line in lines[1:]:
                    if line.startswith("----"):
                        break
                    else:
                        header.append(line)
                lines = lines[len(header) + 2 :]

            # Get title and short title, finally create MetaData object.
            meta = YAML().load("\n".join(header))
            title = meta.get("title", "Untitled")
            short_title = meta.get("short_title", title or "untitled")
            readme = "\n".join([line.rstrip() for line in lines]).lstrip()

            return MetaData(id=identifier, title=title, short_title=short_title, readme=readme)


def load_data(data_source, identifier):
    """Load the data from the data_source URL and identifier."""
    metadata = load_metadata(data_source, identifier)

    logger.info("Loading anndata for %s from %s", data_source, identifier)
    with download_file(data_source, identifier, "data.h5ad") as path_anndata:
        ad = anndata.read_h5ad(path_anndata)
        coords = {}
        for k in ["X_tsne", "X_umap"]:
            if k in ad.obsm.keys():
                coords[k] = pd.DataFrame(
                    ad.obsm[k],
                    index=ad.obs.index,
                    columns=[k[2:].upper() + str(n + 1) for n in range(ad.obsm[k].shape[1])],
                )
        meta = pd.concat(coords.values(), axis=1).join(ad.obs)
        DGE = ad.to_df().T
        # Separate numerical and categorical columns for later.
        numerical_meta = []
        categorical_meta = []
        for col in meta.columns:
            if pd.api.types.is_numeric_dtype(meta[col]):
                numerical_meta.append(col)
            else:
                categorical_meta.append(col)
        genes = DGE.index
        cells = DGE.columns

    logger.info("Loading markers for %s from %s", data_source, identifier)
    with download_file(data_source, identifier, "markers.csv") as path_markers:
        if path_markers:
            logger.info("Reading markers from %s", path_markers)
            with open(path_markers, "rt") as markerf:
                markers = pd.read_csv(markerf, header=0, index_col=0)
            if "gene" not in markers.columns:
                logger.warn('No "gene" column in %s!', path_markers)
                markers = None
        else:
            markers = None

    return Data(
        metadata=metadata,
        ad=ad,
        coords=coords,
        meta=meta,
        DGE=DGE,
        genes=genes,
        cells=cells,
        markers=markers,
        numerical_meta=numerical_meta,
        categorical_meta=categorical_meta,
    )


def fake_data():
    """Create fake ``Data`` to make Dash validation happy."""
    # TODO: actually use and enable callback traceback again!
    return Data(
        metadata=MetaData(
            id="placeholder", title="placeholder", short_title="placeholder", readme="placeholder"
        ),
        ad=None,
        coords=None,
        meta=pd.DataFrame(data={"shape": [0]}),
        DGE=None,
        genes=None,
        cells=None,
        markers=None,
        numerical_meta=None,
        categorical_meta=None,
    )
