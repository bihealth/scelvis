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
from urllib.parse import parse_qs

import anndata
import re
import attr
import fs.path
import fs.tools
from fs.sshfs import SSHFS as SSHFS
from fs.ftpfs import FTPFS
from fs.osfs import OSFS
import s3fs
from fs.tempfs import TempFS
from irods.session import iRODSSession
from irods.ticket import Ticket
from logzero import logger
import numpy as np
import scipy.sparse
import pandas as pd
import requests

from . import settings
from .exceptions import ScelVisException

#: Identifier for fake data
FAKE_DATA_ID = "builtin-fake-data"


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
                purpose=ssl.Purpose.CLIENT_AUTH, cafile=None, capath=None, cadata=None
            )
        }
    else:
        ssl_settings = {}

    if "pam" in url.scheme:
        irods_authentication_scheme = "pam"
    else:
        irods_authentication_scheme = "native"

    logger.info("Creating iRODS session")
    session = iRODSSession(
        host=url.hostname,
        port=(url.port or 1247),
        user=url.username,
        password=url.password or "",
        zone=(url.path or "/").split("/")[1],
        irods_authentication_scheme=irods_authentication_scheme,
        irods_client_server_negotiation=settings.IRODS_CLIENT_SERVER_NEGOTIATION,
        irods_client_server_policy=settings.IRODS_CLIENT_SERVER_POLICY,
        irods_ssl_verify_server=settings.IRODS_SSL_VERIFY_SERVER,
        irods_encryption_algorithm=settings.IRODS_ENCRYPTION_ALGORITHM,
        irods_encryption_key_size=settings.IRODS_ENCRYPTION_KEY_SIZE,
        irods_encryption_num_hash_rounds=settings.IRODS_ENCRYPTION_NUM_HASH_ROUNDS,
        irods_encryption_salt_size=settings.IRODS_ENCRYPTION_SALT_SIZE,
        **ssl_settings
    )
    query = parse_qs(url.query)
    if "ticket" in query:
        logger.info("Setting ticket for session")
        ticket = Ticket(session, query["ticket"][0])
        ticket.supply()
    return session


@attr.s(auto_attribs=True)
class MetaData:
    """Class to bundle the data loaded for SCelVis."""

    #: folder name of the dataset
    name: str
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

    # The raw ad file content.
    ad: anndata.AnnData
    #: The genes data
    genes: pd.Index
    #: The cells data
    cells: pd.Index
    #: The markers data
    markers: pd.DataFrame
    #: coordinates
    coords: object
    #: Numerical data
    numerical_meta: object
    #: Categorical data
    categorical_meta: object


@contextlib.contextmanager
def download_file(url, path=None, *more_components):
    """Download file (if necessary) and yield filename if necessary."""
    if path:
        path = fs.path.join(path, *more_components)
        full_path = fs.path.join(url.path, path)
    else:
        full_path = url.path
        path = fs.path.basename(url.path)
        url = url._replace(path=fs.path.dirname(url.path))
    if url.scheme == "file":
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
                    with open(tmpfs.getospath(basename), "wb") as tmpf:
                        src_fs.download(path, tmpf)
                    logger.info("Download complete.")
                    yield tmpfs.getospath(basename)
                    logger.info("Releasing %s" % tmpfs)
        elif url.scheme == "s3":
            logger.info("Connecting via S3...")
            anon = url.username is None and url.password is None
            s3 = s3fs.S3FileSystem(anon=anon, key=url.username, secret=url.password)
            with TempFS() as tmpfs:
                logger.info("Downloading file %s from %s" % (path, redacted_urlunparse(url)))
                with s3.open("%s/%s" % (url.hostname, path), "rb") as inputf:
                    with open(tmpfs.getospath(basename), "wb") as outputf:
                        shutil.copyfileobj(inputf, outputf, settings.MAX_UPLOAD_DATA_SIZE)
                logger.info("Download complete.")
                yield tmpfs.getospath(basename)
                logger.info("Releasing %s" % tmpfs)
        elif url.scheme.startswith("http"):
            logger.info("Downloading via HTTP(S)...")
            with TempFS() as tmpfs:
                logger.info("Downloading file %s from %s" % (path, redacted_urlunparse(url)))
                r = requests.get(_urlunparse(url._replace(path=full_path)), allow_redirects=True)
                r.raise_for_status()
                with open(tmpfs.getospath(basename), "wb") as outputf:
                    outputf.write(r.content)
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
                    for data_object in collection.data_objects:
                        if data_object.name == name:
                            with data_object.open() as inputf:
                                with open(tmpfs.getospath(basename), "wb") as outputf:
                                    shutil.copyfileobj(
                                        inputf, outputf, settings.MAX_UPLOAD_DATA_SIZE
                                    )
                                    break
                    else:
                        raise ScelVisException(
                            "Could not find %s in %s" % (full_path, redacted_urlunparse(url))
                        )
                    logger.info("Download complete.")
                    yield path_tmp
                    logger.info("Releasing %s" % tmpfs)
        else:
            raise ScelVisException("Invalid URL scheme: %s" % url.scheme)


def load_data(data_source, identifier):
    """Load the data from the data_source URL and identifier."""
    if identifier == FAKE_DATA_ID:
        return fake_data()

    logger.info("Loading anndata for %s from %s", identifier, redacted_urlunparse(data_source))
    with download_file(data_source) as path_anndata:
        ad = anndata.read_h5ad(path_anndata)
        # Extract the meta data
        metadata = MetaData(
            name=identifier,
            title="".join(ad.uns["about_title"]) if "about_title" in ad.uns_keys() else identifier,
            short_title="".join(ad.uns["about_short_title"])
            if "about_short_title" in ad.uns_keys()
            else identifier,
            readme="".join(ad.uns["about_readme"]) if "about_readme" in ad.uns_keys() else None,
        )
        # Separate numerical and categorical columns for later.
        numerical_meta = []
        categorical_meta = []
        for col in ad.obs_keys():
            if pd.api.types.is_numeric_dtype(ad.obs[col]):
                numerical_meta.append(col)
            else:
                categorical_meta.append(col)
                # use only the top 40 most frequent categories to avoid legend blow-up
                if len(ad.obs[col].unique()) >= 40:
                    logger.warn("keeping only the top 40 most frequent items in column %s", col)
                    keep = ad.obs[col].value_counts().index[:39].astype(str)
                    if "other" in keep:
                        logger.warn("value 'other' already present in column %s", col)
                    tmp = ad.obs[col].astype(str)
                    tmp[~np.isin(tmp, keep)] = "other"
                    ad.obs[col] = tmp.astype("category")
        # add coordinates to obs
        coords = {}
        for k in ad.obsm.keys():
            kk = re.sub("^X_|_cell_embeddings$", "", k).upper()
            ndim = min(ad.obsm[k].shape[1], 3)
            coords[kk] = pd.DataFrame(
                ad.obsm[k][:, :ndim],
                index=ad.obs.index,
                columns=[kk + str(n + 1) for n in range(ndim)],
            )
        if len(coords) > 0:
            coords = pd.concat(coords.values(), axis=1)
            ad.obs = coords.join(ad.obs)
        # check if ad contains raw
        if ad.raw:
            logger.info("found ad.raw in object, taking expression values from there")
            ad.X = ad.raw.X
            ad.var = ad.raw.var
        if not scipy.sparse.issparse(ad.X):
            logger.warn("expression matrix is not sparse! this could affect performance!")
        genes = ad.var_names
        cells = ad.obs_names
        markers = {}
        for c in ad.uns_keys():
            if c.startswith("marker"):
                markers[c.replace("marker_", "")] = ad.uns[c]
        if "gene" in markers:
            markers = pd.DataFrame(markers)
        elif len(markers) > 0:
            logger.warn('No "gene" column in h5 file!')
            markers = None
        elif "rank_genes_groups" in ad.uns_keys():
            logger.info("using markers from scanpy")
            markers = pd.DataFrame(
                dict(
                    (k, pd.DataFrame(ad.uns["rank_genes_groups"][k]).stack().astype(t))
                    for k, t in [
                        ("scores", float),
                        ("names", str),
                        ("logfoldchanges", float),
                        ("pvals_adj", float),
                    ]
                )
            )
            markers.columns = ["score", "gene", "LFC", "padj"]
            markers["cluster"] = markers.index.get_level_values(1)
            markers = (
                markers[markers["gene"].isin(ad.var_names)]
                .sort_values(["cluster", "padj"])
                .groupby("cluster")
                .head(10)
            )
        else:
            logger.warn("No markers in h5 file!")
            markers = None

    return Data(
        metadata=metadata,
        ad=ad,
        genes=genes,
        cells=cells,
        markers=markers,
        coords=list(coords.columns) if len(coords) > 0 else [],
        numerical_meta=numerical_meta,
        categorical_meta=categorical_meta,
    )


def fake_data(seed=42):
    """Create fake ``Data`` to make Dash validation happy."""
    np.random.seed(seed)

    ngenes = 100
    ncells = 50
    genes = ["gene_{0}".format(i + 1) for i in range(ngenes)]
    cells = ["cell_{0}".format(i + 1) for i in range(ncells)]
    dge = np.random.negative_binomial(1, 0.5, size=(ncells, ngenes))
    meta = pd.DataFrame(
        {
            "TSNE_1": np.random.random(size=ncells),
            "TSNE_2": np.random.random(size=ncells),
            "cluster": pd.Categorical(
                ["cluster_{0}".format("ABCD"[i]) for i in np.random.randint(4, size=ncells)]
            ),
            "sample": pd.Categorical(
                ["sample_{0}".format("ABC"[i]) for i in np.random.randint(3, size=ncells)]
            ),
            "n_genes": (dge > 0).sum(axis=1),
            "n_counts": dge.sum(axis=1),
        },
        index=cells,
    )
    ad = anndata.AnnData(dge, obs=meta, var=pd.DataFrame([], index=genes))
    coords = ["TSNE_1", "TSNE_2"]
    numerical_meta = ["n_genes", "n_counts"]
    categorical_meta = ["cluster", "sample"]
    markers = pd.DataFrame(
        {
            "cluster": np.random.choice(meta["cluster"], 12),
            "gene": np.random.choice(genes, 12),
            "padj": np.random.rand(12),
        }
    )

    return Data(
        metadata=MetaData(
            name=FAKE_DATA_ID, title="fake data", short_title="fake", readme="fake data"
        ),
        ad=ad,
        genes=genes,
        cells=cells,
        markers=markers,
        coords=coords,
        numerical_meta=numerical_meta,
        categorical_meta=categorical_meta,
    )
