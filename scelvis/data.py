"""Code for the definition of the data and meta data and loading them from disk.

The types of this modules are used commonly in the scelvis app but the actual ``load_`` functions
should only be accessed through the ``.store`` module such that they can be cached by the app.
This module is unaware of the Dash app.
"""

import fs
from fs.tempfs import TempFS
import typing

import anndata
import attr
from logzero import logger
import pandas as pd
from ruamel.yaml import YAML

from . import settings


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


def load_metadata(file_system, path):
    """Load metadata from a dataset directory."""
    logger.info("Loading %s from %s", path, file_system)
    with file_system.open(fs.path.join(path, settings.ABOUT_FILENAME), "rt") as inputf:
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

        return MetaData(
            id=fs.path.basename(path), title=title, short_title=short_title, readme=readme
        )


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


def load_data(file_system, datadir):
    """Load the data from the given directory."""
    logger.info("Loading metadata from %s %s", file_system, datadir)
    metadata = load_metadata(file_system, datadir)

    ad_file = fs.path.join(datadir, "data.h5ad")
    logger.info("Reading all data from %s", ad_file)
    with TempFS() as tmpfs:
        logger.info("Downloading file %s", ad_file)
        fs.copy.copy_file(file_system, ad_file, tmpfs, "tmp.h5ad")
        ad = anndata.read_h5ad(tmpfs.getospath("tmp.h5ad"))
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

    marker_file = fs.path.join(datadir, "markers.csv")
    if file_system.exists(marker_file):
        logger.info("Reading markers from %s", marker_file)
        with file_system.open(marker_file) as markerf:
            markers = pd.read_csv(markerf, header=0, index_col=0)
        if "gene" not in markers.columns:
            logger.warn('No "gene" column in %s!', marker_file)
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
