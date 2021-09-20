"""Conversion from single-cell pipeline output to SCelVis HDF5 files.
"""

import contextlib
import logging
import os
import typing
import json

import anndata
import attr
import scanpy as sc
import scipy.sparse
import re
import scipy.io
import pandas as pd
from logzero import logger

try:
    from ruamel_yaml import YAML
except ModuleNotFoundError:
    from ruamel.yaml import YAML

from .exceptions import ScelVisException


@attr.s(auto_attribs=True)
class About:
    """Class to bundle the information from the ``about.md`` file."""

    #: Title of the data set
    title: str
    #: Short title of the data set
    short_title: str
    #: String with Markdown-formatted README.
    readme: str


def load_about_md(path):
    """Load contents of ``about.md`` file with title and subtitle from ``path``."""
    logger.info("Loading information from about.md at %s", path)
    with open(path, "rt") as inputf:
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

        # Get title and short title, finally create About object.
        meta = YAML().load("\n".join(header)) or {}
        title = meta.get("title", "Untitled")
        short_title = meta.get("short_title", title or "untitled")
        readme = "\n".join([line.rstrip() for line in lines]).lstrip()

        return About(title=title, short_title=short_title, readme=readme)


@contextlib.contextmanager
def with_log_level(logger, level):
    """Set log level to warning temporarily."""
    old_level = logger.level
    logger.setLevel(level)
    yield
    logger.setLevel(old_level)


class CellRangerConverter:
    """Conversion of cell ranger data to SCelVis HDF5 file."""

    def __init__(self, args):
        #: Arguments as parsed by ``argparse.ArgumentParser``.`
        self.args = args

    def run(self):
        """Perform the conversion."""
        logger.info("Starting conversion from CellRanger output to SCelVis HDF5")
        logger.info("Settings: %s", vars(self.args))
        tsne = self._load_tsne()
        clustering = self._load_clustering()
        diffexp = self._load_diffexp()
        ad = self._load_expression(clustering, tsne, diffexp)
        about = self._load_about()

        if self.args.split_species:
            ad = CellRangerConverter._split_species(ad)
        if self.args.split_samples:
            ad = self._split_samples(ad)
        if not self.args.use_raw:
            ad = CellRangerConverter._normalize_filter_dge(ad)

        # Write about meta data into anndata.
        ad.uns["about_title"] = about.title
        ad.uns["about_short_title"] = about.short_title
        ad.uns["about_readme"] = about.readme

        self._write_output(ad)
        logger.info("All done. Have a nice day!")

    def _load_about(self):
        """Load the about.md file (if any)."""
        if self.args.about_md:
            return load_about_md(self.args.about_md)
        elif os.path.exists(os.path.join(self.args.indir, "about.md")):
            return load_about_md(os.path.join(self.args.indir, "about.md"))
        else:
            filename = os.path.basename(self.args.indir)
            return About(title=filename, short_title=filename, readme=None)

    def _load_tsne(self):
        tsne_file = os.path.join(
            self.args.indir, "analysis", "tsne", "2_components", "projection.csv"
        )
        if not os.path.isfile(tsne_file):
            raise ScelVisException("cannot find tSNE output at %s" % tsne_file)  # pragma nocover
        else:
            logger.info("Reading tSNE output from %s", tsne_file)
            return pd.read_csv(tsne_file, header=0, index_col=0)

    def _load_clustering(self):
        clustering_file = os.path.join(
            self.args.indir, "analysis", "clustering", "graphclust", "clusters.csv"
        )
        if not os.path.isfile(clustering_file):
            raise ScelVisException(
                "cannot find clustering output at %s " % clustering_file
            )  # pragma nocover
        else:
            logger.info("Reading clustering output from %s", clustering_file)
            clustering = pd.read_csv(clustering_file, header=0)
            nc = pd.DataFrame({"nclust": clustering["Cluster"].value_counts()})
            clustering = clustering.set_index("Cluster")
            clustering["nclust"] = nc
            clustering["Cluster"] = "Cluster_" + clustering.index.astype(str)
            clustering = clustering.set_index("Barcode")["Cluster"]
            return clustering.astype("category")

    def _load_diffexp(self):
        diffexp_file = os.path.join(
            self.args.indir, "analysis", "diffexp", "graphclust", "differential_expression.csv"
        )
        if not os.path.isfile(diffexp_file):
            raise ScelVisException(
                "cannot find differential expression output at " + diffexp_file
            )  # pragma nocover
        else:
            logger.info("Reading differential expression output from %s", diffexp_file)
            diffexp = pd.read_csv(diffexp_file, header=0, index_col=[0, 1])
            diffexp.columns = (
                diffexp.columns.str.replace("Cluster ", "Cluster_")
                .str.split(n=1, expand=True)
                .rename(["Cluster", "Obs"])
            )
            diffexp = diffexp.stack(level=0).reset_index()
            diffexp.columns = ["GeneID", "gene", "Cluster", "p_adj", "log2_fc", "mean_counts"]
            return diffexp

    def _load_expression(self, clustering, tsne, diffexp):
        expression_file_v3 = os.path.join(self.args.indir, "filtered_feature_bc_matrix.h5")
        expression_file_v2 = os.path.join(self.args.indir, "filtered_gene_bc_matrices_h5.h5")
        if os.path.isfile(expression_file_v3):
            expression_file = expression_file_v3
        elif os.path.isfile(expression_file_v2):
            expression_file = expression_file_v2
        else:
            raise ScelVisException("cannot find expression file at %s" % self.args.indir)
        logger.info("Reading gene expression from %s", expression_file)
        with with_log_level(anndata.utils.logger, logging.WARN):
            ad = sc.read_10x_h5(expression_file)
        ad.var_names_make_unique()
        logger.info("Combining meta data")
        ad.obs["cluster"] = clustering
        ad.obs["n_counts"] = ad.X.sum(1).A1
        ad.obs["n_genes"] = (ad.X > 0).sum(1).A1
        logger.info("Adding coordinates")
        ad.obsm["X_tsne"] = tsne.values
        logger.info("Saving top %d markers per cluster", self.args.nmarkers)
        markers = (
            diffexp[(diffexp["p_adj"] < 0.05) & (diffexp["log2_fc"] > 0)]
            .drop("GeneID", axis=1)
            .sort_values(["Cluster", "p_adj"])
            .groupby("Cluster")
            .head(self.args.nmarkers)
        )
        for col in markers.columns:
            ad.uns["marker_" + col] = markers[col].values

        return ad

    def _split_species(ad):
        logger.info("Determining species mixing")
        species = ad.var_names.str.split("_", n=1).str[0]
        for sp in species.unique():
            ad.obs["nUMI_" + sp] = ad.X[:, species == sp].sum(1).A1
        cols = "nUMI_" + species.unique()
        ratio = ad.obs[cols].divide(ad.obs[cols].sum(axis=1), axis=0)
        ad.obs["species"] = (ratio > 0.95).idxmax(axis=1).str.split("_").str[1]
        ad.obs["species"][ratio.max(axis=1) < 0.95] = "doublet"
        return ad

    def _split_samples(self, ad):
        logger.info("splitting samples")
        samples = ad.obs_names.str.rsplit("-", n=1).str[1]
        summary_json = os.path.join(self.args.indir, "summary.json")
        if not os.path.isfile(summary_json):
            raise ScelVisException("could not find summary.json")
        with open(summary_json) as json_file:
            batches = json.load(json_file)["batches"]
        sample_table = pd.Series(batches, index=map(str, range(1, len(batches) + 1)))
        ad.obs["sample"] = sample_table.loc[samples].values
        return ad

    def _normalize_filter_dge(ad):
        logger.info("Normalizing and filtering DGE")
        sc.pp.filter_cells(ad, min_genes=100)
        sc.pp.filter_genes(ad, min_cells=5)
        sc.pp.normalize_per_cell(ad, counts_per_cell_after=1.0e4)
        sc.pp.log1p(ad)
        return ad

    def _write_output(self, ad):
        if self.args.ncells:
            logger.info("Sampling %s cells from anndata object", self.args.ncells)
            sc.pp.subsample(ad, n_obs=self.args.ncells)
        logger.info("Saving anndata object to %s", self.args.out_file)
        ad.write(self.args.out_file)


class TextConverter:
    """Conversion of text data to SCelVis HDF5 file."""

    def __init__(self, args):
        #: Arguments as parsed by ``argparse.ArgumentParser``.`
        self.args = args

    def run(self):
        """Perform the conversion."""
        logger.info("Starting conversion from text input to SCelVis HDF5")
        logger.info("Settings: %s", vars(self.args))
        coords = self._load_coords()
        annotation = self._load_annotation()
        markers = self._load_markers()
        ad = self._load_expression(coords, annotation, markers)
        about = self._load_about()

        # Write about meta data into anndata.
        ad.uns["about_title"] = about.title
        ad.uns["about_short_title"] = about.short_title
        ad.uns["about_readme"] = about.readme

        self._write_output(ad)
        logger.info("All done. Have a nice day!")

    def _load_about(self):
        """Load the about.md file (if any)."""
        if self.args.about_md:
            return load_about_md(self.args.about_md)
        elif os.path.exists(os.path.join(self.args.indir, "about.md")):
            return load_about_md(os.path.join(self.args.indir, "about.md"))
        else:
            filename = os.path.basename(self.args.indir)
            return About(title=filename, short_title=filename, readme=None)

    def _load_coords(self):
        coords_file = os.path.join(self.args.indir, "coords.tsv")
        if not os.path.isfile(coords_file):
            raise ScelVisException("cannot find coords file at %s" % coords_file)
        else:
            logger.info("Reading coords from %s", coords_file)
            return pd.read_csv(coords_file, header=0, index_col=0, sep="\t")

    def _load_annotation(self):
        annotation_file = os.path.join(self.args.indir, "annotation.tsv")
        if not os.path.isfile(annotation_file):
            raise ScelVisException("cannot find cell annotation at %s " % annotation_file)
        else:
            logger.info("Reading cell annotation from %s", annotation_file)
            return pd.read_csv(annotation_file, header=0, index_col=0, sep="\t")

    def _load_markers(self):
        if self.args.markers:
            logger.info("Reading markers from %s", self.args.markers)
            return pd.read_csv(self.args.markers, header=0, sep="\t")
        marker_file = os.path.join(self.args.indir, "markers.tsv")
        if not os.path.isfile(marker_file):
            logger.info("no markers in %s!", self.args.indir)
            return pd.DataFrame()
        else:
            logger.info("Reading markers from %s", marker_file)
            return pd.read_csv(marker_file, header=0, sep="\t")

    def _load_expression(self, coords, annotation, markers):
        expression_tsv = os.path.join(self.args.indir, "expression.tsv.gz")
        expression_mtx = os.path.join(self.args.indir, "expression.mtx")
        if os.path.isfile(expression_tsv):
            logger.info("Reading gene expression from %s", expression_tsv)
            DGE = pd.read_csv(expression_tsv, header=0, index_col=0, sep="\t")
            X = scipy.sparse.csr_matrix(DGE.values.T)
            cells = DGE.columns
            genes = DGE.index
        elif os.path.isfile(expression_mtx):
            cells_tsv = os.path.join(self.args.indir, "barcodes.tsv")
            genes_tsv = os.path.join(self.args.indir, "genes.tsv")
            if not (os.path.isfile(cells_tsv) and os.path.isfile(genes_tsv)):
                raise ScelVisException(
                    "expression.mtx present at %s, but barcodes.tsv or genes.tsv is missing"
                    % expression_tsv
                )
            logger.info("Reading gene expression from %s", expression_mtx)
            X = scipy.io.mmread(expression_mtx).T.tocsr()
            cells = pd.read_csv(cells_tsv, header=None).squeeze().values
            genes = pd.read_csv(genes_tsv, header=None).squeeze().values
        else:
            raise ScelVisException(
                "cannot find expression data at %s or %s" % (expression_tsv, expression_mtx)
            )
        logger.info("Combining data")
        ad = sc.AnnData(X=X, obs=annotation.loc[cells], var=pd.DataFrame([], index=genes))
        coords_types = set([re.sub("[-_][0-9]*$", "", c).upper() for c in coords.columns])
        for ct in coords_types:
            ct_cols = [c for c in coords.columns if ct in c.upper()]
            ad.obsm["X_" + ct] = coords.loc[cells, ct_cols].values
        for col in markers.columns:
            ad.uns["marker_" + col] = markers[col].values

        return ad

    def _write_output(self, ad):
        if self.args.ncells:
            logger.info("Sampling %s cells from anndata object", self.args.ncells)
            sc.pp.subsample(ad, n_obs=self.args.ncells)
        logger.info("Saving anndata object to %s", self.args.out_file)
        ad.write(self.args.out_file)


class LoomConverter:
    """Conversion of loom files to SCelVis HDF5 file."""

    def __init__(self, args):
        #: Arguments as parsed by ``argparse.ArgumentParser``.`
        self.args = args

    def run(self):
        """Perform the conversion."""
        logger.info("Starting conversion from loom input to SCelVis HDF5")
        logger.info("Settings: %s", vars(self.args))
        markers = self._load_markers()
        ad = self._load_loom(markers)
        about = self._load_about()

        ad.uns["about_title"] = about.title
        ad.uns["about_short_title"] = about.short_title
        ad.uns["about_readme"] = about.readme

        self._write_output(ad)
        logger.info("All done. Have a nice day!")

    def _load_about(self):
        """Load the about.md file (if any)."""
        if self.args.about_md:
            return load_about_md(self.args.about_md)
        elif os.path.exists(os.path.join(self.args.indir, "about.md")):
            return load_about_md(os.path.join(self.args.indir, "about.md"))
        else:
            filename = os.path.basename(self.args.indir)
            return About(title=filename, short_title=filename, readme=None)

    def _load_markers(self):
        if self.args.markers:
            logger.info("Reading markers from %s", self.args.markers)
            return pd.read_csv(self.args.markers, header=0, sep="\t")
        marker_file = os.path.join(self.args.indir, "markers.tsv")
        if not os.path.isfile(marker_file):
            logger.info("no markers in %s!", self.args.indir)
            return pd.DataFrame()
        else:
            logger.info("Reading markers from %s", marker_file)
            return pd.read_csv(marker_file, header=0, sep="\t")

    def _load_loom(self, markers):
        if self.args.indir.endswith(".loom"):
            loom_file = self.args.indir
        else:
            loom_file = os.path.join(self.args.indir, "data.loom")
        if not os.path.isfile(loom_file):
            raise ScelVisException("cannot find loom file at %s" % loom_file)
        logger.info("Reading data from %s", loom_file)
        ad = sc.read_loom(loom_file, X_name="spliced" if self.args.use_raw else "norm_data")
        for layer in list(ad.layers.keys()):
            logger.info("Removing unused layer %s" % layer)
            del ad.layers[layer]
        for col in markers.columns:
            ad.uns["marker_" + col] = markers[col].values

        return ad

    def _write_output(self, ad):
        if self.args.ncells:
            logger.info("Sampling %s cells from anndata object", self.args.ncells)
            sc.pp.subsample(ad, n_obs=self.args.ncells)
        logger.info("Saving anndata object to %s", self.args.out_file)
        ad.write(self.args.out_file)


@attr.s(auto_attribs=True)
class Config:
    """Configuration for the converter."""

    #: Input directory.
    indir: str
    #: Output file.
    out_file: str
    #: Path to about.md file
    about_md: typing.Optional[str] = None
    #: Top N markers to save.
    nmarkers: int = 10
    #: Whether to split by species.
    split_species: bool = False
    #: Whether to split by sample
    split_samples: bool = False
    #: The format to use for conversion.
    format: str = "auto"
    #: Use raw signal.
    use_raw: bool = False
    #: path to markers.tsv file
    markers: typing.Optional[str] = None


def run(args, _parser=None):
    """Main entry point after argument parsing."""
    format_ = args.format
    if format_ == "auto":
        if os.path.exists(os.path.join(args.indir, "coords.tsv")):
            format_ = "text"
        elif os.path.exists(os.path.join(args.indir, "analysis")):
            format_ = "cell-ranger"
        elif os.path.exists(os.path.join(args.indir, "data.loom")):
            format_ = "loom"
        else:
            raise ScelVisException("Could not auto detect the input format")

    converters = {"text": TextConverter, "cell-ranger": CellRangerConverter, "loom": LoomConverter}
    return converters[format_](args).run()


def setup_argparse(parser):
    """Setup argparse sub parser."""
    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        dest="indir",
        help="path to input/pipeline output directory",
    )
    parser.add_argument(
        "-a", "--about-md", help="Path to about.md file to embed in the resulting .h5ad file"
    )
    parser.add_argument(
        "-m", "--markers", help="Path to markers.tsv file to embed in the resulting .h5ad file"
    )
    parser.add_argument(
        "-o", "--output", required=True, dest="out_file", help="Path to the .h5ad file to write to"
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=("auto", "text", "cell-ranger", "loom"),
        default="auto",
        dest="format",
        help="input format",
    )
    parser.add_argument(
        "--use-raw",
        dest="use_raw",
        default=False,
        action="store_true",
        help="Do not normalize expression values (use raw counts)",
    )
    parser.add_argument(
        "--split-species",
        dest="split_species",
        default=False,
        action="store_true",
        help="Split species",
    )
    parser.add_argument(
        "--split-samples",
        dest="split_samples",
        default=False,
        action="store_true",
        help="""split samples according to summary.json file produced by cellranger aggr""",
    )
    parser.add_argument(
        "--nmarkers",
        dest="nmarkers",
        default=10,
        type=int,
        help="Save top n markers per cluster [10]",
    )
    parser.add_argument(
        "--ncells", dest="ncells", type=int, help="sample ncells cells from object [take all]"
    )
    parser.add_argument(
        "--verbose", default=False, action="store_true", help="Enable verbose output"
    )
