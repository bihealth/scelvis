"""Conversion from single-cell pipeline output to SCelVis HDF5 files.

Currently, only conversion from the Cell-Ranger is supported.
"""

import contextlib
import logging
import os

import anndata
import attr
import scanpy.api as sc
import pandas as pd
from logzero import logger

from .exceptions import ScelVisException


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
        pca = self._load_pca()
        clustering = self._load_clustering()
        diffexp = self._load_diffexp()
        ad = self._load_expression(clustering, tsne, pca)

        if self.args.split_species:
            ad = self._split_species(ad)
        if not self.args.use_raw:
            ad = self._normalize_filter_dge(ad)

        self._write_output(ad, diffexp)
        logger.info("All done. Have a nice day!")

    def _load_tsne(self):
        tsne_file = os.path.join(self.args.indir, "tsne", "2_components", "projection.csv")
        if not os.path.isfile(tsne_file):
            raise ScelVisException("cannot find tSNE output at %s" % tsne_file)
        else:
            logging.info("Reading tSNE output from %s", tsne_file)
            return pd.read_csv(tsne_file, header=0, index_col=0)

    def _load_pca(self):
        pca_file = os.path.join(self.args.indir, "pca", "10_components", "projection.csv")
        if not os.path.isfile(pca_file):
            raise ScelVisException("cannot find PCA output at %s" % pca_file)
        else:
            logger.info("Reading PCA output from %s", pca_file)
            return pd.read_csv(pca_file, header=0, index_col=0)

    def _load_clustering(self):
        clustering_file = os.path.join(self.args.indir, "clustering", "graphclust", "clusters.csv")
        if not os.path.isfile(clustering_file):
            raise ScelVisException("cannot find clustering output at %s " % clustering_file)
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
            self.args.indir, "diffexp", "graphclust", "differential_expression.csv"
        )
        if not os.path.isfile(diffexp_file):
            raise ScelVisException("cannot find differential expression output at " + diffexp_file)
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

    def _load_expression(self, clustering, tsne, pca):
        expression_file = os.path.join(self.args.indir, "filtered_feature_bc_matrix.h5")
        if not os.path.isfile(expression_file):
            raise ScelVisException("cannot find expression file at %s" % expression_file)
        else:
            logger.info("Reading gene expression from %s", expression_file)
            with with_log_level(anndata.utils.logger, logging.WARN):
                ad = sc.read_10x_h5(expression_file)
            ad.var_names_make_unique()
            logger.info("Combining meta data")
            # TODO: do we need to make variable names unique here or can we suppress the warning
            # TODO: with ``with_log_level(anndata.utils.logger, logging.WARN)``?
            ad.obs["cluster"] = clustering
            ad.obs["n_counts"] = ad.X.sum(1).A1
            ad.obs["n_genes"] = (ad.X > 0).sum(1).A1
            logger.info("Adding coordinates")
            ad.obsm["X_tsne"] = tsne.values
            ad.obsm["X_pca"] = pca.values
            return ad

    def _split_species(self, ad):
        logger.info("Determining species mixing")
        species = ad.var_names.str.split("_", n=1).str[0]
        for sp in species.unique():
            ad.obs["nUMI_" + sp] = ad.X[:, species == sp].sum(1).A1
        cols = "nUMI_" + species.unique()
        ratio = ad.obs[cols].divide(ad.obs[cols].sum(axis=1), axis=0)
        ad.obs["species"] = (ratio > 0.95).idxmax(axis=1).str.split("_").str[1]
        ad.obs["species"][ratio.max(axis=1) < 0.95] = "other"
        return ad

    def _normalize_filter_dge(self, ad):
        logger.info("Normalizing and filtering DGE")
        # TODO: do we need to make variable names unique here or can we suppress the warning
        # TODO: with ``with_log_level(anndata.utils.logger, logging.WARN)``?
        sc.pp.filter_cells(ad, min_genes=100)
        sc.pp.filter_genes(ad, min_cells=5)
        sc.pp.normalize_per_cell(ad, counts_per_cell_after=1.0e4)
        sc.pp.log1p(ad)
        return ad

    def _write_output(self, ad, diffexp):
        out_file = os.path.join(self.args.outdir, "data.h5ad")
        marker_file = os.path.join(self.args.outdir, "markers.csv")
        logger.info("Saving anndata object to %s", out_file)
        # TODO: explicitely set column types to get rid of warning?
        ad.write(out_file)

        logger.info("Saving top %d markers per cluster to %s", self.args.nmarkers, marker_file)
        diffexp[(diffexp["p_adj"] < 0.05) & (diffexp["log2_fc"] > 0)].drop(
            "GeneID", axis=1
        ).sort_values(["Cluster", "p_adj"]).groupby("Cluster").head(self.args.nmarkers).to_csv(
            marker_file, header=True, index=True
        )


@attr.s(auto_attribs=True)
class Config:
    """Configuration for the converter."""

    #: Input directory.
    indir: str
    #: Output directory.
    outdir: str
    #: Top N markers to save.
    nmarkers: int = 10
    #: Whether to split by species.
    split_species: bool = False
    #: Use raw signal.
    use_raw: bool = False


def run(args, _parser=None):
    """Main entry point after argument parsing."""
    # TODO: detect pipeline output
    return CellRangerConverter(args).run()


def setup_argparse(parser):
    """Setup argparse sub parser."""
    parser.add_argument(
        "-i",
        "--input-dir",
        required=True,
        dest="indir",
        help="path to input/pipeline output directory",
    )
    parser.add_argument("-o", "--output-dir", required=True, dest="outdir", help="output directory")
    parser.add_argument(
        "--use_raw",
        dest="use_raw",
        default=False,
        action="store_true",
        help="Do not normalize DGE (use raw counts)",
    )
    parser.add_argument(
        "--split_species",
        dest="split_species",
        default=False,
        action="store_true",
        help="Split species",
    )
    parser.add_argument(
        "--nmarkers",
        dest="nmarkers",
        default=10,
        type=int,
        help="Save top n markers per cluster in markers.file [10]",
    )
    parser.add_argument(
        "--verbose", default=False, action="store_true", help="Enable verbose output"
    )
