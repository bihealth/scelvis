"""Test for the file conversion."""

import os.path

import pytest

from scelvis import convert


@pytest.fixture(scope="function")
def config(tmp_path):
    outdir = tmp_path / "output"
    outdir.mkdir()
    return convert.Config(
        indir=os.path.join(os.path.dirname(__file__), "..", "examples/hgmm_1k.raw"),
        outdir=str(outdir),
        split_species=True,
        use_raw=False,
    )


def test_run_converter_directly(config):
    convert.CellRangerConverter(config).run()
    assert os.path.exists(os.path.join(config.outdir, "data.h5ad"))


def test_run_main(config):
    assert convert.run(config) is None
    assert os.path.exists(os.path.join(config.outdir, "data.h5ad"))
