"""Test for the file conversion."""

import os.path
import textwrap

import pytest

from scelvis import convert


@pytest.fixture(scope="function")
def config_no_about_md(tmp_path):
    outdir = tmp_path / "output"
    outdir.mkdir()
    return convert.Config(
        indir=os.path.join(os.path.dirname(__file__), "..", "examples/hgmm_1k.raw"),
        out_file=str(outdir / "out.h5ad"),
        about_md=None,
        split_species=True,
        split_samples=False,
        use_raw=False,
    )


@pytest.fixture(scope="function")
def config_with_about_md(tmp_path):
    outdir = tmp_path / "output"
    outdir.mkdir()
    with open(str(tmp_path / "about.md"), "wt") as aboutf:
        about = "\n".join(
            (
                "----",
                "title: This is some title",
                "short_title: Short Title",
                "----",
                "This is the readme contents.",
                "",
                "With several lines.",
            )
        )
        print(textwrap.dedent(about).lstrip(), file=aboutf)
    return convert.Config(
        indir=os.path.join(os.path.dirname(__file__), "..", "examples/hgmm_1k.raw"),
        out_file=str(outdir / "out.h5ad"),
        about_md=str(tmp_path / "about.md"),
        split_species=True,
        split_samples=False,
        use_raw=False,
    )


def test_run_converter_directly_no_about_md(config_no_about_md):
    config = config_no_about_md
    convert.CellRangerConverter(config).run()
    assert os.path.exists(config.out_file)


def test_run_main_no_about_md(config_no_about_md):
    config = config_no_about_md
    assert convert.run(config) is None
    assert os.path.exists(config.out_file)


def test_run_converter_directly_with_about_md(config_with_about_md):
    config = config_with_about_md
    convert.CellRangerConverter(config).run()
    assert os.path.exists(config.out_file)


def test_run_main_with_about_md(config_with_about_md):
    config = config_with_about_md
    assert convert.run(config) is None
    assert os.path.exists(config.out_file)
