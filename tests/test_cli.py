"""Minimal test for CLI"""

import pytest

from scelvis import cli


def test_run_help():
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        cli.main(["--help"])
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
