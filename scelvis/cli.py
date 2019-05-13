"""Command line interface for scelvis.

This is the main program entry point for the ``scelvis`` executable and its sub commands.  The
actual implementation is in the modules ``convert`` and ``webserver``.
"""

import argparse
import logging
import warnings

import logzero

from scelvis import __version__
from .convert import run as run_convert
from .convert import setup_argparse as setup_argparse_convert
from .webserver import run as run_webserver
from .webserver import setup_argparse as setup_argparse_webserver


def run_nocmd(_, parser):
    """No command given, print help and ``exit(1)``."""
    parser.print_help()
    parser.exit(1)


def main(argv=None):
    """Main entry point before parsing command line arguments."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--verbose", action="store_true", default=False, help="Increase verbosity.")
    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)

    subparsers = parser.add_subparsers(dest="cmd")

    setup_argparse_convert(
        subparsers.add_parser("convert", help="Convert pipeline output to SCelVis input.")
    )
    setup_argparse_webserver(subparsers.add_parser("run", help="Run the SCelVis web server."))

    args = parser.parse_args(argv)

    # Setup logging verbosity.
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logzero.loglevel(level=level)

    # Handle the actual command line.
    cmds = {None: run_nocmd, "convert": run_convert, "run": run_webserver}

    # Disable duplicated crypto warnings from paramiko, triggered by fs.sshfs.
    warnings.filterwarnings(
        "once", module="paramiko.ecdsakey", message=".*unsafe construction of public numbers.*"
    )

    return cmds[args.cmd](args, parser)
