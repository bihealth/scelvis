# -*- coding: utf-8 -*-
"""The Dash visualization app web server for visualizing the data."""

import os

from logzero import logger

from . import settings


def run(parser, args):
    """Main entry point after argument parsing."""
    if not args.data_dir:
        parser.error(
            "You either have to specify --data-dir or set environment variable SVIZ_DATADIR"
        )
    elif not os.path.exists(args.data_dir):
        parser.error("The data directory %s does not exit!" % args.data_dir)
    else:
        # Configure the Dash app through ``settings`` module (see module docstring for more info).
        logger.info("Configuring settings from arguments %s", args)
        settings.DATA_DIR = args.data_dir
        # Launch the Dash server.
        logger.info("Starting Dash web server on %s:%d", args.host, args.port)
        from .app import app  # noqa

        app.run_server(host=args.host, port=args.port, debug=args.debug)
        logger.info("Web server stopped. Have a nice day!")


def setup_argparse(parser):
    """Setup argparse sub parser."""
    parser.add_argument("--debug", default=False, action="store_true", help="Enable debug mode")
    parser.add_argument(
        "--host", help="Server host", default=os.environ.get("SVIZ_HOST", "0.0.0.0")
    )
    parser.add_argument(
        "--port", type=int, help="Server port", default=int(os.environ.get("SVIZ_HOST", 8050))
    )
    parser.add_argument(
        "--data-dir", default=os.environ.get("SVIZ_DATADIR"), help="Path to data directory"
    )
