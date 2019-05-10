# -*- coding: utf-8 -*-
"""The Dash visualization app web server for visualizing the data."""

import os
import tempfile

from logzero import logger

from . import settings


def run_server(args):
    """Actually run the Dash server."""
    from .app import app  # noqa

    app.run_server(host=args.host, port=args.port, debug=args.debug)


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
        settings.CACHE_DEFAULT_TIMEOUT = args.cache_default_timeout
        if args.cache_redis_url:
            settings.CACHE_TYPE = "redis"
            settings.CACHE_REDIS_URL = args.cache_redis_url
        elif args.cache_dir:
            settings.CACHE_DIR = args.cache_dir
        # Launch the Dash server.
        logger.info("Starting Dash web server on %s:%d", args.host, args.port)
        if settings.CACHE_TYPE == "filesystem" and not settings.CACHE_DIR:
            with tempfile.TemporaryDirectory(prefix="scelvis.") as tmpdir:
                logger.info("Using cache directory %s", tmpdir)
                settings.CACHE_DIR = tmpdir
                run_server(args)
        else:
            run_server(args)
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
    parser.add_argument(
        "--cache-dir",
        default=os.environ.get("SVIZ_CACHE_DIR"),
        help="Path to cache directory, default is to autocreate one.",
    )
    parser.add_argument(
        "--cache-redis-url",
        default=os.environ.get("SVIZ_CACHE_REDIS_URL"),
        help="Redis URL to use for caching, enables Redis cache",
    )
    parser.add_argument(
        "--cache-default-timeout",
        default=os.environ.get("SVIZ_CACHE_DEFAULT_TIMEOUT", 600),
        help="Default timeout for cache",
    )
