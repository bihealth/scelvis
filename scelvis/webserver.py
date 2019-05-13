# -*- coding: utf-8 -*-
"""The Dash visualization app web server for visualizing the data."""

import os
import tempfile
import urllib

import fs
from logzero import logger

from . import settings


def run_server(args):
    """Actually run the Dash server."""
    from .app import app  # noqa

    app.run_server(host=args.host, port=args.port, debug=args.debug)


def run_temp_dir(args):
    """Setup temporary directory and run server."""
    with tempfile.TemporaryDirectory(prefix="scelvis.tmp") as tmpdir:
        logger.info("Creating temporary directory %s", tmpdir)
        settings.TEMP_DIR = tmpdir
        run_server(args)


def run_upload_dir(args):
    """Setup upload directory if necessary and continue with temp dir setup."""
    if not settings.UPLOAD_ENABLED:
        logger.info("Note that uploads are disabled")
        run_temp_dir(args)
    elif settings.UPLOAD_DIR:
        logger.info("Upload directory is %s", settings.UPLOAD_DIR)
        run_temp_dir(args)
    else:
        with tempfile.TemporaryDirectory(prefix="scelvis.upload") as tmpdir:
            logger.info("Creating upload directory %s", tmpdir)
            settings.UPLOAD_DIR = tmpdir
            run_temp_dir(args)


def run_cache_dir(args):
    """Setup cache directory if necessary and run server."""
    # Launch the Dash server.
    logger.info("Starting Dash web server on %s:%d", args.host, args.port)
    if settings.CACHE_TYPE == "filesystem" and not settings.CACHE_DIR:
        with tempfile.TemporaryDirectory(prefix="scelvis.cache.") as tmpdir:
            logger.info("Using cache directory %s", tmpdir)
            settings.CACHE_DIR = tmpdir
            run_upload_dir(args)
    else:
        run_upload_dir(args)
    logger.info("Web server stopped. Have a nice day!")


def run(args, parser):
    """Main entry point after argument parsing."""
    data_dirs = []
    if os.environ.get("SCELVIS_DATA_DIRS"):
        data_dirs += os.environ.get("SCELVIS_DATA_DIRS").split(";")
    if args.data_dirs:
        data_dirs += args.data_dirs
    if not data_dirs:
        parser.error(
            "You either have to specify --data-dir or set environment variable SCELVIS_DATA_DIRS"
        )
    for data_dir in args.data_dirs:
        url = urllib.parse.urlparse(data_dir)
        the_fs = fs.open_fs("%s://%s/" % (url.scheme or "file", url.netloc))
        if not the_fs.exists(url.path):
            parser.error("The data directory %s [%s] does not exist!" % (the_fs, url.path))
    # Configure the Dash app through ``settings`` module (see module docstring for more info).
    logger.info("Configuring settings from arguments %s", args)
    settings.DATA_DIRS = data_dirs
    settings.CACHE_DEFAULT_TIMEOUT = args.cache_default_timeout
    if args.cache_redis_url:
        settings.CACHE_TYPE = "redis"
        settings.CACHE_REDIS_URL = args.cache_redis_url
    elif args.cache_dir:
        settings.CACHE_DIR = args.cache_dir
    settings.UPLOAD_ENABLED = not args.upload_disabled
    settings.UPLOAD_DIR = args.upload_dir
    settings.CONVERSION_ENABLED = not args.conversion_disabled
    run_cache_dir(args)


def setup_argparse(parser):
    """Setup argparse sub parser."""
    parser.add_argument("--debug", default=False, action="store_true", help="Enable debug mode")
    parser.add_argument(
        "--host", help="Server host", default=os.environ.get("SCELVIS_HOST", "0.0.0.0")
    )
    parser.add_argument(
        "--port", type=int, help="Server port", default=int(os.environ.get("SCELVIS_HOST", 8050))
    )
    parser.add_argument(
        "--data-dir",
        dest="data_dirs",
        default=[],
        action="append",
        help="Path to data directory/ies",
    )
    parser.add_argument(
        "--cache-dir",
        default=os.environ.get("SCELVIS_CACHE_DIR"),
        help="Path to cache directory, default is to autocreate one.",
    )
    parser.add_argument(
        "--cache-redis-url",
        default=os.environ.get("SCELVIS_CACHE_REDIS_URL"),
        help="Redis URL to use for caching, enables Redis cache",
    )
    parser.add_argument(
        "--cache-default-timeout",
        default=os.environ.get("SCELVIS_CACHE_DEFAULT_TIMEOUT", 600),
        help="Default timeout for cache",
    )

    parser.add_argument(
        "--upload-dir",
        default=os.environ.get("SCELVIS_UPLOAD_DIR"),
        help="Directory for visualization uploads, default is to create temporary directory",
    )
    parser.add_argument(
        "--disable-upload",
        default=os.environ.get("SCELVIS_UPLOAD_DISABLED", False),
        dest="upload_disabled",
        action="store_true",
        help="Whether or not to disable visualization uploads",
    )

    parser.add_argument(
        "--disable-conversion",
        default=os.environ.get("SCELVIS_CONVERSION_DISABLED", False),
        dest="conversion_disabled",
        action="store_true",
        help="Directory for visualization uploads, default is to create temporary directory",
    )
