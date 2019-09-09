# -*- coding: utf-8 -*-
"""The Dash visualization app web server for visualizing the data."""

import os
import re
import tempfile
import urllib.parse

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


def parse_url(url_str):
    """Parse URL and fix scheme for ``file://`` URLs."""
    url = urllib.parse.urlparse(url_str)
    if not url.scheme:
        return url._replace(scheme="file", path=os.path.realpath(url.path))
    else:
        return url


def run(args, parser):
    """Main entry point after argument parsing."""
    data_sources = []
    if os.environ.get("SCELVIS_DATA_SOURCES"):
        data_sources += os.environ.get("SCELVIS_DATA_SOURCES").split(";")
    for data_source in args.data_sources:
        data_sources += data_source.split(";")

    data_sources = list(map(parse_url, data_sources))
    if not data_sources and not args.fake_data:
        parser.error(
            "You either have to specify --data-sources or set environment variable SCELVIS_DATA_SOURCES "
            "(or use --fake-data)"
        )
    args.data_sources = list(map(str, args.data_sources))

    # Configure the Dash app through ``settings`` module (see module docstring for more info).
    logger.info("Configuring settings from arguments %s", args)
    settings.PUBLIC_URL_PREFIX = re.sub(r"/+$", "", args.public_url_prefix or "")
    settings.FAKE_DATA = args.fake_data
    settings.DATA_SOURCES = data_sources
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
        "--port", type=int, help="Server port", default=int(os.environ.get("SCELVIS_PORT", 8050))
    )
    parser.add_argument(
        "--fake-data",
        default=False,
        action="store_true",
        help="Enable display of fake data set (for demo purposes).",
    )
    parser.add_argument(
        "--data-source",
        dest="data_sources",
        default=[],
        action="append",
        help="Path to data source(s)",
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
        default=os.environ.get("SCELVIS_UPLOAD_DISABLED", "0") not in ("", "0", "N", "n"),
        dest="upload_disabled",
        action="store_true",
        help="Whether or not to disable visualization uploads",
    )

    parser.add_argument(
        "--public-url-prefix",
        default=os.environ.get("SCELVIS_URL_PREFIX", ""),
        help="The prefix that this app will be served under (e.g., if behind a reverse proxy.)",
    )

    parser.add_argument(
        "--disable-conversion",
        default=os.environ.get("SCELVIS_CONVERSION_DISABLED", "0") not in ("", "0", "N", "n"),
        dest="conversion_disabled",
        action="store_true",
        help="Directory for visualization uploads, default is to create temporary directory",
    )
