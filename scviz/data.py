"""Allows storing the representation of data sets."""

import os

import attr
from ruamel.yaml import YAML


@attr.s(auto_attribs=True)
class MetaData:
    """Class to bundle the data loaded for SCViz."""

    #: ID (= folder name) of the dataset
    id: str
    #: Title of the data set
    title: str
    #: Short title of the data set
    short_title: str
    #: String with Markdown-formatted README.
    readme: str


def load_metadata(path):
    """Load metadata from a dataset directory."""
    with open(os.path.join(path, "about.md"), "rt") as inputf:
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

        # Get title and short title, finally create MetaData object.
        meta = YAML().load("\n".join(header))
        title = meta.get("title", "Untitled")
        short_title = meta.get("title", title or "untitled")
        readme = "\n".join([line.rstrip() for line in lines])

        return MetaData(
            id=os.path.basename(path), title=title, short_title=short_title, readme=readme
        )
