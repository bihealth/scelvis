#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path

from setuptools import setup, find_packages

import versioneer


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

test_requirements = parse_requirements("requirements/test.txt")
# Use special requirements for RTD ("dash" instead of "dash[testing]").
if os.environ.get("READTHEDOCS", False):
    install_requirements = parse_requirements("requirements/rtd.txt")
else:
    install_requirements = parse_requirements("requirements/base.txt")

setup(
    author="Benedikt Obermayer, Manuel Holtgrewe",
    author_email="benedikt.obermayer@bihealth.de, manuel.holtgrewe@bihealth.de",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    entry_points={"console_scripts": ("scelvis = scelvis.cli:main",)},
    description="Single-cell RNA-seq data visualization",
    install_requires=install_requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="single-cell, visualization, bioinformatics",
    name="scelvis",
    packages=find_packages(include=["scelvis"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/bihealth/scelvis",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)
