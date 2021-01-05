#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This file is part of the
#   SCORPIO Project (https://github.com/josegit88/SCORPIO).
# Copyright (c) 2020, Jose Benavides
# License: MIT
#   Full Text: https://github.com/josegit88/SCORPIO/blob/master/LICENSE


# =============================================================================
# DOCS
# =============================================================================

"""This file is for distribute and install SCORPIO
"""


# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib

from ez_setup import use_setuptools

use_setuptools()
from setuptools import setup  # noqa


# =============================================================================
# CONSTANTS
# =============================================================================


REQUIREMENTS = [
    "numpy",
    "astropy",
    "astroquery",
    "matplotlib",
    "attrs",
    "seaborn",
    "retrying",
    "scipy",
]

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

with open(PATH / "README.md") as fp:
    LONG_DESCRIPTION = fp.read()

with open(PATH / "scorpio.py") as fp:
    for idx in fp.readlines():
        if idx.startswith("__version__ = "):
            VERSION = idx.split("=", 1)[-1].replace('"', "").strip()
            break


DESCRIPTION = """SCORPIO is a tool for analyze galaxy pairs."""


# =============================================================================
# FUNCTIONS
# =============================================================================


def do_setup():
    setup(
        name="scorpio_gp",
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type="text/markdown",
        author=["Jose Benavides"],
        author_email="jose.astroph@gmail.com",
        url="https://github.com/josegit88/SCORPIO",
        license="MIT",
        keywords=["scorpio", "survey", "images", "stacking"],
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering",
        ],
        py_modules=["scorpio", "ez_setup"],
        install_requires=REQUIREMENTS,
    )


if __name__ == "__main__":
    do_setup()
