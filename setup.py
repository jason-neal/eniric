import codecs
import os
import sys
from os.path import join
from typing import List

from setuptools import find_packages, setup

if sys.version < "3.6":
    sys.exit(
        "Error: Python 3.6 or greater required for Eniric (using {})".format(
            sys.version
        )
    )


here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with codecs.open(join(here, "README.md")) as f:
    long_description = f.read()

with codecs.open(join(here, "requirements.txt")) as f:
    requirements = f.read().splitlines()

# Conditional pytest-runner install
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner: List[str] = ["pytest-runner>=4"] if needs_pytest else []
setup_requires: List[str] = []

config = {
    "name": "eniric",
    "description": "Eniric: Extended NIR Information Content",
    "long_description": long_description,
    "long_description_content_type": "text/markdown",
    "author": "Jason Neal",
    "url": "https://github.com/jason-neal/eniric.git",
    "download_url": "https://github.com/jason-neal/eniric.git",
    "author_email": "jason.neal@astro.up.pt",
    "version": "1.0rc3",
    "license": "MIT Licence",
    "setup_requires": setup_requires + pytest_runner,
    "tests_require": ["pytest", "hypothesis"],
    "install_requires": requirements,
    "extras_require": {
        "dev": ["check-manifest"],
        "test": ["coverage", "pytest", "pytest-cov", "python-coveralls", "hypothesis"],
        "docs": ["sphinx >= 1.4", "sphinx_rtd_theme", "rstcheck"],
        "ci": [
            "codacy-coverage==1.3.11",
            "codeclimate-test-reporter>=0.2.3",
            "python-coveralls>=2.9.1",
        ],
    },  # $ pip install -e .[dev,test, docs]
    "packages": find_packages(),
    "scripts": [
        "scripts/phoenix_precision.py",
        "scripts/split_atmmodel.py",
        "scripts/barycenter_broaden_atmmodel.py",
        "scripts/untar_here.py",
        "scripts/download/download_eniric_data.sh",
        "scripts/download/ps_download_eniric_data.ps1",
        "scripts/download/download_test_aces.py",
    ],
    "include_package_data": True,
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    "classifiers": [
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
    ],
    # What does your project relate to?
    "keywords": [
        "Astronomy",
        "Radial velocity",
        "Near-infrared spectroscopy",
        "M-dwarfs",
    ],
}

setup(**config)
