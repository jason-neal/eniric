import codecs
import os
import sys

if sys.version < "3.6":
    sys.exit(
        "Error: Python 3.6 or greater required for Eniric (using {})".format(
            sys.version
        )
    )

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with codecs.open(os.path.join(here, "README.md")) as f:
    long_description = f.read()


with codecs.open(os.path.join(here, "requirements.txt")) as f:
    requirements = f.read().splitlines()


config = {
    "name": "eniric",
    "description": "Eniric: Extended NIR Information Content",
    "long_description": long_description,
    "long_description_content_type": "text/markdown",
    "author": "Jason Neal",
    "url": "https://github.com/jason-neal/eniric.git",
    "download_url": "https://github.com/jason-neal/eniric.git",
    "author_email": "jason.neal@astro.up.pt",
    "version": "1.0rc1",
    "license": "MIT Licence",
    "setup_requires": ["pytest-runner"],
    "tests_require": ["pytest", "hypothesis"],
    "install_requires": requirements,
    "extras_require": {
        "dev": ["check-manifest"],
        "test": ["coverage", "pytest", "pytest-cov", "python-coveralls", "hypothesis"],
    },  # $ pip install -e .[dev,test]
    "packages": ["eniric", "eniric_scripts", "eniric.obsolete"],
    "scripts": [
        "eniric_scripts/phoenix_precision.py",
        "eniric_scripts/split_atmmodel.py",
        "eniric_scripts/bary_shift_atmmodel.py",
        "eniric_scripts/untar_here.py",
        "eniric_scripts/download_eniric_data.sh",
        "tests/download_test_PHOENIX_spec.sh",
        "eniric/obsolete/make_test_data.py",
        "eniric/obsolete/nIR_run.py",
        "eniric/obsolete/nIR_precision.py",
        "eniric/obsolete/prepare_data.py",
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
