# [ENIRIC](https://github.com/jason-neal/eniric) - Extended Near InfraRed Information Content

[![Documentation Status](https://readthedocs.org/projects/eniric/badge/?version=latest)](https://eniric.readthedocs.io/en/latest/?badge=latest)
[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/eniric&utm_campaign=badger)
[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jason-neal/eniric&amp;utm_campaign=Badge_Coverage)
[![Build Status](https://travis-ci.org/jason-neal/eniric.svg?branch=master)](https://travis-ci.org/jason-neal/eniric)
[![Coverage Status](https://coveralls.io/repos/github/jason-neal/eniric/badge.svg?branch=master)](https://coveralls.io/github/jason-neal/eniric?branch=master)
[![Updates](https://pyup.io/repos/github/jason-neal/eniric/shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)
[![Python 3](https://pyup.io/repos/github/jason-neal/eniric/python-3-shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)
[![PyPI version](https://badge.fury.io/py/eniric.svg)](https://badge.fury.io/py/eniric)

`Eniric` is a Python 3 software to compute the theoretical Radial Velocity (RV) precision of stellar spectra.
`Eniric` is an overhaul and extension to the code used in [Figueria et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900) to analysis the precision of M-dwarf stars.
Extending the performance and usability, it is able to be used on any synthetic spectra from the [PHOENIX-ACES](http://phoenix.astro.physik.uni-goettingen.de) and [BT-Settl](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/) (CIFIST2001-2015) libraries.

Checkout the documentation on [Read the Docs](https://eniric.readthedocs.io/en/latest/)!

## Features:
`Eniric` contains a number of features to transform and prepare the spectra (observed and synthetic).

- [Spectral broadening](https://eniric.readthedocs.io/en/latest/broadening.html)

     Allows for Rotational and Instrumental broadening of synthetic spectra given a rotation speed ``vsini`` and resolution ``R``.

- [Atmospheric transmission masking](https://eniric.readthedocs.io/en/latest/telluric_corection.html)

   Analyzing the RV precision attainable under the different masking conditions presented in `Figueira et al. 2016`_.

   The three conditions specifically treated are:
   * No contamination or treatment of atmospheric transmission
   * Masking all regions affected by atmospheric absorption of a given depth % over the course of the year.
   * Assuming perfect telluric correction in which the variance of the measured flux is impacted.

- Relative RV precision

  The RV precision can be calculated relative to a specified SNR per pixel in the center of a spectroscopic band.
    The default as used in the Figueira et al. 2016 is a SNR of 100 at the center of the J-band.

- Spectral Resampling

   Allows for resampling of synthetic spectra to ``N`` pixels per FWHM.

- SNR normalization.

   Normalize spectral flux to a defined SNR level.

- Band selection

  Analysis splitable into individual photometric bands ``Z``\ , ``Y``\ , ``J``\ , ``H``\ , ``K``.
  User definable.

- Theoretical RV precision

   Compute spectral RV precision and spectral quality.

- Incremental quality & precision

    Determine the RV precision and spectral quality on narrow wavelength slices across the entire spectrum, similar to that present in Figure 1 of `Artigau et al. 2018 <http://adsabs.harvard.edu/abs/2018AJ....155..198A>`_.

* Analyse relative precision of synthetic libraries

    The RV precision of are present relative to a specified SNR per pixel in the center of a photometric band.
    The default as used in `Figueira et al. 2016`_ is a SNR of 100 at the center of the J-band.
    - Available through [Starfish]'s() grid_tools.
       - [PHOENIX-ACES](http://phoenix.astro.physik.uni-goettingen.de)
       - [BT-Settl](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/)


## Contents

- [Installation](https://eniric.readthedocs.io/en/latest/installation.html)
- [Configuration](https://eniric.readthedocs.io/en/latest/configuration.html)
- [Basic Usage](https://eniric.readthedocs.io/en/latest/basic_usage.html)
- [Broadening](https://eniric.readthedocs.io/en/latest/broadening.html)
- [Atmospheric Transmission](https://eniric.readthedocs.io/en/latest/telluric_corection.html)
- [Normalization](https://eniric.readthedocs.io/en/latest/normalization.html)
- [Resampling](https://eniric.readthedocs.io/en/latest/resampling.html)
- [Theoretical Precision of Synthetic Spectra](https://eniric.readthedocs.io/en/latest/theoretical_precision.html)
- [Scripts](https://eniric.readthedocs.io/en/latest/scripts.html)
- [Example Notebooks](https://eniric.readthedocs.io/en/latest/examples.html)
- [Utilities](https://eniric.readthedocs.io/en/latest/utilities.html)


## Background
The origin of this code was used in [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900).

    P. Figueira, V. Zh. Adibekyan, M. Oshagh, J. J. Neal, B. Rojas-Ayala, C. Lovis, C. Melo, F. Pepe, N. C. Santos, M. Tsantaki, 2016,
    Radial velocity information content of M dwarf spectra in the near-infrared,
    Astronomy and Astrophysics, 586, A101

It had a number of efficiency issues with convolution which were improved upon

To reproduce the updated results for [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900) run

    phoenix_precision.py -t 3900 3500 2800 2600 -l 4.5 -m 0.5 -r 60000 80000 100000 -v 1.0 5.0 10.0 -b Z Y J H K

after installation and configuration.
