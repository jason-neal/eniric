# ENIRIC - Extended Near InfraRed Information Content

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/eniric&utm_campaign=badger)
[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jason-neal/eniric&amp;utm_campaign=Badge_Coverage)
[![Build Status](https://travis-ci.org/jason-neal/eniric.svg?branch=master)](https://travis-ci.org/jason-neal/eniric)
[![Coverage Status](https://coveralls.io/repos/github/jason-neal/eniric/badge.svg?branch=master)](https://coveralls.io/github/jason-neal/eniric?branch=master)
[![Updates](https://pyup.io/repos/github/jason-neal/eniric/shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)
[![Python 3](https://pyup.io/repos/github/jason-neal/eniric/python-3-shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)
[![PyPI version](https://badge.fury.io/py/eniric.svg)](https://badge.fury.io/py/eniric)

Eniric is a Python 3 software written to access the Radial Velocity precision of Near-InfraRed (NIR) Spectra.
Eniric is built upon the code initially used in [Figueria et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900) to analysis the precision of M-dwarf stars, extending the ability to use any synethetic spectra from the PHOENIX-ACES and BT-Settl libraries and making it easier to use.

Checkout the [wiki here](https://github.com/jason-neal/eniric/wiki)!

## Features:
- [Spectral broadening](https://github.com/jason-neal/eniric/wiki/Broadening)
   - Rotational
   - Instrumental
- [Atmospheric transmission masking](https://github.com/jason-neal/eniric/wiki/Atmospheric-Transmission)
- Relative RV precision
   The RV precision can be calculated relative to a specified SNR per pixel in the center of a spectroscopic band.
    The default as used in the Figueira et al. 2016 is a SNR of 100 at the center of the J-band.
- Spectral re-sampling
   - n pixels per FWHM
- Band selection
  - Analysis in individual spectroscopic bands.  
- Incremental quality & precision
- Synthetic libraries
    - PHOENIX-ACES
    - BT-Settl


# Installation
Installation instructions can be found [here](https://github.com/jason-neal/eniric/wiki/Installation).


## Usage
You can now calculate the theoretical RV precision for any PHOENIX-ACES model.
You will need to configure the path to the phoenix models in ´config.yaml´

e.g.

    phoenix_precision.py -t 3900 -l 4.5, -m 0.5 -r 100000 -v 1.0 -b J K

Will calculate the RV precision in the J- and K-band of the PHOENIX-ACES spectra with parameters \[Teff=3900K, logg=4.5, \[Fe/H\]=0.5\] observed at a resolution of 100,000 and rotating with 1.0 km/s.  
For more details on the command line arguments to use see the wiki or type

    phoenix_precision.py -h


## Outline

The code works in two main stages, "spectral preparation" and "precision calculation".

#### Spectrum preparation

`eniric/nIRanalysis.py`

This stage takes in the raw PHOENIX-ACES spectral models and transforms them, saving the results of this computation as .txt files.

It includes:
- Conversion from flux to photon counts.
- Resolution convolution
- Re-sampling

Some scripts are given in `eniric_scripts` to run this preparation over all desired parameters automatically. You will have to modify the paths to things.


#### Precision Calculations

`python eniric_scripts/nIR_precision.py`

This takes in the processed spectra and performs the precision calculations for all 3 conditions outlined in the original paper.
- Cond1. Total information
- Cond2. +/-30km/s telluric line > 2% masking
- Cond3. Perfect telluric correction with variance correction

It also *scales the flux level* to a desired SNR level in a desired band, see below, as this affects the RV precision calculated. By default this is a SNR of 100 in the J band.


## Band SNR Scaling.
By default, in accordance with the initial paper, each spectra band is normalized to 100 SNR in the center of the J band.

This now does this automatically by measuring the SNR in 1 pixel resolution (3 points) in the center of the band. And scales accordingly. This adds a spectral model dependent factor on the RV precision.
To get around you can manually specify the SNR level to normalize to and which specific band to normalize to. (it can be itself for instance).



## Instructions

Create an empty dir to hold your analysis.
Create data dir with re-sampled, results, phoenix_dat
Copy config.yaml and adjust the paths relative to what you created and to the raw phoenix spectra.

eniric_scripts/prepare_spectra.py - This opens the phoenix flux spectra, add wavelength axis in microns and converts flux to photon counts. It saves this in the phoenix_dat dir. (The copy of wavelengths does waste space.)

eniric_scripts/nIR_run.py  - Perform the resolution and rotational convolution on the prepared spectra.

This also does the re-sampling.

e.g. python ../Codes/eniric/eniric_scripts/nIR_run.py -s M0 M3 M6 M9 -b Y J H K -v 1.0 5.0 10.0 -R 60000 80000 100000 --sample_rate 3


## Background
The origin of this code was used in [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900).

    P. Figueira, V. Zh. Adibekyan, M. Oshagh, J. J. Neal, B. Rojas-Ayala, C. Lovis, C. Melo, F. Pepe, N. C. Santos, M. Tsantaki, 2016,
    Radial velocity information content of M dwarf spectra in the near-infrared,
    Astronomy and Astrophysics, 586, A101

It had a number of efficiency issues with convolution which were improved upon

To reproduce the updated results for [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900) run

    phoenix_precision.py -t 3900 3500 2800 2600 -l 4.5, -m 0.5 -r 60000 80000 100000 -v 1.0 5.0 10.0 -b Z Y J H K
