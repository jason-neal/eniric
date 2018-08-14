# ENIRIC - Extended Near InfraRed Information Content

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/eniric&utm_campaign=badger)[![Codacy Badge](https://api.codacy.com/project/badge/Coverage/24d3d525a79d4ae493de8c527540edef)](https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jason-neal/eniric&amp;utm_campaign=Badge_Coverage)
[![Build Status](https://travis-ci.org/jason-neal/eniric.svg?branch=master)](https://travis-ci.org/jason-neal/eniric)[![Coverage Status](https://coveralls.io/repos/github/jason-neal/eniric/badge.svg?branch=master)](https://coveralls.io/github/jason-neal/eniric?branch=master)[![Updates](https://pyup.io/repos/github/jason-neal/eniric/shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)[![Code Climate](https://codeclimate.com/github/jason-neal/eniric/badges/gpa.svg)](https://codeclimate.com/github/jason-neal/eniric)[![Issue Count](https://codeclimate.com/github/jason-neal/eniric/badges/issue_count.svg)](https://codeclimate.com/github/jason-neal/eniric)[![Test Coverage](https://codeclimate.com/github/jason-neal/eniric/badges/coverage.svg)](https://codeclimate.com/github/jason-neal/eniric/coverage)[![Python 3](https://pyup.io/repos/github/jason-neal/eniric/python-3-shield.svg)](https://pyup.io/repos/github/jason-neal/eniric/)

Eniric is a Python 3 software written to access the Radial Velocity precision of Near-InfraRed (NIR) Spectra. 
Eniric is built upon the code initially used in [Figueria et al. 2016]() to analysis the precision of M-dwarf stars, extending the ability to use any synethetic spectra from the PHOENIX-ACES library and making it easier to use.

## Features:
- Spectral broadening
   Allows for Rotational and Instrumental broadening of synthetic spectra given rotation speed `vsini` and resolution `R`.
- Atmospheric transmission
 Analysis RV precision attainable under 3 conditions presented in [Figueria et al. 2016]().
    - No treatment of atmospheric transmission
    - Masking all regions affected by atmospheric absorption of a given % over the course of the year.
    - Assuming perfect telluric correction in which the variance of the measured flux is impacted.
- Relative RV precision 
   The RV precision are present relative to a specified SNR per pixel in the center of a photometric band. The default as usesd in the Figueira et al. 2016 is a SNR of 100 at the center of the J-band.
- Re-sampling
   Allows for re-sampling of spectra to n pixels per FWHM. Default=3.
- Band selection
  Analysis the RV precision attainable by the individual photometric bands `Z`, `Y`, `J`, `H`, `K` in the NIR.   
- Incremental quality/precision
    As well as the photometric band precision you can determine the incremental spectral quality or RV precision on narrow wavelenght slices across the entire spectrum, similar to that present in Figure 1 of [Artigau et al. 2018](http://adsabs.harvard.edu/abs/2018AJ....155..198A).

## Usage
You can now calculate the theoretical RV precision for any PHOENIX-ACES model.
You will need to configure the path to the phoenix models in ´config.yaml´

e.g.

    aces_precision.py -t 3900 -l 4.5, -m 0.5 -r 100000 -v 1.0 -b J K
    
Will calculate the RV precision in the J- and K-band of the PHOENIX-ACES spectra with parameters \[Teff=3900K, logg=4.5, \[Fe/H\]=0.5\] observed at a resolution of 100,000 and rotating with 1.0 km/s.  
For more details on the command line arguments to use see the wiki or type

    aces_precision.py -h 


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
The origin of this code was used in [Figueira et al. 2016](https://arxiv.org/abs/1511.07468).

    P. Figueira, V. Zh. Adibekyan, M. Oshagh, J. J. Neal, B. Rojas-Ayala, C. Lovis, C. Melo, F. Pepe, N. C. Santos, M. Tsantaki, 2016,
    Radial velocity information content of M dwarf spectra in the near-infrared,
    Astronomy and Astrophysics, 586, A101

It had a number of efficiency issues with convolution which were improved upo

## Bugs:
A number of bugs were found when adapting this code. Mainly affecting the condition involving telluric line masking. 
This alters the RV in this column, sometimes significantly. This, however, does **NOT** alter the conclusions in the published paper. A corrected table will be included in an upcomming publication Neal et. al. 2018 (in prep.).

n. The largest were
    1) Use numpy mapping slicing instead of comprehension lists.  (~>250 times faster)
    2) Use `Joblib` to parallelize the convolutions.
    3) Caching results to avoid repeating the same convolutions.

In addressing 1) the convolution speed for a particular test case went from `1hour 22 minutes` down to `27 seconds`.
