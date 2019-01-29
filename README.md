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
- [Spectral broadening](https://github.com/jason-neal/eniric/wiki/Broadening)
   - Rotational
   - Instrumental
- [Atmospheric transmission masking](https://github.com/jason-neal/eniric/wiki/Atmospheric-Transmission)
- Relative RV precision
  - The RV precision can be calculated relative to a specified SNR per pixel in the center of a spectroscopic band.
    The default as used in the Figueira et al. 2016 is a SNR of 100 at the center of the J-band.
- Spectral re-sampling
   - n pixels per FWHM
- Band selection
  - Analysis in individual spectroscopic bands.
- Incremental quality & precision
- Synthetic libraries available
    - Available through [Starfish]'s() grid_tools.
       - [PHOENIX-ACES](http://phoenix.astro.physik.uni-goettingen.de)
       - [BT-Settl](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/)


# [Installation](https://github.com/jason-neal/eniric/wiki/Installation)

You can install ``Eniric`` by cloning the current ``master`` branch.

It is recommended to use a ``conda`` or ``virtualenv`` environment.
To use the most up-to-date packages install the pinned requirements from the ``requirements_dev.txt`` file

    git clone https://github.com/jason-neal/eniric
    cd eniric
    pip install -r requirements_dev.txt
    python setup.py develop


Installation from the github repository should also be possible with `pip`.

``` 
  pip install https://github.com/jason-neal/eniric/archive/develop.zip#egg=eniric
```


If there are issues with installing Starfish see [https://github.com/iancze/Starfish](https://github.com/iancze/Starfish)

If you are not going to use `eniric` to analyze PHOENIX-ACES or BT-Settl synthetic spectral models then you may
get away with not installing it (some tests will just xfail).

##### Requirements for `Eniric` :
These should be automatically installed (if not present) when installing `eniric`.
- astropy
- joblib>=0.12.3
- matplotlib
- multiprocess
- numpy>=0.15.4
- pandas
- pyyaml
- scipy
- Starfish
- tqdm

## Preparation
#### Data download
To download the data for eniric, an atmopsheric transmission spectra and some test Spectra run the
following from the main repo directory

Linux:
`download_eniric_data.sh`

Windows:
`... `

This should place the data in `data/atmmodel` and `data/testdata` where it can be found for testing.

### Configuration
`Eniric` uses a `config.yaml` file which is required in directory where you are running `Eniric`. (i.e. the current directory)
to specify some paths, such as the location the the synthetic spectral library.

```
paths:
   phoenix-raw: path/to/phoenix/aces/spectra
   btsettl-raw: ["path", "to", "btsettl" ,"spectra"]
   ...
```
The paths can either be a string or a list of strings to pass to `os.path.join` (os independant).

You can use the `config.yaml` to specify custom wavelength ranges to use
```
bands:
  all: [..., myband]  # add myband to all list

custom_bands:
    myband: [1.5, 1.6] # micron
```

You can then pass `myband` to the `band` arguments in `Eniric` scripts/functions.

This based off `Starfish` and although many keywords are needed to be present
for `Starfish` to run they are not used for `Eniric`'s usage of `Starfish` and are fine left blank.


#### Atmospheric data:
To perform telluric masking and account for the transmission of Earth's atmosphere a telluric spectra is required.
`Eniric` includes the telluric spectra uses in Figueira et al. 2016, averaged over 2014.
To automatically prepare the telluric masks, splitting into bands and applying the barycentric expansion run the following scripts:
- `split_atmmodel.py`
- `barycenter_broaden_atmmodel.py`

These will split the large telluirc spectra into the bands specified in the `config.yaml` so that the
 opening and slicing of the large telluric spectrum is not performed continually.

To change the telluric line cutoff depth you to 4% can pass (default = 2%) you can pass it like this

    `split_atmmodel.py --cutoff-depth 4`

You can specify your own telluric mask instead.
By keeping it in the same format and setting atmmodel parameters in `config.yaml` you can make use of the
these scripts and the`Atmosphere` class which can perform the mask cutoff and doppler shifting.

Or you can manually apply your own masking function as the mask parameter to the `rv_precision` function.

### Convolutions
The most computational component in `Eniric` is the convolutions. To help with this we use parallel prcessing and caching.

- *Caching*:
The convolution results are cached using Joblib to avoid repeating the convoutions. This can be disabled by
setting the `location: None` in the `config.yaml`.

- *Parallel Processing*:
The default number of processors used is one less then the total number of cores (N-1).
You can change this by specifying the `num_procs`.
Setting `num_procs = 0` or `1` disables parallel processing.

## Usage
You can now calculate the theoretical RV precision for any PHOENIX-ACES model.
You will need to configure the path to the phoenix models in ´config.yaml´

e.g.

    phoenix_precision.py -t 3900 -l 4.5, -m 0.5 -r 100000 -v 1.0 -b J K

Will calculate the RV precision in the `J` and `K`-band of the PHOENIX-ACES spectra with parameters \[Teff=3900K, logg=4.5, \[Fe/H\]=0.5\] observed at a resolution of 100,000 and rotating with 1.0 km/s.
For more details on the command line arguments to use see the [wiki](https://github.com/jason-neal/eniric/wiki) or type

    phoenix_precision.py -h


# The Readme below this point needs amended....

## Outline

The code works in two main stages, "spectral preparation" and "precision calculation".

#### Spectrum preparation

`eniric/nIRanalysis.py`

This stage takes in the raw PHOENIX-ACES spectral models and transforms them, saving the results of this computation as .dat files.

It includes:
- Conversion from flux to photon counts.
- Resolution convolution
- Re-sampling

Some scripts are given in `scripts` to run this preparation over all desired parameters automatically. You will have to modify the paths to things.


#### Precision Calculations

`python scripts/nIR_precision.py`

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

scripts/prepare_spectra.py - This opens the phoenix flux spectra, add wavelength axis in microns and converts flux to photon counts. It saves this in the phoenix_dat dir. (The copy of wavelengths does waste space.)

scripts/nIR_run.py - Perform the resolution and rotational convolution on the prepared spectra.

This also does the re-sampling.

e.g. python ../Codes/eniric/scripts/nIR_run.py -s M0 M3 M6 M9 -b Y J H K -v 1.0 5.0 10.0 -R 60000 80000 100000 --sample_rate 3


## Background
The origin of this code was used in [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900).

    P. Figueira, V. Zh. Adibekyan, M. Oshagh, J. J. Neal, B. Rojas-Ayala, C. Lovis, C. Melo, F. Pepe, N. C. Santos, M. Tsantaki, 2016,
    Radial velocity information content of M dwarf spectra in the near-infrared,
    Astronomy and Astrophysics, 586, A101

It had a number of efficiency issues with convolution which were improved upon

To reproduce the updated results for [Figueira et al. 2016](http://dx.doi.org/10.1051/0004-6361/201526900) run

    phoenix_precision.py -t 3900 3500 2800 2600 -l 4.5, -m 0.5 -r 60000 80000 100000 -v 1.0 5.0 10.0 -b Z Y J H K
