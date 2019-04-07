# ENIRIC CHANGELOG

To get a list of commit messages since last version to help write change log try
    `git log YOUR_LAST_VERSION_TAG..HEAD --no-merges --format=%B`

### v1.0.0
----------
- Add JOSS paper
- User configurable wavelength bands in `config.yaml`.
- Refactor atmospheric transmission into `Atmosphere` class.
- Updated `split_atmmodel.py` and `bary_shift_atmmodel.py` scripts for atmmodel preparation.
- Add precision of `BT-Settl` library models in `phoenix_precision.py`.
- Functions for incremental quality and RV.
- Updated plotting script, cumulative rv.
- Doppler shift in `phoenix_precision.py`.
- Update Example Notebooks.
- Add `Contributing.md` and `Contributors.md`.
- Zip required testing data, Delete unneeded data.
- Refactor `eniric.Qcalculator` module to `eniric.precision`.
- Depreciate `--model phoenix` from `phoenix_precision.py`.

Other:
- `Blacken` source code.
- Add `pre-commit` hooks.
- Add `Python 3.7` testing.
- Drop `Python 3.5` support.
- Refactoring for 1.0 release.
- Relax pining in `requirements.txt`.
- Use `requirements.txt` in `setup.py`
- move download scripts to scripts/download.
- General documentation overhaul.


### v0.5.1
----------
- Adjust setup/readme and \_\_init\_\_ files.


### v0.5
--------
Modified for the NIRPS ETC calculations
- New precision script
    - eniric_scripts/any_spectral_quality.py
    - Precision for any PHOENIX-ACES spectra
    - Adds Starfish dependency with GridTools
- Calculate spectral quality Q.
- Drop python 2.7 compatibility and remove some support code.


### v0.4
--------
- Upgrades for SPIRou ETC calculations
- Relative SNR in any band

...

### v0.0.0
----------
- Aquire Figueira et al. 2016 original code base.
