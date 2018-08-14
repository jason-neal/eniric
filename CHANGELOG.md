# ENIRIC CHANGELOG

To get a list of commit messages since last version to help write change log try
    `git log YOUR_LAST_VERSION_TAG..HEAD --no-merges --format=%B`

### Upcomming
-------
- Blacken source code
- Add pre-commit hooks
- Use configurable Bands in config.yaml
- Refactor atmosphere handling code into `Atmosphere` class
- Updated plotting script, cumulatve rv
- Updated `split_atmmodel.py` and `bary_shift_atmmodel.py` scripts for atmmodel preparation.
- Precision of `BT-Settl` models in aces_precision.py
- Test on Python 3.7
- Add Contributing.md and Contributors.md.

### v0.5
--------
Modified for the NIRPS ETC calculations
- New precision script
    - eniric_scripts/any_spectral_quality.py
    - Precision for any PHOENIX-ACES spectra
    - Adds Starfish dependency with GridTools
- Calculate spectral quality Q.
- Drop python 2.7 compatibility and remove some support code.
-
-

### v0.4
----
- Upgrades for SPIRou ETC calculations
- Relative SNR in any band
-
-

### v0.3
----
...
