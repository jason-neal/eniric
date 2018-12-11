.. Eniric documentation master file, created by
   sphinx-quickstart on Tue Dec 11 23:58:47 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Eniric's documentation!
==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation

`Eniric` is a Python 3 software to compute the theoretical Radial Velocity (RV) precision of stellar spectra.
`Eniric` is an overhaul and extension to the code used in `Figueria et al. 2016 <http://dx.doi.org/10.1051/0004-6361/201526900>`_ to analysis the precision of M-dwarf stars.
Extending the performance and usability, it is able to be used on any synthetic spectra from the `PHOENIX-ACES <http://phoenix.astro.physik.uni-goettingen.de>`_ and `BT-Settl <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/>`_ (CIFIST2001-2015) libraries.

Checkout the `wiki here <https://github.com/jason-neal/eniric/wiki>`_!

Features
--------
`Eniric` contains a number of features to transform and prepare the spectra (observed and synthetic).
- `Spectral broadening <https://github.com/jason-neal/eniric/wiki/Broadening>`_
   - Rotational
   - Instrumental
- `Atmospheric transmission masking <https://github.com/jason-neal/eniric/wiki/Atmospheric-Transmission>`_
- Relative RV precision
  - The RV precision can be calculated relative to a specified SNR per pixel in the center of a spectroscopic band.
    The default as used in the Figueira et al. 2016 is a SNR of 100 at the center of the J-band.
- Spectral re-sampling
   - n pixels per FWHM
- Band selection
  - Analysis in individual spectroscopic bands.
- Incremental quality & precision
- Synthetic libraries available
    - Available through Starfish's `grid_tools <https://iancze.github.io/Starfish/current/grid_tools.html>`_.
       - `PHOENIX-ACES <http://phoenix.astro.physik.uni-goettingen.de>`_
       - `BT-Settl <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/>`_


Support
-------

If you are having issues, please let us know.
Submit an issue on `Github <https://github.com/jason-neal/eniric/issues>`_.

License
-------

The project is licensed under the MIT License.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
