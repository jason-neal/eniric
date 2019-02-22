.. Eniric documentation master file, created by
   sphinx-quickstart on Tue Dec 11 23:58:47 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://img.shields.io/badge/python-3.6+-blue.svg
   :target: https://www.python.org/downloads/release/python-360/
   :alt: Python 3.6+


.. image:: https://travis-ci.org/jason-neal/eniric.svg?branch=master
   :target: https://travis-ci.org/jason-neal/eniric
   :alt: Build Status


.. image:: https://coveralls.io/repos/github/jason-neal/eniric/badge.svg?branch=master
   :target: https://coveralls.io/github/jason-neal/eniric?branch=master
   :alt: Coverage Status


.. image:: https://api.codacy.com/project/badge/Grade/24d3d525a79d4ae493de8c527540edef
   :target: https://www.codacy.com/app/jason-neal/eniric?utm_source=github.com&utm_medium=referral&utm_content=jason-neal/eniric&utm_campaign=badger
   :alt: Codacy Badge


.. image:: https://badge.fury.io/py/eniric.svg
   :target: https://badge.fury.io/py/eniric
   :alt: PyPI version


Welcome to Eniric's documentation!
==================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Installation
   Configuration
   Basic-usage
   Broadening
   Atmospheric-Transmission
   Theoretical-Precision-of-PHONEIX-ACES-Spectra


`Eniric` is a Python 3.6+ software written to calculate the theoretical Radial Velocity (RV) precision of Near-InfraRed (NIR) stellar spectra.
`Eniric` is an overhaul and extension to the code initially used in `Figueira et al. 2016`_ to analysis the precision of M-dwarf stars in the NIR, extending its ability to use any synthetic spectra from the PHOENIX-ACES and BT-Settl libraries, or user provided spectra, and making it easier/faster to use.
Extending the performance and usability, it is able to be used on any synthetic spectra from the `PHOENIX-ACES`_ and `BT-Settl`_ (CIFIST2001-2015) libraries.

To get started see the :doc:`Installation <./Installation>`, or the :doc:`Basic Usage <./Basic-usage>`.

Features
--------

`Eniric` contains a number of features to transform and prepare the spectra; both observed and synthetic.

* `Spectral broadening <https://github.com/jason-neal/eniric/wiki/Broadening>`_\ :
   Allows for Rotational and Instrumental broadening of synthetic spectra given rotation speed ``vsini`` and resolution ``R``.
* `Atmospheric transmission <https://github.com/jason-neal/eniric/wiki/Atmospheric-Transmission>`_ masking:
    Analyzing the RV precision attainable under the different masking conditions presented in `Figueira et al. 2016`_.

  * No treatment of atmospheric transmission
  * Masking all regions affected by atmospheric absorption of a given % over the course of the year.
  * Assuming perfect telluric correction in which the variance of the measured flux is impacted.

* Relative RV precision
   The RV precision are present relative to a specified SNR per pixel in the center of a photometric band.
   The default as used in the `Figueira et al. 2016`_ is a SNR of 100 at the center of the J-band.
* Spectral Re-sampling
   Allows for re-sampling of spectra to `N` pixels per FWHM. Default is 3.
* Photometric Band selection

  * Analysis the RV precision attainable by individual photometric bands ``Z``\ , ``Y``\ , ``J``\ , ``H``\ , ``K``.
  * User definable

* Incremental quality/precision
    Determine the incremental spectral quality or RV precision on narrow wavelength slices across the entire spectrum, similar to that present in Figure 1 of `Artigau et al. 2018 <http://adsabs.harvard.edu/abs/2018AJ....155..198A>`_.
* Analyse precision of synthetic libraries
    - Available through Starfish's `grid_tools <https://iancze.github.io/Starfish/current/grid_tools.html>`_.
       - `PHOENIX-ACES`_
       - `BT-Settl`_


.. figure:: _static/precisions.png
   :target: https://github.com/jason-neal/eniric/blob/master/paper/precisions.png
   :alt: Example relative precisions
   :align: center

   Precision achieved with *eniric* as a function of spectral band for stars with a rotational velocity of vsini=1.0\ kms$^{-1}$ and temperatures 3900 K, 3500 K, 2800 K, 2600 K, corresponding to spectral types M0, M3, M6, and M9 respectively.
   The dashed line represents the theoretical limits imposed by condition 1, and the filled area represents the values within the limits set by conditions 2 (circles) and 3 (triangles); blue, green, and red represent the results obtained for resolutions of 60000, 80000, and 100000, respectively.
   The spectra were normalized to have a SNR of 100 per resolution element as measured at the center of the J-band.
   This is similar to Figure 1 from `Figueira et al. 2016`_ but with updated precision values.


Support
-------

If you are having issues, please let us know.
Submit an issue on `Github <https://github.com/jason-neal/eniric/issues>`_.

License
-------

The project is licensed under the `MIT License <https://github.com/jason-neal/eniric/LICENCE>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`Figueira et al. 2016`: http://dx.doi.org/10.1051/0004-6361/201526900
.. _`PHOENIX-ACES`: http://phoenix.astro.physik.uni-goettingen.de
.. _`BT-Settl`: https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/
