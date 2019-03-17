==========================================
Theoretical Precision of Synthetic Spectra
==========================================

The theoretical precision of M-dwarf synthetic spectra is extensively explored in `Figueira et al. (2016)`_.
``Eniric`` extends this work to any spectra in the `PHOENIX-ACES <http://phoenix.astro.physik.uni-goettingen.de/>`_ or the `BT-Settl-CIFIST2011_2015 <https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/>`_ spectral libraries.

The script ``phoenix_precision.py`` is provided to easily generate synthetic RV precision values, similarly to and beyond the work of `Figueira et al. (2016)`_.
The precision of a spectrum can be obtained by providing its library parameters Teff, logg, Fe/H and, alpha. This has been mainly used on M-dwarf spectra with temperatures < 4000K, but higher temperatures also work.

For the library selection and loading of the spectra `Starfish's Grid Tools <https://iancze.github.io/Starfish/current/grid_tools.html>`_ is used.
The spectral libraries need to be pre-downloaded, and the ``raw_path`` to their location configured in the ``config.yaml`` file in the working directory.
For more information see the Starfish documentation `here <https://iancze.github.io/Starfish/current/grid_tools.html#downloading-model-spectra>`_.


Example Usage
-------------

``phoenix_precision.py`` is called by passing the stellar and other parameters on the command line.
Several flags can take multiple values separated by a space. All combinations of the parameters supplied will be computed.

The available inputs parameters are:

.. code-block:: text

    $ phoenix_precision.py -h

    usage: phoenix_precision.py [-h] [-t TEMP [TEMP ...]] [-l LOGG [LOGG ...]]
                                [-m METAL [METAL ...]] [-a ALPHA [ALPHA ...]]
                                [-s SAMPLING [SAMPLING ...]]
                                [-r RESOLUTION [RESOLUTION ...]]
                                [-v VSINI [VSINI ...]]
                                [-b {K,H,J,Y,Z,TEST} [{K,H,J,Y,Z,TEST} ...]]
                                [--model {aces,btsettl,phoenix}] [--snr SNR]
                                [--ref_band {SELF,K,H,J,Y,Z,TEST}]
                                [--num_procs NUM_PROCS] [-o OUTPUT]
                                [--rv RV [RV ...]] [--add_rv] [--air] [--correct]
                                [-V] [--disable_normalization]

    Calculate precision and quality for synthetic spectra.

    optional arguments:
      -h, --help            show this help message and exit
      -t TEMP [TEMP ...], --temp TEMP [TEMP ...]
                            Temperature, default=[3900].
      -l LOGG [LOGG ...], --logg LOGG [LOGG ...]
                            Logg, default = [4.5].
      -m METAL [METAL ...], --metal METAL [METAL ...]
                            Metallicity, default=[0.0].
      -a ALPHA [ALPHA ...], --alpha ALPHA [ALPHA ...]
                            Alpha, default=[0.0].
      -s SAMPLING [SAMPLING ...], --sampling SAMPLING [SAMPLING ...]
                            Sampling, default=[3.0].
      -r RESOLUTION [RESOLUTION ...], --resolution RESOLUTION [RESOLUTION ...]
                            Instrumental resolution, default=[50000]
      -v VSINI [VSINI ...], --vsini VSINI [VSINI ...]
                            Rotational Velocity, default = [1.0]
      -b {K,H,J,Y,Z,TEST} [{K,H,J,Y,Z,TEST} ...], --bands {K,H,J,Y,Z,TEST} [{K,H,J,Y,Z,TEST} ...]
                            Wavelength bands to select, default=['J'].
      --model {aces,btsettl,phoenix}
                            Spectral models to use, default='aces'.
      --snr SNR             Mid-band SNR scaling, default=100.
      --ref_band {SELF,K,H,J,Y,Z,TEST}
                            SNR reference band, default='J'.
                            'self' scales each band relative to the SNR itself.
      --num_procs NUM_PROCS
                            Number of processors to use,
                            default = (Total cores - 1)
      -o OUTPUT, --output OUTPUT
                            Filename for result file, default='precisions.csv'.
      --rv RV [RV ...]      Radial velocity value, default=[0.0]
      --add_rv              Include a radial velocity shift.
      --air                 Convert to air wavelengths.
      --correct             Apply Artigau et al. (2018) RV band corrections.
      -V, --verbose         Enable verbose.
      --disable_normalization
                            Disable the convolution normalization.


.. Note:: The parameter space is multiplicative so the runtime increases when increasing the number of values for each parameter.


Output File
-----------

The script returns a table with input parameters and the calculated precisions for all three tellruic conditions presented in Figueira et al. (2016) for the parameter combinations requested.
E.g:

::

   # precisions.csv
   temp,logg,feh,alpha,band,resolution,vsini,sampling,correctflag,quality,cond1,cond2,cond3
   3900,4.5,0.0,0.0,K,100k,1.0,3.0,0,4916,7.4,33.6,8.0

The different columns of the output file are given in the table below.

+--------+--------------+----------------------------------------------------------------------+
| Col.   | Name         | Description                                                          |
+========+==============+======================================================================+
| 1      | temp         | Library stellar effective temperature (K).                           |
+--------+--------------+----------------------------------------------------------------------+
| 2      | logg         | Library stellar surface gravity.                                     |
+--------+--------------+----------------------------------------------------------------------+
| 3      | feh          | Library stellar metallicity.                                         |
+--------+--------------+----------------------------------------------------------------------+
| 4      | alpha        | Library stellar alpha ratio.                                         |
+--------+--------------+----------------------------------------------------------------------+
| 5      | band         | Wavelength band letter.                                              |
+--------+--------------+----------------------------------------------------------------------+
| 6      | resolution   | Instrumental resolution.                                             |
+--------+--------------+----------------------------------------------------------------------+
| 7      | vsini        | Stellar rotation (km/s).                                             |
+--------+--------------+----------------------------------------------------------------------+
| 8      | sampling     | Spectral sampling - N points per resolution element.                 |
+--------+--------------+----------------------------------------------------------------------+
| 9      | correctflag  | Indicate if `Artigau et al. (2018)`_ precision correction is applied.|
+--------+--------------+----------------------------------------------------------------------+
| 10     | quality      | Theoretical spectral quality.                                        |
+--------+--------------+----------------------------------------------------------------------+
| 11     | cond1        | RV precision with no masking (m/s). (Condition 1)                    |
+--------+--------------+----------------------------------------------------------------------+
| 12     | cond2        | RV precision with binary masking (m/s). (Condition 2)                |
+--------+--------------+----------------------------------------------------------------------+
| 13     | cond3        | RV precision with transmission masking (m/s). (Condition 3)          |
+--------+--------------+----------------------------------------------------------------------+

The first 9 columns uniquely identify a set of input parameter values, this is used to avoid repeating an identical computaion.
In this way ``precsions.csv`` can be appended to with new values, while keeping the other values, if desired.


Calculating Precisions
----------------------
Below are some specific examples of using ``phoenix_precision.py``.

This script has been used to generate nIR RV precision values across the M-dwarf temperature range.
These were requested by the NIRPS and SPIRou consortia for use as into their respective Exposure Time Calculators.
The commands to use to generate these datasets are provided below.


Figueira et al. 2016
++++++++++++++++++++

To reproduce the calculations of `Figueira et al. (2016)`_ you can use the

.. code-block:: bash

   phoenix_precision.py -t 2600, 2900, 3500, 3900 -m 0.0, -l 4.5 --snr 100 -b Z Y J H K --ref_band J


NIRPS
+++++
For the NIRPS ETC precisions were calculated for the whole M-dwarf range between 2500 and 4000 K.
These were centred on the H-band centering with a SNR of 100. This also included R=75000 tailored to the NIRPS instrument.

.. code-block:: bash

   phoenix_precision.py -t 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    -m 0.0, -l 5.0 --snr 100 -b Z Y J H K --ref_band H -r 60000 75000 80000 100000 -v 1.0 5.0 10.0


SPIRou
++++++
For the SPIRou ETC the parameter combinations are the same as `Figueira et al. (2016)`_ but
calculated relative to a SNR of 100 in each respective bands.

.. code-block:: bash

   phoenix_precision.py -t 2600, 2900, 3500, 3900 -m 0.0, -l 4.5 --snr 100 -b Z Y J H K --ref_band self


BT-SETTL
--------
To use the BT-Settl sectral library  use the ``--model`` flag.

.. code-block:: bash

   phoenix_precision.py -t 2600, 2900, 3500, 3900 -b Z Y J H K --model btsettl


.. _`Figueira et al. (2016)`: http://dx.doi.org/10.1051/0004-6361/201526900
.. _`Artigau et al. (2018)`: http://adsabs.harvard.edu/abs/2018AJ....155..
