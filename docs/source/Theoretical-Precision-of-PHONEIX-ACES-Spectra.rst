
Theoretical Precision of Synthetic Spectra
------------------------------------------

Eniric provides a script ``phoenix_precision.py`` to generate RV precision values for any spectra in the PHOENIX-ACES spectral library (Available at `http://phoenix.astro.physik.uni-goettingen.de/ <http://phoenix.astro.physik.uni-goettingen.de/>`_\ ). The precision of a spectra can be obtained by providing its library parameters Teff, logg, Fe/H and, alpha. This has been mainly used on M-dwarf spectra (< 4000K), but it higher temperatures should work to.

For the library selection and loading of the spectra we use `Starfish's Grid Tools <https://iancze.github.io/Starfish/current/grid_tools.html>`_. You need to have the PHOENIX-ACES spectra downloaded, and the ``raw_path`` configured in the ``config.yaml`` file in the working directory. For more information see the Starfish documentation `here <https://iancze.github.io/Starfish/current/grid_tools.html#downloading-model-spectra>`_.

After downloading the PHOENIX-ACES spectra spectra you need to, set the ``raw_path`` to the downloaded files in the ``config.yaml`` file in the working directory.   

The script returns a table with input parameters and calculated spectral quality and precisions for all three conditions presented in Figueira et al. 2016 for all parameter combinations requested.   

::

   # quality.csv
   Temp, logg, [Fe/H], Alpha, Band, Resolution, vsini, Sampling, Quality, Cond. 1, Cond. 2, Cond. 3, correct flag
   3900, 4.5,     0.0,   0.0,    J,       100k,   1.0,      3.0,    3501,     7.6,    20.2,     8.0,            0

Example Usage
^^^^^^^^^^^^^

You use this script by passing it stellar parameters for the spectral library as well as broadening and wavelength parameters.
Multiple parameters can be passed to each flag, separated by a space.
To see all the available inputs parameters run:

.. code-block:: bash

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

    Calculate Calculate precision and quality for synthetic spectra.

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


This script has been used to generate nIR RV precision values across the M-dwarf temperature range. These were requested by the NIRPS and SPIRou consortia for use as into their respective Exposure Time Calculators. The commands to use to generate these datasets are provided below. 

NIRPS
"""""

H-band centering, SNR 100, whole M-dwarf range 3500-4000 K

.. code-block:: bash

   phoenix_precision.py -t 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000
    -m 0.0, -l 5.0 --snr 100 -b Z Y J H K --ref H -r 60000 75000 80000 100000 -v 1.0 5.0 10.0


SPIRou
""""""

Same as `Figueira et al. 2016`_ but all spectra relative to a SNR of 100 in each respective bands.

.. code-block:: bash

   phoenix_precision.py -t 2600, 2900, 3500, 3900 -m 0.0, -l 4.5 --snr 100 -b Z Y J H K --ref self


Remember: The parameter space is multiplicative so the runtime increases when increasing the number of values for each parameter.

Future
~~~~~~

Using this script as a starting point it should be fairly straight forward to do the same for your own favorite spectra library. (Especially if it already has a Starfish interface). For example, to derive RV precisions with the BT-Settl  spectral library.


.. _`Figueira et al. 2016`: http://dx.doi.org/10.1051/0004-6361/201526900
