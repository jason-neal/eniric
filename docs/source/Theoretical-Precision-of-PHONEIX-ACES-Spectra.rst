
Theoretical Precision of PHONEIX-ACES Spectra
---------------------------------------------

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

You use this script by passing it stellar parameters for the spectral library as well as broadening and wavelength parameters. Multiple parameters can be passed to each flag, separated by a space.
To see all the available inputs parameters run:
.. code-block:: bash

     phoenix_precision.py -h

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