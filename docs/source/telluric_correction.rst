
Atmospheric Transmission
========================

For ground-based observations it is important to understand the how the atmospheric absorption affects the Radial Velocity precision. ``Eniric`` does this in two scenarios; masking out all wavelength regions affected by telluric lines over adjusting for barycentric motion, or assuming perfect atmospheric correction, in which the variance of the photon noise is adjusted.

For this ``eniric`` contains a average telluric transmission model. It is an average model was created from weekly transmission spectra in 2014 simulated with the `TAPAS`_ software (\ `Bertaux et al 2014`_\) as detailed in (`Figueira et al. 2016`_).
The main parmeters are


* La Silla Observatory
* airmass = 1.2 (z = 33.5 degrees)
* R = 100,000
* sampling = 10

This average model that is used is available at ``data/atmmodel/Average_TAPAS_2014.txt``. It has 4 columns, wavelength (in nano-meters), transmission (between 0-1), std, mask (0 or 1). The mask is 0 for telluric lines deeper than 2%, 1 elsewhere. The std is the standard deviation of the transmission across the year but is not used in deriving precision.

You may supply your own transmission atmospheric models if you wish., e.g. ``my_telluric_model.txt``
To make individual precision calculations faster this model can be split into the separate spectroscopic bands using the scrip ``scripts/split_atmmodels.py``\ , making smaller files to load later, but taking up more disk space.

e.g.
Use ``split_atmmodels.py -h`` to get the description of flags to use.

.. code-block:: bash

   split_atmmodels.py My_telluric_model.txt -b H K

Will create the subsets of ``my_telluric_model_H.txt`` and ``my_telluric_model_K.txt`` from ``my_teluric_model.txt``.
If the the band subset cannot be found the ``my_telluric_model.txt`` will be used as a fallback.

You can configure the location of and atmmodel to use in config.yaml.

Barycentric Motion
==================

To exclude wavelength regions that will be affected by telluric lines at some point during the year you can extend out the telluric maskk.

.. code-block:: python

   new_mask = barycentric_broaden(wav, transmission, mask)

This extends the regions that are masked out, you can check that the mask continues to mask out deep lines like so...

.. code-block:: python

   assert np.all(transmission[mask] > 0.98)     # check old mask works
   assert np.all(transmission[new_mask] > 0.98) # Check new mask
   # More points should be zero in extended mask
   assert np.sum(new_mask) < np.sum(mask)

The fraction ``(np.sum(new_mask) - np.sum(mask)) /np.sum(mask)`` can also indicate the increase in masked out wavelengths.


Module
------
.. automodule:: eniric.atmosphere


.. _TAPAS: https://cds-espri.ipsl.upmc.fr/tapas/
.. _`Bertaux et al 2014`: http://adsabs.harvard.edu/abs/2014A%26A...564A..46B
.. _`Figueira et al. 2016`: http://dx.doi.org/10.1051/0004-6361/201526900
