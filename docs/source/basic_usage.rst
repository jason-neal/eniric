===========
Basic Usage
===========
To calculate the RV precision given a spectrum with ``wavelength`` and ``flux`` you can use the ``rv_precision`` function.

.. code-block:: python

   from eniric.precision import rv_precision
   rv = rv_precision(wavelength, flux)

or for the flux independent spectral ``quality``, Q,

.. code-block:: python

   from eniric.precision import quality
   q = quality(wavelength, flux)


.. autofunction:: eniric.precision.rv_precision
.. autofunction:: eniric.precision.quality


.. _atmospheremasking_ref:

Masking
-------
It is possible to include custom pixel masks both :func:`rv_precision` and :func:`quality`, to selectively select, exclude, or weight the spectra.
These are passed as the 3rd argument to either function.

.. code-block:: python

   rv = rv_precision(wavelength, flux, mask=my_mask)
   q = quality(wavelength, flux, mask=my_mask)

Not supplying a mask or setting ``mask=None`` is equivalent to a mask of all ``1``\ 's which will have no impact on the precision (i.e. all pixels are weighted equally).

To use the default (supplied) telluric line model to apply spectral masking you can use the Atmosphere class:

.. code-block:: python

   from eniric.atmosphere import Atmosphere
   # Assuming K is the correct band, and you want to mask for seasonal variation.
   atm = Atmosphere.from_band("K", bary=True)

   # Obtain closet telluric model values at the wavelength values (this telluric mask is super sampled).
   atm = atm.at(wavelength)

   # Boolean telluric line mask at 2% deep.
   rv2 = rv_precision(wavelength, flux, mask=atm.mask)

   # Perfect telluric correction mask.
   rv3 = rv_precision(wavelength, flux, mask=atm.transmission**2)

For more details about about the masks used here see :ref:`Masking_ref`.
