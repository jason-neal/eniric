===========
Basic Usage
===========

To calculate the RV precision given a spectrum with ``wavelength`` and ``flux`` you can use the ``rv_precsion`` function.

.. autofunction eniric.precision :: rv_precision

.. code-block:: python

   from eniric.precision import rv_precision
   rv = rv_precision(wavelength, flux)

or for the flux independent spectral ``quality``, Q,

.. code-block:: python

   from eniric.precision import quality
   q = quality(wavelength, flux)

To apply a mask to the optimal pixel weights pass it as the 3rd argument to either function.

.. code-block:: python

   rv = rv_precision(wavelength, flux, mask=my_mask)
   q = quality(wavelength, flux, mask=my_mask)

Not supplying a mask or setting ``mask=None`` is equivalent to a mask of all ``1``\ 's which will have no impact on the precision.

To use the default (supplied) telluric line model to apply spectral masking you can use the Atmopshere class:

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

For more details about about the mask parameter see :ref:`Masking_ref`.
