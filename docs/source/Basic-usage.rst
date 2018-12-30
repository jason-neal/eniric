
Basic Usage
===========

To calculate the RV precision given a spectrum with ``wavelength`` and ``flux`` you can use the ``rv_precsion`` function

.. code-block::

   from eniric.Qcalculator import rv_precision
   rv = rv_precision(wavelength, flux)

or for the flux independent spectral ``quality``

.. code-block::

   from eniric.Qcalculator import quality
   q = quality(wavelength, flux)

To apply a mask to the optimal pixel weights pass it as the 3rd argument.

.. code-block::

   rv = rv_precision(wavelength, flux, my_mask)

or 

.. code-block::

   rv = rv_precision(wavelength, flux, mask=my_mask)

Not suppling a mask or setting ``mask=None`` is equivalent to a mask of all ``1``\ 's which will have no impact on the precision.

To use the default telluric line model to apply masking:

.. code-block::

   from eniric.atmosphere import Atmosphere
   # Assuming K is the correct band, and you want to mask seasonal variation.
   atm = Atmosphere.from_band("K, bary=True)

   # Obtain closet telluric model values at the wavelength values (telluric mask is super sampled).
   atm = atm.at(wavelength)  

   # Boolean telluric line mask at 2% deep.
   rv2 = rv_precision(wavelength, flux, mask=atm.mask)

   # Perfect telluric correction mask.
   rv3 = rv_precision(wavelength, flux, mask=atm.transmission**2)

...

As the presence of absorption lines in Earth's atmosphere blocks affects the spectral precision.
There are two options to handle this. Complete masking of telluric lines within 30 km/s, reduction in spectral variance by the transmission spectrum.

Weight Masking
~~~~~~~~~~~~~~

Three cases: 


* all 1's
* 1s and 0s
* Transmission spectrum
  The masking function is squared which only affects the 3rd option, If you want to alter the weigths in a specific way you will need to account for this when deriving the mask.
