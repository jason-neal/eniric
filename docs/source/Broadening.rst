
Broadening
==========

Two mechanisms are used to broaden the spectra:


* Rotational broadening

  * Convolution of the spectra with a rotational broadening function with a given velocity ``vsini`` (in km/s).
    The default limb darkening coefficient is ´epsilon= 0.6´.

* Instrumental broadening

  * Convolution by a Gaussian with a FWHM equivalent to the resolving power ``R`` at each wavelength. The convolution is extended either side to a ``fwhm_lim`` of 5 by default.

When analyzing the spectral libraries, rotational broadening is preformed first, followed by the instrumental broadening.

Our convolution functions use wavelength dependent kernels and do not require uniform spacing between points, unlike `PyAstronomy <https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/broadening.html>`_.
This means our convolutions are slower but are more precise. We compare the convolutions in the `Convolution_speeds.ipynb <https://github.com/jason-neal/eniric/blob/develop/docs/Notebooks/Convolution_speeds.ipynb>`_ notebook.

Caching of the convolution stages is performed to avoid re-computation of this slow component when possible using `\ ``joblib.Memory``  <https://joblib.readthedocs.io/en/latest/memory.html>`_. The caching directory defaults to ``~/.joblib`` but can be changed in the configuration file ``config.yaml``.

Caching can be disabled by setting ``location=None``


.. automodule:: eniric.broaden

