************
Installation
************

You can install ``Eniric`` by cloning the current ``master`` branch.

It is recommended to use a ``conda`` or ``virtualenv`` environment.
To use the most up-to-date packages install the pinned requirements from the ``requirements_dev.txt`` file

.. code-block:: bash

   git clone https://github.com/jason-neal/eniric
   cd eniric
   pip install -r requirements_dev.txt
   python setup.py develop

Installation from the github repository should also be possible.

.. code-block:: bash

    pip install https://github.com/jason-neal/eniric/archive/develop.zip#egg=eniric


Requirements
------------

The requirements for ``eniric`` are:

* astropy
* joblib
* matplotlib
* multiprocess
* numpy
* pandas
* pyyaml
* scipy
* `Starfish`__
* tqdm

There are some specific version requirements given in ``requirements.txt`` but the most recent versions are pinned in ``requirements_dev.txt``
to install these run

.. code-block:: bash

     pip install -r requirements_dev.txt

from the cloned ``eniric`` repo.


Starfish
^^^^^^^^

Eniric makes use of Starfish's ``GridTools`` module to load the synthetic libraries: PHOENIX-ACES and BT-Settl models.

``Starfish`` is currently going through API changes and not yet available on pypi.

A custom fixed branch of Starfish can be used on both Linux and Windows.
``Starfish`` should be automatically installed during the installation of ``eniric``, but if not you can install it with.

.. code-block:: bash

    pip install https://github.com/jason-neal/Starfish/archive/eniric_suitable.zip#egg=Starfish

If issues arise regarding ``Starfish`` see `github.com/iancze/Starfish <Starfishgithub_>`_:,

Other requirements required with Starfish are:

*   corner
*   cycler
*   cython
*   emcee
*   h5py
*   kiwisolver
*   pyparsing
*   python-dateutil
*   pyyaml


OS
--

Eniric has been tested to work on both  **Linux** and **Windows**.

.. _Starfishgithub: https://github.com/iancze/Starfish.git

__ Starfishgithub_
