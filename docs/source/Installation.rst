
Installation
^^^^^^^^^^^^

You can install ``Eniric`` by cloning the current ``master`` branch.

It is recommended to use a ``conda`` or ``virtualenv`` environment.
To use the most up-to-date packages install the pinned requirements from the ``requirements_dev.txt`` file

.. code-block:: bash

   git clone https://github.com/jason-neal/eniric
   cd eniric
   pip install -r requirements_dev.txt
   python setup.py develop

Installation from the github repository should also be possible.

.. code-block: bash

    pip install https://github.com/jason-neal/eniric/archive/develop.zip#egg=eniric


Requirements
~~~~~~~~~~~~

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
to install these you will need to run

.. code-block:: bash

     pip install -r requirements_dev.txt

from the cloned ``eniric`` repo.


Starfish
~~~~~~~~~~

Eniric makes use of Starfish's ``GridTools`` module to load the synthetic libraries: PHOENIX-ACES and BT-Settl models.
``Starfish`` is only a requirement if using the ``eniric.utilities.load_aces_spectrum``\ , ``eniric.utilities.load_btsettl_spectrum`` functions or the ``phoenix_precision.py`` script.

``Starfish`` is currently going through major changes.
``Starfish`` should be automatically installed during installation of ``eniric`` with this custom fixed branch which is suitable for both Windows and Linux installation.

.. code-block:: bash
    pip install https://github.com/jason-neal/Starfish/archive/eniric_suitable.zip#egg=Starfish

If issues arise regarding ``Starfish`` see `github.com/iancze/Starfish <Starfishgithub_>`_:,
If you need to install Starfish manually first you may need to remove the Starfish https link from ``requirements.txt``.


OS
~~

Eniric has been tested to work on both  **Linux** and **Windows**.

.. _Starfishgithub: https://github.com/iancze/Starfish.git

__ Starfishgithub_