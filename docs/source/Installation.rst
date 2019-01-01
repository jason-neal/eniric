
Installation
^^^^^^^^^^^^

You can install ``Eniric`` by cloning the current ``master`` branch.

It is recommended to use a ``conda`` or ``virtualenv`` environment.
To use the most up-to-date packages install requirements from the ``requirements.txt`` file

.. code-block::

   git clone https://github.com/jason-neal/eniric
   cd eniric
   pip install -r requirements.txt
   python setup.py develop



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
* `Starfish <https://github.com/iancze/Starfish.git>`_
* tqdm

There are some specific version requirements given in ``requirements.txt`` but the most recent versions are pinned in ``requirements_dev.txt``
to install these you will need to run
     ``pip install -r requirements_dev.txt``
from the cloned ``eniric`` repo.


Starfish
~~~~~~~~

Eniric makes use of Starfish's ``GridTools`` module to load the synthetic libraries: PHOENIX-ACES and BT-Settl models.
``Starfish`` is only a requirement if using the ``eniric.utilities.load_aces_spectrum``\ , ``eniric.utilities.load_btsettl_spectrum`` functions or the ``phoenix_precision.py`` script.

``Starfish`` should be automatically installed during installation of ``eniric``. If issues arise regarding ``Starfish`` see `github.com/iance/Starfish <https://github.com/iancze/Starfish>`_\ :

.. code-block::

    pip install https://github.com/iancze/Starfish/archive/develop.zip#egg=Starfish


If are using Windows you will need to clone a fork of ``Starfish`` and checkout the ``windows`` branch

.. code-block::

   git https://github.com/jason-neal/Starfish
   cd Starfish
   git checkout windows
   python setup.py build_ext --inplace
   python setup.py develop


This is due to compiling the extensions on Windows.

.. note::

   This information on Windows is out of date.

OS
~~

Eniric has been tested to work on both  **Linux** and **Windows**.

