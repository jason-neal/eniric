.. _install_ref:

============
Installation
============

You can install ``eniric`` by cloning the current ``master`` branch.

It is recommended to use a ``conda`` or ``virtualenv`` environment.
To use the most up-to-date packages install the pinned requirements from the ``requirements_dev.txt`` file

.. code-block:: bash

   git clone https://github.com/jason-neal/eniric
   cd eniric
   pip install -r requirements_dev.txt
   python setup.py develop

Installation from the github repository should also be possible.

.. code-block:: bash

    pip install https://github.com/jason-neal/eniric/archive/master.zip#egg=eniric

or for the develop branch

.. code-block:: bash

    pip install https://github.com/jason-neal/eniric/archive/develop.zip#egg=eniric


Requirements
------------

The requirements for ``eniric`` are:

* astropy
* joblib>=0.12.3
* matplotlib
* multiprocess
* numpy>=0.15.4
* pandas
* pyyaml
* oyaml
* scipy
* `Starfish`__
* tqdm

There are some specific version requirements given in ``requirements.txt`` but the most recent versions are pinned in ``requirements_dev.txt``
to install these run

.. code-block:: bash

     pip install -r requirements_dev.txt

from the cloned ``eniric`` repository.


Starfish
^^^^^^^^

``Eniric`` makes use of Starfish's ``GridTools`` module to load the synthetic libraries: PHOENIX-ACES and BT-Settl models.

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

``Eniric`` has been tested to work on both  **Linux** and **Windows**.

.. _Starfishgithub: https://github.com/iancze/Starfish.git

__ Starfishgithub_


Eniric Data
-----------
There are some data files not included in the ``eniric`` repository that are *necessary for testing*.
These can be easily downloaded using the provided scripts.

.. code-block:: bash

    $ download_eniric_data.sh

.. Note:: If you have an issue connecting to dropbox you can also try
        adding the ``--no-check-certificate`` flag.

or on **Windows** in a PowerShell

.. code-block:: bash

    ps_download_eniric_data.ps1

This includes an atmospheric transmission spectrum, located at ``data/atmos/Average_TAPAS_2014.dat``, which can be used for spectral masking.

.. Note:: This is to keep the size of the git repository small.


Testing
-------
To test ``eniric`` is installed try

.. code-block:: bash

    python -c "import eniric"


To run the test suite run ``pytest`` from the root directory of the repository (requires pytest).
This will result in an output similar to:

.. code-block:: text

    $ pytest

    ============================= test session starts ==============================
    platform linux -- Python 3.6.7, pytest-4.3.0, py-1.7.0, pluggy-0.8.0
    hypothesis profile 'default' -> database=DirectoryBasedExampleDatabase('/home/travis/build/jason-neal/eniric/.hypothesis/examples')
    rootdir: /home/travis/build/jason-neal/eniric, inifile: setup.cfg
    plugins: cov-2.6.1, hypothesis-4.7.17
    collected 718 items

    ...

    ======= 610 passed, 84 xfailed, 24 xpassed, 2 warnings in 33.00 seconds ========
   The command "pytest" exited with 0.


The requirements for the test suite can be installed from the root of the repository using

.. code-block:: bash

   $ python setup.py install .[test]

.. Note:: A users copied ``config.yaml`` file in the repository home directory may interfere with the test results, causing some failures.
