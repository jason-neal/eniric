
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


Hopefully there will be a pip installable version soon...

Starfish
~~~~~~~~

Eniric makes use of ``Starfish``\ s' ``GridTools`` to load the synthetic library PHOENIX-ACES/BT-Settl models.
``Starfish`` is only a requirement if using the ``eniric.utilities.load_aces_spectrum``\ , ``eniric.utilities.load_btsettl_spectrum`` functions or the ``phoenix_precision.py`` script.

``Starfish`` can be downloaded and installed using the following commands, more info can be found `here <https://github.com/iancze/Starfish>`_\ :

.. code-block::

   git clone https://github.com/iancze/Starfish.git
   cd Starfish
   python setup.py build_ext --inplace
   python setup.py develop


If are using Windows you will need to clone a fork of ``Starfish`` and checkout the ``windows`` branch

.. code-block::

   git https://github.com/jason-neal/Starfish
   cd Starfish
   git checkout windows
   python setup.py build_ext --inplace
   python setup.py develop


This is due to compiling the extensions on Windows.

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
* Starfish*
* tqdm

Other requirements required for Starfish are:


* corner
* cython
* emcee
* h5py
* scikit-learn

There are no specific version requirements but the most recent version are pinned in ``requirements.txt``
to install these you will need to run
     ``pip install -r requirements.txt`` 
from the cloned ``eniric`` repo.


OS
~~

Eniric has been tested to work on both  **Linux** and **Windows**.

