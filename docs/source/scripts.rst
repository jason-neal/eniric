=======
Scripts
=======
Scripts provided with eniric.


.. automodule:: scripts.phoenix_precision
    :members:


.. automodule:: scripts.barycenter_broaden_atmmodel
    :members:


.. automodule:: scripts.split_atmmodel
    :members:


.. automodule:: scripts.precision_four_panel
    :members:


.. figure:: _static/precisions.png
   :target: https://github.com/jason-neal/eniric/blob/master/paper/precisions.png
   :alt: Example relative precisions from ``precision_four_panel.py``
   :align: center

.. automodule:: scripts.csv2tsv
    :members:


.. automodule:: scripts.untar_here
    :members:


Download Scripts
----------------
Scripts to download the eniric data and phoenix data.
They can be run from the command line.

.. code-block:: bash

    scripts/download/download_eniric_data.sh

This is also available as a powershell script.

.. code-block:: bash

    scripts/download/ps_download_eniric_data.ps1

The test data from the PHOENIX-ACES library is specifically
downloaded using Starfish utilities in

.. code-block:: bash

    scripts/download/download_test_aces.py.
