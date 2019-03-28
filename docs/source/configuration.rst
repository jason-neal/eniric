.. _config_ref:

=============
Configuration
=============

.. py:module:: eniric._config
    :synopsis: Handle eniric configuration.

Configuration is preformed using a ``config.yaml``  placed in the working directory you wish to run ``eniric`` from.
It configures the path locations to the spectral libraries and allows for user defined spectral bands.

If a ``config.yaml``  file does not exist then the default is used, located at ``eniric/config.yaml``.

You can copy the default configuration to a directory located at ``<path_to_dir>`` using ``config.copy_file(<path_to_dir>)``

For example, to copy it to the current directory from the command line use:

.. code-block:: bash

    $ python -c "from eniric import config; config.copy_file('.')"

The configuration values can be changed in python and saved back to the config file.
For example

.. code-block:: python

    from eniric import config
    config.paths["precision"] = ["new", "path", "to", "precision"]  # or "new/path/to/precision"
    config.update()

will update precision path in the configuration file.


.. note:: The default ``config.yaml`` file cannot be overwritten.


Eniric configuration
--------------------
The default configuration file is shown below. The comments explain the keywords needed.

.. literalinclude:: ../../eniric/config.yaml
    :linenos:
    :language: yaml


Config Class
------------
.. autoclass:: eniric._config.Config
    :members:

    .. automethod:: copy_file(directory=os.getcwd(), switch=True)
