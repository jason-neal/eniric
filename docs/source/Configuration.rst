
Configuration
=============

Configuration is preformed using a config.yaml placed in the directory you wish to run eniric from.
It configures the path locations to the spectral libraries, as well as user spectral band configuration.

config.yaml

::

   # YAML configuration script

   name: eniric_default

   # The parameters defining the spectral libraries live here.
   paths:
     phoenix_raw: "../../data/PHOENIX-ALL/PHOENIX"
     phoenix_dat: "data/PHOENIX-ACES_spectra"
     results: "data/results"
     resampled: "data/resampled"
     test_data: "data/test_data/"
     atmmodel: "data/atmmodel"
     precision_results: "precision/"

   bands:
     all: ["VIS", "K", "H", "J", "Y", "Z", "CONT", "NIR", "GAP"]


   # Keywords needed for Starfish to be used
   outdir : results/

   plotdir : plots/

   # Starfish parameters defining your raw spectral library live here.
   grid:
       raw_path: "../../data/PHOENIX-ALL/PHOENIX"
       parname: ["temp", "logg", "Z"]
       hdf5_path: ""
       # raw_path: "../libraries/raw/CIFIST/"
       key_name: "t{0:.0f}g{1:.1f}" # Specifies how the params are stored
       parname: ["temp", "logg", "Z"]
       parrange: [[2300, 3700], [4.0, 5.5], [-1.5, 1.5]]
       wl_range: [6300, 6360]
       buffer: 50. # AA

   data:
       grid_name: "Eniric"
       instruments: []
