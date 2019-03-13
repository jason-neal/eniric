#!/usr/bin/env bash

# Inspired by Starfish
# Invoke this script from the top directory.

DIRECTORY="tests/data/phoenix-raw/"
DIRECTORY2="tests/data/btsettl-raw/"

# Check to see if libraries/raw/ directory exists, if not, make it.
if [ ! -d "$DIRECTORY" ]; then
  echo $DIRECTORY does not exist, creating.
  mkdir -p $DIRECTORY
fi
if [ ! -d "$DIRECTORY2" ]; then
  echo $DIRECTORY2 does not exist, creating.
  mkdir -p $DIRECTORY2
fi

(
  cd $DIRECTORY

  # Download the wavelength file
  wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits

  # Create the Z=0.0, directories, and download the appropriate spectra
  Z0="Z-0.0"

  mkdir $Z0
  cd $Z0
  wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
  wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
  wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
  wget ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
)

(
  cd $DIRECTORY2
  # Download the BT-Settl
  wget https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz
)
