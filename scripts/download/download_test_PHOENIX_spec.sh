#!/usr/bin/env bash

# Inspired by Starfish
# Invoke this script from the top directory.


# Check for an --no-check-certificates argument
function display_usage(){
   echo "
   Usage:
       $0 [--no-check-certificates]

   Downloads phoenix test spectra from dropbox and the phoenix aces library."
}

# check whether user had supplied -h or --help . If yes display usage
if [[ ( "$@" == "--help" ) ||  ( "$@" == "-h" ) ]]
    then
		display_usage;
		exit 0
fi

# Allow one or zero args only
if [ $# -gt 1 ];
    then
    display_usage;
    exit 1
fi

# Check if $1 equals --no-check-certificates,
# otherwise an unknown flag was used:
if [ $# -eq 0 ];
    then
    certcheck=""
elif [ "$1" = "--no-check-certificates" ];
    then
    certcheck=$1
else
    display_usage;
    exit 1
 fi


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

#(
#  cd $DIRECTORY

  # Download the wavelength file
#  wget $certcheck http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits

  # Create the Z=0.0, directories, and download the appropriate spectra
#  Z0="Z-0.0"

#  mkdir $Z0
#  cd $Z0
#  wget $certcheck http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
#  wget $certcheck http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
#  wget $certcheck http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
#  wget $certcheck http://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits
#)
download_test_aces.py

(
  cd $DIRECTORY2
  # Download the BT-Settl
  wget $certcheck https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz
)
