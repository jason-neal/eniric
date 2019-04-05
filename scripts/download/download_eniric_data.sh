#!/usr/bin/env bash

# Inspired by Starfish
# Invoke this script from the top eniric directory.

# Check for an --no-check-certificate argument

function display_usage(){
   echo "
   Usage:
       $0 [--no-check-certificate]

   Downloads data for eniric from dropbox and the phoenix aces library."
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

# Check if $1 equals --no-check-certificate,
# otherwise an unknown flag was used:
if [ $# -eq 0 ];
    then
    certcheck=""
elif [ "$1" = "--no-check-certificate" ];
    then
    certcheck=$1
else
    display_usage;
    exit 1
 fi


DIRECTORY1="data/"
DIRECTORY2="data/atmmodel/"
DIRECTORY3="tests/data/phoenix-raw/"
DIRECTORY4="tests/data/btsettl-raw/"


# Precision data
if [ ! -d "$DIRECTORY1" ]; then
  echo $DIRECTORY1 does not exist, creating.
  mkdir -p $DIRECTORY1
fi
(cd $DIRECTORY1
    wget $certcheck "https://www.dropbox.com/s/raw/i4cxjrhcbx6e37x/precision.tar.gz"
   ../scripts/untar_here.py precision.tar.gz
    rm precision.tar.gz
)

# Atmmodel data
if [ ! -d "$DIRECTORY2" ]; then
  echo $DIRECTORY2 does not exist, creating.
  mkdir -p $DIRECTORY2
fi
(cd $DIRECTORY2
    wget $certcheck "https://www.dropbox.com/s/raw/uab283lmkaptsib/Average_TAPAS_2014.dat.tar.gz"
    ../../scripts/untar_here.py Average_TAPAS_2014.dat.tar.gz
)

# ACES test data
if [ ! -d "$DIRECTORY3" ]; then
  echo $DIRECTORY3 does not exist, creating.
  mkdir -p $DIRECTORY3
fi
# This should download the PHOENIX ACES spectra
download_test_aces.py -o $DIRECTORY3


# BTsettl test data
if [ ! -d "$DIRECTORY4" ]; then
  echo $DIRECTORY4 does not exist, creating.
  mkdir -p $DIRECTORY4
fi
(
  cd $DIRECTORY4
  wget $certcheck https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz
)
