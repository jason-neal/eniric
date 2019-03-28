#!/usr/bin/env bash

# Inspired by Starfish
# Invoke this script from the top directory.


# Check for an --no-check-certificates argument

function display_usage(){
   echo "
   Usage:
       $0 [--no-check-certificates]

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

### DOWNLOAD Transmission spectrum
DIRECTORY="data/atmmodel/"

# Check to see if libraries/raw/ directory exists, if not, make it.
if [ ! -d "$DIRECTORY" ]; then
  echo $DIRECTORY does not exist, creating.
  mkdir -p $DIRECTORY
fi
(cd $DIRECTORY
    wget $certcheck "https://www.dropbox.com/s/raw/uab283lmkaptsib/Average_TAPAS_2014.dat.tar.gz"
    ../../scripts/untar_here.py Average_TAPAS_2014.dat.tar.gz
)


DIRECTORY2="tests/data/"
# Check to see if libraries/raw/ directory exists, if not, make it.
if [ ! -d "$DIRECTORY2" ]; then
  echo $DIRECTORY2 does not exist, creating.
  mkdir -p $DIRECTORY2
fi


DIRECTORY3="data/"

# Check to see if libraries/raw/ directory exists, if not, make it.
if [ ! -d "$DIRECTORY3" ]; then
  echo $DIRECTORY2 does not exist, creating.
  mkdir -p $DIRECTORY3
fi

(cd $DIRECTORY3
    wget $certcheck "https://www.dropbox.com/s/raw/i4cxjrhcbx6e37x/precision.tar.gz"
   ../scripts/untar_here.py precision.tar.gz
    rm precision.tar.gz
)


# This should download the PHOENIX spectra
download_test_PHOENIX_spec.sh

# FTP on travis gives several delays when getting phoenix data.
# So putting Phoenix data in dropbox now also.
#(cd $DIRECTORY2
#    wget $certcheck "https://www.dropbox.com/s/raw/skg8zwi7vnxgesj/data_raw.tar.gz"
#    ../../scripts/untar_here.py data_raw.tar.gz
#    rm data_raw.tar.gz
#)
