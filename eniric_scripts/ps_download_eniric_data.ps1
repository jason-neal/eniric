# Script to Download the data needed for Eniric on Windows.

$HOMEDIR = Convert-Path "."

# Atmmodel data

$DIRECTORY = ".\data\atmmodel\"

Write-Host $DIRECTORY

if ( -not (Test-Path $DIRECTORY -PathType Container))
{
Write-Host "$DIRECTORY is not present. Creating"
mkdir $DIRECTORY
}

wget "https://www.dropbox.com/s/raw/uab283lmkaptsib/Average_TAPAS_2014.dat.tar.gz" -UseBasicParsing -OutFile ".\data\atmmodel\Average_TAPAS_2014.dat.tar.gz"


# Precision data 
cd ".\data\"
wget "https://www.dropbox.com/s/raw/i4cxjrhcbx6e37x/precision.tar.gz" -UseBasicParsing -OutFile ".\precision.tar.gz"
dir
python ..\eniric_scripts\untar_here.py "precision.tar.gz"
cd $HOMEDIR


# Test data

$DIRECTORY = ".\data\test_data\"

Write-Host $DIRECTORY

if ( -not (Test-Path $DIRECTORY -PathType Container))
{
Write-Host "$DIRECTORY is not present. Creating"
mkdir $DIRECTORY
}


cd $DIRECTORY
dir
# obsolete data
wget "https://www.dropbox.com/s/raw/oq2x7dsjeuxrf7t/obsolete.tar.gz" -UseBasicParsing -OutFile "obsolete.tar.gz"
python ..\..\eniric_scripts\untar_here.py ".\obsolete.tar.gz"
del "obsolete.tar.gz"

# Phoenix-raw data  and BT-settl raw data
wget "https://www.dropbox.com/s/raw/skg8zwi7vnxgesj/data_raw.tar.gz" -UseBasicParsing -OutFile "phoenix-raw.tar.gz"
dir
python ..\..\eniric_scripts\untar_here.py ".\data_raw.tar.gz"
del "data_raw.tar.gz"

cd $HOMEDIR

# # Phoenix-raw data
# wget "https://www.dropbox.com/s/raw/09fqard9l5s87zo/phoenix-raw.tar.gz" -UseBasicParsing -OutFile "phoenix-raw.tar.gz"
# dir
# python ..\..\eniric_scripts\untar_here.py ".\phoenix-raw.tar.gz"
# del "phoenix-raw.tar.gz"
#
# cd $HOMEDIR
#
# # BTsettl
#
# $DIRECTORY = ".\data\test_data\btsettl-raw\"
#
# Write-Host $DIRECTORY
#
# if ( -not (Test-Path $DIRECTORY -PathType Container))
# {
# Write-Host "$DIRECTORY is not present. Creating"
# mkdir $DIRECTORY
# }
# cd $DIRECTORY
# wget https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz -UseBasicParsing -OutFile "lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz"
# cd $HOMEDIR