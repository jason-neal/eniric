# Script to Download the data needed for Eniric on Windows.

# Invoke this script from the top eniric directory.

$HOMEDIR = Convert-Path "."

$DIRECTORY1 = ".\data\"
$DIRECTORY2 = ".\data\atmmodel\"
$DIRECTORY3 = ".\tests\data\phoenix-raw\"
$DIRECTORY4 = ".\tests\data\btsettl-raw\"


# Precision data
Write-Host $DIRECTORY1
if ( -not (Test-Path $DIRECTORY1 -PathType Container))
{
Write-Host "$DIRECTORY1 is not present. Creating"
mkdir $DIRECTORY1
}
cd $DIRECTORY1
wget "https://www.dropbox.com/s/raw/i4cxjrhcbx6e37x/precision.tar.gz" -UseBasicParsing -OutFile ".\precision.tar.gz"
dir
python ..\scripts\untar_here.py "precision.tar.gz"
cd $HOMEDIR


# Atmmodel data
Write-Host $DIRECTORY2
if ( -not (Test-Path $DIRECTORY2 -PathType Container))
{
Write-Host "$DIRECTORY2 is not present. Creating"
mkdir $DIRECTORY2
}
wget "https://www.dropbox.com/s/raw/uab283lmkaptsib/Average_TAPAS_2014.dat.tar.gz" -UseBasicParsing -OutFile "$DIRECTORY2\Average_TAPAS_2014.dat.tar.gz"
cd $DIRECTORY2
python ..\..\scripts\untar_here.py "Average_TAPAS_2014.dat.tar.gz"
cd $HOMEDIR


# ACES test data
Write-Host $DIRECTORY3
if ( -not (Test-Path $DIRECTORY3 -PathType Container))
{
Write-Host "$DIRECTORY3 is not present. Creating"
mkdir $DIRECTORY3
}
python .\scripts\download\download_test_aces.py -o $DIRECTORY3


# BTsettl test data
Write-Host $DIRECTORY4
if ( -not (Test-Path $DIRECTORY4 -PathType Container))
{
Write-Host "$DIRECTORY4 is not present. Creating"
mkdir $DIRECTORY4
}
cd $DIRECTORY4
wget https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz -UseBasicParsing -OutFile "lte039.0-4.5-0.0a+0.0.BT-Settl.spec.fits.gz"
cd $HOMEDIR
