"""
Created on Sun Dec 14 15:43:13 2014

@author: pfigueira
"""
from __future__ import print_function, division

import numpy as np
from os import listdir
from os.path import isfile, join
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib import rc
#set stuff for latex usage
rc('text', usetex=True)

from eniric.original_code.IOmodule import read_2col, read_3col

from eniric.original_code.Qcalculator import *

data_rep = "../data/nIRmodels/"
results_dir = "../data/original_code/results/"
resampled_dir = "../data/original_code/resampled/"

#models form PHOENIX-ACES
M0_ACES = data_rep+"PHOENIX-ACES/PHOENIX-ACES-AGSS-COND-2011-HiRes/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M3_ACES = data_rep+"PHOENIX-ACES/PHOENIX-ACES-AGSS-COND-2011-HiRes/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M6_ACES = data_rep+"PHOENIX-ACES/PHOENIX-ACES-AGSS-COND-2011-HiRes/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M9_ACES = data_rep+"PHOENIX-ACES/PHOENIX-ACES-AGSS-COND-2011-HiRes/lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

def read_spectrum(spec_name):
    """
    function that reads spectra from the database
    """
    wav, flux = read_2col(spec_name)
    wav = np.array(wav, dtype="float64")*1.0e-4 # conversion to microns
    flux = np.array(flux, dtype="float64")

    flux_photons = flux * wav

    return [wav[:], flux_photons[:]]

def band_selector(wav, flux, band):
    if(band in ["ALL", "all", ""]):
        return [wav, flux]
    elif(band == "Y"):
        bandmin = 1.0
        bandmax = 1.1
    elif(band == "J"):
        bandmin = 1.17
        bandmax = 1.33
    elif(band == "H"):
        bandmin = 1.5
        bandmax = 1.75
    elif(band == "K"):
        bandmin = 2.07
        bandmax = 2.35
    else:
        print("Unrecognized band tag.")
        exit()

    #select values form the band
    wav_band, flux_band = wav_selector(wav, flux, bandmin, bandmax)

    return [wav_band, flux_band]

def plotter(spectrum, band, vsini=0, R=0):
    """
    Reads and plots the selected spectrum in a given band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)

    plt.figure(1)
    plt.xlabel(r"wavelength [ $\mu$m ])")
    plt.ylabel(r"flux [counts] ")
    plt.plot(wav_band, flux_band, color ='k', marker="o", linestyle="-")
    plt.show()
    plt.close()

def convolution(spectrum, band, vsini, R, epsilon=0.6, FWHM_lim=5.0, plot=True):

    """
    function that convolves a given spectra to a resolution of R
    R = 60 000 , R = 80 000, R = 100 000
    """

    print("Reading the data...")
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)
    print("Done.")

    #we need to calculate the FWHM at this value in order to set the starting point for the convolution
    FWHM_min = wav_band[0]/R    #FWHM at the extremes of vector
    FWHM_max = wav_band[-1]/R

    #performing convolution with rotation kernel
    print("Starting the Rotation convolution for vsini=%.2f..." % (vsini),
          "At", dt.datetime.now())
    delta_lambda_min = wav_band[0]*vsini/3.0e5
    delta_lambda_max = wav_band[-1]*vsini/3.0e5
    #widest wavelength bin for the rotation convolution
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, wav_band[0]-delta_lambda_min-FWHM_lim*FWHM_min, wav_band[-1]+delta_lambda_max+FWHM_lim*FWHM_max)
    wav_ext_rotation = np.array(wav_ext_rotation, dtype="float64")
    flux_ext_rotation = np.array(flux_ext_rotation, dtype="float64")

    #wide wavelength bin for the resolution_convolution
    wav_extended, flux_extended = wav_selector(wav, flux, wav_band[0]-FWHM_lim*FWHM_min, wav_band[-1]+FWHM_lim*FWHM_max)
    wav_extended = np.array(wav_extended, dtype="float64")
    flux_extended = np.array(flux_extended, dtype="float64")

    flux_conv_rot = []
    counter=0
    for wav in wav_extended:
    # select all values such that they are within the FWHM limits
        delta_lambda_L = wav*vsini/3.0e5
        indexes = [ i for i in range(len(wav_ext_rotation)) if ((wav - delta_lambda_L) < wav_ext_rotation[i] < (wav + delta_lambda_L))]
        flux_2convolve = flux_ext_rotation[indexes[0]:indexes[-1]+1]
        rotation_profile = rotation_kernel(wav_ext_rotation[indexes[0]:indexes[-1]+1]-wav, delta_lambda_L, vsini, epsilon)
        flux_conv_rot.append(np.sum(rotation_profile*flux_2convolve))
        if(len(flux_conv_rot)%(len(wav_extended)/100 ) == 0):
            counter = counter+1
            print("Rotation Convolution at %d%%..." % (counter))

    print("Done At", dt.datetime.now(), "\n")
    flux_conv_rot = np.array(flux_conv_rot, dtype="float64")

    print("Starting the Resolution convolution... At", dt.datetime.now())

    flux_conv_res = []
    counter=0
    for wav in wav_band:
        # select all values such that they are within the FWHM limits
        FWHM = wav/R
        indexes = [ i for i in range(len(wav_extended)) if ((wav - FWHM_lim*FWHM) < wav_extended[i] < (wav + FWHM_lim*FWHM))]
        flux_2convolve = flux_conv_rot[indexes[0]:indexes[-1]+1]
        IP = unitary_Gauss(wav_extended[indexes[0]:indexes[-1]+1], wav, FWHM)
        flux_conv_res.append(np.sum(IP*flux_2convolve))
        if(len(flux_conv_res)%(len(wav_band)/100 ) == 0):
            counter = counter+1
            print("Resolution Convolution at %d%%..." % (counter))
    flux_conv_res = np.array(flux_conv_res, dtype="float64")
    print("Done. At ", dt.datetime.now(), "\n")

    print("Saving results... At", dt.datetime.now())

    #Note: difference in sampling at 1.0 and 1.5 microns makes jumps in the beginning of Y and H bands

    name_model = name_assignment(spectrum)
    filename = results_dir+"Spectrum_"+name_model+"_"+band+"band_vsini"+str(vsini)+"_R"+str(int(R/1000))+"k.txt"
    write_3col(filename, wav_band, flux_band, flux_conv_res)
    print("Done.")

    if(plot):
        fig=plt.figure(1)
        plt.xlabel(r"wavelength [ $\mu$m ])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_band, flux_band/np.max(flux_band), color ='k', linestyle="-", label="Original spectra")
        plt.plot(wav_band, flux_conv_res/np.max(flux_conv_res), color ='b', linestyle="-", label="%s spectrum observed at vsini=%.2f and R=%d ." % (name_model, vsini, R))
        plt.legend(loc='best')
        plt.show()

        fig.savefig(filename[:-3]+"pdf", facecolor='w', format='pdf',bbox_inches='tight')
        plt.close()

    return [wav_band, flux_conv_res]


###############################################################################
#
#   Auxiliary functions
#
###############################################################################

def wav_selector(wav, flux, wav_min, wav_max):
    """
    function that returns wavelength and flux withn a giving range
    """
    wav_sel = np.array([value for value in wav if(wav_min < value < wav_max)], dtype="float64")
    flux_sel = np.array([value[1] for value in zip(wav,flux) if(wav_min < value[0] < wav_max)], dtype="float64")

    return [wav_sel, flux_sel]

def unitary_Gauss(x, center, FWHM):
    """
    Gaussian_function of area=1

    p[0] = A;
    p[1] = mean;
    p[2] = FWHM;
    """

    sigma = np.abs(FWHM) /( 2 * np.sqrt(2 * np.log(2)) );
    Amp = 1.0 / (sigma*np.sqrt(2*np.pi))
    tau = -((x - center)**2) / (2*(sigma**2))
    result = Amp * np.exp( tau );

    return result;

def rotation_kernel(delta_lambdas, delta_lambda_L, vsini, epsilon):
    """
    calculates the rotation kernel for a given wavelength
    """
    c1 = 2.0*(1.0-epsilon) / (np.pi*vsini*(1.0-epsilon/3.0))
    c2 = 0.5*np.pi*epsilon / (np.pi*vsini*(1.0-epsilon/3.0))

    return (c1*np.sqrt(1.0-(delta_lambdas/delta_lambda_L)**2.0) + c2*(1.0-(delta_lambdas/delta_lambda_L)**2.0) )


###############################################################################

def resample_allfiles(folder=results_dir):
    """
    reample all files inside folder
    """
    #getting a list of all the files
    onlyfiles = [ f for f in listdir(folder) if isfile(join(folder,f)) ]

    [resampler(results_dir+spectrum_file) for spectrum_file in onlyfiles if spectrum_file[-4:]==".txt"]
def resampler(spectrum_name = "results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R60k.txt", sampling=3.0, plottest=False):
    """
    resamples a spectrum by interpolation onto a grid with a sampling of 3 pixels per resolution element
    """
    wavelength, theoretical_spectrum, spectrum  = read_3col(spectrum_name)

    wavelength_start = wavelength[1] #because fo border effects
    wavelength_end = wavelength[-1]
    resolution_string = spectrum_name[-8:-5]

    if(resolution_string[0]=="R"):
        resolution = int(resolution_string[1:])*1000
    else:
        resolution = int(resolution_string)*1000

    wav_grid = [wavelength_start]
    while(wav_grid[-1]< wavelength_end):
         wav_grid.append(wav_grid[-1]*(1.0+1.0/(sampling*resolution)))
    wav_grid = np.array(wav_grid)

    interpolated_flux = np.interp(wav_grid, wavelength, spectrum)
    filetowrite = resampled_dir + spectrum_name[8:-4]+"_res"+str(int(sampling))+".txt"
    write_2col(filetowrite, wav_grid, interpolated_flux)

    if(plottest):
        plt.figure(1)
        plt.xlabel(r"wavelength [ $\mu$m ])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wavelength, spectrum, color ='k', linestyle="-", label="Original spectrum")
        plt.plot(wav_grid, interpolated_flux, color ='b', linestyle="-", label="Interpolated spectrum")
        plt.legend(loc='best')
        plt.show()

        plt.close()

###############################################################################

def name_assignment(spectrum):
    """
    assigns a name to the filename in which the spectrum is going to be saved
    """
    if (spectrum == M0_ACES):
        name = "M0-PHOENIX-ACES"
    elif(spectrum == M3_ACES):
        name = "M3-PHOENIX-ACES"
    elif(spectrum == M6_ACES):
        name = "M6-PHOENIX-ACES"
    elif(spectrum == M9_ACES):
        name = "M9-PHOENIX-ACES"
    else:
        print("Name not found!")
        #exit()
        name = "Unknown_name"
    return name

def write_2col(filename, data1, data2):
# Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
 #       f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t%e\t\t%e\n" % (data1[i], data2[i]))

    f.close();

def write_3col(filename, data1, data2, data3):
# Writes data in 3 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
 #       f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t%e\t\t%e\t\t%e\n" % (data1[i], data2[i], data3[i]))

    f.close();
###############################################################################

def list_creator(spectrum, band):
    """
    creates a list of potential lines from a brute-force analysis of the band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)

    line_centers = []
    print(band+" band list:")
    for i in range(2, len(wav_band)-2):
        if(flux_band[i-2]> flux_band[i-1] > flux_band[i] and flux_band[i] < flux_band[i+1] < flux_band[i+2]):
            line_centers.append(wav_band[i])
            print("\t ", wav_band[i]*1.0e4)
    print("In a spectrum with %d points, %d lines were found." % (len(wav_band), len(line_centers)))
    return line_centers
