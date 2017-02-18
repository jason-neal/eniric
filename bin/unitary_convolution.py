
"""Convolutions on unity for rotation and resolution.

Perform rotational and resolution convolutions on a vector of ones.
This was used to normalize for the effect of convolution in the original paper.

Was given some updated code of Pedros in which he does this in
nIRanalysis_CONT, saves the result to file and then divides the
results in the resampler.

"""
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from eniric.nIRanalysis import convolution, write_2col, write_3col, \
                               read_spectrum, save_convolution_results

from eniric.nIRanalysis import rotational_convolution, resolution_convolution
import eniric.IOmodule as io
# New code from PEDRO
from eniric.updated_code.nIRanalysis_CONT import convolution_CONT

base_dir = "../data/nIRmodels/"   # relative to script location in bin

create_dat = True
unitary_name = base_dir + "spectrum_of_ones.dat"
if create_dat:
        # Load the Spectrum to get the full wavelength range
    filename = base_dir + "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
    wav = fits.getdata(filename)
    wav = np.asarray(wav, dtype=np.float64)
    wav_header = fits.getheader(filename)

    # Convert to microns
    # wav *= 1e-4 - this is done in read_spectrum
    ones_photons = np.ones_like(wav)
    # read_spectrum converts to photons by multiplying flux by wav.
    # add divide here to counter it
    write_2col(unitary_name, wav, ones_photons / (wav*1e-4))  # convert wav to microns

else:
    wav, ones_photons = read_spectrum(unitary_name)
    # wav, ones_photons = io.pdread_2col(unitary_name, noheader=True)
# Now that we have a spectrum file we can perform the convolutions
name_model = "UNITARY"

bands = ["VIS", "GAP", "Z", "Y", "J", "H", "K"]
vsini = [1.0, 5.0, 10.0]
R = [60000, 80000, 100000]

sampling = ["3"]

bands = ["K"]    # single run with a known band (in current develop version)
vsini = [1.0]
R = [100000]
"""Applying convolution stage of nIRanalysis for all bands vsini and
resolution of paper.

"""
do_convolutions = False
# Applying Rotational and Resolution convolution to unitary spectum
results_dir = base_dir + "../unitary_convolution/"
for band in bands:
    for vel in vsini:
        for Res in R:
                # Provide own path as this model does not convolve to the phoenix
                # name assignment

            filename = ("Spectrum_" +
                        name_model + "_" + band + "band_vsini" + str(vel) +
                        "_R" + str(int(Res/1000)) + "k.txt")

            filename_unnorm = ("Spectrum_" + name_model + "_" + band +
                             "band_vsini" + str(vel) + "_R" +
                             str(int(Res/1000)) + "k_unnormalized.txt")
            if do_convolutions:

                wav, flux = read_spectrum(unitary_name)
                # without normalization
                wav1, flux1, flux1_conv = convolution(wav, flux, vel, Res,
                                                      band=band)
                save_convolution_results(results_dir + filename, wav1, flux1,
                                         flux1_conv)

                # With normalization (should be 1s only)
                wav2, flux2, flux2_conv = convolution(wav, flux, vel, Res,
                                                      band=band, normalize=False)
                save_convolution_results(results_dir + filename_unnorm, wav2,
                                         flux2, flux2_conv)

                cont_wav, cont_flux = convolution_CONT(unitary_name, band, vel,
                                                       Res, plot=False,
                                                       return_only=True)
                write_2col(results_dir + "result_from_convolution_CONT.txt",
                           cont_wav, cont_flux)

            else:
    # Just load data from files
                wav1, flux1, flux1_conv = io.pdread_3col(results_dir + filename, noheader=True)
                wav2, flux2, flux2_conv = io.pdread_3col(results_dir + filename_unnorm, noheader=True)
                cont_wav, cont_flux = io.pdread_2col(results_dir + "result_from_convolution_CONT.txt", noheader=True)

            plt.plot(wav1, flux1_conv, label="Normalized")
            plt.plot(wav2, flux2_conv, "o-", label="Unnormalized")
            plt.plot(cont_wav, cont_flux, label="Pedros CONT")
            plt.legend()
            plt.show()
            # assert np.all(new_wav == cont_wav)
            # assert np.all(new_flux == cont_flux)

"""
Unitary value over full Phoenix Spectrum wavelength range.

Performing rotational and resolution convolutions independantly.
There will be edge effects but these are well outside the normal bands
we are interested in.

"""
# Rotational convolution on ones.

for vel in vsini:
    flux_conv_rot = rotational_convolution(wav, wav, ones_photons, vel,
                                           epsilon=0.6)
    filename = (base_dir + "../unitary_convolution/" + "Spectrum_" +
                name_model + "_vsini" + str(vel) + ".txt")

    write_3col(filename, wav, ones_photons, flux_conv_rot)

# Resolution convolution on ones.
for Res in R:
    flux_conv_res = resolution_convolution(wav, wav, ones_photons, Res,
                                           fwhm_lim=5.0)
    filename = (base_dir + "../unitary_convolution/" + "Spectrum_" +
                name_model + "_R" + str(int(Res/1000)) + "k.txt")

    write_3col(filename, wav, ones_photons, flux_conv_res)
