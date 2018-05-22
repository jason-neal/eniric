#!/usr/bin/env python
# Testing script for nIR analysis
# Run new and old code to test output.S
# Jason Neal
# December 2016
import eniric.IOmodule as io
from eniric.Qcalculator import RVprec_calc

spectrum_name = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

data_rep = "../data/PHOENIX-ACES_spectra/"
results_dir = "../data/results/"
resampled_dir = "../data/resampled/"

spectrum_path = data_rep + spectrum_name
# Some test parameters
bands = ["GAP", "Y", "J", "H", "K"]  # VIS takes 24 minutes to compute.
#                                    # Compared to 24 seconds for the rest???
band = "Y"
R = 100000
vsini = 10
epsilon = 0.6
fwhm_lim = 5
plot = False
num_procs = 4

for band in bands:
    pass
if False:
    print("Starting band", band)
    # New version
    wav, flux = read_spectrum(spectrum_path)
    # raise NotImplementedError("This has broken due to changes in dir structure.")
    wav_band, flux_band, flux_conv = convolution(wav, flux, vsini, R, epsilon,
                                                 fwhm_lim, band=band,
                                                 num_procs=num_procs,
                                                 results_dir="../data/results/unnorm/",
                                                 normalize=False)

    resample_allfiles(results_dir="../data/results/unnorm/",
                      resampled_dir="../data/resampled/unnorm/")

    wav_band_norm, flux_band_norm, flux_conv_norm = convolution(wav, flux, vsini, R,
                                                                epsilon, fwhm_lim, band=band,
                                                                num_procs=num_procs,
                                                                results_dir="../data/results/norm/",
                                                                normalize=True)

    resample_allfiles(results_dir="../data/results/norm/",
                      resampled_dir="../data/resampled/norm/")

    if plot:
        # Plot results together
        plt.plot(wav_band, flux_conv, label='Unnormalized')
        plt.plot(wav_band_norm, flux_conv_norm, label='Normalized (res only)')
        plt.legend(loc=0)
        plt.show()

# Calculate RV Precision
print("Radial velocity Precision, vsini = {0}, R = {1}".format(vsini, R))
for band in bands:
    # Argument unpacking with *
    norm_wav, norm_flux = io.pdread_2col("../data/resampled/norm/"
                                         "Spectrum_M0-PHOENIX-ACES_"
                                         "{0}band_vsini{1}_R{2}k_res3.txt".format(band, int(vsini),
                                                                                  int(R / 1000)), noheader=True)

    unnorm_wav, unnorm_flux = io.pdread_2col("../data/resampled/unnorm/"
                                             "Spectrum_M0-PHOENIX-ACES_"
                                             "{0}band_vsini{1}_R{2}k_res3.txt".format(band, int(vsini),
                                                                                      int(R / 1000)), noheader=True)
    norm_prec = RVprec_calc(norm_wav, norm_flux)

    unnorm_prec = RVprec_calc(unnorm_wav, unnorm_flux)

    # Normalization with Pedro's mysterious factor / ((1.634e4)**2.0)
    norm_factor = 1.634e4 ** 2.0
    pedros_norm_prec = RVprec_calc(norm_wav, norm_flux / norm_factor)
    pedros_unnorm_prec = RVprec_calc(unnorm_wav, unnorm_flux / norm_factor)

    print("Unnormalized RV_Precision                           \t{0} band \t= {1:0.4f}".format(band, unnorm_prec))
    print("Unnormalized RV_Precision with Pedro's /(1.634e4)**2,\t{0} band \t= {1:6.4f}".format(band,
                                                                                                pedros_unnorm_prec))
    print("Normalized RV_Precision                             \t{0} band \t= {1:6.4f}".format(band, norm_prec))
    print("Normalized RV_Precision with Pedro's /(1.634e4)**2,  \t{0} band \t= {1:6.4f}".format(band,
                                                                                                pedros_norm_prec))
