"""
Created on Fri Feb  6 15:42:03 2015

@author: pfigueira

Updated for eniric/python3 - Janurary 2017
@author: Jason Neal
"""

import numpy as np
from sys import exit
import matplotlib.pyplot as plt

# to remove labels in one tick
from matplotlib.ticker import MaxNLocator

import eniric.IOmodule as IOmodule
import eniric.Qcalculator as Qcalculator

from eniric.utilities import band_selector

from matplotlib import rc
# set stuff for latex usage
rc('text', usetex=True)

atmmodel = "../data/atmmodel/Average_TAPAS_2014.txt"
resampled_dir = "resampled_cont/"
resampled_dir_OLD = "resampled/"

spectral_types = ["M0", "M3", "M6", "M9"]
bands = ["Z", "Y", "J", "H", "K"]
vsini = ["1.0", "5.0", "10.0"]
R = ["60k", "80k", "100k"]
sampling = ["3"]


def read_contfit():
    """
    function that reads continuum fitting
    """
    M0_contfit = "PHOENIX_ACES_spectra/M0_nIRcont.txt"
    M3_contfit = "PHOENIX_ACES_spectra/M3_nIRcont.txt"
    M6_contfit = "PHOENIX_ACES_spectra/M6_nIRcont.txt"
    M9_contfit = "PHOENIX_ACES_spectra/M9_nIRcont.txt"

    wav_M0, flux_M0 = IOmodule.pdread_2col(M0_contfit)
    wav_M3, flux_M3 = IOmodule.pdread_2col(M3_contfit)
    wav_M6, flux_M6 = IOmodule.pdread_2col(M6_contfit)
    wav_M9, flux_M9 = IOmodule.pdread_2col(M9_contfit)

    wav_M0 = np.array(wav_M0, dtype="float64")*1.0e-4  # conversion to microns
    flux_M0 = np.array(flux_M0, dtype="float64")*wav_M0
    wav_M3 = np.array(wav_M3, dtype="float64")*1.0e-4  # conversion to microns
    flux_M3 = np.array(flux_M3, dtype="float64")*wav_M3
    wav_M6 = np.array(wav_M6, dtype="float64")*1.0e-4  # conversion to microns
    flux_M6 = np.array(flux_M6, dtype="float64")*wav_M6
    wav_M9 = np.array(wav_M9, dtype="float64")*1.0e-4  # conversion to microns
    flux_M9 = np.array(flux_M9, dtype="float64")*wav_M9

    return [wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9]


def read_nIRspectra():
    """
    function that reads nIR spectra
    """
    M0_contfit = "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M3_contfit = "PHOENIX_ACES_spectra/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M6_contfit = "PHOENIX_ACES_spectra/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M9_contfit = "PHOENIX_ACES_spectra/lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"

    print("Reading PHOENIX original spectrum...")
    wav_M0, flux_M0 = IOmodule.pdread_2col(M0_contfit)
    wav_M3, flux_M3 = IOmodule.pdread_2col(M3_contfit)
    wav_M6, flux_M6 = IOmodule.pdread_2col(M6_contfit)
    wav_M9, flux_M9 = IOmodule.pdread_2col(M9_contfit)
    print("Done.")

    wav_M0 = np.array(wav_M0, dtype="float64")*1.0e-4  # conversion to microns
    flux_M0 = np.array(flux_M0, dtype="float64")*wav_M0
    wav_M3 = np.array(wav_M3, dtype="float64")*1.0e-4  # conversion to microns
    flux_M3 = np.array(flux_M3, dtype="float64")*wav_M3
    wav_M6 = np.array(wav_M6, dtype="float64")*1.0e-4  # conversion to microns
    flux_M6 = np.array(flux_M6, dtype="float64")*wav_M6
    wav_M9 = np.array(wav_M9, dtype="float64")*1.0e-4  # conversion to microns
    flux_M9 = np.array(flux_M9, dtype="float64")*wav_M9

    wav_M0 = wav_M0[1000:-1000]
    wav_M3 = wav_M3[1000:-1000]
    wav_M6 = wav_M6[1000:-1000]
    wav_M9 = wav_M9[1000:-1000]

    flux_M0 = moving_average(flux_M0, 2000)[1000:-1000]
    flux_M3 = moving_average(flux_M3, 2000)[1000:-1000]
    flux_M6 = moving_average(flux_M6, 2000)[1000:-1000]
    flux_M9 = moving_average(flux_M9, 2000)[1000:-1000]

    return [wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9]


def ratios_calc(wav_bin):
    """
    funtion that calculates ratios as set in the continuum to apply to the spectrum
    """
    wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9 = read_contfit()

    index_M0_l = np.searchsorted(wav_M0, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M3_l = np.searchsorted(wav_M3, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M6_l = np.searchsorted(wav_M6, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M9_l = np.searchsorted(wav_M9, [0.83, 1.0, 1.17, 1.5, 2.07])

    index_M0_r = np.searchsorted(wav_M0, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M3_r = np.searchsorted(wav_M3, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M6_r = np.searchsorted(wav_M6, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M9_r = np.searchsorted(wav_M9, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])

    return [flux_bin(flux_M0, index_M0_l, index_M0_r), flux_bin(flux_M3, index_M3_l, index_M3_r), flux_bin(flux_M6, index_M6_l, index_M6_r), flux_bin(flux_M9, index_M9_l, index_M9_r)]


def flux_bin(flux, index_left, index_right):
    fluxes = []
    for ind_l, ind_r in zip(index_left, index_right):
        fluxes.append(np.average(flux[ind_l: ind_r]))
    return fluxes


def prepare_atmopshere():
    """ Read in atmopheric model and prepare. """
    wav_atm, flux_atm, std_flux_atm, mask_atm = IOmodule.pdread_4col(atmmodel)
    # pandas lready returns numpy arrays
    wav_atm = wav_atm / 1000.0  # conversion from nanometers to micrometers
    mask_atm = np.array(mask_atm, dtype=bool)
    return wav_atm, flux_atm, std_flux_atm, mask_atm


def Barycenter_shift(wav_atm, mask_atm, offset_RV=0.0):
    """ Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentic motion of the earth.
    """

    pixels_start = np.sum(mask_atm)
    pixels_total = len(mask_atm)

    mask_atm_30kms = []
    for value in zip(wav_atm, mask_atm):
        if (value[1] is False and offset_RV == 666.0):    # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(value[1])

        else:

            delta_lambda = value[0] * 3.0e4/Qcalculator.c
            starting_lambda = value[0] * offset_RV*1.0e3/Qcalculator.c
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda+value[0]-delta_lambda, starting_lambda+value[0]+delta_lambda])
            indexes_30kmslice = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in indexes_30kmslice]
            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]], dtype=bool)    # selecting only the slice in question

            # if(False in mask_atm_30kmslice):
            #    mask_atm_30kms.append(False)
            # else:
            #    mask_atm_30kms.append(True)

            mask_atm_30kmslice_reversed = [not i for i in mask_atm_30kmslice]

            clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice_reversed))[0]+1)[::2]

            tester = True
            for block in clump:
                if len(clump) >= 3:
                    tester = False
                    break

            mask_atm_30kms.append(tester)

    mask_atm = np.array(mask_atm_30kms, dtype=bool)
    pixels_end = np.sum(mask_atm)
    print(("Barycentric impact masks out {:04.1}\% more of the atmospheric"
          " spectrum").format((pixels_end-pixels_start)/pixels_total))
    return mask_atm


def calculate_prec(plot_atm=False, plot_ste=False, plot_flux=True, paper_plots=True, offset_RV=0.0):

    print("Reading atmospheric model...")
    wav_atm, flux_atm, std_flux_atm, mask_atm = prepare_atmopshere()
    print(("There were {:d} unmasked pixels out of {:d}., or {:.1%}."
          "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)))
    print("The model ranges from {:4.2f} to {:4.2f} micron.".format(wav_atm[0], wav_atm[-1]))
    print("Done.")

    print("Calculating impact of Barycentric movement on mask...")
    mask_atm = Barycenter_shift(wav_atm, mask_atm, offset_RV=offset_RV)  # Extend masked regions
    print(("There were {:d} unmasked pixels out of {:d}, or {:.1%}."
           "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) /
                      len(mask_atm)))


    # calculating the number of pixels inside the mask
    wav_Z, mask_Z = band_selector(wav_atm, mask_atm, "Z")
    wav_Y, mask_Y = band_selector(wav_atm, mask_atm, "Y")
    wav_J, mask_J = band_selector(wav_atm, mask_atm, "J")
    wav_H, mask_H = band_selector(wav_atm, mask_atm, "H")
    wav_K, mask_K = band_selector(wav_atm, mask_atm, "K")

    bands_masked = np.concatenate((mask_Z, mask_Y, mask_J, mask_H, mask_K))

    print(("Inside the bands, there were {:.0f} unmasked pixels out of {:d}"
           ", or {:.1%}.").format(np.sum(bands_masked), len(bands_masked),
            np.sum(bands_masked) / len(bands_masked)))

    if(plot_atm):

        # identify non-masked pixels
        selected_transmission = wav_atm[mask_atm]
        dummy_vector = [1.0 for value in selected_transmission]
        """
        plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"Transmission [ ] ")
        plt.plot(wav_atm, flux_atm, color='k')
        plt.vlines(selected_transmission, 0.8, 1.0, colors="b")
        plt.xlim(wav_atm[0], wav_atm[-1])
        plt.ylim(0.0, 1.0)
        # plt.legend(loc='best')
        plt.show()
        plt.close()
        """
        # do the same per band
        wav_atm_Z, flux_atm_Z = band_selector(wav_atm, flux_atm, "Z")
        wav_atm_Y, flux_atm_Y = band_selector(wav_atm, flux_atm, "Y")
        wav_atm_J, flux_atm_J = band_selector(wav_atm, flux_atm, "J")
        wav_atm_H, flux_atm_H = band_selector(wav_atm, flux_atm, "H")
        wav_atm_K, flux_atm_K = band_selector(wav_atm, flux_atm, "K")

        wav_mask_Z, dummy = band_selector(selected_transmission, dummy_vector, "Z")
        wav_mask_Y, dummy = band_selector(selected_transmission, dummy_vector, "Y")
        wav_mask_J, dummy = band_selector(selected_transmission, dummy_vector, "J")
        wav_mask_H, dummy = band_selector(selected_transmission, dummy_vector, "H")
        wav_mask_K, dummy = band_selector(selected_transmission, dummy_vector, "K")

        print("Z:", len(wav_mask_Z)/float(len(wav_atm_Z)), np.average(flux_atm_Z), np.median(flux_atm_Z))
        print("Y:", len(wav_mask_Y)/float(len(wav_atm_Y)), np.average(flux_atm_Y), np.median(flux_atm_Y))
        print("J:", len(wav_mask_J)/float(len(wav_atm_J)), np.average(flux_atm_J), np.median(flux_atm_J))
        print("H:", len(wav_mask_H)/float(len(wav_atm_H)), np.average(flux_atm_H), np.median(flux_atm_H))
        print("K:", len(wav_mask_K)/float(len(wav_atm_K)), np.average(flux_atm_K), np.median(flux_atm_K))

        f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=False, sharey=False)

        ax1.plot(wav_atm_Z, flux_atm_Z, color='k')
        # ax1.vlines(selected_transmission, 0.0, 0.2, colors="b")
        ax1.set_xlim([wav_atm_Z[0], wav_atm_Z[-1]])
        ax1.yaxis.set_major_locator(MaxNLocator(nbins=len(ax1.get_yticklabels()), prune='lower'))
        ax1.text(0.9, 0.8, ' Z band ', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12, color='black', bbox=dict(facecolor='white'))

        ax2.plot(wav_atm_Y, flux_atm_Y, color='k')
        # ax2.vlines(selected_transmission, 0.0, 0.2, colors="b")
        ax2.set_xlim([wav_atm_Y[0], wav_atm_Y[-1]])
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='lower'))
        ax2.text(0.9, 0.8, ' Y band ', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12, color='black', bbox=dict(facecolor='white'))

        ax3.plot(wav_atm_J, flux_atm_J, color='k')
        # ax3.vlines(selected_transmission, 0.0, 0.2, colors="b")
        ax3.set_xlim([wav_atm_J[0], wav_atm_J[-1]])
        ax3.yaxis.set_major_locator(MaxNLocator(nbins=len(ax3.get_yticklabels()), prune='upper'))
        ax3.text(0.9, 0.8, ' J band', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12, color='black', bbox=dict(facecolor='white'))

        ax4.plot(wav_atm_H, flux_atm_H, color='k')
        # ax4.vlines(selected_transmission, 0.0, 0.2, colors="b")
        ax4.set_xlim([wav_atm_H[0], wav_atm_H[-1]])
        ax4.yaxis.set_major_locator(MaxNLocator(nbins=len(ax4.get_yticklabels()), prune='lower'))
        ax4.text(0.9, 0.8, ' H band ', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12, color='black', bbox=dict(facecolor='white'))

        ax5.plot(wav_atm_K, flux_atm_K, color='k')
        # ax5.vlines(selected_transmission, 0.0, 0.2, colors="b")
        ax5.set_xlim([wav_atm_K[0], wav_atm_K[-1]])
        ax5.yaxis.set_major_locator(MaxNLocator(nbins=len(ax5.get_yticklabels()), prune='upper'))
        ax5.text(0.9, 0.8, ' K band ', horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes, size=12, color='black', bbox=dict(facecolor='white'))

        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r'Atmospheric Transmission [ ]', ha='center', va='center', rotation='vertical')
        plt.show()
        # f.savefig('AtmosphericTransmission.pdf', format='pdf')

        plt.close()
        exit()

    # theoretical ratios calculation
    wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9 = read_nIRspectra()

    results = {}    # creating empty dictionary for the results
    wav_plot_M0 = []   # creating empty lists for the plots
    flux_plot_M0 = []
    wav_plot_M3 = []
    flux_plot_M3 = []
    wav_plot_M6 = []
    flux_plot_M6 = []
    wav_plot_M9 = []
    flux_plot_M9 = []
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        # print("Working on "+file_to_read+".")
                        wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir+file_to_read)
                        # removing boundary effects
                        wav_stellar = wav_stellar[2:-2]
                        flux_stellar = flux_stellar[2:-2]

                        id_string = star+"-"+band+"-"+vel+"-"+resolution    # sample was left aside because only one value existed

                        # gettin the wav, flux and mask from the atm model closes to the stellar

                        index_atm = np.searchsorted(wav_atm, wav_stellar)
                        # replace indexes outside the array, at the very end, by the value at the very end
                        index_atm = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in index_atm]

                        wav_atm_selected = wav_atm[index_atm]
                        flux_atm_selected = flux_atm[index_atm]
                        mask_atm_selected = mask_atm[index_atm]

                        if("M0" in id_string):
                            norm_constant = 1607

                        elif("M3" in id_string):
                            norm_constant = 1373

                        elif("M6" in id_string):
                            if("1.0" in id_string):
                                norm_constant = 933
                            elif("5.0" in id_string):
                                norm_constant = 967
                            else:
                                norm_constant = 989

                        elif("M9" in id_string):
                            if("1.0" in id_string):
                                norm_constant = 810
                            elif("5.0" in id_string):
                                norm_constant = 853
                            else:
                                norm_constant = 879
                        else:
                            print("Constant not defined. Aborting...")
                            exit()

                        flux_stellar = flux_stellar / ((norm_constant/100.0)**2.0)

                        if(id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k", "M6-J-1.0-100k", "M9-J-1.0-100k"]):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the %s reference model was of %4.2f." % (id_string, SN_estimate))
                        elif("J" in id_string):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the %s non-reference model was of %4.2f." % (id_string, SN_estimate))
                        if(plot_ste or plot_ste == id_string):
                            # Plot the stellar spectr as considered

                            selected_transmission_stellar = wav_atm_selected[mask_atm_selected]

                            plt.figure(1)
                            plt.xlabel(r"wavelength [$\mu$m])")
                            plt.ylabel(r"Flux_stellar [ ] ")
                            plt.plot(wav_stellar, flux_stellar, color='k')
                            plt.vlines(selected_transmission_stellar, np.min(flux_stellar), 0.3*np.max(flux_stellar), colors="b")
                            plt.xlim(wav_stellar[0], wav_stellar[-1])
                            plt.ylim(np.min(flux_stellar)-0.1*(np.max(flux_stellar)-np.min(flux_stellar)), np.max(flux_stellar)+0.1*(np.max(flux_stellar)-np.min(flux_stellar)))
                            # plt.legend(loc='best')
                            plt.show()
                            plt.close()

                        if(plot_flux and id_string in ["M0-Z-1.0-100k", "M0-Y-1.0-100k", "M0-J-1.0-100k", "M0-H-1.0-100k", "M0-K-1.0-100k"]):
                            wav_plot_M0.append(wav_stellar)
                            flux_plot_M0.append(flux_stellar)
                        if(plot_flux and id_string in ["M3-Z-1.0-100k", "M3-Y-1.0-100k", "M3-J-1.0-100k", "M3-H-1.0-100k", "M3-K-1.0-100k"]):
                            wav_plot_M3.append(wav_stellar)
                            flux_plot_M3.append(flux_stellar)
                        if(plot_flux and id_string in ["M6-Z-1.0-100k", "M6-Y-1.0-100k", "M6-J-1.0-100k", "M6-H-1.0-100k", "M6-K-1.0-100k"]):
                            wav_plot_M6.append(wav_stellar)
                            flux_plot_M6.append(flux_stellar)
                        if(plot_flux and id_string in ["M9-Z-1.0-100k", "M9-Y-1.0-100k", "M9-J-1.0-100k", "M9-H-1.0-100k", "M9-K-1.0-100k"]):
                            wav_plot_M9.append(wav_stellar)
                            flux_plot_M9.append(flux_stellar)

                        # precision given by the first method:
                        print("Performing analysis for: ", id_string)
                        prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

                        # precision as given by the second_method
                        """
                        Example Joao
                        a = np.array([1, 5, 6, 8, 16, 34, 5, 7, 10, 83, 12, 6, 17, 18])
                        b = np.array([1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1], dtype=bool)

                        # this will give you a list of numpy arrays
                        c = np.array_split(a, np.where(np.diff(b))[0]+1)[::2]

                        # this will give you a list of lists
                        d = [list(cc) for cc in c]
                        print(d)
                        >>> [[1, 5], [16, 34, 5], [83, 12], [17, 18]]
                        """

                        wav_stellar_chunks_unformated = np.array_split(wav_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        wav_stellar_chunks = [list(chunk) for chunk in wav_stellar_chunks_unformated]

                        """
                        # test section
                        print("check that lengths are the same", len(wav_stellar), len(mask_atm_selected))
                        print("size of spectra %d vs number of chunks %d" % (len(wav_stellar), len(wav_stellar_chunks)))
                        print("number of true elements in all chunks: %d" % (len(mask_atm_selected[mask_atm_selected])))
                        """

                        flux_stellar_chunks_unformated = np.array_split(flux_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        flux_stellar_chunks = [list(chunk) for chunk in flux_stellar_chunks_unformated]

                        """
                        # histogram checking
                        lengths = [len(chunk) for chunk in flux_stellar_chunks_unformated]
                        n, bins, patches = plt.hist(lengths, 500, range=[0.5, 500.5], histtype='stepfilled')
                        plt.title(id_string)
                        plt.show()
                        """

                        prec_2 = Qcalculator.RVprec_calc_chunks(wav_stellar_chunks, flux_stellar_chunks)

                        # precision as given by the third_method
                        prec_3 = Qcalculator.RV_prec_calc_Trans(wav_stellar, flux_stellar, flux_atm_selected)

                        # adding result to the dictionary
                        results[id_string] = [prec_1, prec_2, prec_3]
    if(plot_flux):

        # print the plot for flux comparison
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        ax1.plot(wav_plot_M0[0], np.array(flux_plot_M0[0]), color='0.1')
        ax1.plot(wav_plot_M0[1], np.array(flux_plot_M0[1]), color='0.1')
        ax1.plot(wav_plot_M0[2], np.array(flux_plot_M0[2]), color='0.1')
        ax1.plot(wav_plot_M0[3], np.array(flux_plot_M0[3]), color='0.1')
        ax1.plot(wav_plot_M0[4], np.array(flux_plot_M0[4]), color='0.1')
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)
        ax1.plot(wav_M0, (flux_M0/flux_M0[0])*flux_plot_M0[0][1000], color='0.1', linestyle='--')

        ax2.plot(wav_plot_M3[0], np.array(flux_plot_M3[0]), color='0.3')
        ax2.plot(wav_plot_M3[1], np.array(flux_plot_M3[1]), color='0.3')
        ax2.plot(wav_plot_M3[2], np.array(flux_plot_M3[2]), color='0.3')
        ax2.plot(wav_plot_M3[3], np.array(flux_plot_M3[3]), color='0.3')
        ax2.plot(wav_plot_M3[4], np.array(flux_plot_M3[4]), color='0.3')
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.plot(wav_M3, (flux_M3/flux_M3[0])*flux_plot_M3[0][1000], color='0.1', linestyle='--')

        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))   # remove upper element from yaxis
        # ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[0], np.array(flux_plot_M6[0]), color='0.4')
        ax3.plot(wav_plot_M6[1], np.array(flux_plot_M6[1]), color='0.4')
        ax3.plot(wav_plot_M6[2], np.array(flux_plot_M6[2]), color='0.4')
        ax3.plot(wav_plot_M6[3], np.array(flux_plot_M6[3]), color='0.4')
        ax3.plot(wav_plot_M6[4], np.array(flux_plot_M6[4]), color='0.4')
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.plot(wav_M6, (flux_M6/flux_M6[0])*flux_plot_M6[0][1000]*1.2, color='0.1', linestyle='--')

        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[0], np.array(flux_plot_M9[0]), color='0.6')
        ax4.plot(wav_plot_M9[1], np.array(flux_plot_M9[1]), color='0.6')
        ax4.plot(wav_plot_M9[2], np.array(flux_plot_M9[2]), color='0.6')
        ax4.plot(wav_plot_M9[3], np.array(flux_plot_M9[3]), color='0.6')
        ax4.plot(wav_plot_M9[4], np.array(flux_plot_M9[4]), color='0.6')
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.plot(wav_M9, (flux_M9/flux_M9[0])*flux_plot_M9[0][1000]*1.4, color='0.1', linestyle='--')

        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.text(0.06, 0.5, r'Flux [ ]', ha='center', va='center', rotation='vertical')
        plt.xlabel(r"wavelength [$\mu$m]")

        plt.show()
        plt.close()

        # print the Z band spectra
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r"Flux in Z band [ ]", ha='center', va='center', rotation='vertical')

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.set_xlim(np.min(wav_plot_M0[0]), np.max(wav_plot_M0[0]))
        ax1.set_ylim(0.0, 5.0e3)

        ax1.plot(wav_plot_M0[0], np.array(flux_plot_M0[0]), color='0.1', label="M0")
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)
        ax1.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

        ax2.plot(wav_plot_M3[0], np.array(flux_plot_M3[0]), color='0.3', label="M3")
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[0], np.array(flux_plot_M6[0]), color='0.4', label="M6")
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis
        ax3.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

        ax4.plot(wav_plot_M9[0], np.array(flux_plot_M9[0]), color='0.6', label="M9")
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis
        ax4.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        plt.close()

        # print the Y band spectra
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r"Flux in Y band [ ]", ha='center', va='center', rotation='vertical')

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        ax1.plot(wav_plot_M0[1], np.array(flux_plot_M0[1]), color='0.1', label="M0")
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

        ax2.plot(wav_plot_M3[1], np.array(flux_plot_M3[1]), color='0.3', label="M3")
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[1], np.array(flux_plot_M6[1]), color='0.4', label="M6")
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[1], np.array(flux_plot_M9[1]), color='0.6', label="M9")
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        plt.close()

        # print the J band spectra
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r"Flux in J band [ ]", ha='center', va='center', rotation='vertical')

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        ax1.set_xlim(np.min(wav_plot_M0[2]), np.max(wav_plot_M0[2]))

        ax1.plot(wav_plot_M0[2], np.array(flux_plot_M0[2]), color='0.1', label="M0")
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

        ax2.plot(wav_plot_M3[2], np.array(flux_plot_M3[2]), color='0.3', label="M3")
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[2], np.array(flux_plot_M6[2]), color='0.4', label="M6")
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[2], np.array(flux_plot_M9[2]), color='0.6', label="M9")
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        plt.close()

        # print the H band spectra
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r"Flux in H band [ ]", ha='center', va='center', rotation='vertical')

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        ax1.plot(wav_plot_M0[3], np.array(flux_plot_M0[3]), color='0.1', label="M0")
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

        ax2.plot(wav_plot_M3[3], np.array(flux_plot_M3[3]), color='0.3', label="M3")
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[3], np.array(flux_plot_M6[3]), color='0.4', label="M6")
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[3], np.array(flux_plot_M9[3]), color='0.6', label="M9")
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        plt.show()
        plt.close()

    # print the K band spectra
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
        plt.xlabel(r"wavelength [$\mu$m]")
        f.text(0.06, 0.5, r"Flux in K band [ ]", ha='center', va='center', rotation='vertical')

        ax1.set_xlim([np.min(wav_plot_M0[4]), np.max(wav_plot_M0[4])])
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        ax1.plot(wav_plot_M0[4], np.array(flux_plot_M0[4]), color='0.1', label="M0")
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

        ax2.plot(wav_plot_M3[4], np.array(flux_plot_M3[4]), color='0.3', label="M3")
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[4], np.array(flux_plot_M6[4]), color='0.4', label="M6")
        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[4], np.array(flux_plot_M9[4]), color='0.6', label="M9")
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

        plt.show()
        plt.close()

    if(paper_plots):
        # print the paper plots
        """
        In every plot we will assume sample=3 and vsini=1
        y = RVprec between prec3 and prec2
        x = different bands

        different panels for different spectral types

        different colors will represent the different resolutions
        """

        print("Results for vsini of 1.0 km/s")
        # preparation of data: plot1
        y1_60k = [results["M0-Z-1.0-60k"], results["M0-Y-1.0-60k"], results["M0-J-1.0-60k"], results["M0-H-1.0-60k"], results["M0-K-1.0-60k"]]
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]

        y1_60k_lim = [y[0] for y in y1_60k]

        y1_80k = [results["M0-Z-1.0-80k"], results["M0-Y-1.0-80k"], results["M0-J-1.0-80k"], results["M0-H-1.0-80k"], results["M0-K-1.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]

        y1_80k_lim = [y[0] for y in y1_80k]

        y1_100k = [results["M0-Z-1.0-100k"], results["M0-Y-1.0-100k"], results["M0-J-1.0-100k"], results["M0-H-1.0-100k"], results["M0-K-1.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]

        y1_100k_lim = [y[0] for y in y1_100k]

        # preparation of data: plot2
        y2_60k = [results["M3-Z-1.0-60k"], results["M3-Y-1.0-60k"], results["M3-J-1.0-60k"], results["M3-H-1.0-60k"], results["M3-K-1.0-60k"]]
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]

        y2_60k_lim = [y[0] for y in y2_60k]

        y2_80k = [results["M3-Z-1.0-80k"], results["M3-Y-1.0-80k"], results["M3-J-1.0-80k"], results["M3-H-1.0-80k"], results["M3-K-1.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]

        y2_80k_lim = [y[0] for y in y2_80k]

        y2_100k = [results["M3-Z-1.0-100k"], results["M3-Y-1.0-100k"], results["M3-J-1.0-100k"], results["M3-H-1.0-100k"], results["M3-K-1.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k]

        y2_100k_lim = [y[0] for y in y2_100k]

        # preparation of data: plot3
        y3_60k = [results["M6-Z-1.0-60k"], results["M6-Y-1.0-60k"], results["M6-J-1.0-60k"], results["M6-H-1.0-60k"], results["M6-K-1.0-60k"]]
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]

        y3_60k_lim = [y[0] for y in y3_60k]

        y3_80k = [results["M6-Z-1.0-80k"], results["M6-Y-1.0-80k"], results["M6-J-1.0-80k"], results["M6-H-1.0-80k"], results["M6-K-1.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]

        y3_80k_lim = [y[0] for y in y3_80k]

        y3_100k = [results["M6-Z-1.0-100k"], results["M6-Y-1.0-100k"], results["M6-J-1.0-100k"], results["M6-H-1.0-100k"], results["M6-K-1.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k]

        y3_100k_lim = [y[0] for y in y3_100k]

        # preparation of data: plot4
        y4_60k = [results["M9-Z-1.0-60k"], results["M9-Y-1.0-60k"], results["M9-J-1.0-60k"], results["M9-H-1.0-60k"], results["M9-K-1.0-60k"]]
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]

        y4_60k_lim = [y[0] for y in y4_60k]

        y4_80k = [results["M9-Z-1.0-80k"], results["M9-Y-1.0-80k"], results["M9-J-1.0-80k"], results["M9-H-1.0-80k"], results["M9-K-1.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]

        y4_80k_lim = [y[0] for y in y4_80k]

        y4_100k = [results["M9-Z-1.0-100k"], results["M9-Y-1.0-100k"], results["M9-J-1.0-100k"], results["M9-H-1.0-100k"], results["M9-K-1.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]

        y4_100k_lim = [y[0] for y in y4_100k]

        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])

        # data for correction plot
        y1_60k_vsini1 = (np.array(y1_60k_top) - np.array(y1_60k_bottom)) / np.array(y1_60k_bottom)*100.0
        y1_80k_vsini1 = (np.array(y1_80k_top) - np.array(y1_80k_bottom)) / np.array(y1_80k_bottom)*100.0
        y1_100k_vsini1 = (np.array(y1_100k_top) - np.array(y1_100k_bottom)) / np.array(y1_100k_bottom)*100.0

        y2_60k_vsini1 = (np.array(y2_60k_top) - np.array(y2_60k_bottom)) / np.array(y2_60k_bottom)*100.0
        y2_80k_vsini1 = (np.array(y2_80k_top) - np.array(y2_80k_bottom)) / np.array(y2_80k_bottom)*100.0
        y2_100k_vsini1 = (np.array(y2_100k_top) - np.array(y2_100k_bottom)) / np.array(y2_100k_bottom)*100.0

        y3_60k_vsini1 = (np.array(y3_60k_top) - np.array(y3_60k_bottom)) / np.array(y3_60k_bottom)*100.0
        y3_80k_vsini1 = (np.array(y3_80k_top) - np.array(y3_80k_bottom)) / np.array(y3_80k_bottom)*100.0
        y3_100k_vsini1 = (np.array(y3_100k_top) - np.array(y3_100k_bottom)) / np.array(y3_100k_bottom)*100.0

        y4_60k_vsini1 = (np.array(y4_60k_top) - np.array(y4_60k_bottom)) / np.array(y4_60k_bottom)*100.0
        y4_80k_vsini1 = (np.array(y4_80k_top) - np.array(y4_80k_bottom)) / np.array(y4_80k_bottom)*100.0
        y4_100k_vsini1 = (np.array(y4_100k_top) - np.array(y4_100k_bottom)) / np.array(y4_100k_bottom)*100.0

        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)

        ax1.fill_between(range(1, len(bands) + 1), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)

        ax1.plot(range(1, len(bands) + 1), y1_60k_lim, color="b", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_80k_lim, color="g", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_100k_lim, color="r", linestyle="--")

        ax1.scatter(range(1, len(bands) + 1), y1_60k_bottom, marker='^', color="b", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_60k_top, marker='o', color="b", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_80k_bottom, marker='^', color="g", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_top, marker='o', color="g", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_100k_bottom, marker='^', color="r", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_top, marker='o', color="r", alpha=0.4)

        ax1.text(1.0, positiony_max, "M0", size=14)

        # ticks and labels
        # ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(bands)+0.5)
        ax1.set_xticks(range(1, len(bands) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax1.tick_params(axis='both', which='major', labelsize=15)

        ax2 = fig.add_subplot(222)

        ax2.fill_between(range(1, len(bands) + 1), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.plot(range(1, len(bands) + 1), y2_60k_lim, color="b", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_80k_lim, color="g", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_100k_lim, color="r", linestyle="--")

        ax2.scatter(range(1, len(bands) + 1), y2_60k_bottom, marker='^', color="b", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_60k_top, marker='o', color="b", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_80k_bottom, marker='^', color="g", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_top, marker='o', color="g", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_100k_bottom, marker='^', color="r", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_top, marker='o', color="r", alpha=0.4)

        ax2.text(1.0, positiony_max, "M3", size=14)

        # ticks and labels
        ax2.set_xlim(0.5, len(bands)+0.5)
        ax2.set_xticks(range(1, len(bands) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax2.tick_params(axis='both', which='major', labelsize=15)

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1, len(bands) + 1), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.plot(range(1, len(bands) + 1), y3_60k_lim, color="b", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_80k_lim, color="g", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_100k_lim, color="r", linestyle="--")

        ax3.scatter(range(1, len(bands) + 1), y3_60k_bottom, marker='^', color="b", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_60k_top, marker='o', color="b", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_80k_bottom, marker='^', color="g", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_top, marker='o', color="g", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_100k_bottom, marker='^', color="r", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_top, marker='o', color="r", alpha=0.4)

        ax3.text(1.0, positiony_max, "M6", size=14)

        # ticks and labels
        # ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands)+0.5)
        ax3.set_xticks(range(1, len(bands) + 1))
        ax3.set_xticklabels(bands)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3.tick_params(axis='both', which='major', labelsize=15)

        ax4 = fig.add_subplot(224)

        ax4.fill_between(range(1, len(bands) + 1), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.plot(range(1, len(bands) + 1), y4_60k_lim, color="b", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_80k_lim, color="g", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_100k_lim, color="r", linestyle="--")

        ax4.scatter(range(1, len(bands) + 1), y4_60k_bottom, marker='^', color="b", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_60k_top, marker='o', color="b", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_80k_bottom, marker='^', color="g", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_top, marker='o', color="g", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_100k_bottom, marker='^', color="r", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_top, marker='o', color="r", alpha=0.4)

        ax4.text(1.0, positiony_max, "M9", size=14)

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands)+0.5)
        ax4.set_xticks(range(1, len(bands) + 1))
        ax4.set_xticklabels(bands)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4.tick_params(axis='both', which='major', labelsize=15)

        fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        """
        same plot for total precision
        """

        positiony_max = np.max([np.max(RV_cumulative(y1_60k_top)), np.max(RV_cumulative(y1_60k_bottom)), np.max(RV_cumulative(y1_80k_top)), np.max(RV_cumulative(y1_80k_bottom)), np.max(RV_cumulative(y1_100k_top)), np.max(RV_cumulative(y1_100k_bottom)),\
                                   np.max(RV_cumulative(y2_60k_top)), np.max(RV_cumulative(y2_60k_bottom)), np.max(RV_cumulative(y2_80k_top)), np.max(RV_cumulative(y2_80k_bottom)), np.max(RV_cumulative(y2_100k_top)), np.max(RV_cumulative(y2_100k_bottom)),\
                                   np.max(RV_cumulative(y3_60k_top)), np.max(RV_cumulative(y3_60k_bottom)), np.max(RV_cumulative(y3_80k_top)), np.max(RV_cumulative(y3_80k_bottom)), np.max(RV_cumulative(y3_100k_top)), np.max(RV_cumulative(y3_100k_bottom)),\
                                   np.max(RV_cumulative(y4_60k_top)), np.max(RV_cumulative(y4_60k_bottom)), np.max(RV_cumulative(y4_80k_top)), np.max(RV_cumulative(y4_80k_bottom)), np.max(RV_cumulative(y4_100k_top)), np.max(RV_cumulative(y4_100k_bottom))])

        positiony_min = np.min([np.min(RV_cumulative(y1_60k_top)), np.min(RV_cumulative(y1_60k_bottom)), np.min(RV_cumulative(y1_80k_top)), np.min(RV_cumulative(y1_80k_bottom)), np.min(RV_cumulative(y1_100k_top)), np.min(RV_cumulative(y1_100k_bottom)),\
                                   np.min(RV_cumulative(y2_60k_top)), np.min(RV_cumulative(y2_60k_bottom)), np.min(RV_cumulative(y2_80k_top)), np.min(RV_cumulative(y2_80k_bottom)), np.min(RV_cumulative(y2_100k_top)), np.min(RV_cumulative(y2_100k_bottom)),\
                                   np.min(RV_cumulative(y3_60k_top)), np.min(RV_cumulative(y3_60k_bottom)), np.min(RV_cumulative(y3_80k_top)), np.min(RV_cumulative(y3_80k_bottom)), np.min(RV_cumulative(y3_100k_top)), np.min(RV_cumulative(y3_100k_bottom)),\
                                   np.min(RV_cumulative(y4_60k_top)), np.min(RV_cumulative(y4_60k_bottom)), np.min(RV_cumulative(y4_80k_top)), np.min(RV_cumulative(y4_80k_bottom)), np.min(RV_cumulative(y4_100k_top)), np.min(RV_cumulative(y4_100k_bottom))])

        bands_total = ["ZY", "ZYJ", "ZYJH", "ZYJHK"]
        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)

        ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_bottom), RV_cumulative(y1_60k_top), color="b", alpha=0.2)
        ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_bottom), RV_cumulative(y1_80k_top), color="g", alpha=0.2)
        ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_bottom), RV_cumulative(y1_100k_top), color="r", alpha=0.2)

        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_bottom), marker='^', color="b", alpha=0.4)
        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_top), marker='o', color="b", alpha=0.4)

        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_bottom), marker='^', color="g", alpha=0.4)
        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_top), marker='o', color="g", alpha=0.4)

        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_bottom), marker='^', color="r", alpha=0.4)
        ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_top), marker='o', color="r", alpha=0.4)

        ax1.text(1.0, positiony_max, "M0", size=14)

        # ticks and labels
        # ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(bands_total)+0.5)
        ax1.set_xticks(range(1, len(bands_total) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax1.tick_params(axis='both', which='major', labelsize=15)

        ax2 = fig.add_subplot(222)

        ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_bottom), RV_cumulative(y2_60k_top), color="b", alpha=0.2)
        ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_bottom), RV_cumulative(y2_80k_top), color="g", alpha=0.2)
        ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_bottom), RV_cumulative(y2_100k_top), color="r", alpha=0.2)

        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_bottom), marker='^', color="b", alpha=0.4)
        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_top), marker='o', color="b", alpha=0.4)

        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_bottom), marker='^', color="g", alpha=0.4)
        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_top), marker='o', color="g", alpha=0.4)

        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_bottom), marker='^', color="r", alpha=0.4)
        ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_top), marker='o', color="r", alpha=0.4)

        ax2.text(1.0, positiony_max, "M3", size=14)

        # ticks and labels
        ax2.set_xlim(0.5, len(bands_total)+0.5)
        ax2.set_xticks(range(1, len(bands_total) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax2.tick_params(axis='both', which='major', labelsize=15)

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_60k_bottom), RV_cumulative(y3_60k_top), color="b", alpha=0.2)
        ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_80k_bottom), RV_cumulative(y3_80k_top), color="g", alpha=0.2)
        ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_100k_bottom), RV_cumulative(y3_100k_top), color="r", alpha=0.2)

        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_60k_bottom), marker='^', color="b", alpha=0.4)
        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_60k_top), marker='o', color="b", alpha=0.4)

        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_80k_bottom), marker='^', color="g", alpha=0.4)
        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_80k_top), marker='o', color="g", alpha=0.4)

        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_100k_bottom), marker='^', color="r", alpha=0.4)
        ax3.scatter(range(1, len(bands_total) + 1), RV_cumulative(y3_100k_top), marker='o', color="r", alpha=0.4)

        ax3.text(1.0, positiony_max, "M6", size=14)

        # ticks and labels
        # ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands_total)+0.5)
        ax3.set_xticks(range(1, len(bands_total) + 1))
        ax3.set_xticklabels(bands_total)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3.tick_params(axis='both', which='major', labelsize=15)

        ax4 = fig.add_subplot(224)

        ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_bottom), RV_cumulative(y4_60k_top), color="b", alpha=0.2)
        ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_bottom), RV_cumulative(y4_80k_top), color="g", alpha=0.2)
        ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_bottom), RV_cumulative(y4_100k_top), color="r", alpha=0.2)

        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_bottom), marker='^', color="b", alpha=0.4)
        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_top), marker='o', color="b", alpha=0.4)

        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_bottom), marker='^', color="g", alpha=0.4)
        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_top), marker='o', color="g", alpha=0.4)

        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_bottom), marker='^', color="r", alpha=0.4)
        ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_top), marker='o', color="r", alpha=0.4)

        ax4.text(1.0, positiony_max, "M9", size=14)

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands_total)+0.5)
        ax4.set_xticks(range(1, len(bands_total) + 1))
        ax4.set_xticklabels(bands_total)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4.tick_params(axis='both', which='major', labelsize=15)

        fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        """
        VSINI=5km/s
        """
        print("Results for vsini of 5.0 km/s")

        # preparation of data: plot1
        y1_60k = [results["M0-Z-5.0-60k"], results["M0-Y-5.0-60k"], results["M0-J-5.0-60k"], results["M0-H-5.0-60k"], results["M0-K-5.0-60k"]]
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]

        y1_80k = [results["M0-Z-5.0-80k"], results["M0-Y-5.0-80k"], results["M0-J-5.0-80k"], results["M0-H-5.0-80k"], results["M0-K-5.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]

        y1_100k = [results["M0-Z-5.0-100k"], results["M0-Y-5.0-100k"], results["M0-J-5.0-100k"], results["M0-H-5.0-100k"], results["M0-K-5.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]

        y1_60k_lim = [y[0] for y in y1_60k]

        y1_80k_lim = [y[0] for y in y1_80k]

        y1_100k_lim = [y[0] for y in y1_100k]

        # preparation of data: plot2
        y2_60k = [results["M3-Z-5.0-60k"], results["M3-Y-5.0-60k"], results["M3-J-5.0-60k"], results["M3-H-5.0-60k"], results["M3-K-5.0-60k"]]
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]

        y2_80k = [results["M3-Z-5.0-80k"], results["M3-Y-5.0-80k"], results["M3-J-5.0-80k"], results["M3-H-5.0-80k"], results["M3-K-5.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]

        y2_100k = [results["M3-Z-5.0-100k"], results["M3-Y-5.0-100k"], results["M3-J-5.0-100k"], results["M3-H-5.0-100k"], results["M3-K-5.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k]

        y2_60k_lim = [y[0] for y in y2_60k]

        y2_80k_lim = [y[0] for y in y2_80k]

        y2_100k_lim = [y[0] for y in y2_100k]

        # preparation of data: plot3
        y3_60k = [results["M6-Z-5.0-60k"], results["M6-Y-5.0-60k"], results["M6-J-5.0-60k"], results["M6-H-5.0-60k"], results["M6-K-5.0-60k"]]
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]

        y3_80k = [results["M6-Z-5.0-80k"], results["M6-Y-5.0-80k"], results["M6-J-5.0-80k"], results["M6-H-5.0-80k"], results["M6-K-5.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]

        y3_100k = [results["M6-Z-5.0-100k"], results["M6-Y-5.0-100k"], results["M6-J-5.0-100k"], results["M6-H-5.0-100k"], results["M6-K-5.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k]

        y3_60k_lim = [y[0] for y in y3_60k]

        y3_80k_lim = [y[0] for y in y3_80k]

        y3_100k_lim = [y[0] for y in y3_100k]

        # preparation of data: plot4
        y4_60k = [results["M9-Z-5.0-60k"], results["M9-Y-5.0-60k"], results["M9-J-5.0-60k"], results["M9-H-5.0-60k"], results["M9-K-5.0-60k"]]
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]

        y4_80k = [results["M9-Z-5.0-80k"], results["M9-Y-5.0-80k"], results["M9-J-5.0-80k"], results["M9-H-5.0-80k"], results["M9-K-5.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]

        y4_100k = [results["M9-Z-5.0-100k"], results["M9-Y-5.0-100k"], results["M9-J-5.0-100k"], results["M9-H-5.0-100k"], results["M9-K-5.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]

        y4_60k_lim = [y[0] for y in y4_60k]

        y4_80k_lim = [y[0] for y in y4_80k]

        y4_100k_lim = [y[0] for y in y4_100k]

        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])

        """
        # data for correction plot
        y1_60k_vsini5 = (np.array(y1_60k_top) - np.array(y1_60k_bottom)) / np.array(y1_60k_bottom)*100.0
        y1_80k_vsini5 = (np.array(y1_80k_top) - np.array(y1_80k_bottom)) / np.array(y1_80k_bottom)*100.0
        y1_100k_vsini5 = (np.array(y1_100k_top) - np.array(y1_100k_bottom)) / np.array(y1_100k_bottom)*100.0

        y2_60k_vsini5 = (np.array(y2_60k_top) - np.array(y2_60k_bottom)) / np.array(y2_60k_bottom)*100.0
        y2_80k_vsini5 = (np.array(y2_80k_top) - np.array(y2_80k_bottom)) / np.array(y2_80k_bottom)*100.0
        y2_100k_vsini5 = (np.array(y2_100k_top) - np.array(y2_100k_bottom)) / np.array(y2_100k_bottom)*100.0

        y3_60k_vsini5 = (np.array(y3_60k_top) - np.array(y3_60k_bottom)) / np.array(y3_60k_bottom)*100.0
        y3_80k_vsini5 = (np.array(y3_80k_top) - np.array(y3_80k_bottom)) / np.array(y3_80k_bottom)*100.0
        y3_100k_vsini5 = (np.array(y3_100k_top) - np.array(y3_100k_bottom)) / np.array(y3_100k_bottom)*100.0

        y4_60k_vsini5 = (np.array(y4_60k_top) - np.array(y4_60k_bottom)) / np.array(y4_60k_bottom)*100.0
        y4_80k_vsini5 = (np.array(y4_80k_top) - np.array(y4_80k_bottom)) / np.array(y4_80k_bottom)*100.0
        y4_100k_vsini5 = (np.array(y4_100k_top) - np.array(y4_100k_bottom)) / np.array(y4_100k_bottom)*100.0
        """

        fig = plt.figure(2)
        ax1 = fig.add_subplot(221)

        ax1.fill_between(range(1, len(bands) + 1), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)

        ax1.plot(range(1, len(bands) + 1), y1_60k_lim, color="b", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_80k_lim, color="g", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_100k_lim, color="r", linestyle="--")

        ax1.scatter(range(1, len(bands) + 1), y1_60k_bottom, marker='^', color="b", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_60k_top, marker='o', color="b", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_80k_bottom, marker='^', color="g", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_top, marker='o', color="g", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_100k_bottom, marker='^', color="r", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_top, marker='o', color="r", alpha=0.4)

        ax1.text(1.0, positiony_max, "M0", size=14)

        # ticks and labels
        # ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(bands)+0.5)
        ax1.set_xticks(range(1, len(bands) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax1.tick_params(axis='both', which='major', labelsize=15)

        ax2 = fig.add_subplot(222)

        ax2.fill_between(range(1, len(bands) + 1), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.plot(range(1, len(bands) + 1), y2_60k_lim, color="b", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_80k_lim, color="g", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_100k_lim, color="r", linestyle="--")

        ax2.scatter(range(1, len(bands) + 1), y2_60k_bottom, marker='^', color="b", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_60k_top, marker='o', color="b", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_80k_bottom, marker='^', color="g", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_top, marker='o', color="g", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_100k_bottom, marker='^', color="r", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_top, marker='o', color="r", alpha=0.4)

        ax2.text(1.0, positiony_max, "M3", size=14)

        # ticks and labels
        ax2.set_xlim(0.5, len(bands)+0.5)
        ax2.set_xticks(range(1, len(bands) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax2.tick_params(axis='both', which='major', labelsize=15)

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1, len(bands) + 1), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.plot(range(1, len(bands) + 1), y3_60k_lim, color="b", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_80k_lim, color="g", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_100k_lim, color="r", linestyle="--")

        ax3.scatter(range(1, len(bands) + 1), y3_60k_bottom, marker='^', color="b", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_60k_top, marker='o', color="b", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_80k_bottom, marker='^', color="g", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_top, marker='o', color="g", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_100k_bottom, marker='^', color="r", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_top, marker='o', color="r", alpha=0.4)

        ax3.text(1.0, positiony_max, "M6", size=14)

        # ticks and labels
        # ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands)+0.5)
        ax3.set_xticks(range(1, len(bands) + 1))
        ax3.set_xticklabels(bands)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3.tick_params(axis='both', which='major', labelsize=15)

        ax4 = fig.add_subplot(224)

        ax4.fill_between(range(1, len(bands) + 1), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.plot(range(1, len(bands) + 1), y4_60k_lim, color="b", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_80k_lim, color="g", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_100k_lim, color="r", linestyle="--")

        ax4.scatter(range(1, len(bands) + 1), y4_60k_bottom, marker='^', color="b", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_60k_top, marker='o', color="b", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_80k_bottom, marker='^', color="g", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_top, marker='o', color="g", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_100k_bottom, marker='^', color="r", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_top, marker='o', color="r", alpha=0.4)

        ax4.text(1.0, positiony_max, "M9", size=14)

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands)+0.5)
        ax4.set_xticks(range(1, len(bands) + 1))
        ax4.set_xticklabels(bands)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4.tick_params(axis='both', which='major', labelsize=15)

        fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        """
        VSINI=10km/s
        """
        print("Results for vsini of 10.0 km/s")

        # preparation of data: plot1
        y1_60k = [results["M0-Z-10.0-60k"], results["M0-Y-10.0-60k"], results["M0-J-10.0-60k"], results["M0-H-10.0-60k"], results["M0-K-10.0-60k"]]
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]

        y1_80k = [results["M0-Z-10.0-80k"], results["M0-Y-10.0-80k"], results["M0-J-10.0-80k"], results["M0-H-10.0-80k"], results["M0-K-10.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]

        y1_100k = [results["M0-Z-10.0-100k"], results["M0-Y-10.0-100k"], results["M0-J-10.0-100k"], results["M0-H-10.0-100k"], results["M0-K-10.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]

        y1_60k_lim = [y[0] for y in y1_60k]

        y1_80k_lim = [y[0] for y in y1_80k]

        y1_100k_lim = [y[0] for y in y1_100k]

        # preparation of data: plot2
        y2_60k = [results["M3-Z-10.0-60k"], results["M3-Y-10.0-60k"], results["M3-J-10.0-60k"], results["M3-H-10.0-60k"], results["M3-K-10.0-60k"]]
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]

        y2_80k = [results["M3-Z-10.0-80k"], results["M3-Y-10.0-80k"], results["M3-J-10.0-80k"], results["M3-H-10.0-80k"], results["M3-K-10.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]

        y2_100k = [results["M3-Z-10.0-100k"], results["M3-Y-10.0-100k"], results["M3-J-10.0-100k"], results["M3-H-10.0-100k"], results["M3-K-10.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k]

        y2_60k_lim = [y[0] for y in y2_60k]

        y2_80k_lim = [y[0] for y in y2_80k]

        y2_100k_lim = [y[0] for y in y2_100k]

        # preparation of data: plot3
        y3_60k = [results["M6-Z-10.0-60k"], results["M6-Y-10.0-60k"], results["M6-J-10.0-60k"], results["M6-H-10.0-60k"], results["M6-K-10.0-60k"]]
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]

        y3_80k = [results["M6-Z-10.0-80k"], results["M6-Y-10.0-80k"], results["M6-J-10.0-80k"], results["M6-H-10.0-80k"], results["M6-K-10.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]

        y3_100k = [results["M6-Z-10.0-100k"], results["M6-Y-10.0-100k"], results["M6-J-10.0-100k"], results["M6-H-10.0-100k"], results["M6-K-10.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k]

        y3_60k_lim = [y[0] for y in y3_60k]

        y3_80k_lim = [y[0] for y in y3_80k]

        y3_100k_lim = [y[0] for y in y3_100k]

        # preparation of data: plot4
        y4_60k = [results["M9-Z-10.0-60k"], results["M9-Y-10.0-60k"], results["M9-J-10.0-60k"], results["M9-H-10.0-60k"], results["M9-K-10.0-60k"]]
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]

        y4_80k = [results["M9-Z-10.0-80k"], results["M9-Y-10.0-80k"], results["M9-J-10.0-80k"], results["M9-H-10.0-80k"], results["M9-K-10.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]

        y4_100k = [results["M9-Z-10.0-100k"], results["M9-Y-10.0-100k"], results["M9-J-10.0-100k"], results["M9-H-10.0-100k"], results["M9-K-10.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]

        y4_60k_lim = [y[0] for y in y4_60k]

        y4_80k_lim = [y[0] for y in y4_80k]

        y4_100k_lim = [y[0] for y in y4_100k]

        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])

        # data for correction plot
        y1_60k_vsini10 = (np.array(y1_60k_top) - np.array(y1_60k_bottom)) / np.array(y1_60k_bottom)*100.0
        y1_80k_vsini10 = (np.array(y1_80k_top) - np.array(y1_80k_bottom)) / np.array(y1_80k_bottom)*100.0
        y1_100k_vsini10 = (np.array(y1_100k_top) - np.array(y1_100k_bottom)) / np.array(y1_100k_bottom)*100.0

        y2_60k_vsini10 = (np.array(y2_60k_top) - np.array(y2_60k_bottom)) / np.array(y2_60k_bottom)*100.0
        y2_80k_vsini10 = (np.array(y2_80k_top) - np.array(y2_80k_bottom)) / np.array(y2_80k_bottom)*100.0
        y2_100k_vsini10 = (np.array(y2_100k_top) - np.array(y2_100k_bottom)) / np.array(y2_100k_bottom)*100.0

        y3_60k_vsini10 = (np.array(y3_60k_top) - np.array(y3_60k_bottom)) / np.array(y3_60k_bottom)*100.0
        y3_80k_vsini10 = (np.array(y3_80k_top) - np.array(y3_80k_bottom)) / np.array(y3_80k_bottom)*100.0
        y3_100k_vsini10 = (np.array(y3_100k_top) - np.array(y3_100k_bottom)) / np.array(y3_100k_bottom)*100.0

        y4_60k_vsini10 = (np.array(y4_60k_top) - np.array(y4_60k_bottom)) / np.array(y4_60k_bottom)*100.0
        y4_80k_vsini10 = (np.array(y4_80k_top) - np.array(y4_80k_bottom)) / np.array(y4_80k_bottom)*100.0
        y4_100k_vsini10 = (np.array(y4_100k_top) - np.array(y4_100k_bottom)) / np.array(y4_100k_bottom)*100.0

        fig = plt.figure(3)
        ax1 = fig.add_subplot(221)

        ax1.fill_between(range(1, len(bands) + 1), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)

        ax1.plot(range(1, len(bands) + 1), y1_60k_lim, color="b", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_80k_lim, color="g", linestyle="--")
        ax1.plot(range(1, len(bands) + 1), y1_100k_lim, color="r", linestyle="--")

        ax1.scatter(range(1, len(bands) + 1), y1_60k_bottom, marker='^', color="b", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_60k_top, marker='o', color="b", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_80k_bottom, marker='^', color="g", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_top, marker='o', color="g", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_100k_bottom, marker='^', color="r", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_top, marker='o', color="r", alpha=0.4)

        ax1.text(1.0, positiony_max, "M0", size=14)

        # ticks and labels
        # ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(bands)+0.5)
        ax1.set_xticks(range(1, len(bands) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax1.tick_params(axis='both', which='major', labelsize=15)

        ax2 = fig.add_subplot(222)

        ax2.fill_between(range(1, len(bands) + 1), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.plot(range(1, len(bands) + 1), y2_60k_lim, color="b", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_80k_lim, color="g", linestyle="--")
        ax2.plot(range(1, len(bands) + 1), y2_100k_lim, color="r", linestyle="--")

        ax2.scatter(range(1, len(bands) + 1), y2_60k_bottom, marker='^', color="b", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_60k_top, marker='o', color="b", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_80k_bottom, marker='^', color="g", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_top, marker='o', color="g", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_100k_bottom, marker='^', color="r", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_top, marker='o', color="r", alpha=0.4)

        ax2.text(1.0, positiony_max, "M3", size=14)

        # ticks and labels
        ax2.set_xlim(0.5, len(bands)+0.5)
        ax2.set_xticks(range(1, len(bands) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax2.tick_params(axis='both', which='major', labelsize=15)

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1, len(bands) + 1), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.plot(range(1, len(bands) + 1), y3_60k_lim, color="b", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_80k_lim, color="g", linestyle="--")
        ax3.plot(range(1, len(bands) + 1), y3_100k_lim, color="r", linestyle="--")

        ax3.scatter(range(1, len(bands) + 1), y3_60k_bottom, marker='^', color="b", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_60k_top, marker='o', color="b", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_80k_bottom, marker='^', color="g", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_top, marker='o', color="g", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_100k_bottom, marker='^', color="r", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_top, marker='o', color="r", alpha=0.4)

        ax3.text(1.0, positiony_max, "M6", size=14)

        # ticks and labels
       # ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands)+0.5)
        ax3.set_xticks(range(1, len(bands) + 1))
        ax3.set_xticklabels(bands)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3.tick_params(axis='both', which='major', labelsize=15)

        ax4 = fig.add_subplot(224)

        ax4.fill_between(range(1, len(bands) + 1), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.plot(range(1, len(bands) + 1), y4_60k_lim, color="b", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_80k_lim, color="g", linestyle="--")
        ax4.plot(range(1, len(bands) + 1), y4_100k_lim, color="r", linestyle="--")

        ax4.scatter(range(1, len(bands) + 1), y4_60k_bottom, marker='^', color="b", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_60k_top, marker='o', color="b", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_80k_bottom, marker='^', color="g", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_top, marker='o', color="g", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_100k_bottom, marker='^', color="r", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_top, marker='o', color="r", alpha=0.4)

        ax4.text(1.0, positiony_max, "M9", size=14)

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands)+0.5)
        ax4.set_xticks(range(1, len(bands) + 1))
        ax4.set_xticklabels(bands)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4.tick_params(axis='both', which='major', labelsize=15)

        fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        # correction figure
        positiony_max = np.max([np.max(y1_60k_vsini1), np.max(y1_60k_vsini10), np.max(y1_80k_vsini1), np.max(y1_80k_vsini10), np.max(y1_100k_vsini1), np.max(y1_100k_vsini10),\
                                   np.max(y2_60k_vsini1), np.max(y2_60k_vsini10), np.max(y2_80k_vsini1), np.max(y2_80k_vsini10), np.max(y2_100k_vsini1), np.max(y2_100k_vsini10),\
                                   np.max(y3_60k_vsini1), np.max(y3_60k_vsini10), np.max(y3_80k_vsini1), np.max(y3_80k_vsini10), np.max(y3_100k_vsini1), np.max(y3_100k_vsini10),\
                                   np.max(y4_60k_vsini1), np.max(y4_60k_vsini10), np.max(y4_80k_vsini1), np.max(y4_80k_vsini10), np.max(y4_100k_vsini1), np.max(y4_100k_vsini10)])

        positiony_min = np.min([np.min(y1_60k_vsini1), np.min(y1_60k_vsini10), np.min(y1_80k_vsini1), np.min(y1_80k_vsini10), np.min(y1_100k_vsini1), np.min(y1_100k_vsini10),\
                                   np.min(y2_60k_vsini1), np.min(y2_60k_vsini10), np.min(y2_80k_vsini1), np.min(y2_80k_vsini10), np.min(y2_100k_vsini1), np.min(y2_100k_vsini10),\
                                   np.min(y3_60k_vsini1), np.min(y3_60k_vsini10), np.min(y3_80k_vsini1), np.min(y3_80k_vsini10), np.min(y3_100k_vsini1), np.min(y3_100k_vsini10),\
                                   np.min(y4_60k_vsini1), np.min(y4_60k_vsini10), np.min(y4_80k_vsini1), np.min(y4_80k_vsini10), np.min(y4_100k_vsini1), np.min(y4_100k_vsini10)])

        rc('xtick', labelsize=15)
        rc('ytick', labelsize=15)

        fig = plt.figure(4)
        ax1 = fig.add_subplot(221)

        ax1.scatter(range(1, len(bands) + 1), y1_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

        ax1.scatter(range(1, len(bands) + 1), y1_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

        ax1.text(1.0, positiony_max*0.90, "M0", size=14)

        # ticks and labels
        ax1.set_xlim(0.5, len(bands)+0.5)
        ax1.set_xticks(range(1, len(bands) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.05*(positiony_max-positiony_min), positiony_max+0.05*(positiony_max-positiony_min))

        ax2 = fig.add_subplot(222)

        ax2.scatter(range(1, len(bands) + 1), y2_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

        ax2.scatter(range(1, len(bands) + 1), y2_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

        ax2.text(1.0, positiony_max*0.90, "M3", size=14)

        # ticks and labels
        ax2.set_xlim(0.5, len(bands)+0.5)
        ax2.set_xticks(range(1, len(bands) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.05*(positiony_max-positiony_min), positiony_max+0.05*(positiony_max-positiony_min))

        ax3 = fig.add_subplot(223)

        ax3.scatter(range(1, len(bands) + 1), y3_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

        ax3.scatter(range(1, len(bands) + 1), y3_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

        ax3.text(1.0, positiony_max*0.90, "M6", size=14)

        # ticks and labels
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands)+0.5)
        ax3.set_xticks(range(1, len(bands) + 1))
        ax3.set_xticklabels(bands)
        ax3.set_ylim(positiony_min-0.05*(positiony_max-positiony_min), positiony_max+0.05*(positiony_max-positiony_min))

        ax4 = fig.add_subplot(224)

        ax4.scatter(range(1, len(bands) + 1), y4_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

        ax4.scatter(range(1, len(bands) + 1), y4_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

        ax4.text(1.0, positiony_max*0.90, "M9", size=14)

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands)+0.5)
        ax4.set_xticks(range(1, len(bands) + 1))
        ax4.set_xticklabels(bands)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.05*(positiony_max-positiony_min), positiony_max+0.05*(positiony_max-positiony_min))

        fig.text(0.06, 0.5, r'Loss in RV precision [\%]', ha='center', va='center', rotation='vertical', size=14)
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        # if the paper_plots option is on, then print the results
        print("\n")
        for star in spectral_types:
            for band in bands:
                for vel in vsini:
                    for resolution in R:
                        for smpl in sampling:
                            id_string = star+"-"+band+"-"+vel+"-"+resolution    # sample was left aside because only one value existed
                            precision = results[id_string]
                            print("%s: & %.1f\t & %.1f\t & %.1f \\\\" % (id_string, precision[0], precision[1], precision[2]))
    else:
        return results


###############################################################################
#
#           OUTDATED FUNCTION!
#
###############################################################################
def calculate_prec_VIS(plot_atm=False, plot_ste=False, plot_flux=True, paper_plots=True):

    print("Reading atmospheric model...")
    wav_atm, flux_atm, std_flux_atm, mask_atm = prepare_atmopshere()
    print(("There were {:d} unmasked pixels out of {:d}., or {:.1%}."
          "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)))
    print("The model ranges from {:4.2f} to {:4.2f} micron.".format(wav_atm[0], wav_atm[-1]))
    print("Done.")

    print("Calculating impact of Barycentric movement on mask...")
    mask_atm_30kms = []
    for value in zip(wav_atm, mask_atm):
        if (value[1] is False):
            mask_atm_30kms.append(value[1])
        else:
            delta_lambda = value[0] * 3.0e4/Qcalculator.c
            indexes_30kmslice = np.searchsorted(wav_atm, [value[0]-delta_lambda, value[0]+delta_lambda])
            indexes_30kmslice = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in indexes_30kmslice]
            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]], dtype=bool)    # selecting only the slice in question
            mask_atm_30kms.append(all(mask_atm_30kmslice))

    mask_atm = np.array(mask_atm_30kms, dtype=bool)
    print(("There were {:d} unmasked pixels out of {:d}, or {:.1%}."
           "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) /
                      len(mask_atm)))

    if(plot_atm):

        # identify non-masked pixels
        selected_transmission = wav_atm[mask_atm]

        plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"Transmission [ ] ")
        plt.plot(wav_atm, flux_atm, color='k')
        plt.vlines(selected_transmission, 0.8, 1.0, colors="b")
        plt.xlim(wav_atm[0], wav_atm[-1])
        plt.ylim(0.0, 1.0)
        # plt.legend(loc='best')
        plt.show()
        plt.close()

    bands = ["VIS", "Z", "Y", "J", "H", "K"]

    results = {}    # creating empty dictionary for the results
    wav_plot_M0 = []   # creating empty lists for the plots
    flux_plot_M0 = []
    wav_plot_M3 = []
    flux_plot_M3 = []
    wav_plot_M6 = []
    flux_plot_M6 = []
    wav_plot_M9 = []
    flux_plot_M9 = []
    for star in spectral_types:
        for band in bands:
            for vel in ["1.0"]:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        # print("Working on "+file_to_read+".")
                        wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir+file_to_read)
                        # removing boundary effects
                        wav_stellar = wav_stellar[2:-2]
                        flux_stellar = flux_stellar[2:-2]  # / ((1.634e4)**2.0)

                        id_string = star+"-"+band+"-"+vel+"-"+resolution    # sample was left aside because only one value existed

                        # gettin the wav, flux and mask from the atm model closes to the stellar

                        index_atm = np.searchsorted(wav_atm, wav_stellar)
                        # replace indexes outside the array, at the very end, by the value at the very end
                        index_atm = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in index_atm]

                        wav_atm_selected = wav_atm[index_atm]
                        flux_atm_selected = flux_atm[index_atm]
                        mask_atm_selected = mask_atm[index_atm]

                        if("VIS" in id_string):
                            index_0p5 = np.searchsorted(wav_stellar, [0.5])
                            ratio_0p5 = flux_stellar[index_0p5+2] / flux_stellar[index_0p5-2]
                            flux_stellar[index_0p5:] = np.array(flux_stellar[index_0p5:])/ratio_0p5
                            flux_0p5 = flux_stellar[index_0p5]
                        elif("Z" in id_string):
                            flux_stellar = np.array(flux_stellar)/ratio_0p5

                        if("M0" in id_string):
                            # correcting for the different jumps:
                            if("Y" in id_string or "J" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M0[1]/ratios_M0[0])*(flux_0p5/flux_stellar[1])
                            elif("H" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M0[2]/ratios_M0[0])*(flux_0p5/flux_stellar[1])
                            elif("K" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M0[3]/ratios_M0[0])*(flux_0p5/flux_stellar[1])

                            norm_constant = 4712425.8785

                        elif("M3" in id_string):
                            if("Y" in id_string or "J" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M3[1]/ratios_M3[0])*(flux_0p5/flux_stellar[1])
                            elif("H" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M3[2]/ratios_M3[0])*(flux_0p5/flux_stellar[1])
                            elif("K" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M3[3]/ratios_M3[0])*(flux_0p5/flux_stellar[1])

                            norm_constant = 3873172.8673

                        elif("M6" in id_string):
                            if("Y" in id_string or "J" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M6[1]/ratios_M6[0])*(flux_0p5/flux_stellar[1])
                            elif("H" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M6[2]/ratios_M6[0])*(flux_0p5/flux_stellar[1])
                            elif("K" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M6[3]/ratios_M6[0])*(flux_0p5/flux_stellar[1])

                            if("1.0" in id_string):
                                norm_constant = 1713987.3947
                            elif("5.0" in id_string):
                                norm_constant = 1713987.3947
                            else:
                                norm_constant = 1713987.3947

                        elif("M9" in id_string):
                            if("Y" in id_string or "J" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M9[1]/ratios_M9[0])*(flux_0p5/flux_stellar[1])
                            elif("H" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M9[2]/ratios_M9[0])*(flux_0p5/flux_stellar[1])
                            elif("K" in id_string):
                                flux_stellar = np.array(flux_stellar)*(ratios_M9[3]/ratios_M9[0])*(flux_0p5/flux_stellar[1])

                            if("1.0" in id_string):
                                norm_constant = 1189097.2706
                            elif("5.0" in id_string):
                                norm_constant = 1189097.2706
                            else:
                                norm_constant = 1189097.2706

                        else:
                            print("Constant not defined. Aborting...")
                            exit()
                        flux_stellar = flux_stellar / ((norm_constant/100.0)**2.0)

                        if(id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k", "M6-J-1.0-100k", "M9-J-1.0-100k"]):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the %s reference model was of %4.2f." % (id_string, SN_estimate))
                        elif("J" in id_string):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the %s non-reference model was of %4.2f." % (id_string, SN_estimate))
                        if(plot_ste or plot_ste == id_string):
                            # Plot the stellar spectr as considered

                            selected_transmission_stellar = wav_atm_selected[mask_atm_selected]

                            plt.figure(1)
                            plt.xlabel(r"wavelength [$\mu$m])")
                            plt.ylabel(r"Flux_stellar [ ] ")
                            plt.plot(wav_stellar, flux_stellar, color='k')
                            plt.vlines(selected_transmission_stellar, np.min(flux_stellar), 0.3*np.max(flux_stellar), colors="b")
                            plt.xlim(wav_stellar[0], wav_stellar[-1])
                            plt.ylim(np.min(flux_stellar)-0.1*(np.max(flux_stellar)-np.min(flux_stellar)), np.max(flux_stellar)+0.1*(np.max(flux_stellar)-np.min(flux_stellar)))
                            # plt.legend(loc='best')
                            plt.show()
                            plt.close()

                        if(plot_flux and id_string in ["M0-VIS-1.0-100k", "M0-Z-1.0-100k", "M0-Y-1.0-100k", "M0-J-1.0-100k", "M0-H-1.0-100k", "M0-K-1.0-100k"]):
                            wav_plot_M0.append(wav_stellar)
                            flux_plot_M0.append(flux_stellar)
                        if(plot_flux and id_string in ["M3-VIS-1.0-100k", "M3-Z-1.0-100k", "M3-Y-1.0-100k", "M3-J-1.0-100k", "M3-H-1.0-100k", "M3-K-1.0-100k"]):
                            wav_plot_M3.append(wav_stellar)
                            flux_plot_M3.append(flux_stellar)
                        if(plot_flux and id_string in ["M6-VIS-1.0-100k", "M6-Z-1.0-100k", "M6-Y-1.0-100k", "M6-J-1.0-100k", "M6-H-1.0-100k", "M6-K-1.0-100k"]):
                            wav_plot_M6.append(wav_stellar)
                            flux_plot_M6.append(flux_stellar)
                        if(plot_flux and id_string in ["M9-VIS-1.0-100k", "M9-Z-1.0-100k", "M9-Y-1.0-100k", "M9-J-1.0-100k", "M9-H-1.0-100k", "M9-K-1.0-100k"]):
                            wav_plot_M9.append(wav_stellar)
                            flux_plot_M9.append(flux_stellar)

                        # precision given by the first method:
                        prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

                        # precision as given by the second_method
                        """
                        Example Joao
                        a = np.array([1, 5, 6, 8, 16, 34, 5, 7, 10, 83, 12, 6, 17, 18])
                        b = np.array([1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1], dtype=bool)

                        # this will give you a list of numpy arrays
                        c = np.array_split(a, np.where(np.diff(b))[0]+1)[::2]

                        # this will give you a list of lists
                        d = [list(cc) for cc in c]
                        print(d)
                        >>> [[1, 5], [16, 34, 5], [83, 12], [17, 18]]
                        """

                        wav_stellar_chunks_unformated = np.array_split(wav_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        wav_stellar_chunks = [list(chunk) for chunk in wav_stellar_chunks_unformated]

                        """
                        # test section
                        print("check that lengths are the same", len(wav_stellar), len(mask_atm_selected))
                        print("size of spectra %d vs number of chunks %d" % (len(wav_stellar), len(wav_stellar_chunks)))
                        print("number of true elements in all chunks: %d" % (len(mask_atm_selected[mask_atm_selected])))
                        """

                        flux_stellar_chunks_unformated = np.array_split(flux_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        flux_stellar_chunks = [list(chunk) for chunk in flux_stellar_chunks_unformated]

                        prec_2 = Qcalculator.RVprec_calc_chunks(wav_stellar_chunks, flux_stellar_chunks)

                        # precision as given by the third_method
                        prec_3 = Qcalculator.RV_prec_calc_Trans(wav_stellar, flux_stellar, flux_atm_selected)

                        # adding result to the dictionary
                        results[id_string] = [prec_1, prec_2, prec_3]
    if(plot_flux):
        # print the plot for flux comparison
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)

        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        ax1.plot(wav_plot_M0[0], np.array(flux_plot_M0[0]), color='0.1')
        ax1.plot(wav_plot_M0[1], np.array(flux_plot_M0[1]), color='0.1')
        ax1.plot(wav_plot_M0[2], np.array(flux_plot_M0[2]), color='0.1')
        ax1.plot(wav_plot_M0[3], np.array(flux_plot_M0[3]), color='0.1')
        ax1.plot(wav_plot_M0[4], np.array(flux_plot_M0[4]), color='0.1')
        ax1.plot(wav_plot_M0[5], np.array(flux_plot_M0[5]), color='0.1')
        ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)

        ax2.plot(wav_plot_M3[0], np.array(flux_plot_M3[0]), color='0.3')
        ax2.plot(wav_plot_M3[1], np.array(flux_plot_M3[1]), color='0.3')
        ax2.plot(wav_plot_M3[2], np.array(flux_plot_M3[2]), color='0.3')
        ax2.plot(wav_plot_M3[3], np.array(flux_plot_M3[3]), color='0.3')
        ax2.plot(wav_plot_M3[4], np.array(flux_plot_M3[4]), color='0.3')
        ax2.plot(wav_plot_M3[5], np.array(flux_plot_M3[5]), color='0.3')
        ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
        ax2.get_yaxis().get_offset_text().set_visible(False)            # remove offset from the yaxis
        ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='upper'))   # remove upper element from yaxis

        ax3.plot(wav_plot_M6[0], np.array(flux_plot_M6[0]), color='0.4')
        ax3.plot(wav_plot_M6[1], np.array(flux_plot_M6[1]), color='0.4')
        ax3.plot(wav_plot_M6[2], np.array(flux_plot_M6[2]), color='0.4')
        ax3.plot(wav_plot_M6[3], np.array(flux_plot_M6[3]), color='0.4')
        ax3.plot(wav_plot_M6[4], np.array(flux_plot_M6[4]), color='0.4')
        ax3.plot(wav_plot_M6[5], np.array(flux_plot_M6[5]), color='0.4')

        ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes)
        ax3.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        ax4.plot(wav_plot_M9[0], np.array(flux_plot_M9[0]), color='0.6')
        ax4.plot(wav_plot_M9[1], np.array(flux_plot_M9[1]), color='0.6')
        ax4.plot(wav_plot_M9[2], np.array(flux_plot_M9[2]), color='0.6')
        ax4.plot(wav_plot_M9[3], np.array(flux_plot_M9[3]), color='0.6')
        ax4.plot(wav_plot_M9[4], np.array(flux_plot_M9[4]), color='0.6')
        ax4.plot(wav_plot_M9[5], np.array(flux_plot_M9[5]), color='0.6')
        ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes)
        ax4.get_yaxis().get_offset_text().set_visible(False)             # remove offset from the yaxis

        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        f.text(0.06, 0.5, r'Flux [ ]', ha='center', va='center', rotation='vertical')
        plt.xlabel(r"wavelength [$\mu$m])")

        plt.show()
        plt.close()

    if(paper_plots):
        # print the paper plots
        """
        In every plot we will assume sample=3 and vsini=1
        y = RVprec between prec3 and prec2
        x = different bands

        different panels for different spectral types

        different colors will represent the different resolutions
        """

        print("Results for vsini of 1.0 km/s")
        # preparation of data: plot1
        y1_60k = [results["M0-VIS-1.0-60k"], results["M0-Z-1.0-60k"], results["M0-Y-1.0-60k"], results["M0-J-1.0-60k"], results["M0-H-1.0-60k"], results["M0-K-1.0-60k"]]
        y1_60k_top = [y[1] for y in y1_60k]
        y1_60k_bottom = [y[2] for y in y1_60k]

        y1_80k = [results["M0-VIS-1.0-80k"], results["M0-Z-1.0-80k"], results["M0-Y-1.0-80k"], results["M0-J-1.0-80k"], results["M0-H-1.0-80k"], results["M0-K-1.0-80k"]]
        y1_80k_top = [y[1] for y in y1_80k]
        y1_80k_bottom = [y[2] for y in y1_80k]

        y1_100k = [results["M0-VIS-1.0-100k"], results["M0-Z-1.0-100k"], results["M0-Y-1.0-100k"], results["M0-J-1.0-100k"], results["M0-H-1.0-100k"], results["M0-K-1.0-100k"]]
        y1_100k_top = [y[1] for y in y1_100k]
        y1_100k_bottom = [y[2] for y in y1_100k]

        # preparation of data: plot2
        y2_60k = [results["M3-VIS-1.0-60k"], results["M3-Z-1.0-60k"], results["M3-Y-1.0-60k"], results["M3-J-1.0-60k"], results["M3-H-1.0-60k"], results["M3-K-1.0-60k"]]
        y2_60k_top = [y[1] for y in y2_60k]
        y2_60k_bottom = [y[2] for y in y2_60k]

        y2_80k = [results["M3-VIS-1.0-80k"], results["M3-Z-1.0-80k"], results["M3-Y-1.0-80k"], results["M3-J-1.0-80k"], results["M3-H-1.0-80k"], results["M3-K-1.0-80k"]]
        y2_80k_top = [y[1] for y in y2_80k]
        y2_80k_bottom = [y[2] for y in y2_80k]

        y2_100k = [results["M3-VIS-1.0-100k"], results["M3-Z-1.0-100k"], results["M3-Y-1.0-100k"], results["M3-J-1.0-100k"], results["M3-H-1.0-100k"], results["M3-K-1.0-100k"]]
        y2_100k_top = [y[1] for y in y2_100k]
        y2_100k_bottom = [y[2] for y in y2_100k]

        # preparation of data: plot3
        y3_60k = [results["M6-VIS-1.0-60k"], results["M6-Z-1.0-60k"], results["M6-Y-1.0-60k"], results["M6-J-1.0-60k"], results["M6-H-1.0-60k"], results["M6-K-1.0-60k"]]
        y3_60k_top = [y[1] for y in y3_60k]
        y3_60k_bottom = [y[2] for y in y3_60k]

        y3_80k = [results["M6-VIS-1.0-80k"], results["M6-Z-1.0-80k"], results["M6-Y-1.0-80k"], results["M6-J-1.0-80k"], results["M6-H-1.0-80k"], results["M6-K-1.0-80k"]]
        y3_80k_top = [y[1] for y in y3_80k]
        y3_80k_bottom = [y[2] for y in y3_80k]

        y3_100k = [results["M6-VIS-1.0-100k"], results["M6-Z-1.0-100k"], results["M6-Y-1.0-100k"], results["M6-J-1.0-100k"], results["M6-H-1.0-100k"], results["M6-K-1.0-100k"]]
        y3_100k_top = [y[1] for y in y3_100k]
        y3_100k_bottom = [y[2] for y in y3_100k]

        # preparation of data: plot4
        y4_60k = [results["M9-VIS-1.0-60k"], results["M9-Z-1.0-60k"], results["M9-Y-1.0-60k"], results["M9-J-1.0-60k"], results["M9-H-1.0-60k"], results["M9-K-1.0-60k"]]
        y4_60k_top = [y[1] for y in y4_60k]
        y4_60k_bottom = [y[2] for y in y4_60k]

        y4_80k = [results["M9-VIS-1.0-80k"], results["M9-Z-1.0-80k"], results["M9-Y-1.0-80k"], results["M9-J-1.0-80k"], results["M9-H-1.0-80k"], results["M9-K-1.0-80k"]]
        y4_80k_top = [y[1] for y in y4_80k]
        y4_80k_bottom = [y[2] for y in y4_80k]

        y4_100k = [results["M9-VIS-1.0-100k"], results["M9-Z-1.0-100k"], results["M9-Y-1.0-100k"], results["M9-J-1.0-100k"], results["M9-H-1.0-100k"], results["M9-K-1.0-100k"]]
        y4_100k_top = [y[1] for y in y4_100k]
        y4_100k_bottom = [y[2] for y in y4_100k]

        positiony_max = np.max([np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top), np.max(y1_100k_bottom),\
                                   np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top), np.max(y2_100k_bottom),\
                                   np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top), np.max(y3_100k_bottom),\
                                   np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top), np.max(y4_100k_bottom)])

        positiony_min = np.min([np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top), np.min(y1_100k_bottom),\
                                   np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top), np.min(y2_100k_bottom),\
                                   np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top), np.min(y3_100k_bottom),\
                                   np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top), np.min(y4_100k_bottom)])

        fig = plt.figure(1)
        ax1 = fig.add_subplot(221)

        ax1.fill_between(range(1, len(bands) + 1), y1_60k_bottom, y1_60k_top, color="b", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_80k_bottom, y1_80k_top, color="g", alpha=0.2)
        ax1.fill_between(range(1, len(bands) + 1), y1_100k_bottom, y1_100k_top, color="r", alpha=0.2)

        ax1.scatter(range(1, len(bands) + 1), y1_60k_bottom, marker='^', color="b", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_60k_top, marker='o', color="b", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_80k_bottom, marker='^', color="g", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_80k_top, marker='o', color="g", alpha=0.4)

        ax1.scatter(range(1, len(bands) + 1), y1_100k_bottom, marker='^', color="r", alpha=0.4)
        ax1.scatter(range(1, len(bands) + 1), y1_100k_top, marker='o', color="r", alpha=0.4)

        ax1.text(1.0, positiony_max, "M0")

        # ticks and labels
        # ax1.set_ylabel('Precision [m/s]')
        ax1.set_xlim(0.5, len(bands)+0.5)
        ax1.set_xticks(range(1, len(bands) + 1))
        ax1.set_xticklabels([])
        ax1.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax2 = fig.add_subplot(222)

        ax2.fill_between(range(1, len(bands) + 1), y2_60k_bottom, y2_60k_top, color="b", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_80k_bottom, y2_80k_top, color="g", alpha=0.2)
        ax2.fill_between(range(1, len(bands) + 1), y2_100k_bottom, y2_100k_top, color="r", alpha=0.2)

        ax2.scatter(range(1, len(bands) + 1), y2_60k_bottom, marker='^', color="b", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_60k_top, marker='o', color="b", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_80k_bottom, marker='^', color="g", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_80k_top, marker='o', color="g", alpha=0.4)

        ax2.scatter(range(1, len(bands) + 1), y2_100k_bottom, marker='^', color="r", alpha=0.4)
        ax2.scatter(range(1, len(bands) + 1), y2_100k_top, marker='o', color="r", alpha=0.4)

        ax2.text(1.0, positiony_max, "M3")

        # ticks and labels
        ax2.set_xlim(0.5, len(bands)+0.5)
        ax2.set_xticks(range(1, len(bands) + 1))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax3 = fig.add_subplot(223)

        ax3.fill_between(range(1, len(bands) + 1), y3_60k_bottom, y3_60k_top, color="b", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_80k_bottom, y3_80k_top, color="g", alpha=0.2)
        ax3.fill_between(range(1, len(bands) + 1), y3_100k_bottom, y3_100k_top, color="r", alpha=0.2)

        ax3.scatter(range(1, len(bands) + 1), y3_60k_bottom, marker='^', color="b", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_60k_top, marker='o', color="b", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_80k_bottom, marker='^', color="g", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_80k_top, marker='o', color="g", alpha=0.4)

        ax3.scatter(range(1, len(bands) + 1), y3_100k_bottom, marker='^', color="r", alpha=0.4)
        ax3.scatter(range(1, len(bands) + 1), y3_100k_top, marker='o', color="r", alpha=0.4)

        ax3.text(1.0, positiony_max, "M6")

        # ticks and labels
        # ax3.set_ylabel('Precision [m/s]')
        ax3.set_xlabel('Bands')
        ax3.set_xlim(0.5, len(bands)+0.5)
        ax3.set_xticks(range(1, len(bands) + 1))
        ax3.set_xticklabels(bands)
        ax3.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        ax4 = fig.add_subplot(224)

        ax4.fill_between(range(1, len(bands) + 1), y4_60k_bottom, y4_60k_top, color="b", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_80k_bottom, y4_80k_top, color="g", alpha=0.2)
        ax4.fill_between(range(1, len(bands) + 1), y4_100k_bottom, y4_100k_top, color="r", alpha=0.2)

        ax4.scatter(range(1, len(bands) + 1), y4_60k_bottom, marker='^', color="b", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_60k_top, marker='o', color="b", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_80k_bottom, marker='^', color="g", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_80k_top, marker='o', color="g", alpha=0.4)

        ax4.scatter(range(1, len(bands) + 1), y4_100k_bottom, marker='^', color="r", alpha=0.4)
        ax4.scatter(range(1, len(bands) + 1), y4_100k_top, marker='o', color="r", alpha=0.4)

        ax4.text(1.0, positiony_max, "M9")

        # ticks and labels
        ax4.set_xlabel('Bands')
        ax4.set_xlim(0.5, len(bands)+0.5)
        ax4.set_xticks(range(1, len(bands) + 1))
        ax4.set_xticklabels(bands)
        ax4.set_yticklabels([])
        ax4.set_ylim(positiony_min-0.1*(positiony_max-positiony_min), positiony_max+0.1*(positiony_max-positiony_min))

        fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical')
        fig.subplots_adjust(hspace=0, wspace=0)

        plt.show()
        plt.close()

        # if the paper_plots option is on, then print the results
        print("\n")
        for star in spectral_types:
            for band in bands:
                for vel in ["1.0"]:
                    for resolution in R:
                        for smpl in sampling:
                            id_string = star+"-"+band+"-"+vel+"-"+resolution    # sample was left aside because only one value existed
                            precision = results[id_string]
                            print("%s: & %.1f\t & %.1f\t & %.1f \\\\" % (id_string, precision[0], precision[1], precision[2]))
    else:
        return results


###############################################################################
def compare_runs():
    """
    Function that compares spectra as resampled in the two versions of the code
    """
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        # print "Working on "+file_to_read+"."
                        wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir+file_to_read)
                        wav_stellar = wav_stellar
                        flux_stellar = flux_stellar / ((1.634e4)**2.0)

                        wav_stellar_OLD, flux_stellar_OLD = IOmodule.pdread_2col(resampled_dir_OLD+file_to_read)
                        flux_stellar_OLD = np.array(flux_stellar_OLD) / ((1.634e4)**2.0)

                        plt.figure(1)
                        plt.xlabel(r"wavelength [$\mu$m])")
                        plt.ylabel(r"Flux_stellar difference [ ] ")
                        plt.plot(wav_stellar, flux_stellar-flux_stellar_OLD, color='k')
                        plt.show()
                        plt.close()


def compare_output():
    """
    function that compares a spectrum prior to convolution, after, and after resampling
    """

    pre_convolution = "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    pre_wav, pre_flux = IOmodule.pdread_2col(pre_convolution)
    pre_wav = np.array(pre_wav, dtype="float64")*1.0e-4  # conversion to microns
    pre_flux = np.array(pre_flux, dtype="float64")*pre_wav

    convolved = "results_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k.txt"
    sampled = "resampled_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt"

    conv_wav, theor_flux, conv_flux = IOmodule.pdread_3col(convolved)
    sampled_wav, sampled_flux = IOmodule.pdread_2col(sampled)

    theor_flux = np.array(theor_flux)
    conv_flux = np.array(conv_flux)

    ratio_flux = moving_average(conv_flux, 300) / moving_average(theor_flux, 300)
    ratio_flux = ratio_flux/ratio_flux[0]

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux[ ] ")
    plt.plot(conv_wav, np.array(theor_flux)/theor_flux[0], color='k')
    plt.plot(conv_wav, np.array(conv_flux)/conv_flux[0], color='b')
    plt.plot(conv_wav, ratio_flux, color='g', linestyle='--')
    plt.show()
    plt.close()

    conv_flux_corrected = conv_flux / ratio_flux

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux corrected[ ] ")
    plt.plot(conv_wav, np.array(theor_flux)/theor_flux[0], color='k')
    plt.plot(conv_wav, np.array(conv_flux_corrected)/conv_flux_corrected[0], color='b')
    plt.show()
    plt.close()


def RV_cumulative(RV_vector):
    """
    funtion that calculates the cumulative RV vector weighted_error
    """

    return[weighted_error(RV_vector[:2]), weighted_error(RV_vector[:3]), weighted_error(RV_vector[:4]), weighted_error(RV_vector)]


def weighted_error(RV_vector):
    """
    function that calculates the average weighted error from a vector of errors
    """

    RV_vector = np.array(RV_vector)
    RV_value = 1.0/(np.sqrt(np.sum((1.0/RV_vector)**2.0)))

    return RV_value


def moving_average(x, window_size):
    """
    moving average
    """
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(x, window, 'same')


###############################################################################

if __name__ == "__main__":
    calculate_prec()
