# Placed here to keep it around.


# Functions from rv precision calculations that were tried for first paper but not used in results
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


# Also removed from nIR_precision
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
