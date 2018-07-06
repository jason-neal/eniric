"""
Plotting functions for nIR_precision.

"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
# to remove labels in one tick
from matplotlib.ticker import MaxNLocator

from eniric.utilities import band_selector

# set stuff for latex usage
rc('text', usetex=True)


def plot_atmosphere_model(wav_atm, flux_atm, mask_atm):
    """Plot atmospheric transmission model for each band.

    This is run when plot_atom=True in calculate_prec()"""
    # identify non-masked pixels
    selected_transmission = wav_atm[mask_atm]
    dummy_vector = np.ones_like(selected_transmission)
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

    wav_mask_Z, __ = band_selector(selected_transmission, dummy_vector, "Z")
    wav_mask_Y, __ = band_selector(selected_transmission, dummy_vector, "Y")
    wav_mask_J, __ = band_selector(selected_transmission, dummy_vector, "J")
    wav_mask_H, __ = band_selector(selected_transmission, dummy_vector, "H")
    wav_mask_K, __ = band_selector(selected_transmission, dummy_vector, "K")

    print("Z:", len(wav_mask_Z) / float(len(wav_atm_Z)), np.average(flux_atm_Z), np.median(flux_atm_Z))
    print("Y:", len(wav_mask_Y) / float(len(wav_atm_Y)), np.average(flux_atm_Y), np.median(flux_atm_Y))
    print("J:", len(wav_mask_J) / float(len(wav_atm_J)), np.average(flux_atm_J), np.median(flux_atm_J))
    print("H:", len(wav_mask_H) / float(len(wav_atm_H)), np.average(flux_atm_H), np.median(flux_atm_H))
    print("K:", len(wav_mask_K) / float(len(wav_atm_K)), np.average(flux_atm_K), np.median(flux_atm_K))

    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=False, sharey=False)

    ax1.plot(wav_atm_Z, flux_atm_Z, color='k')
    # ax1.vlines(selected_transmission, 0.0, 0.2, colors="b")
    ax1.set_xlim([wav_atm_Z[0], wav_atm_Z[-1]])
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=len(ax1.get_yticklabels()), prune='lower'))
    ax1.text(0.9, 0.8, ' Z band ', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,
             size=12, color='black', bbox=dict(facecolor='white'))

    ax2.plot(wav_atm_Y, flux_atm_Y, color='k')
    # ax2.vlines(selected_transmission, 0.0, 0.2, colors="b")
    ax2.set_xlim([wav_atm_Y[0], wav_atm_Y[-1]])
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='lower'))
    ax2.text(0.9, 0.8, ' Y band ', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,
             size=12, color='black', bbox=dict(facecolor='white'))

    ax3.plot(wav_atm_J, flux_atm_J, color='k')
    # ax3.vlines(selected_transmission, 0.0, 0.2, colors="b")
    ax3.set_xlim([wav_atm_J[0], wav_atm_J[-1]])
    ax3.yaxis.set_major_locator(MaxNLocator(nbins=len(ax3.get_yticklabels()), prune='upper'))
    ax3.text(0.9, 0.8, ' J band', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes,
             size=12, color='black', bbox=dict(facecolor='white'))

    ax4.plot(wav_atm_H, flux_atm_H, color='k')
    # ax4.vlines(selected_transmission, 0.0, 0.2, colors="b")
    ax4.set_xlim([wav_atm_H[0], wav_atm_H[-1]])
    ax4.yaxis.set_major_locator(MaxNLocator(nbins=len(ax4.get_yticklabels()), prune='lower'))
    ax4.text(0.9, 0.8, ' H band ', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes,
             size=12, color='black', bbox=dict(facecolor='white'))

    ax5.plot(wav_atm_K, flux_atm_K, color='k')
    # ax5.vlines(selected_transmission, 0.0, 0.2, colors="b")
    ax5.set_xlim([wav_atm_K[0], wav_atm_K[-1]])
    ax5.yaxis.set_major_locator(MaxNLocator(nbins=len(ax5.get_yticklabels()), prune='upper'))
    ax5.text(0.9, 0.8, ' K band ', horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes,
             size=12, color='black', bbox=dict(facecolor='white'))

    plt.xlabel(r"wavelength [$\mu$m]")
    f.text(0.06, 0.5, r'Atmospheric Transmission [ ]', ha='center', va='center', rotation='vertical')
    plt.show()
    # f.savefig('AtmosphericTransmission.pdf', format='pdf')

    plt.close()
    sys.exit(0)


def plot_stellar_spectum(wav_stellar, flux_stellar, wav_atm_selected, mask_atm_selected):
    """Plot the stellar spectrum as considered."""
    selected_transmission_stellar = wav_atm_selected[mask_atm_selected]

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux_stellar [ ] ")
    plt.plot(wav_stellar, flux_stellar, color='k')
    plt.vlines(selected_transmission_stellar, np.min(flux_stellar),
               0.3 * np.max(flux_stellar), colors="b")
    plt.xlim(wav_stellar[0], wav_stellar[-1])
    plt.ylim(np.min(flux_stellar) - 0.1 * (np.max(flux_stellar) - np.min(flux_stellar)),
             np.max(flux_stellar) + 0.1 * (np.max(flux_stellar) - np.min(flux_stellar)))
    # plt.legend(loc='best')
    plt.show()
    plt.close()


def plot_nIR_flux():
    """Plot_flux code from nIR_precision"""
    # print the plot for flux comparison
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)

    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    ax1.plot(wav_plot_m0[0], np.array(flux_plot_m0[0]), color='0.1')
    ax1.plot(wav_plot_m0[1], np.array(flux_plot_m0[1]), color='0.1')
    ax1.plot(wav_plot_m0[2], np.array(flux_plot_m0[2]), color='0.1')
    ax1.plot(wav_plot_m0[3], np.array(flux_plot_m0[3]), color='0.1')
    ax1.plot(wav_plot_m0[4], np.array(flux_plot_m0[4]), color='0.1')
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)
    ax1.plot(wav_m0, (flux_m0 / flux_m0[0]) * flux_plot_m0[0][1000], color='0.1', linestyle='--')

    ax2.plot(wav_plot_m3[0], np.array(flux_plot_m3[0]), color='0.3')
    ax2.plot(wav_plot_m3[1], np.array(flux_plot_m3[1]), color='0.3')
    ax2.plot(wav_plot_m3[2], np.array(flux_plot_m3[2]), color='0.3')
    ax2.plot(wav_plot_m3[3], np.array(flux_plot_m3[3]), color='0.3')
    ax2.plot(wav_plot_m3[4], np.array(flux_plot_m3[4]), color='0.3')
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.plot(wav_m3, (flux_m3 / flux_m3[0]) * flux_plot_m3[0][1000], color='0.1', linestyle='--')

    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))  # remove upper element from yaxis
    # ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='upper'))   # remove upper element from yaxis

    ax3.plot(wav_plot_m6[0], np.array(flux_plot_m6[0]), color='0.4')
    ax3.plot(wav_plot_m6[1], np.array(flux_plot_m6[1]), color='0.4')
    ax3.plot(wav_plot_m6[2], np.array(flux_plot_m6[2]), color='0.4')
    ax3.plot(wav_plot_m6[3], np.array(flux_plot_m6[3]), color='0.4')
    ax3.plot(wav_plot_m6[4], np.array(flux_plot_m6[4]), color='0.4')
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.plot(wav_m6, (flux_m6 / flux_m6[0]) * flux_plot_m6[0][1000] * 1.2, color='0.1', linestyle='--')

    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    ax4.plot(wav_plot_m9[0], np.array(flux_plot_m9[0]), color='0.6')
    ax4.plot(wav_plot_m9[1], np.array(flux_plot_m9[1]), color='0.6')
    ax4.plot(wav_plot_m9[2], np.array(flux_plot_m9[2]), color='0.6')
    ax4.plot(wav_plot_m9[3], np.array(flux_plot_m9[3]), color='0.6')
    ax4.plot(wav_plot_m9[4], np.array(flux_plot_m9[4]), color='0.6')
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.plot(wav_m9, (flux_m9 / flux_m9[0]) * flux_plot_m9[0][1000] * 1.4, color='0.1', linestyle='--')

    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

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
    ax1.set_xlim(np.min(wav_plot_m0[0]), np.max(wav_plot_m0[0]))
    ax1.set_ylim(0.0, 5.0e3)

    ax1.plot(wav_plot_m0[0], np.array(flux_plot_m0[0]), color='0.1', label="M0")
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

    ax2.plot(wav_plot_m3[0], np.array(flux_plot_m3[0]), color='0.3', label="M3")
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))  # remove upper element from yaxis

    ax3.plot(wav_plot_m6[0], np.array(flux_plot_m6[0]), color='0.4', label="M6")
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax3.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))

    ax4.plot(wav_plot_m9[0], np.array(flux_plot_m9[0]), color='0.6', label="M9")
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
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

    ax1.plot(wav_plot_m0[1], np.array(flux_plot_m0[1]), color='0.1', label="M0")
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

    ax2.plot(wav_plot_m3[1], np.array(flux_plot_m3[1]), color='0.3', label="M3")
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))  # remove upper element from yaxis

    ax3.plot(wav_plot_m6[1], np.array(flux_plot_m6[1]), color='0.4', label="M6")
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    ax4.plot(wav_plot_m9[1], np.array(flux_plot_m9[1]), color='0.6', label="M9")
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()
    plt.close()

    # print the J band spectra
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
    plt.xlabel(r"wavelength [$\mu$m]")
    f.text(0.06, 0.5, r"Flux in J band [ ]", ha='center', va='center', rotation='vertical')

    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax1.set_xlim(np.min(wav_plot_m0[2]), np.max(wav_plot_m0[2]))

    ax1.plot(wav_plot_m0[2], np.array(flux_plot_m0[2]), color='0.1', label="M0")
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

    ax2.plot(wav_plot_m3[2], np.array(flux_plot_m3[2]), color='0.3', label="M3")
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))  # remove upper element from yaxis

    ax3.plot(wav_plot_m6[2], np.array(flux_plot_m6[2]), color='0.4', label="M6")
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    ax4.plot(wav_plot_m9[2], np.array(flux_plot_m9[2]), color='0.6', label="M9")
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()
    plt.close()

    # print the H band spectra
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
    plt.xlabel(r"wavelength [$\mu$m]")
    f.text(0.06, 0.5, r"Flux in H band [ ]", ha='center', va='center', rotation='vertical')

    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    ax1.plot(wav_plot_m0[3], np.array(flux_plot_m0[3]), color='0.1', label="M0")
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

    ax2.plot(wav_plot_m3[3], np.array(flux_plot_m3[3]), color='0.3', label="M3")
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))  # remove upper element from yaxis

    ax3.plot(wav_plot_m6[3], np.array(flux_plot_m6[3]), color='0.4', label="M6")
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    ax4.plot(wav_plot_m9[3], np.array(flux_plot_m9[3]), color='0.6', label="M9")
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()
    plt.close()

    # print the K band spectra
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True, sharey=True)
    plt.xlabel(r"wavelength [$\mu$m]")
    f.text(0.06, 0.5, r"Flux in K band [ ]", ha='center', va='center', rotation='vertical')

    ax1.set_xlim([np.min(wav_plot_m0[4]), np.max(wav_plot_m0[4])])
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    ax1.plot(wav_plot_m0[4], np.array(flux_plot_m0[4]), color='0.1', label="M0")
    ax1.text(0.9, 0.8, 'M0', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, size=12)

    ax2.plot(wav_plot_m3[4], np.array(flux_plot_m3[4]), color='0.3', label="M3")
    ax2.text(0.9, 0.8, 'M3', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, size=12)
    ax2.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))  # remove upper element from yaxis

    ax3.plot(wav_plot_m6[4], np.array(flux_plot_m6[4]), color='0.4', label="M6")
    ax3.text(0.9, 0.8, 'M6', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, size=12)
    ax3.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    ax4.plot(wav_plot_m9[4], np.array(flux_plot_m9[4]), color='0.6', label="M9")
    ax4.text(0.9, 0.8, 'M9', horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, size=12)
    ax4.get_yaxis().get_offset_text().set_visible(False)  # remove offset from the yaxis

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    plt.show()
    plt.close()


def plot_paper_plots():
    """Print the paper plots.

    In every plot we will assume sample=3 and vsini=1
    y = RVprec between prec3 and prec2
    x = different bands

    different panels for different spectral types

    different colors will represent the different resolutions
    """

    print("Results for vsini of 1.0 km/s")
    # preparation of data: plot1
    y1_60k = [results["M0-Z-1.0-60k"], results["M0-Y-1.0-60k"], results["M0-J-1.0-60k"], results["M0-H-1.0-60k"],
              results["M0-K-1.0-60k"]]
    y1_60k_top = [y[1] for y in y1_60k]
    y1_60k_bottom = [y[2] for y in y1_60k]

    y1_60k_lim = [y[0] for y in y1_60k]

    y1_80k = [results["M0-Z-1.0-80k"], results["M0-Y-1.0-80k"], results["M0-J-1.0-80k"], results["M0-H-1.0-80k"],
              results["M0-K-1.0-80k"]]
    y1_80k_top = [y[1] for y in y1_80k]
    y1_80k_bottom = [y[2] for y in y1_80k]

    y1_80k_lim = [y[0] for y in y1_80k]

    y1_100k = [results["M0-Z-1.0-100k"], results["M0-Y-1.0-100k"], results["M0-J-1.0-100k"], results["M0-H-1.0-100k"],
               results["M0-K-1.0-100k"]]
    y1_100k_top = [y[1] for y in y1_100k]
    y1_100k_bottom = [y[2] for y in y1_100k]

    y1_100k_lim = [y[0] for y in y1_100k]

    # preparation of data: plot2
    y2_60k = [results["M3-Z-1.0-60k"], results["M3-Y-1.0-60k"], results["M3-J-1.0-60k"], results["M3-H-1.0-60k"],
              results["M3-K-1.0-60k"]]
    y2_60k_top = [y[1] for y in y2_60k]
    y2_60k_bottom = [y[2] for y in y2_60k]

    y2_60k_lim = [y[0] for y in y2_60k]

    y2_80k = [results["M3-Z-1.0-80k"], results["M3-Y-1.0-80k"], results["M3-J-1.0-80k"], results["M3-H-1.0-80k"],
              results["M3-K-1.0-80k"]]
    y2_80k_top = [y[1] for y in y2_80k]
    y2_80k_bottom = [y[2] for y in y2_80k]

    y2_80k_lim = [y[0] for y in y2_80k]

    y2_100k = [results["M3-Z-1.0-100k"], results["M3-Y-1.0-100k"], results["M3-J-1.0-100k"], results["M3-H-1.0-100k"],
               results["M3-K-1.0-100k"]]
    y2_100k_top = [y[1] for y in y2_100k]
    y2_100k_bottom = [y[2] for y in y2_100k]

    y2_100k_lim = [y[0] for y in y2_100k]

    # preparation of data: plot3
    y3_60k = [results["M6-Z-1.0-60k"], results["M6-Y-1.0-60k"], results["M6-J-1.0-60k"], results["M6-H-1.0-60k"],
              results["M6-K-1.0-60k"]]
    y3_60k_top = [y[1] for y in y3_60k]
    y3_60k_bottom = [y[2] for y in y3_60k]

    y3_60k_lim = [y[0] for y in y3_60k]

    y3_80k = [results["M6-Z-1.0-80k"], results["M6-Y-1.0-80k"], results["M6-J-1.0-80k"], results["M6-H-1.0-80k"],
              results["M6-K-1.0-80k"]]
    y3_80k_top = [y[1] for y in y3_80k]
    y3_80k_bottom = [y[2] for y in y3_80k]

    y3_80k_lim = [y[0] for y in y3_80k]

    y3_100k = [results["M6-Z-1.0-100k"], results["M6-Y-1.0-100k"], results["M6-J-1.0-100k"], results["M6-H-1.0-100k"],
               results["M6-K-1.0-100k"]]
    y3_100k_top = [y[1] for y in y3_100k]
    y3_100k_bottom = [y[2] for y in y3_100k]

    y3_100k_lim = [y[0] for y in y3_100k]

    # preparation of data: plot4
    y4_60k = [results["M9-Z-1.0-60k"], results["M9-Y-1.0-60k"], results["M9-J-1.0-60k"], results["M9-H-1.0-60k"],
              results["M9-K-1.0-60k"]]
    y4_60k_top = [y[1] for y in y4_60k]
    y4_60k_bottom = [y[2] for y in y4_60k]

    y4_60k_lim = [y[0] for y in y4_60k]

    y4_80k = [results["M9-Z-1.0-80k"], results["M9-Y-1.0-80k"], results["M9-J-1.0-80k"], results["M9-H-1.0-80k"],
              results["M9-K-1.0-80k"]]
    y4_80k_top = [y[1] for y in y4_80k]
    y4_80k_bottom = [y[2] for y in y4_80k]

    y4_80k_lim = [y[0] for y in y4_80k]

    y4_100k = [results["M9-Z-1.0-100k"], results["M9-Y-1.0-100k"], results["M9-J-1.0-100k"], results["M9-H-1.0-100k"],
               results["M9-K-1.0-100k"]]
    y4_100k_top = [y[1] for y in y4_100k]
    y4_100k_bottom = [y[2] for y in y4_100k]

    y4_100k_lim = [y[0] for y in y4_100k]

    positiony_max = np.max(
        [np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top),
         np.max(y1_100k_bottom),
         np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top),
         np.max(y2_100k_bottom),
         np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top),
         np.max(y3_100k_bottom),
         np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top),
         np.max(y4_100k_bottom)])

    positiony_min = np.min(
        [np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top),
         np.min(y1_100k_bottom),
         np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top),
         np.min(y2_100k_bottom),
         np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top),
         np.min(y3_100k_bottom),
         np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top),
         np.min(y4_100k_bottom)])

    # data for correction plot
    y1_60k_vsini1 = (np.array(y1_60k_top) - np.array(y1_60k_bottom)) / np.array(y1_60k_bottom) * 100.0
    y1_80k_vsini1 = (np.array(y1_80k_top) - np.array(y1_80k_bottom)) / np.array(y1_80k_bottom) * 100.0
    y1_100k_vsini1 = (np.array(y1_100k_top) - np.array(y1_100k_bottom)) / np.array(y1_100k_bottom) * 100.0

    y2_60k_vsini1 = (np.array(y2_60k_top) - np.array(y2_60k_bottom)) / np.array(y2_60k_bottom) * 100.0
    y2_80k_vsini1 = (np.array(y2_80k_top) - np.array(y2_80k_bottom)) / np.array(y2_80k_bottom) * 100.0
    y2_100k_vsini1 = (np.array(y2_100k_top) - np.array(y2_100k_bottom)) / np.array(y2_100k_bottom) * 100.0

    y3_60k_vsini1 = (np.array(y3_60k_top) - np.array(y3_60k_bottom)) / np.array(y3_60k_bottom) * 100.0
    y3_80k_vsini1 = (np.array(y3_80k_top) - np.array(y3_80k_bottom)) / np.array(y3_80k_bottom) * 100.0
    y3_100k_vsini1 = (np.array(y3_100k_top) - np.array(y3_100k_bottom)) / np.array(y3_100k_bottom) * 100.0

    y4_60k_vsini1 = (np.array(y4_60k_top) - np.array(y4_60k_bottom)) / np.array(y4_60k_bottom) * 100.0
    y4_80k_vsini1 = (np.array(y4_80k_top) - np.array(y4_80k_bottom)) / np.array(y4_80k_bottom) * 100.0
    y4_100k_vsini1 = (np.array(y4_100k_top) - np.array(y4_100k_bottom)) / np.array(y4_100k_bottom) * 100.0

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
    ax1.set_xlim(0.5, len(bands) + 0.5)
    ax1.set_xticks(range(1, len(bands) + 1))
    ax1.set_xticklabels([])
    ax1.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax2.set_xlim(0.5, len(bands) + 0.5)
    ax2.set_xticks(range(1, len(bands) + 1))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax3.set_xlim(0.5, len(bands) + 0.5)
    ax3.set_xticks(range(1, len(bands) + 1))
    ax3.set_xticklabels(bands)
    ax3.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax4.set_xlim(0.5, len(bands) + 0.5)
    ax4.set_xticks(range(1, len(bands) + 1))
    ax4.set_xticklabels(bands)
    ax4.set_yticklabels([])
    ax4.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

    ax4.tick_params(axis='both', which='major', labelsize=15)

    fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
    fig.subplots_adjust(hspace=0, wspace=0)

    plt.show()
    plt.close()

    """
    same plot for total precision
    """

    positiony_max = np.max(
        [np.max(RV_cumulative(y1_60k_top)), np.max(RV_cumulative(y1_60k_bottom)), np.max(RV_cumulative(y1_80k_top)),
         np.max(RV_cumulative(y1_80k_bottom)), np.max(RV_cumulative(y1_100k_top)),
         np.max(RV_cumulative(y1_100k_bottom)),
         np.max(RV_cumulative(y2_60k_top)), np.max(RV_cumulative(y2_60k_bottom)), np.max(RV_cumulative(y2_80k_top)),
         np.max(RV_cumulative(y2_80k_bottom)), np.max(RV_cumulative(y2_100k_top)),
         np.max(RV_cumulative(y2_100k_bottom)),
         np.max(RV_cumulative(y3_60k_top)), np.max(RV_cumulative(y3_60k_bottom)), np.max(RV_cumulative(y3_80k_top)),
         np.max(RV_cumulative(y3_80k_bottom)), np.max(RV_cumulative(y3_100k_top)),
         np.max(RV_cumulative(y3_100k_bottom)),
         np.max(RV_cumulative(y4_60k_top)), np.max(RV_cumulative(y4_60k_bottom)), np.max(RV_cumulative(y4_80k_top)),
         np.max(RV_cumulative(y4_80k_bottom)), np.max(RV_cumulative(y4_100k_top)),
         np.max(RV_cumulative(y4_100k_bottom))])

    positiony_min = np.min(
        [np.min(RV_cumulative(y1_60k_top)), np.min(RV_cumulative(y1_60k_bottom)), np.min(RV_cumulative(y1_80k_top)),
         np.min(RV_cumulative(y1_80k_bottom)), np.min(RV_cumulative(y1_100k_top)),
         np.min(RV_cumulative(y1_100k_bottom)),
         np.min(RV_cumulative(y2_60k_top)), np.min(RV_cumulative(y2_60k_bottom)), np.min(RV_cumulative(y2_80k_top)),
         np.min(RV_cumulative(y2_80k_bottom)), np.min(RV_cumulative(y2_100k_top)),
         np.min(RV_cumulative(y2_100k_bottom)),
         np.min(RV_cumulative(y3_60k_top)), np.min(RV_cumulative(y3_60k_bottom)), np.min(RV_cumulative(y3_80k_top)),
         np.min(RV_cumulative(y3_80k_bottom)), np.min(RV_cumulative(y3_100k_top)),
         np.min(RV_cumulative(y3_100k_bottom)),
         np.min(RV_cumulative(y4_60k_top)), np.min(RV_cumulative(y4_60k_bottom)), np.min(RV_cumulative(y4_80k_top)),
         np.min(RV_cumulative(y4_80k_bottom)), np.min(RV_cumulative(y4_100k_top)),
         np.min(RV_cumulative(y4_100k_bottom))])

    bands_total = ["ZY", "ZYJ", "ZYJH", "ZYJHK"]
    fig = plt.figure(1)
    ax1 = fig.add_subplot(221)

    ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_bottom), RV_cumulative(y1_60k_top), color="b",
                     alpha=0.2)
    ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_bottom), RV_cumulative(y1_80k_top), color="g",
                     alpha=0.2)
    ax1.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_bottom), RV_cumulative(y1_100k_top),
                     color="r", alpha=0.2)

    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_bottom), marker='^', color="b", alpha=0.4)
    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_60k_top), marker='o', color="b", alpha=0.4)

    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_bottom), marker='^', color="g", alpha=0.4)
    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_80k_top), marker='o', color="g", alpha=0.4)

    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_bottom), marker='^', color="r", alpha=0.4)
    ax1.scatter(range(1, len(bands_total) + 1), RV_cumulative(y1_100k_top), marker='o', color="r", alpha=0.4)

    ax1.text(1.0, positiony_max, "M0", size=14)

    # ticks and labels
    # ax1.set_ylabel('Precision [m/s]')
    ax1.set_xlim(0.5, len(bands_total) + 0.5)
    ax1.set_xticks(range(1, len(bands_total) + 1))
    ax1.set_xticklabels([])
    ax1.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

    ax1.tick_params(axis='both', which='major', labelsize=15)

    ax2 = fig.add_subplot(222)

    ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_bottom), RV_cumulative(y2_60k_top), color="b",
                     alpha=0.2)
    ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_bottom), RV_cumulative(y2_80k_top), color="g",
                     alpha=0.2)
    ax2.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_bottom), RV_cumulative(y2_100k_top),
                     color="r", alpha=0.2)

    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_bottom), marker='^', color="b", alpha=0.4)
    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_60k_top), marker='o', color="b", alpha=0.4)

    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_bottom), marker='^', color="g", alpha=0.4)
    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_80k_top), marker='o', color="g", alpha=0.4)

    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_bottom), marker='^', color="r", alpha=0.4)
    ax2.scatter(range(1, len(bands_total) + 1), RV_cumulative(y2_100k_top), marker='o', color="r", alpha=0.4)

    ax2.text(1.0, positiony_max, "M3", size=14)

    # ticks and labels
    ax2.set_xlim(0.5, len(bands_total) + 0.5)
    ax2.set_xticks(range(1, len(bands_total) + 1))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

    ax2.tick_params(axis='both', which='major', labelsize=15)

    ax3 = fig.add_subplot(223)

    ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_60k_bottom), RV_cumulative(y3_60k_top), color="b",
                     alpha=0.2)
    ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_80k_bottom), RV_cumulative(y3_80k_top), color="g",
                     alpha=0.2)
    ax3.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y3_100k_bottom), RV_cumulative(y3_100k_top),
                     color="r", alpha=0.2)

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
    ax3.set_xlim(0.5, len(bands_total) + 0.5)
    ax3.set_xticks(range(1, len(bands_total) + 1))
    ax3.set_xticklabels(bands_total)
    ax3.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

    ax3.tick_params(axis='both', which='major', labelsize=15)

    ax4 = fig.add_subplot(224)

    ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_bottom), RV_cumulative(y4_60k_top), color="b",
                     alpha=0.2)
    ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_bottom), RV_cumulative(y4_80k_top), color="g",
                     alpha=0.2)
    ax4.fill_between(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_bottom), RV_cumulative(y4_100k_top),
                     color="r", alpha=0.2)

    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_bottom), marker='^', color="b", alpha=0.4)
    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_60k_top), marker='o', color="b", alpha=0.4)

    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_bottom), marker='^', color="g", alpha=0.4)
    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_80k_top), marker='o', color="g", alpha=0.4)

    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_bottom), marker='^', color="r", alpha=0.4)
    ax4.scatter(range(1, len(bands_total) + 1), RV_cumulative(y4_100k_top), marker='o', color="r", alpha=0.4)

    ax4.text(1.0, positiony_max, "M9", size=14)

    # ticks and labels
    ax4.set_xlabel('Bands')
    ax4.set_xlim(0.5, len(bands_total) + 0.5)
    ax4.set_xticks(range(1, len(bands_total) + 1))
    ax4.set_xticklabels(bands_total)
    ax4.set_yticklabels([])
    ax4.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    y1_60k = [results["M0-Z-5.0-60k"], results["M0-Y-5.0-60k"], results["M0-J-5.0-60k"], results["M0-H-5.0-60k"],
              results["M0-K-5.0-60k"]]
    y1_60k_top = [y[1] for y in y1_60k]
    y1_60k_bottom = [y[2] for y in y1_60k]

    y1_80k = [results["M0-Z-5.0-80k"], results["M0-Y-5.0-80k"], results["M0-J-5.0-80k"], results["M0-H-5.0-80k"],
              results["M0-K-5.0-80k"]]
    y1_80k_top = [y[1] for y in y1_80k]
    y1_80k_bottom = [y[2] for y in y1_80k]

    y1_100k = [results["M0-Z-5.0-100k"], results["M0-Y-5.0-100k"], results["M0-J-5.0-100k"], results["M0-H-5.0-100k"],
               results["M0-K-5.0-100k"]]
    y1_100k_top = [y[1] for y in y1_100k]
    y1_100k_bottom = [y[2] for y in y1_100k]

    y1_60k_lim = [y[0] for y in y1_60k]

    y1_80k_lim = [y[0] for y in y1_80k]

    y1_100k_lim = [y[0] for y in y1_100k]

    # preparation of data: plot2
    y2_60k = [results["M3-Z-5.0-60k"], results["M3-Y-5.0-60k"], results["M3-J-5.0-60k"], results["M3-H-5.0-60k"],
              results["M3-K-5.0-60k"]]
    y2_60k_top = [y[1] for y in y2_60k]
    y2_60k_bottom = [y[2] for y in y2_60k]

    y2_80k = [results["M3-Z-5.0-80k"], results["M3-Y-5.0-80k"], results["M3-J-5.0-80k"], results["M3-H-5.0-80k"],
              results["M3-K-5.0-80k"]]
    y2_80k_top = [y[1] for y in y2_80k]
    y2_80k_bottom = [y[2] for y in y2_80k]

    y2_100k = [results["M3-Z-5.0-100k"], results["M3-Y-5.0-100k"], results["M3-J-5.0-100k"], results["M3-H-5.0-100k"],
               results["M3-K-5.0-100k"]]
    y2_100k_top = [y[1] for y in y2_100k]
    y2_100k_bottom = [y[2] for y in y2_100k]

    y2_60k_lim = [y[0] for y in y2_60k]

    y2_80k_lim = [y[0] for y in y2_80k]

    y2_100k_lim = [y[0] for y in y2_100k]

    # preparation of data: plot3
    y3_60k = [results["M6-Z-5.0-60k"], results["M6-Y-5.0-60k"], results["M6-J-5.0-60k"], results["M6-H-5.0-60k"],
              results["M6-K-5.0-60k"]]
    y3_60k_top = [y[1] for y in y3_60k]
    y3_60k_bottom = [y[2] for y in y3_60k]

    y3_80k = [results["M6-Z-5.0-80k"], results["M6-Y-5.0-80k"], results["M6-J-5.0-80k"], results["M6-H-5.0-80k"],
              results["M6-K-5.0-80k"]]
    y3_80k_top = [y[1] for y in y3_80k]
    y3_80k_bottom = [y[2] for y in y3_80k]

    y3_100k = [results["M6-Z-5.0-100k"], results["M6-Y-5.0-100k"], results["M6-J-5.0-100k"], results["M6-H-5.0-100k"],
               results["M6-K-5.0-100k"]]
    y3_100k_top = [y[1] for y in y3_100k]
    y3_100k_bottom = [y[2] for y in y3_100k]

    y3_60k_lim = [y[0] for y in y3_60k]

    y3_80k_lim = [y[0] for y in y3_80k]

    y3_100k_lim = [y[0] for y in y3_100k]

    # preparation of data: plot4
    y4_60k = [results["M9-Z-5.0-60k"], results["M9-Y-5.0-60k"], results["M9-J-5.0-60k"], results["M9-H-5.0-60k"],
              results["M9-K-5.0-60k"]]
    y4_60k_top = [y[1] for y in y4_60k]
    y4_60k_bottom = [y[2] for y in y4_60k]

    y4_80k = [results["M9-Z-5.0-80k"], results["M9-Y-5.0-80k"], results["M9-J-5.0-80k"], results["M9-H-5.0-80k"],
              results["M9-K-5.0-80k"]]
    y4_80k_top = [y[1] for y in y4_80k]
    y4_80k_bottom = [y[2] for y in y4_80k]

    y4_100k = [results["M9-Z-5.0-100k"], results["M9-Y-5.0-100k"], results["M9-J-5.0-100k"], results["M9-H-5.0-100k"],
               results["M9-K-5.0-100k"]]
    y4_100k_top = [y[1] for y in y4_100k]
    y4_100k_bottom = [y[2] for y in y4_100k]

    y4_60k_lim = [y[0] for y in y4_60k]

    y4_80k_lim = [y[0] for y in y4_80k]

    y4_100k_lim = [y[0] for y in y4_100k]

    positiony_max = np.max(
        [np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top),
         np.max(y1_100k_bottom),
         np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top),
         np.max(y2_100k_bottom),
         np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top),
         np.max(y3_100k_bottom),
         np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top),
         np.max(y4_100k_bottom)])

    positiony_min = np.min(
        [np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top),
         np.min(y1_100k_bottom),
         np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top),
         np.min(y2_100k_bottom),
         np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top),
         np.min(y3_100k_bottom),
         np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top),
         np.min(y4_100k_bottom)])

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
    ax1.set_xlim(0.5, len(bands) + 0.5)
    ax1.set_xticks(range(1, len(bands) + 1))
    ax1.set_xticklabels([])
    ax1.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax2.set_xlim(0.5, len(bands) + 0.5)
    ax2.set_xticks(range(1, len(bands) + 1))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax3.set_xlim(0.5, len(bands) + 0.5)
    ax3.set_xticks(range(1, len(bands) + 1))
    ax3.set_xticklabels(bands)
    ax3.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax4.set_xlim(0.5, len(bands) + 0.5)
    ax4.set_xticks(range(1, len(bands) + 1))
    ax4.set_xticklabels(bands)
    ax4.set_yticklabels([])
    ax4.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    y1_60k = [results["M0-Z-10.0-60k"], results["M0-Y-10.0-60k"], results["M0-J-10.0-60k"], results["M0-H-10.0-60k"],
              results["M0-K-10.0-60k"]]
    y1_60k_top = [y[1] for y in y1_60k]
    y1_60k_bottom = [y[2] for y in y1_60k]

    y1_80k = [results["M0-Z-10.0-80k"], results["M0-Y-10.0-80k"], results["M0-J-10.0-80k"], results["M0-H-10.0-80k"],
              results["M0-K-10.0-80k"]]
    y1_80k_top = [y[1] for y in y1_80k]
    y1_80k_bottom = [y[2] for y in y1_80k]

    y1_100k = [results["M0-Z-10.0-100k"], results["M0-Y-10.0-100k"], results["M0-J-10.0-100k"],
               results["M0-H-10.0-100k"], results["M0-K-10.0-100k"]]
    y1_100k_top = [y[1] for y in y1_100k]
    y1_100k_bottom = [y[2] for y in y1_100k]

    y1_60k_lim = [y[0] for y in y1_60k]

    y1_80k_lim = [y[0] for y in y1_80k]

    y1_100k_lim = [y[0] for y in y1_100k]

    # preparation of data: plot2
    y2_60k = [results["M3-Z-10.0-60k"], results["M3-Y-10.0-60k"], results["M3-J-10.0-60k"], results["M3-H-10.0-60k"],
              results["M3-K-10.0-60k"]]
    y2_60k_top = [y[1] for y in y2_60k]
    y2_60k_bottom = [y[2] for y in y2_60k]

    y2_80k = [results["M3-Z-10.0-80k"], results["M3-Y-10.0-80k"], results["M3-J-10.0-80k"], results["M3-H-10.0-80k"],
              results["M3-K-10.0-80k"]]
    y2_80k_top = [y[1] for y in y2_80k]
    y2_80k_bottom = [y[2] for y in y2_80k]

    y2_100k = [results["M3-Z-10.0-100k"], results["M3-Y-10.0-100k"], results["M3-J-10.0-100k"],
               results["M3-H-10.0-100k"], results["M3-K-10.0-100k"]]
    y2_100k_top = [y[1] for y in y2_100k]
    y2_100k_bottom = [y[2] for y in y2_100k]

    y2_60k_lim = [y[0] for y in y2_60k]

    y2_80k_lim = [y[0] for y in y2_80k]

    y2_100k_lim = [y[0] for y in y2_100k]

    # preparation of data: plot3
    y3_60k = [results["M6-Z-10.0-60k"], results["M6-Y-10.0-60k"], results["M6-J-10.0-60k"], results["M6-H-10.0-60k"],
              results["M6-K-10.0-60k"]]
    y3_60k_top = [y[1] for y in y3_60k]
    y3_60k_bottom = [y[2] for y in y3_60k]

    y3_80k = [results["M6-Z-10.0-80k"], results["M6-Y-10.0-80k"], results["M6-J-10.0-80k"], results["M6-H-10.0-80k"],
              results["M6-K-10.0-80k"]]
    y3_80k_top = [y[1] for y in y3_80k]
    y3_80k_bottom = [y[2] for y in y3_80k]

    y3_100k = [results["M6-Z-10.0-100k"], results["M6-Y-10.0-100k"], results["M6-J-10.0-100k"],
               results["M6-H-10.0-100k"], results["M6-K-10.0-100k"]]
    y3_100k_top = [y[1] for y in y3_100k]
    y3_100k_bottom = [y[2] for y in y3_100k]

    y3_60k_lim = [y[0] for y in y3_60k]

    y3_80k_lim = [y[0] for y in y3_80k]

    y3_100k_lim = [y[0] for y in y3_100k]

    # preparation of data: plot4
    y4_60k = [results["M9-Z-10.0-60k"], results["M9-Y-10.0-60k"], results["M9-J-10.0-60k"], results["M9-H-10.0-60k"],
              results["M9-K-10.0-60k"]]
    y4_60k_top = [y[1] for y in y4_60k]
    y4_60k_bottom = [y[2] for y in y4_60k]

    y4_80k = [results["M9-Z-10.0-80k"], results["M9-Y-10.0-80k"], results["M9-J-10.0-80k"], results["M9-H-10.0-80k"],
              results["M9-K-10.0-80k"]]
    y4_80k_top = [y[1] for y in y4_80k]
    y4_80k_bottom = [y[2] for y in y4_80k]

    y4_100k = [results["M9-Z-10.0-100k"], results["M9-Y-10.0-100k"], results["M9-J-10.0-100k"],
               results["M9-H-10.0-100k"], results["M9-K-10.0-100k"]]
    y4_100k_top = [y[1] for y in y4_100k]
    y4_100k_bottom = [y[2] for y in y4_100k]

    y4_60k_lim = [y[0] for y in y4_60k]

    y4_80k_lim = [y[0] for y in y4_80k]

    y4_100k_lim = [y[0] for y in y4_100k]

    positiony_max = np.max(
        [np.max(y1_60k_top), np.max(y1_60k_bottom), np.max(y1_80k_top), np.max(y1_80k_bottom), np.max(y1_100k_top),
         np.max(y1_100k_bottom),
         np.max(y2_60k_top), np.max(y2_60k_bottom), np.max(y2_80k_top), np.max(y2_80k_bottom), np.max(y2_100k_top),
         np.max(y2_100k_bottom),
         np.max(y3_60k_top), np.max(y3_60k_bottom), np.max(y3_80k_top), np.max(y3_80k_bottom), np.max(y3_100k_top),
         np.max(y3_100k_bottom),
         np.max(y4_60k_top), np.max(y4_60k_bottom), np.max(y4_80k_top), np.max(y4_80k_bottom), np.max(y4_100k_top),
         np.max(y4_100k_bottom)])

    positiony_min = np.min(
        [np.min(y1_60k_top), np.min(y1_60k_bottom), np.min(y1_80k_top), np.min(y1_80k_bottom), np.min(y1_100k_top),
         np.min(y1_100k_bottom),
         np.min(y2_60k_top), np.min(y2_60k_bottom), np.min(y2_80k_top), np.min(y2_80k_bottom), np.min(y2_100k_top),
         np.min(y2_100k_bottom),
         np.min(y3_60k_top), np.min(y3_60k_bottom), np.min(y3_80k_top), np.min(y3_80k_bottom), np.min(y3_100k_top),
         np.min(y3_100k_bottom),
         np.min(y4_60k_top), np.min(y4_60k_bottom), np.min(y4_80k_top), np.min(y4_80k_bottom), np.min(y4_100k_top),
         np.min(y4_100k_bottom)])

    # data for correction plot
    y1_60k_vsini10 = (np.array(y1_60k_top) - np.array(y1_60k_bottom)) / np.array(y1_60k_bottom) * 100.0
    y1_80k_vsini10 = (np.array(y1_80k_top) - np.array(y1_80k_bottom)) / np.array(y1_80k_bottom) * 100.0
    y1_100k_vsini10 = (np.array(y1_100k_top) - np.array(y1_100k_bottom)) / np.array(y1_100k_bottom) * 100.0

    y2_60k_vsini10 = (np.array(y2_60k_top) - np.array(y2_60k_bottom)) / np.array(y2_60k_bottom) * 100.0
    y2_80k_vsini10 = (np.array(y2_80k_top) - np.array(y2_80k_bottom)) / np.array(y2_80k_bottom) * 100.0
    y2_100k_vsini10 = (np.array(y2_100k_top) - np.array(y2_100k_bottom)) / np.array(y2_100k_bottom) * 100.0

    y3_60k_vsini10 = (np.array(y3_60k_top) - np.array(y3_60k_bottom)) / np.array(y3_60k_bottom) * 100.0
    y3_80k_vsini10 = (np.array(y3_80k_top) - np.array(y3_80k_bottom)) / np.array(y3_80k_bottom) * 100.0
    y3_100k_vsini10 = (np.array(y3_100k_top) - np.array(y3_100k_bottom)) / np.array(y3_100k_bottom) * 100.0

    y4_60k_vsini10 = (np.array(y4_60k_top) - np.array(y4_60k_bottom)) / np.array(y4_60k_bottom) * 100.0
    y4_80k_vsini10 = (np.array(y4_80k_top) - np.array(y4_80k_bottom)) / np.array(y4_80k_bottom) * 100.0
    y4_100k_vsini10 = (np.array(y4_100k_top) - np.array(y4_100k_bottom)) / np.array(y4_100k_bottom) * 100.0

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
    ax1.set_xlim(0.5, len(bands) + 0.5)
    ax1.set_xticks(range(1, len(bands) + 1))
    ax1.set_xticklabels([])
    ax1.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax2.set_xlim(0.5, len(bands) + 0.5)
    ax2.set_xticks(range(1, len(bands) + 1))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax3.set_xlim(0.5, len(bands) + 0.5)
    ax3.set_xticks(range(1, len(bands) + 1))
    ax3.set_xticklabels(bands)
    ax3.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

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
    ax4.set_xlim(0.5, len(bands) + 0.5)
    ax4.set_xticks(range(1, len(bands) + 1))
    ax4.set_xticklabels(bands)
    ax4.set_yticklabels([])
    ax4.set_ylim(positiony_min - 0.1 * (positiony_max - positiony_min),
                 positiony_max + 0.1 * (positiony_max - positiony_min))

    ax4.tick_params(axis='both', which='major', labelsize=15)

    fig.text(0.06, 0.5, r'Precision [m/s]', ha='center', va='center', rotation='vertical', size=14)
    fig.subplots_adjust(hspace=0, wspace=0)

    plt.show()
    plt.close()

    # correction figure
    positiony_max = np.max(
        [np.max(y1_60k_vsini1), np.max(y1_60k_vsini10), np.max(y1_80k_vsini1), np.max(y1_80k_vsini10),
         np.max(y1_100k_vsini1), np.max(y1_100k_vsini10),
         np.max(y2_60k_vsini1), np.max(y2_60k_vsini10), np.max(y2_80k_vsini1), np.max(y2_80k_vsini10),
         np.max(y2_100k_vsini1), np.max(y2_100k_vsini10),
         np.max(y3_60k_vsini1), np.max(y3_60k_vsini10), np.max(y3_80k_vsini1), np.max(y3_80k_vsini10),
         np.max(y3_100k_vsini1), np.max(y3_100k_vsini10),
         np.max(y4_60k_vsini1), np.max(y4_60k_vsini10), np.max(y4_80k_vsini1), np.max(y4_80k_vsini10),
         np.max(y4_100k_vsini1), np.max(y4_100k_vsini10)])

    positiony_min = np.min(
        [np.min(y1_60k_vsini1), np.min(y1_60k_vsini10), np.min(y1_80k_vsini1), np.min(y1_80k_vsini10),
         np.min(y1_100k_vsini1), np.min(y1_100k_vsini10),
         np.min(y2_60k_vsini1), np.min(y2_60k_vsini10), np.min(y2_80k_vsini1), np.min(y2_80k_vsini10),
         np.min(y2_100k_vsini1), np.min(y2_100k_vsini10),
         np.min(y3_60k_vsini1), np.min(y3_60k_vsini10), np.min(y3_80k_vsini1), np.min(y3_80k_vsini10),
         np.min(y3_100k_vsini1), np.min(y3_100k_vsini10),
         np.min(y4_60k_vsini1), np.min(y4_60k_vsini10), np.min(y4_80k_vsini1), np.min(y4_80k_vsini10),
         np.min(y4_100k_vsini1), np.min(y4_100k_vsini10)])

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

    ax1.text(1.0, positiony_max * 0.90, "M0", size=14)

    # ticks and labels
    ax1.set_xlim(0.5, len(bands) + 0.5)
    ax1.set_xticks(range(1, len(bands) + 1))
    ax1.set_xticklabels([])
    ax1.set_ylim(positiony_min - 0.05 * (positiony_max - positiony_min),
                 positiony_max + 0.05 * (positiony_max - positiony_min))

    ax2 = fig.add_subplot(222)

    ax2.scatter(range(1, len(bands) + 1), y2_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
    ax2.scatter(range(1, len(bands) + 1), y2_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
    ax2.scatter(range(1, len(bands) + 1), y2_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

    ax2.scatter(range(1, len(bands) + 1), y2_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
    ax2.scatter(range(1, len(bands) + 1), y2_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
    ax2.scatter(range(1, len(bands) + 1), y2_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

    ax2.text(1.0, positiony_max * 0.90, "M3", size=14)

    # ticks and labels
    ax2.set_xlim(0.5, len(bands) + 0.5)
    ax2.set_xticks(range(1, len(bands) + 1))
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.set_ylim(positiony_min - 0.05 * (positiony_max - positiony_min),
                 positiony_max + 0.05 * (positiony_max - positiony_min))

    ax3 = fig.add_subplot(223)

    ax3.scatter(range(1, len(bands) + 1), y3_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
    ax3.scatter(range(1, len(bands) + 1), y3_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
    ax3.scatter(range(1, len(bands) + 1), y3_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

    ax3.scatter(range(1, len(bands) + 1), y3_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
    ax3.scatter(range(1, len(bands) + 1), y3_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
    ax3.scatter(range(1, len(bands) + 1), y3_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

    ax3.text(1.0, positiony_max * 0.90, "M6", size=14)

    # ticks and labels
    ax3.set_xlabel('Bands')
    ax3.set_xlim(0.5, len(bands) + 0.5)
    ax3.set_xticks(range(1, len(bands) + 1))
    ax3.set_xticklabels(bands)
    ax3.set_ylim(positiony_min - 0.05 * (positiony_max - positiony_min),
                 positiony_max + 0.05 * (positiony_max - positiony_min))

    ax4 = fig.add_subplot(224)

    ax4.scatter(range(1, len(bands) + 1), y4_60k_vsini1, marker="<", s=100.0, color="b", alpha=0.2)
    ax4.scatter(range(1, len(bands) + 1), y4_80k_vsini1, marker="<", s=100.0, color="g", alpha=0.2)
    ax4.scatter(range(1, len(bands) + 1), y4_100k_vsini1, marker="<", s=100.0, color="r", alpha=0.2)

    ax4.scatter(range(1, len(bands) + 1), y4_60k_vsini10, marker=">", s=100.0, color="b", alpha=0.2)
    ax4.scatter(range(1, len(bands) + 1), y4_80k_vsini10, marker=">", s=100.0, color="g", alpha=0.2)
    ax4.scatter(range(1, len(bands) + 1), y4_100k_vsini10, marker=">", s=100.0, color="r", alpha=0.2)

    ax4.text(1.0, positiony_max * 0.90, "M9", size=14)

    # ticks and labels
    ax4.set_xlabel('Bands')
    ax4.set_xlim(0.5, len(bands) + 0.5)
    ax4.set_xticks(range(1, len(bands) + 1))
    ax4.set_xticklabels(bands)
    ax4.set_yticklabels([])
    ax4.set_ylim(positiony_min - 0.05 * (positiony_max - positiony_min),
                 positiony_max + 0.05 * (positiony_max - positiony_min))

    fig.text(0.06, 0.5, r'Loss in RV precision [\%]', ha='center', va='center', rotation='vertical', size=14)
    fig.subplots_adjust(hspace=0, wspace=0)

    plt.show()
    plt.close()

    # if the paper_plots option is on, then print the results
    print("\n")
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for res in R:
                    for smpl in sampling:
                        id_string = "{0:s}-{1:s}-{2:.01f}-{3:s}-{4:2.01f}".format(star, band, float(vel),
                                                                                  res, float(smpl))
                        precision = results[id_string]
                        print("{0!s}: & {1:.1f}\t & {2:.1f}\t & {3:.1f} \\\\".format(id_string, precision[0],
                                                                                     precision[1], precision[2]))
