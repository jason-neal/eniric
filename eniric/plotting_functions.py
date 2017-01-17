
"""
Plotting functions for nIR_precision.

"""
import numpy as np
from sys import exit
import matplotlib.pyplot as plt

# to remove labels in one tick
from matplotlib.ticker import MaxNLocator
from eniric.utilities import band_selector

from matplotlib import rc
# set stuff for latex usage
rc('text', usetex=True)


def plot_atmopshere_model(wav_atm, flux_atm, mask_atm):
    """ Plot atmospheric transmission model for each band.

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
    exit(0)

def plot_stellar_spectum(wav_stellar, flux_stellar, wav_atm_selected, mask_atm_selected):
    # Plot the stellar spectrum as considered
    selected_transmission_stellar = wav_atm_selected[mask_atm_selected]

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux_stellar [ ] ")
    plt.plot(wav_stellar, flux_stellar, color='k')
    plt.vlines(selected_transmission_stellar, np.min(flux_stellar),
               0.3*np.max(flux_stellar), colors="b")
    plt.xlim(wav_stellar[0], wav_stellar[-1])
    plt.ylim(np.min(flux_stellar) - 0.1*(np.max(flux_stellar) - np.min(flux_stellar)),
             np.max(flux_stellar) + 0.1*(np.max(flux_stellar) - np.min(flux_stellar)))
    # plt.legend(loc='best')
    plt.show()
    plt.close()



def plot_nIR_flux():
    """Plot_flux code from nIR_precision"""
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

