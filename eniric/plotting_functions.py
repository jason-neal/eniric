
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
