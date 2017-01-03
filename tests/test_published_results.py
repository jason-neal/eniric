
from __future__ import division, print_function
import pytest
import numpy as np
import eniric.original_code.exorunner.nIRanalysis as exonIR
#import eniric.original_code.nIRanalysis as orgnIR
import eniric.original_code.exorunner.Qcalculator as exoQ
import eniric.original_code.Qcalculator as orgQ
import eniric.nIRanalysis as nIR
import eniric.Qcalculator as Q
# To test if the new code produces the same precision values on the published results.

# Test with and without fudge factor

path = "data/Published_Results/resampled/"


def test_RVprec_using_pd_load():
    """Test with data loaded with pandas """
    for vsini in [1, 5, 10]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(int(vsini))

        wav, flux = nIR.pdread_2col(path+name, noheader=True)
        wav = np.array(wav)
        flux = np.array(flux)

        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        assert exo_rv == org_rv
        assert new_rv == org_rv
        assert exo_rv == new_rv


def test_RVprec_using_exo_load():
    """Test  with data loaded with exorunner IOmodule  """
    for vsini in [1, 5, 10]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)

        wav, flux = exonIR.read_2col(path+name)
        wav = np.array(wav)
        flux = np.array(flux)

        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value
        assert exo_rv == org_rv
        assert new_rv == org_rv
        assert exo_rv == new_rv


def test_RVprec_using_scaling():
    """Test RVprec with pedros flux scaling value. """
    for vsini in [1, 5, 10]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)

        wav, flux = exonIR.read_2col(path+name)
        wav = np.array(wav)
        flux = np.array(flux)
        # scaling
        flux = flux / ((1.634e4)**2.0)
        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        assert exo_rv == org_rv
        assert new_rv == org_rv
        assert exo_rv == new_rv

def test_RVprec_equals_published_values_vsini1():
    """Test RVprec with pedros flux scaling value. """
    Published_Results = {1: 3.8, 5: 9.1, 10: 20.7}
    for vsini in [1]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)

        wav, flux = exonIR.read_2col(path+name)
        wav = np.array(wav)
        flux = np.array(flux)
        # scaling
        flux = flux / ((1.634e4)**2.0)
        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        this_result = Published_Results[vsini]
        assert np.round(exo_rv, 1) == this_result
        assert np.round(org_rv, 1) == this_result
        assert np.round(exo_rv, 1) == this_result

def test_RVprec_equals_published_values_vsini5():
    """Test RVprec with pedros flux scaling value. """
    Published_Results = {1: 3.8, 5: 9.1, 10: 20.7}
    for vsini in [5]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)

        wav, flux = exonIR.read_2col(path+name)
        wav = np.array(wav)
        flux = np.array(flux)
        # scaling
        flux = flux / ((1.634e4)**2.0)
        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        this_result = Published_Results[vsini]
        assert np.round(exo_rv, 1) == this_result
        assert np.round(org_rv, 1) == this_result
        assert np.round(exo_rv, 1) == this_result

def test_RVprec_equals_published_values_vsini10():
    """Test RVprec with pedros flux scaling value. """
    Published_Results = {1: 3.8, 5: 9.1, 10: 20.7}
    for vsini in [10]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)

        wav, flux = exonIR.read_2col(path+name)
        wav = np.array(wav)
        flux = np.array(flux)
        # scaling
        flux = flux / ((1.634e4)**2.0)
        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        this_result = Published_Results[vsini]
        assert np.round(exo_rv, 1) == this_result
        assert np.round(org_rv, 1) == this_result
        assert np.round(exo_rv, 1) == this_result
