
from __future__ import division, print_function
import pytest
import numpy as np
import eniric.original_code.exorunner.nIRanalysis as exonIR
# import eniric.original_code.nIRanalysis as orgnIR
import eniric.original_code.exorunner.Qcalculator as exoQ
import eniric.original_code.Qcalculator as orgQ
# import eniric.nIRanalysis as nIR
from eniric.IOmodule import pdread_2col
import eniric.Qcalculator as Q
from bin.prec_1 import calc_prec1
# To test if the new code produces the same precision values on the published results.

# Test with and without fudge factor
# For python2.X compatibility
file_error_to_catch = getattr(__builtins__,'FileNotFoundError', IOError)

path = "data/Published_Results/resampled/"

@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_RVprec_using_pd_load():
    """Test with data loaded with pandas """
    for vsini in [1, 5, 10]:
        name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(int(vsini))

        wav, flux = pdread_2col(path+name, noheader=True)
        wav = np.array(wav)
        flux = np.array(flux)

        exo_rv = exoQ.RVprec_calc(wav, flux)
        org_rv = orgQ.RVprec_calc(wav, flux)
        new_rv = Q.RVprec_calc(wav, flux).value

        assert exo_rv == org_rv
        assert new_rv == org_rv
        assert exo_rv == new_rv

@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
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


@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_presicion_1():
    """ New precision 1 test that works."""
    published_results = {1: 3.8, 5: 9.1, 10: 20.7}
    path = "data/resampled/"
    for vsini in [1, 5, 10]:
        #name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{}.0_R100k_res3.txt".format(vsini)
        __, p1 = calc_prec1("M0", "Y", vsini, "100k", 3, resampled_dir=path)

        assert np.round(p1,1).value == published_results[vsini]
