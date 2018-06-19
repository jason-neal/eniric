
# coding: utf-8

# ## Testing barycentric shift in eniric
# 

# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
# to remove labels in one tick
from matplotlib.ticker import MaxNLocator

import eniric.IOmodule as IOmodule
import eniric.Qcalculator as Qcalculator
from eniric.plotting_functions import (plot_atmosphere_model,
                                       plot_stellar_spectum)
from eniric.utilities import band_selector, wav_selector

# In[ ]:


def prepare_atmosphere():
    """ Read in atmopheric model and prepare. """
    wav_atm, flux_atm, std_flux_atm, mask_atm = IOmodule.pdread_4col(atmmodel)
    # pandas lready returns numpy arrays
    wav_atm = wav_atm / 1000.0  # conversion from nanometers to micrometers
    mask_atm = np.array(mask_atm, dtype=bool)
    return wav_atm, flux_atm, std_flux_atm, mask_atm


def barycenter_shift(wav_atm, mask_atm, offset_RV=0.0):
    """ Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentic motion of the earth.
    
    offset_RV: float
        Radial veloctiy offset in km/s.
    """
    RV_bary = 3.0e4   # 30 km/s barycentric velocity
    offset_RV = offset_RV * 1.0e3     # Trun offset into m/s
    pixels_start = np.sum(mask_atm)
    pixels_total = len(mask_atm)

    mask_atm_30kms = []
    for wav_val, mask_val in zip(wav_atm, mask_atm):
        if (mask_val is False and offset_RV == 666.0):    # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(mask_val)

        else:
            delta_lambda = wav_val * RV_bary / Qcalculator.c.value    # doppler shift amount   delta_lambda = lambda*v/c
            
            starting_lambda = wav_val * offset_RV / Qcalculator.c.value   # offset value
            
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda+wav_val-delta_lambda, starting_lambda+wav_val+delta_lambda])
            
            # Replaces any values outside len(wav_atm) with end index value
            indexes_30kmslice = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in indexes_30kmslice]
            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]], dtype=bool)    # selecting only the slice in question

            # if(False in mask_atm_30kmslice):
            #    mask_atm_30kms.append(False)
            # else:
            #    mask_atm_30kms.append(True)

            mask_atm_30kmslice_reversed = [not i for i in mask_atm_30kmslice]    # comp list of array is bad
            print("mask_atm_30kmslice_reversed", mask_atm_30kmslice_reversed)
            
            clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice_reversed))[0]+1)[::2]
            print("clump", clump)
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


# In[ ]:


def inner_function():
    pass

def new_inner_function():
    pass


# In[ ]:


atmmodel = "../data/atmmodel/Average_TAPAS_2014_KJ.txt"
# w_min = 2.0
# w_max = 2.1

wav_atm_org, flux_atm_org, std_flux_atm_org, mask_atm_org = prepare_atmosphere()

# wav_atm, flux_atm = wav_selector(wav_atm_org, flux_atm_org, w_min, w_max)
# __, std_flux_atm = wav_selector(wav_atm_org, std_flux_atm_org, w_min, w_max)
# __, mask_atm = wav_selector(wav_atm_org, mask_atm_org, w_min, w_max)


# In[ ]:


print(wav_atm_org)
print(mask_atm_org)
print(sum(mask_atm_org))


# In[ ]:


new_mask_atm = barycenter_shift(wav_atm_org, mask_atm_org, offset_RV=0.0)


# In[ ]:


pixels_start = np.sum(mask_atm)
pixels_total = len(mask_atm)
print(pixels_start, pixels_total)

offset_RV=0

mask_atm_30kms = []
for wav_val, mask_val in zip(wav_atm, mask_atm):
        if (mask_val is False and offset_RV == 666.0):    # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(mask_val)

        else:
            delta_lambda = wav_val * 3.0e4 / Qcalculator.c.value
            print("delta_lambda", delta_lambda)
            
            starting_lambda = wav_val * offset_RV * 1.0e3 / Qcalculator.c.value
            print("starting_lambda", starting_lambda)
            
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda+wav_val-delta_lambda, starting_lambda+wav_val+delta_lambda])
            
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
