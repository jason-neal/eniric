# ENIRIC - Extended Near InfraRed Information Content
Analysis of near infrared spectra information content.

## Purpose
To analysis which spectroscopic bands contain the most information for radial velocity measurements. 
Model spectra are used to analyse the information conent of different bands.
They undergo two convolutions, one for rotational broadining and one for instrumental broadening.
The slope of the spectra are used as a proxy for the radial velocity precision attainable.

The purpose of this work is to:
- Extend the previous analysis to a range of different metallicity values.


## Background
The origin of this code was used in [this paper](https://arxiv.org/abs/1511.07468).
    
    P. Figueira, V. Zh. Adibekyan, M. Oshagh, J. J. Neal, B. Rojas-Ayala, C. Lovis, C. Melo, F. Pepe, N. C. Santos, M. Tsantaki, 2016,
    Radial velocity information content of M dwarf spectra in the near-infrared,
    Astronomy and Astrophysics, 586, A101

It has a number of effecincy problems which need to be improved upon before the new analysis is performed.

1) Use numpy mapping slicing instead of comprehension lists.  (~>250 times faster)
2) Use joblib to parallelize the convolutions.


## Runtime results:
Comparing the same calculation perfromed between the old and new code after a series of changes.
...
