# ENIRIC - Extended Near InfraRed Information Content
Analysis of near infrared spectra information content.

[![Build Status](https://travis-ci.org/jason-neal/eniric.svg?branch=master)](https://travis-ci.org/jason-neal/eniric)[![Coverage Status](https://coveralls.io/repos/github/jason-neal/eniric/badge.svg?branch=master)](https://coveralls.io/github/jason-neal/eniric?branch=master)[![Code Climate](https://codeclimate.com/github/jason-neal/eniric/badges/gpa.svg)](https://codeclimate.com/github/jason-neal/eniric)[![Code Issues](https://www.quantifiedcode.com/api/v1/project/3ed0e4393e4c40498f882855a77636c7/badge.svg)](https://www.quantifiedcode.com/app/project/3ed0e4393e4c40498f882855a77636c7)[![Issue Count](https://codeclimate.com/github/jason-neal/eniric/badges/issue_count.svg)](https://codeclimate.com/github/jason-neal/eniric)[![Test Coverage](https://codeclimate.com/github/jason-neal/eniric/badges/coverage.svg)](https://codeclimate.com/github/jason-neal/eniric/coverage)

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
On laptop after replacing the comprehension lists:

    new convolution = 27 seconds
    old convolution = 1hr 22 minutes
Rediculous!


## Bugs:
A number of bugs were found when improving this code. Manily affecting the condition invovling telluric line masking. This alters the RV in this column, sometimes significantly. This, however, does **NOT** alter the conclusions in the published paper. 
