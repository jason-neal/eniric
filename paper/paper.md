---
title: 'Eniric: Extended NIR Information Content'
tags:
  - Python
  - astronomy
  - radial velocity precision
  - near-infrared
authors:
 - name: J.J. Neal
   orcid: 0000-0003-0390-9896
   affiliation: "1, 2" # (Multiple affiliations must be quoted)
 - name: P. Figueira
   orcid: 0000
   affiliation: "3, 1"


affiliations:
 - name: Instituto de Astrofísica e Ciências do Espaço, Universidade do Porto, CAUP, Rua das Estrelas, 4150-762, Porto, Portugal
   index: 1
 - name: Departamento de Física e Astronomia, Faculdade de Ciências, Universidade do Porto, Rua do Campo Alegre, 4169-007, Porto, Portugal
   index: 2
 - name: European Southern Observatory, Alonso de Córdova 3107, Vitacura, Casilla 19001, Santiago 19, Chile
   index: 3

date: XXXX 2018
bibliography: paper.bib
---

Eniric provides a simple way to calcualte the spectral quality and/or precsion of stellar spectra.

Functions are available to calcualte the precision for any PHOENIX-ACES library spectra.
Or can be apllied to  observed spectra.

The fundamental RV precision as provided by photon noise as formulated in [@Connes1985, @bouchy_fundamental_2001] is calculated for model spectra.

2 situations:  quality (flux independant), precision of models, full, with masking of telluric lines, or with perfect telluric correction

SNR scaling...
Accounts for stellar parameters, rotational and instrumental broadening, sampling rate, SNR level, (and doppler shift)

This code has been used to calcualte the RV precision of M-dwarfs in the NIR bands to access [@figueira_radial_2016], and its results have provided for use in the  Exposure Time Calculators of two new NIR spectrographs, NIRPS [@bouchy_nearinfrared_2017] and SPIRou [@artigau_spirou_2014] . 


## Acknowledgements

This work was supported by Funda\c{c}\~ao para a Ci\^encia e a Tecnologia (FCT) (Portugal) research grants through national funds and from FEDER through COMPETE2020 by the following grants: UID/FIS/04434/2013 & POCI--01--0145-FEDER--007672, PTDC/FIS-AST/1526/2014 & POCI--01--0145-FEDER--016886 & PTDC/FIS-AST/7073/2014 & POCI-01-0145-FEDER-016880.
<!--  -->
J.J.N. acknowledges support from FCT though the PhD::Space fellowship PD/BD/52700/2014.

  
