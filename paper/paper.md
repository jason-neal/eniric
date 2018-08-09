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

Eniric provides a simple way to calculate the spectral quality and/or precsion of stellar spectra.

Functions are available to calculate the precision for any PHOENIX-ACES library spectra.
Or can be applied to  observed spectra.

The fundamental RV precision as provided by photon noise as formulated in [@Connes1985, @bouchy_fundamental_2001] is calculated for model spectra.

2 situations:  quality (flux independent), precision of models, full, with masking of telluric lines, or with perfect telluric correction

SNR scaling...
Accounts for stellar parameters, rotational and instrumental broadening, sampling rate, SNR level, (and doppler shift)

This code has been used to calculate the RV precision of M-dwarfs in the NIR bands to access [@figueira_radial_2016], and its results have provided for use in the  Exposure Time Calculators of two new NIR spectrographs, NIRPS [@bouchy_nearinfrared_2017] and SPIRou [@artigau_spirou_2014] . 


![Precision achieved with *eniric* as a function of spectral band for stars with a rotational velocity of vsini=1.0 kms\(^{−1}\) and temperatures 3900 K, 3500 K, 2800 K, 2600 K, corresponding to spectral types M0, M3, M6, and M9 respectively.
The dashed line represents the theoretical limits imposed by condition 1, and the ﬁlled area represents the values within the limits set by conditions 2 (circles) and 3 (triangles); blue, green, and red represent the results obtained for resolutions of 60000, 80000, and 100000, respectively. 
The spectra were normalized to have a S/N of 100 per resolution element as measured at the center of the J-band.
This is similar to Figure 1 from [@figueira_radial_2016] but with updated precision values.](./precisions.png)


## Acknowledgements

This work was supported by Funda\c{c}\~ao para a Ci\^encia e a Tecnologia (FCT) (Portugal) research grants through national funds and from FEDER through COMPETE2020 by the following grants: UID/FIS/04434/2013 & POCI--01--0145-FEDER--007672, PTDC/FIS-AST/1526/2014 & POCI--01--0145-FEDER--016886 & PTDC/FIS-AST/7073/2014 & POCI-01-0145-FEDER-016880.
<!--  -->
J.J.N. acknowledges support from FCT though the PhD::Space fellowship PD/BD/52700/2014.
