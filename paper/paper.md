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

With recent high-precision spectrographs targeting RV precision at the 10cm/s level [@Pepe2014] in the quest to find smallest planets it is important to understand the theoretical precision attainable in the spectra of the host star.
Eniric, writen in Python3, provides a simple way to calculate the theoretical spectral quality and RV precision of stellar spectra.


Eniric calculates the fundamental photon noise RV precision as formulated in [@Connes1985] and [@bouchy_fundamental_2001] under the three conditions explored in [@Figueira2016].

 - Condition 1: Considering the full spectral information within each spectral band.

 - Condition 2: Pixels less that 30 kms\(^{-1}\) from telluric absorption lines deeper than 2% discarded.

 - Condition 3: Considering a perfect correction but with the photon noise contribution of the spectrum being amplified by the telluric correction.

By default eniric uses a telluric transmission spectrum obtained by averaging a year of weekly spectra as seen from La Sillia Observatory at an airmass of 1.2(*z*=35\(^o\)) simulated with TAPAS[@Bertaux2014]. This can however be replaced by the user.

2 situations:  quality (flux independent), precision of models, full, with masking of telluric lines, or with perfect telluric correction

Relative precision
Convolutions
Enable the exploration between synthetic libraries and observed spectra similarly to [@Artigua2018]

SNR scaling...
Accounts for stellar parameters, rotational and instrumental broadening, sampling rate, SNR level, (and doppler shift)


A script is available to calculate the relative precision of any synthetic spectra in the PHOENIX-ACES and BT-Settl (CIFIST2015)

Functions are available to calculate the precision for any PHOENIX-ACES library spectra.
Or can be applied to observed spectra.

Eniric is an improved version of the software used to calculate the RV precision of M-dwarfs in the NIR bands in [@figueira_radial_2016].
Expanded to explore the wider parameter space available in the full spectral library, and to other SNR's, bands, resolutions, and velocities.
It has been recently used to calculate RV precision of M-dwarf spectra for use in the Exposure Time Calculators of two new NIR spectrographs, NIRPS [@bouchy_nearinfrared_2017] and SPIRou [@artigau_spirou_2014]. For SPIRou these were specifically calculated at the instruments resolution (R=75,000) and for all stellar temperatures between 2500-4000~K.


![Precision achieved with *eniric* as a function of spectral band for stars with a rotational velocity of vsini=1.0 kms\(^{−1}\) and temperatures 3900 K, 3500 K, 2800 K, 2600 K, corresponding to spectral types M0, M3, M6, and M9 respectively.
The dashed line represents the theoretical limits imposed by condition 1, and the ﬁlled area represents the values within the limits set by conditions 2 (circles) and 3 (triangles); blue, green, and red represent the results obtained for resolutions of 60000, 80000, and 100000, respectively.
The spectra were normalized to have a S/N of 100 per resolution element as measured at the center of the J-band.
This is similar to Figure 1 from [@figueira_radial_2016] but with updated precision values.](./precisions.png)


## Acknowledgements

This work was supported by Funda\c{c}\~ao para a Ci\^encia e a Tecnologia (FCT) (Portugal) research grants through national funds and from FEDER through COMPETE2020 by the following grants: UID/FIS/04434/2013 & POCI--01--0145-FEDER--007672, PTDC/FIS-AST/1526/2014 & POCI--01--0145-FEDER--016886 & PTDC/FIS-AST/7073/2014 & POCI-01-0145-FEDER-016880.
<!--  -->
J.J.N. acknowledges support from FCT though the PhD::Space fellowship PD/BD/52700/2014.
