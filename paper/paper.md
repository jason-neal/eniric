---
title: 'Eniric: Extended NIR Information Content'

tags:
  - Python
  - Astronomy
  - Radial velocity
  - Near-infrared

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

*Eniric* is an improved version of the software used to calculate the RV precision of M-dwarfs in the NIR bands in [@figueira_radial_2016]. It has been expanded to be user configurable to allow precisions calculated for different SNR's, spectral bands, resolutions, and rotational velocities as well as extending the available spectra to all spectra in the PHOENIX-ACES  [@Husser2013] and BT-Settl [@Barraffe2015] synthetic libraries.

*Eniric* calculates the fundamental photon noise RV precision as formulated in [@Connes1985] and [@bouchy_fundamental_2001].

The RV precision is defined as:

    RV_{rms} = c / Q \sqrt{Ne} = c / {\Sigma(W_i)}

where were c is the speed of light, Q is the spectral quality factor, Ne is the number of photoelectrons counted and \(W_i\) are the optimal pixel weigths given by:

    W_i = \lambda_i^2 (dA/d\lambda_i)^2 / A

Different scenarios can be explored by multiplication of the pixel weigths with a masking function \(M\), which can be defined by the user.
 
    W_i = W_i * M
    
As an example the three scenarios explored in [@Figueira2016] are applied with the following weights:

 - Condition 1: Considering the full spectral information within each spectral band. \(M_1 = 1\)

 - Condition 2: Pixels less that 30 kms\(^{-1}\) from telluric absorption lines deeper than 2% discarded.  \(M_2 = 0 (T < 0.98), 1 (T >= 0.98)\)

 - Condition 3: Considering a perfect correction but with the photon noise contribution of the spectrum being amplified by the telluric correction. \(M_3 = T^2\)

By default *eniric* contains and uses a telluric transmission spectrum \(T\) obtained by averaging a year of weekly spectra as seen from La Sillia Observatory at an airmass of 1.2(*z*=35\(^o\)) simulated with TAPAS[@Bertaux2014]. This can however be replaced by the user if desired.

*Eniric* can preform rotational and instrumental broadening of spectra through convolution with a rotational kernel (Grey ...) and gaussian kernel respectively. 
Both kernels are wavelength dependant and do **not** require a uniformly spaced wavelength vector, unlike the convolution functions given in PyAstronomy.
*Eniric* utilizes the **embarisingly parallel** nature of the convolution (each pixel can be calculated independently of its neighbours) to compute the convolutions in parallel; the convolution results are also cached using [Joblib](https://joblib.readthedocs.io/en/latest/) to avoid recomputation. 
This improves the convolution performance but not to the level achievable by algorithms that require an equal wavelength spacing and use fixed kernels (only valid for small wavelength regions), i.e. the fast convolutions provided in [PyAstronomy](https://github.com/sczesla/PyAstronomy).

**Eniric** enables the relative precision between synthetic spectra by allowing for normalization to a user defined SNR per pixel at a specific wavelength. The default choice is a SNR of 100 at the center of the J-band (1.25 \mu m) as done in [@Figueria2016] but is user definable.

As well as calculating precsion of the spectoscopic bands, the precision on smaller wavelength slices are possible. 
This allows for the exploration of varying precision with wavelength and for the comparision between synthetic libraries and observed spectra similarly to [@Artigua2018].

A script is provided to easily calculate the relative precisions of synthetic spectra with a large choice of parameters,
It provides the results of each parameter combination in a tabular format for further analysis.
The rv precision function can be used outside of this script to also preform the same calculations on other spectra defined by the user. 

*Eniric* has been recently used to provide relative RV precision of M-dwarf spectra for use in the Exposure Time Calculators of NIRPS [@bouchy_nearinfrared_2017] and SPIRou [@artigau_spirou_2014], two new high-resolution NIR spectrographs. For SPIRou these were specifically calculated for the instruments resolution (R=75,000) and for all stellar temperatures between 2500-4000~K. In both instances the SNR normalization was preformed relative to each band, rather than only the J-band. 


![Precision achieved with *eniric* as a function of spectral band for stars with a rotational velocity of vsini=1.0 kms\(^{−1}\) and temperatures 3900 K, 3500 K, 2800 K, 2600 K, corresponding to spectral types M0, M3, M6, and M9 respectively.
The dashed line represents the theoretical limits imposed by condition 1, and the ﬁlled area represents the values within the limits set by conditions 2 (circles) and 3 (triangles); blue, green, and red represent the results obtained for resolutions of 60000, 80000, and 100000, respectively.
The spectra were normalized to have a S/N of 100 per resolution element as measured at the center of the J-band.
This is similar to Figure 1 from [@figueira_radial_2016] but with updated precision values.](./precisions.png)


## Acknowledgements

This work was supported by Funda\c{c}\~ao para a Ci\^encia e a Tecnologia (FCT) (Portugal) research grants through national funds and from FEDER through COMPETE2020 by the following grants: UID/FIS/04434/2013 & POCI--01--0145-FEDER--007672, PTDC/FIS-AST/1526/2014 & POCI--01--0145-FEDER--016886 & PTDC/FIS-AST/7073/2014 & POCI-01-0145-FEDER-016880.
<!--  -->
J.J.N. acknowledges support from FCT though the PhD::Space fellowship PD/BD/52700/2014.
**P.F acknowledges ....**