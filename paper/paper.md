---
title: 'Eniric: Extended NIR Information Content'

tags:
  - Python
  - Astronomy
  - Radial velocity precision
  - Near-infrared
  - Spectral quality

authors:
 - name: J.J. Neal
   orcid: 0000-0003-0390-9896
   affiliation: "1, 2" # (Multiple affiliations must be quoted)
 - name: P. Figueira
   orcid: 0000-0001-8504-283X
   affiliation: "3, 1"

affiliations:
 - name: Instituto de Astrofísica e Ciências do Espaço, Universidade do Porto, CAUP, Rua das Estrelas, 4150-762, Porto, Portugal
   index: 1
 - name: Departamento de Física e Astronomia, Faculdade de Ciências, Universidade do Porto, Rua do Campo Alegre, 4169-007, Porto, Portugal
   index: 2
 - name: European Southern Observatory, Alonso de Córdova 3107, Vitacura, Casilla 19001, Santiago 19, Chile
   index: 3

date: September 2018
bibliography: paper.bib
---

With recent high-precision spectrographs targeting radial velocity (RV) precision at the 10\ cms$^{-1}$ level\ [e.g. @pepe_espresso_2014] in the quest to find smallest planets, it is important to understand the theoretical precision attainable in stellar spectra.
*Eniric* provides a simple way to calculate the theoretical spectral quality
and RV precision (i.e., information content) of synthetic and/or observed stellar spectra given vectors of wavelength and photon flux.

Written in *Python 3*, *Eniric* calculates the fundamental photon noise RV precision as formulated in @connes_absolute_1985 and @bouchy_fundamental_2001.
It is an improved version of the software used in @figueira_radial_2016 for calculating the RV precision of synthetic M-dwarf spectra in the near-infrared (NIR) bands.
The code was refactored, with hard-coded constraints removed, making it faster and simpler to explore a larger combination of parameters (e.g. not limited to M-dwarfs and NIR wavelengths).

*Eniric* contains several independent functions to transform observed and synthetic spectra, such as wavelength selection,
broadening, SNR normalization and to compute RV precisions.

*Eniric* performs rotational and instrumental broadening of spectra through convolution with a rotational kernel [@gray_observation_2005] and gaussian kernel respectively.
Both kernels are wavelength dependant and do **not** require a uniformly spaced wavelength vector, unlike the convolution functions given in [PyAstronomy](https://github.com/sczesla/PyAstronomy).
*Eniric* utilizes the *embarrassingly parallel* nature of the convolutions (each pixel can be calculated independently of its neighbours) to compute the convolutions in parallel; the convolution results are also cached using [Joblib](https://joblib.readthedocs.io/en/latest/) to avoid re-computation.
This improves the convolution performance but not to the level achievable by algorithms that require an equal wavelength spacing and use fixed kernels (only valid for small wavelength regions), e.g. the "fast" convolutions provided in [PyAstronomy](https://github.com/sczesla/PyAstronomy).

*Eniric* enables the relative precision between synthetic spectra by allowing for normalization to a user defined signal-to-noise ratio (SNR) per pixel at a specific wavelength.
 Although user definable the default choice is a SNR of 100 at the center of the J-band (1.25\ $\mu$m) as used in @figueira_radial_2016.

The precision calculations are not limited to the large spectroscopic bands, but can also be performed on narrow wavelength slices along the entire spectrum.
This allows one to explore the RV precision across the entire spectrum and perform comparision between observations and synthetic libraries [e.g. @artigau_optical_2018].

Extraneous information not included in the spectra (i.e., not photon noise nor line content information)
can be included in the precision calculation through the use of a spectral mask. This mask can be used to
indicate which spectral lines are to be included/excluded (via a binary mask) or if some spectral lines should receive more
statistical weight for an external reason.
For example, masks derived from an atmospheric absorption spectrum can be used to explore the treatment
and correction the atmospheric absorption on the RV precision [@figueira_radial_2016].
A typical treatment for the Earths atmosphere is the exclusion (via masking) of any region within 30\ kms$^{-1}$ of an atmospheric absorption line deeper than 2%, to account for Earths barycentric motion.

The script [`phoenix_precision.py`](https://github.com/jason-neal/eniric/blob/develop/scripts/phoenix_precision.py) is provided to easily compute relative RV precisions for synthetic spectra from the
[*PHOENIX-ACES*](http://phoenix.astro.physik.uni-goettingen.de) [@husser_new_2013] and
[*BT-Settl* (CIFIST2011-2015)](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/) [@baraffe_new_2015] libraries.
The [Grid Tools](https://iancze.github.io/Starfish/current/grid_tools.html) module from *Starfish* [@czekala_constructing_2015] is used to load in the library spectra given the identifying stellar parameters i.e, temperature, logg, metalicity and alpha.
The RV precision is computed under the three separate conditions of @figueira_radial_2016.
These are: the full spectrum, the exclusion of regions within 30\ kms$^{-1}$ of deep telluric lines (binary mask), and the full spectrum assuming a ``perfect'' telluric correction in which the spectrum is recovered but with an increased variance dependant on the telluric line depth (weighted mask).
The results are tabulated for all combinations of
the spectral parameters, SNR, instrument resolutions, rotational velocities, pixel sampling, and
wavelength choices provided to the script. An example of relative precision results is shown in Figure 1.

*Eniric* has been recently used to calculate relative RV precision of M-dwarf spectra for use in the Exposure Time Calculators of *NIRPS*\ [@bouchy_nearinfrared_2017] and *SPIRou*\ [@artigau_spirou_2014], two new high-resolution NIR spectrographs.
For *SPIRou* these were specifically calculated for the instruments resolution (R=75,000) and for all stellar temperatures between 2500--4000\ K.
It is also currently being used to compare the theoretical precision between observed *CARMENES* spectra\ [@reiners_carmenes_2018] and their synthetic library counterparts.

![Precision achieved with *Eniric* as a function of spectral band for stars with a rotational velocity of vsini=1.0\ kms$^{-1}$ and temperatures 3900 K, 3500 K, 2800 K, 2600 K, corresponding to spectral types M0, M3, M6, and M9 respectively.
The dashed line represents the theoretical limits imposed by the full spectrum, and the filled area represents the values within the limits set by the complete removal of atmospheric absorption regions (circles) and the perfect correction of atmospheric absorption (triangles); blue, green, and red represent the results obtained for resolutions of 60000, 80000, and 100000, respectively.
The spectra were normalized to have a SNR of 100 per resolution element as measured at the center of the J-band.
This is similar to Figure 1 from [@figueira_radial_2016] but with updated precision values.](./precisions.png)

[*Eniric*](https://github.com/jason-neal/eniric) is available on [*Github*](https://github.com/jason-neal/eniric) with documentation found in the `README.md` and at [ReadtheDocs](https://eniric.readthedocs.io/en/latest/) with usage examples provided as [Jupyter notebooks](https://github.com/jason-neal/eniric/tree/master/docs/Notebooks).
It utilizes packages from the scientific Python stack including [NumPy](http://www.numpy.org/) [@oliphant_guide_2015] and [SciPy](https://www.scipy.org/) [@scipy_scipy.org_2019], [Matplotlib](https://matplotlib.org/) [@hunter_matplotlib_2007], [Pandas](http://pandas.pydata.org/) [@mckinney_data_2010], and [Astropy](http://docs.astropy.org/en/stable/) [@astropy_collaboration_astropy_2013,@astropy_collaboration_astropy_2018]. It also uses [Joblib](https://joblib.readthedocs.io/en/latest/)[@joblib_joblib_2019], [Starfish](https://starfish.readthedocs.io/en/latest/) [@czekala_constructing_2015], and [tqdm](https://tqdm.github.io/) [@tqdm/tqdm].

[Comment]: # ([*PyPI*](https://pypi.org/project/eniric/) and is installable with *pip*.)

## Acknowledgements

This work was supported by Fundação para a Ciência e Tecnologia (FCT) (Portugal) research grants
through national funds and from FEDER (Fundo Europeu de Desenvolvimento Regional) through COMPETE2020 by the following grants:
POCI--01--0145-FEDER--007672, POCI--01--0145-FEDER--016880, POCI--01--0145-FEDER--016886, POCI-01-0145-FEDER-028953, POCI-01-0145-FEDER-032113, PTDC/FIS-AST/1526/2014, PTDC/FIS-AST/7073/2014, & UID/FIS/04434/2013.
J.J.N. acknowledges support from FCT though the PhD::Space fellowship PD/BD/52700/2014.


## References
