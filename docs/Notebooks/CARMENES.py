
# coding: utf-8

# # Table of Contents
# * [Explore CARMENES Spectrum](#Explore-CARMENES-Spectrum)
# 

# # Explore CARMENES Spectrum

# - Explore how to load CARMENES Spectra
# - Splitting/verse correcting Telluric lines
# - Try extract some precisions and compare to Reiners 2017 [1](#cite-reiners_carmenes_2017)
# 
# Using Barnard's star so can also compare to Artigau 2018 [2](#cite-artigau_optical_2018)

# In[ ]:


<!--bibtex
@article{bouchy_fundamental_2001,
  title = {Fundamental Photon Noise Limit to Radial Velocity Measurements},
  volume = {374},
  issn = {0004-6361, 1432-0756},
  doi = {10.1051/0004-6361:20010730},
  number = {2},
  journal = {Astronomy and Astrophysics},
  author = {Bouchy, F. and Pepe, F. and Queloz, D.},
  month = aug,
  year = {2001},
  pages = {733-739},
  file = {/home/jneal/.mozilla/firefox/2rcitvnq.default/zotero/storage/UGSFCJ25/aa1316.pdf},
  note = {00258}
}

@article{artigau_optical_2018,
  archivePrefix = {arXiv},
  eprinttype = {arxiv},
  eprint = {1803.07646},
  title = {Optical and {{Near}}-{{Infrared Radial Velocity Content}} of {{M Dwarfs}}: {{Testing Models}} with {{Barnard}}'s {{Star}}},
  volume = {155},
  issn = {1538-3881},
  shorttitle = {Optical and {{Near}}-{{Infrared Radial Velocity Content}} of {{M Dwarfs}}},
  doi = {10/gddxj4},
  number = {5},
  journal = {The Astronomical Journal},
  author = {Artigau, {\'E}tienne and Malo, Lison and Doyon, Ren{\'e} and Figueira, Pedro and Delfosse, Xavier and Astudillo-Defru, Nicola},
  month = apr,
  year = {2018},
  keywords = {Astrophysics - Instrumentation and Methods for Astrophysics,Astrophysics - Solar and Stellar Astrophysics},
  pages = {198},
  file = {/home/jneal/.mozilla/firefox/2rcitvnq.default/zotero/storage/7FULJUP6/Artigau et al_2018_Optical and Near-Infrared Radial Velocity Content of M Dwarfs.pdf;/home/jneal/.mozilla/firefox/2rcitvnq.default/zotero/storage/VXVUZDVW/1803.html}
}

@article{reiners_carmenes_2017,
  title = {The {{CARMENES}} Search for Exoplanets around {{M}} Dwarfs: {{High}}-Resolution Optical and near-Infrared Spectroscopy of 324 Survey Stars},
  volume = {1711},
  shorttitle = {The {{CARMENES}} Search for Exoplanets around {{M}} Dwarfs},
  journal = {ArXiv e-prints},
  author = {Reiners, A. and Zechmeister, M. and Caballero, J. A. and Ribas, I. and Morales, J. C. and Jeffers, S. V. and Sch{\"o}fer, P. and Tal-Or, L. and Quirrenbach, A. and Amado, P. J. and Kaminski, A. and Seifert, W. and Abril, M. and Aceituno, J. and Alonso-Floriano, F. J. and Ammler-von Eiff, M. and Antona, R. and Anglada-Escud{\'e}, G. and Anwand-Heerwart, H. and Arroyo-Torres, B. and Azzaro, M. and Baroch, D. and Barrado, D. and Bauer, F. F. and Becerril, S. and B{\'e}jar, V. J. S. and Ben{\'\i}tez, D. and Berdi{\~n}as, Z. M. and Bergond, G. and Bl{\"u}mcke, M. and Brinkm{\"o}ller, M. and {del Burgo}, C. and Cano, J. and C{\'a}rdenas V{\'a}zquez, M. C. and Casal, E. and Cifuentes, C. and Claret, A. and Colom{\'e}, J. and Cort{\'e}s-Contreras, M. and Czesla, S. and D{\'\i}ez-Alonso, E. and Dreizler, S. and Feiz, C. and Fern{\'a}ndez, M. and Ferro, I. M. and Fuhrmeister, B. and Galad{\'\i}-Enr{\'\i}quez, D. and Garcia-Piquer, A. and Garc{\'\i}a Vargas, M. L. and Gesa, L. and G{\'o}mez, V. and {Galera} and Gonz{\'a}lez Hern{\'a}ndez, J. I. and Gonz{\'a}lez-Peinado, R. and Gr{\"o}zinger, U. and Grohnert, S. and Gu{\`a}rdia, J. and Guenther, E. W. and Guijarro, A. and {de Guindos}, E. and Guti{\'e}rrez-Soto, J. and Hagen, H.-J. and Hatzes, A. P. and Hauschildt, P. H. and Hedrosa, R. P. and Helmling, J. and Henning, Th. and Hermelo, I. and Hern{\'a}ndez Arab{\'\i}, R. and Hern{\'a}ndez Casta{\~n}o, L. and Hern{\'a}ndez Hernando, F. and Herrero, E. and Huber, A. and Huke, P. and Johnson, E. and {de Juan}, E. and Kim, M. and Klein, R. and Kl{\"u}ter, J. and Klutsch, A. and K{\"u}rster, M. and Lafarga, M. and Lamert, A. and Lamp{\'o}n, M. and Lara, L. M. and Laun, W. and Lemke, U. and Lenzen, R. and Launhardt, R. and {L{\'o}pez del Fresno}, M. and L{\'o}pez-Gonz{\'a}lez, J. and L{\'o}pez-Puertas, M. and L{\'o}pez Salas, J. F. and L{\'o}pez-Santiago, J. and Luque, R. and Mag{\'a}n Madinabeitia, H. and Mall, U. and Mancini, L. and Mandel, H. and Marfil, E. and Mar{\'\i}n Molina, J. A. and Maroto, D. and {Fern{\'a}ndez} and Mart{\'\i}n, E. L. and Mart{\'\i}n-Ruiz, S. and Marvin, C. J. and Mathar, R. J. and Mirabet, E. and Montes, D. and Moreno-Raya, M. E. and Moya, A. and Mundt, R. and Nagel, E. and Naranjo, V. and Nortmann, L. and Nowak, G. and Ofir, A. and Oreiro, R. and Pall{\'e}, E. and Panduro, J. and Pascual, J. and Passegger, V. M. and Pavlov, A. and Pedraz, S. and P{\'e}rez-Calpena, A. and P{\'e}rez Medialdea, D. and Perger, M. and Perryman, M. A. C. and Pluto, M. and Rabaza, O. and Ram{\'o}n, A. and Rebolo, R. and Redondo, P. and Reffert, S. and Reinhart, S. and Rhode, P. and Rix, H.-W. and Rodler, F. and Rodr{\'\i}guez, E. and Rodr{\'\i}guez-L{\'o}pez, C. and Rodr{\'\i}guez Trinidad, A. and Rohloff, R.-R. and Rosich, A. and Sadegi, S. and S{\'a}nchez-Blanco, E. and S{\'a}nchez Carrasco, M. A. and S{\'a}nchez-L{\'o}pez, A. and Sanz-Forcada, J. and Sarkis, P. and Sarmiento, L. F. and Sch{\"a}fer, S. and Schmitt, J. H. M. M. and Schiller, J. and Schweitzer, A. and Solano, E. and Stahl, O. and Strachan, J. B. P. and St{\"u}rmer, J. and Su{\'a}rez, J. C. and Tabernero, H. M. and Tala, M. and Trifonov, T. and Tulloch, S. M. and Ulbrich, R. G. and Veredas, G. and Vico Linares, J. I. and Vilardell, F. and Wagner, K. and Winkler, J. and Wolthoff, V. and Xu, W. and Yan, F. and Zapatero Osorio, M. R.},
  month = nov,
  year = {2017},
  keywords = {Astrophysics - Earth and Planetary Astrophysics,Astrophysics - Solar and Stellar Astrophysics},
  pages = {arXiv:1711.06576}
}

-->


# # References
# 
# <a name="cite-reiners_carmenes_2017"/><sup>[^](#ref-1) </sup>Reiners, A. and Zechmeister, M. and Caballero, J. A. and Ribas, I. and Morales, J. C. and Jeffers, S. V. and Sch&ouml;fer, P. and Tal-Or, L. and Quirrenbach, A. and Amado, P. J. and Kaminski, A. and Seifert, W. and Abril, M. and Aceituno, J. and Alonso-Floriano, F. J. and Ammler-von Eiff, M. and Antona, R. and Anglada-Escud&eacute;, G. and Anwand-Heerwart, H. and Arroyo-Torres, B. and Azzaro, M. and Baroch, D. and Barrado, D. and Bauer, F. F. and Becerril, S. and B&eacute;jar, V. J. S. and Ben\'\itez, D. and Berdi&ntilde;as, Z. M. and Bergond, G. and Bl&uuml;mcke, M. and Brinkm&ouml;ller, M. and del Burgo, C. and Cano, J. and C&aacute;rdenas V&aacute;zquez, M. C. and Casal, E. and Cifuentes, C. and Claret, A. and Colom&eacute;, J. and Cort&eacute;s-Contreras, M. and Czesla, S. and D\'\iez-Alonso, E. and Dreizler, S. and Feiz, C. and Fern&aacute;ndez, M. and Ferro, I. M. and Fuhrmeister, B. and Galad\'\i-Enr\'\iquez, D. and Garcia-Piquer, A. and Garc\'\ia Vargas, M. L. and Gesa, L. and G&oacute;mez, V. and Galera and Gonz&aacute;lez Hern&aacute;ndez, J. I. and Gonz&aacute;lez-Peinado, R. and Gr&ouml;zinger, U. and Grohnert, S. and Gu&agrave;rdia, J. and Guenther, E. W. and Guijarro, A. and de Guindos, E. and Guti&eacute;rrez-Soto, J. and Hagen, H.-J. and Hatzes, A. P. and Hauschildt, P. H. and Hedrosa, R. P. and Helmling, J. and Henning, Th. and Hermelo, I. and Hern&aacute;ndez Arab\'\i, R. and Hern&aacute;ndez Casta&ntilde;o, L. and Hern&aacute;ndez Hernando, F. and Herrero, E. and Huber, A. and Huke, P. and Johnson, E. and de Juan, E. and Kim, M. and Klein, R. and Kl&uuml;ter, J. and Klutsch, A. and K&uuml;rster, M. and Lafarga, M. and Lamert, A. and Lamp&oacute;n, M. and Lara, L. M. and Laun, W. and Lemke, U. and Lenzen, R. and Launhardt, R. and L&oacute;pez del Fresno, M. and L&oacute;pez-Gonz&aacute;lez, J. and L&oacute;pez-Puertas, M. and L&oacute;pez Salas, J. F. and L&oacute;pez-Santiago, J. and Luque, R. and Mag&aacute;n Madinabeitia, H. and Mall, U. and Mancini, L. and Mandel, H. and Marfil, E. and Mar\'\in Molina, J. A. and Maroto, D. and Fern&aacute;ndez and Mart\'\in, E. L. and Mart\'\in-Ruiz, S. and Marvin, C. J. and Mathar, R. J. and Mirabet, E. and Montes, D. and Moreno-Raya, M. E. and Moya, A. and Mundt, R. and Nagel, E. and Naranjo, V. and Nortmann, L. and Nowak, G. and Ofir, A. and Oreiro, R. and Pall&eacute;, E. and Panduro, J. and Pascual, J. and Passegger, V. M. and Pavlov, A. and Pedraz, S. and P&eacute;rez-Calpena, A. and P&eacute;rez Medialdea, D. and Perger, M. and Perryman, M. A. C. and Pluto, M. and Rabaza, O. and Ram&oacute;n, A. and Rebolo, R. and Redondo, P. and Reffert, S. and Reinhart, S. and Rhode, P. and Rix, H.-W. and Rodler, F. and Rodr\'\iguez, E. and Rodr\'\iguez-L&oacute;pez, C. and Rodr\'\iguez Trinidad, A. and Rohloff, R.-R. and Rosich, A. and Sadegi, S. and S&aacute;nchez-Blanco, E. and S&aacute;nchez Carrasco, M. A. and S&aacute;nchez-L&oacute;pez, A. and Sanz-Forcada, J. and Sarkis, P. and Sarmiento, L. F. and Sch&auml;fer, S. and Schmitt, J. H. M. M. and Schiller, J. and Schweitzer, A. and Solano, E. and Stahl, O. and Strachan, J. B. P. and St&uuml;rmer, J. and Su&aacute;rez, J. C. and Tabernero, H. M. and Tala, M. and Trifonov, T. and Tulloch, S. M. and Ulbrich, R. G. and Veredas, G. and Vico Linares, J. I. and Vilardell, F. and Wagner, K. and Winkler, J. and Wolthoff, V. and Xu, W. and Yan, F. and Zapatero Osorio, M. R.. 2017. _The CARMENES Search for Exoplanets around M Dwarfs: High-Resolution Optical and near-Infrared Spectroscopy of 324 Survey Stars_.
# 
# <a name="cite-artigau_optical_2018"/><sup>[^](#ref-2) </sup>Artigau, &Eacute;tienne and Malo, Lison and Doyon, Ren&eacute; and Figueira, Pedro and Delfosse, Xavier and Astudillo-Defru, Nicola. 2018. _Optical and Near-Infrared Radial Velocity Content of M Dwarfs: Testing Models with Barnard's Star_.
# 
# 
