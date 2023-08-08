BlackHawk Version 2.0 (13 October 2021)
-----------------------------------
By Alexandre Arbey (alexandre.arbey@ens-lyon.fr) and Jeremy Auffinger (j.auffinger@ipnl.in2p3.fr)

INTRODUCTION
------------
The most recent version of this program can be obtained from: 
https://blackhawk.hepforge.org/

This C code (in C99 standard) provides a detailed computation of the Hawking emission spectrum
for any mass and spin distribution of Black Holes in diverse metric (Kerr, Reissner-Nordström,
higher-dimensional, polimerized).

The program has been tested on Linux, Mac and Windows (using Cygwin) machines with gcc, clang and icc.

If you use BlackHawk to write a paper, please cite:

A. Arbey and J. Auffinger, Eur. Phys. J. C79 (2019) 693, arXiv:1905.04268v3 [gr-qc]
A. Arbey and J. Auffinger, arXiv:2108.02737 [gr-qc]

as well as PYTHIA, HERWIG or Hazma (respectively) depending on the choice of hadronization code:

T. Sjostrand, S. Ask, J. R. Christiansen, R. Corke, N. Desai, P. Ilten, S. Mrenna, S. Prestel, C. O. Rasmussen, and P. Z.
Skands, Comput. Phys. Commun. 191, 159 (2015), arXiv:1410.3012 [hep-ph]
J. Bellm et al., Eur. Phys. J. C 76, 196 (2016), arXiv:1512.01178 [hep-ph]
A. Coogan, L. Morrison, and S. Profumo, JCAP 01, 056, arXiv:1907.11846 [hep-ph]

Installation and Compilation
----------------------------
- tar xzvf blackhawk_vX.X.tgz
- cd blackhawk_vX.X
- in Makefile, define your C compiler
- compile with: make
- create the executable with: make BlackHawk_*, where * is "tot" or "inst"

Included Files
--------------
- Procedures in src/:
evolution.c general.c primary.c secondary.c spectrum.c technical.c hadro_herwig.c hadro_pythia.c hadro_pythianew.c hadro_hazma.c

- Main programs:
BlackHawk_inst.c: calculation of the instantaneous Hawking spectra
BlackHawk_tot.c: calculation of the time-dependent Hawking spectra

- Headers in src/:
include.h: definitions and prototypes

- src/tables/:
Numerical tables used in the code

- manual/:
A .pdf of the up-to-date manual of the code

- results/:
The folder where run data will be saved. The result are given in CGS units.

- scripts/:
C, C++ and Mathematica scripts used to compute the numerical tables (cosmology_scripts/, greybody_scripts/, herwig_scripts/, pythia_scripts/ and pythianew_scripts/)

- Others:
Makefile(s)
README.txt files


History
-------

v1.0 - 22/05/2019 - First public release.
v1.1 - 14/11/2019 - Various bugs corrected. Addition of a spin distribution.
V2.0 - 06/08/2021 - Non standard BH metrics. Dark matter emission. Corrected PYTHIA tables. Various bug corrections.

LICENSE
-------
    BlackHawk Copyright (C) 2019 A. Arbey, J. Auffinger

    This program is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any 
    later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    See <http://www.gnu.org/licenses/>.  
