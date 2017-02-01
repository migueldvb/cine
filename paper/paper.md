---
title: 'Cine: A solar-pumped fluorescence model for cometary atmospheres'
tags:
  - Python
  - Fluorescence
  - Comets
  - Planetary Atmospheres
authors:
  - name: Miguel de Val-Borro
    orcid: 0000-0002-0455-9384
    affiliation: 1,2
  - name: Martin A. Cordiner
    orcid: 0000-0001-8233-2436
    affiliation: 1,2
  - name: Stefanie N. Milam
    orcid: 0000-0001-7694-4129
    affiliation: 1
  - name: Steven B. Charnley
    affiliation: 1
affiliations:
  - name: NASA Goddard Space Flight Center
    index: 1
  - name: The Catholic University of America
    index: 2
date: 1 Feb 2017
bibliography: paper.bib
---

# Summary

*CINE* is a Python module for calculating infrared pumping efficiencies that can
be applied to the most common molecules found in cometary comae such as water,
hydrogen cyanide or methanol.  Excitation by solar radiation of vibrational
bands followed by radiative decay to the ground vibrational state is one of the
main mechanisms for molecular excitation in comets.  This code calculates the
effective pumping rates for rotational levels in the ground vibrational state
scaled by the heliocentric distance of the comet as detailed in @bensch and @crovisier.
Line transitions are queried from the latest version of the HITRAN
spectroscopic repository [@hitran] using the astroquery affiliated package of
astropy [@astroquery].  Molecular data are obtained from the LAMDA database
[@lamda].

These coefficients are useful for modeling rotational emission lines observed
in cometary spectra at sub-millimeter wavelengths. Combined with computational
methods to solve the radiative transfer equations based, e.g., on the Monte
Carlo algorithm [@lime], this model can retrieve production rates and
rotational temperatures from the observed emission spectrum.

The code and issue tracker of *CINE* are available on GitHub [@cine_github] and
any questions or bug reports can be raised there.  The source code for *CINE*
is also available from the Python Package Index (PyPI).

# References
