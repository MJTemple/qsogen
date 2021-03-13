#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration file with parameters for qsosed.py
Parameters derived from SDSS DR16Q x UKIDSS LAS x unWISE as described in
Temple, Hewett & Banerji (MNRAS submitted).

Idea is to load emission line templates, host galaxy template, and extinction
curve once to minimise i/o.

@author: Matthew Temple
Last edit 2021 March 13
"""
import numpy as np

f1 = 'qsosed_emlines.dat'
emline_template = np.genfromtxt(f1, unpack=True)
# wav, median_emlines, continuum, peaky_line, windy_lines, narrow_lines

f2 = 'S0_template_norm.sed'
galaxy_template = np.genfromtxt(f2, unpack=True)
# S0 galaxy template from SWIRE
# https://ui.adsabs.harvard.edu/abs/2008MNRAS.386..697R/

f3 = 'pl_ext_comp_03.sph'
reddening_curve = np.genfromtxt(f3, unpack=True)

# fit to DR16Q median 2sigma-clipped colours in multi-imag bins
params = dict(plslp1=-0.347,
              plslp2=0.739,
              plstep=-1.0,  # (not fit for)
              plbrk1=3800.,
              tbb=1187.,
              plbrk3=1200.,  # (not fit for)
              bbnorm=4.08,
              scal_emline=-1.,
              emline_type=None,
              scal_halpha=1.,
              scal_lya=1.,
              scal_nlr=1.,
              emline_template=emline_template,
              galaxy_template=galaxy_template,
              reddening_curve=reddening_curve,
              zlum_lumval=np.array([[0.23, 0.34, 0.6, 1.0, 1.4, 1.8, 2.2,
                                     2.6, 3.0, 3.3, 3.7, 4.13, 4.5],
                                    [-21.76, -22.9, -24.1, -25.4, -26.0,
                                     -26.6, -27.1, -27.6, -27.9, -28.1, -28.4,
                                     -28.6, -28.9]]),
              M_i=None,
              beslope=0.189,
              benorm=-27.,  # (not fit for)
              bcnorm=False,
              lyForest=True,
              lylim=912.,  # (not fit for)
              gflag=True,
              fragal=0.247,
              gplind=0.723)
