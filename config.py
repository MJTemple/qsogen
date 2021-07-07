#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration file with parameters for qsosed.py
Parameters derived from SDSS DR16Q x UKIDSS LAS x unWISE as described in
Temple, Hewett & Banerji (MNRAS submitted).

Idea is to load emission line templates, host galaxy template, and extinction
curve once to minimise i/o.

@author: Matthew Temple
Last edit 2021 July 02: v20210625 emline_templates and associated model params.
"""
import numpy as np

f1 = 'qsosed_emlines_20210625.dat'
emline_template = np.genfromtxt(f1, unpack=True)
# wav, median_emlines, continuum, peaky_line, windy_lines, narrow_lines

f2 = 'S0_template_norm.sed'
galaxy_template = np.genfromtxt(f2, unpack=True)
# S0 galaxy template from SWIRE
# https://ui.adsabs.harvard.edu/abs/2008MNRAS.386..697R/

f3 = 'pl_ext_comp_03.sph'
reddening_curve = np.genfromtxt(f3, unpack=True)
# Extinction curve, format: [lambda, E(lambda-V)/E(B-V)]
# Recall flux_reddened(lambda) = flux(lambda)*10^(-A(lambda)/2.5)
# where A(lambda) = E(B-V)*[E(lambda-V)/E(B-V) + R] 
# so taking R=3.1, A(lambda) = E(B-V)*[Col#2 + 3.1]

# fit to DR16Q median 2sigma-clipped colours in multi-imag bins
params = dict(plslp1=-0.349,
              plslp2=0.593,
              plstep=-1.0,    # (not fit for)
              plbrk1=3880.,
              tbb=1243.6,
              plbrk3=1200,   # (not fit for)
              bbnorm=3.961,
              scal_emline=-0.9936,
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
              beslope=0.183,
              benorm=-27.,    # (not fit for)
              bcnorm=False,
              lyForest=True,
              lylim=912,   # (not fit for)
              gflag=True,
              fragal=0.244,
              gplind=0.684)
