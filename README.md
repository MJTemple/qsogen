
# qsogen <a href="https://ascl.net/2205.003"><img src="https://img.shields.io/badge/ascl-2205.003-blue.svg?colorB=262255" alt="ascl:2205.003" /></a>

------------
Introduction
------------

`qsogen` is a collection of Python code which models quasar colours, magnitudes and
SEDs. It implements an empirically-motivated parametric model to efficiently
account for the observed emission-line properties, host-galaxy contribution, dust
reddening, hot dust emission, and IGM suppression in the rest-frame 900-30000A
wavelength range for quasars with a wide range of redshift and luminosity.

The code is packaged with a set of empirically-derived emission-line templates and
an empirically-derived quasar dust extinction curve.

`qsogen` requires `numpy`, `scipy`, and `astropy`; the examples below also use `matplotlib`. 
The code has been written on a RHEL7 Linux machine running `Python 3.6.7`,
with `astropy 4.0`, `matplotlib 3.0.2`, `numpy 1.16.2` and `scipy 1.2.1` from `anaconda3`.

The code was written by [Matthew Temple](https://mjtemple.github.io), with core functionality translated from
earlier FORTRAN code written by Paul Hewett. The structure of the model is
described in
[Temple, Hewett & Banerji, 2021, MNRAS, 508, 737](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508..737T)
[(arXiv:2109.04472)](https://arxiv.org/abs/2109.04472),
which should be read in conjunction with this file.

If you use this code in a scientific publication, please cite
[Temple, Hewett & Banerji (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508..737T).

Comments, questions, suggestions and bug reports are welcomed by the first
author: Matthew.Temple[at]mail.udp.cl

Contents of this README file:
* [Introduction](https://github.com/MJTemple/qsogen#introduction)
* [File descriptions](https://github.com/MJTemple/qsogen#File-descriptions)
* [A note on measures of luminosity](https://github.com/MJTemple/qsogen#A-note-on-measures-of-luminosity)
* [Example use cases](https://github.com/MJTemple/qsogen#Example-use-cases)
* [Exceptions](https://github.com/MJTemple/qsogen#Exceptions)

-----------------
File descriptions
-----------------

The required code files are:

1. `qsosed.py`

Defines a class `Quasar_sed`, which generates an instance of the model SED.
    To print documentation for `Quasar_sed`, including a description of input
    parameters, run `python qsosed.py` from the terminal.
    
2. `config.py`

Contains a dictionary of parameters, `params`, which are passed to `Quasar_sed`.
    All parameters can be overruled by passsing them as `**kwargs` to `Quasar_sed`.
    
3. `model_colours.py`

Defines functions `get_colours` and `get_mags` which return arrays of model
    colours and model magnitudes, respectively, for a given set of redshifts.
    Model parameters can be passed as `**kwargs` to `Quasar_sed` which overwrite
    settings from `config.py`.

The model also requires some additional input files, which are:

4. `qsosed_emlines_20210625.dat` 

the emission line templates, with columns
[wavelength in A, median line template, reference continuum, high-EW template, high-blueshift template, narrow line template],

5. `S0_template_norm.sed`

the S0 galaxy template from SWIRE, with columns [wavelength in A, f_lambda], and

6. `pl_ext_comp_03.sph`  

the Quasar extinction curve, with columns [wavelength in A, E(lambda-V)/E(B-V)].

If you want to generate synthetic photometry using `model_colours.py`,
you also need the relevant filter response files in the same working directory,
which must have the form [wavelength in A, filter response].

Response files are currently provided (in the `filters/` subdirectory) for the following filters:

* DECam ugrizY
* Euclid YJH
* GALEX FUV and NUV
* HSC grizY, i2, r2
* LSST ugrizy
* SDSS ugriz
* UKIDSS ZYJHK
* VISTA ZYJHKs
* WISE W1234

In addition, if you want to compute Vega zero-point magnitudes for any new filters,
you need the Vega spectrum:
* vega_2007.lis

--------------------------------
A note on measures of luminosity
--------------------------------

There are two different measures of luminosity which can be input to the model.

**`LogL3000`** controls the output monochromatic 3000A continuum luminosity in erg/s.
    All this does is re-scale the output model flux appropriately.
    This should be used when you care about the absolute value of the flux units
    e.g. when using `get_mags` to return synthetic magnitudes, or using `get_mags`
    to fit the model to observed magnitudes and returning the luminosity of the
    observed quasar
    (NB. in such cases `ebv` should also be a free parameter).
    The default is to normalise the model such that L3000 = 10^46 erg/s.

**`M_i`** represents the absolute i-band magnitude at z=2, as defined by 
[Richards et al., 2006, AJ, 131, 2766](https://ui.adsabs.harvard.edu/abs/2006AJ....131.2766R/abstract)
    and as adopted by the SDSS quasar catalogues.
    This parameter is used 'under the hood' to control the emission line 
    properties and the relative contribution of the host galaxy component.

The default settings are to use the average `M_i` as a function of redshift
    from the sample of 18.6<i_AB<19.1 quasars in SDSS DR16Q. For exploring the
    properties of significantly brighter or fainter populations it is recommended
    to input your own luminosity-redshift relation using the parameter `zlum_lumval`.
    Alternatively, one can override the default settings and specify `M_i` exactly
    for individual objects.

Predicted values of `M_i` as a function of redshift and apparent i_AB magnitude
    can be approximated by 

`M_i(z, i) = -log10(z)*(0.250*(i/20) + 5.050) - (17.40*(20/i) + 6.82)`

A very crude approximation for converting between these parameters is

`LogLbol = -0.4*M_i + 36 = LogL3000 + 0.7`

-----------------
Example use cases
-----------------


Create and plot a z=2 quasar model using the default parameters in rest frame:
```python
>>> from qsosed import Quasar_sed
>>> import matplotlib.pyplot as plt
>>> Quasar2 = Quasar_sed(z=2)
>>> plt.subplots()
>>> plt.plot(Quasar2.wavlen, Quasar2.flux)
>>> plt.xlabel('Rest Wavelength [A]')
>>> plt.ylabel('Flux density per unit wavelength')
```
Plot the same model in the observed frame:
```python
>>> plt.subplots()
>>> plt.loglog(Quasar2.wavred, Quasar2.wavred*Quasar2.flux, label='z=2 model')
>>> plt.xlabel('Observed Wavelength [A]')
>>> plt.ylabel('Flux density $\lambda F_\lambda$')
```
 Add a z=4 quasar model to the plot:
```python
>>> Quasar4 = Quasar_sed(z=4)
>>> plt.loglog(Quasar4.wavred, Quasar4.wavred*Quasar4.flux, label='z=4 model')
>>> plt.legend()
```
Return intra-band colours for default filters SDSS-UKIDSS-WISE for `z=[1,2,3]`:
```python
>>> from model_colours import get_colours
>>> colours = get_colours([1,2,3])
# colours[0] is (u-g, g-r, r-i, i-z, z-Y, Y-J, J-H, H-K, K-W1, W1-W2) for z=1
# colours[1] is the same for z=2
# colours[2] is the same for z=3
```
`get_colours()` takes the same `**kwargs` as `Quasar_sed`:
```python
>>> redder_colours = get_colours([1,2,3], ebv=0.1)
# redder_colours[0] is an array of colours for z=1 with extinction E(B-V)=0.1
```
If you want to find colours for a different filter set, use a list of filters:
```python
>>> LSST_colours = get_colours([1,2,3],
>>>                            filters=['LSST_u_AB', 'LSST_g_AB', 'LSST_r_AB',
>>>                                     'LSST_i_AB', 'LSST_z_AB', 'LSST_y_AB',
>>>                                     'Euclid_Y_AB', 'Euclid_J_AB',
>>>                                     'Euclid_H_AB', 'WISE_W1_Vega',
>>>                                     'WISE_W2_Vega'])
```
`get_mags()` has same syntax as `get_colours()`, but returns synthetic magnitudes
instead of synthetic colours. See note on measures of luminosity above.

Example: fit model to observed ugriz photometry to estimate L3000 and E(B-V)
          for object of known redshift `z<2.3` (i.e. not u-band dropout)

Note the first argument of `get_mags` is an array of redshifts `[z]`, and it
 returns an array of arrays, so `get_mags([2, 3])[0]` is an array of magnitudes
 for `z=2` and `get_mags([2, 3])[1]` is an array of magnitudes for `z=3`.
```python
>>> from model_colours import get_mags
>>> import numpy as np
>>> from scipy.optimize import curve_fit
>>> observed_mags = np.array([u, g, r, i, z])
>>> observed_mag_errs = np.array([uerr, gerr, rerr, ierr, zerr])
>>> def f(Z, LogL3000, ebv):
>>>     mod_mags = get_mags([Z],
>>>                         LogL3000=LogL3000,
>>>                         ebv=ebv)[0] 
>>>     return mod_mags[:5]  # keep ugriz
>>> popt, pcov = curve_fit(f, Z, observed_mags, p0=[46., 0.], sigma=observed_mag_errs)
>>> popt[0]  # best-fitting LogL3000 of unreddened quasar, i.e. assuming E(B-V)=0
>>> popt[1]  # best-ftting E(B-V)
```
Explore the effect of dust reddening E(B-V), galaxy contribution, and emission
 lines on colour-redshift tracks:
```python
>>> from model_colours import get_colours
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> zarr = np.linspace(0.1, 5., num=100)
>>> base_colours = get_colours(zarr).T  # transpose to get colours as fn of z
>>> red_colours = get_colours(zarr, ebv=0.1).T
>>> redder_colours = get_colours(zarr, ebv=0.2).T
>>> nogal_colours = get_colours(zarr, fragal=0.1).T  # 'low galaxy' model
>>> pky_colours = get_colours(zarr, emline_type=1.).T  # High-EW emission lines
>>> wdy_colours = get_colours(zarr, emline_type=-1.).T  # High-blueshift emission lines
>>> 
>>> names = ['$u-g$', '$g-r$', '$r-i$', '$i-z$', '$z-Y$', 
>>>          '$Y-J$', '$J-H$', '$H-K$', '$K-W1$', '$W1-W2$']
>>> fig, ax = plt.subplots(10, 1, figsize=(6, 12), sharex=True)
>>> axarr = ax.flatten()
>>> for ind in range(10):
>>>     axarr[9-ind].plot(zarr, base_colours[ind], label='Base model')
>>>     axarr[9-ind].plot(zarr, red_colours[ind], label='E(B-V)=0.1 model')
>>>     axarr[9-ind].plot(zarr, redder_colours[ind], label='E(B-V)=0.2 model')
>>>     axarr[9-ind].plot(zarr, nogal_colours[ind], label='Low galaxy model')
>>>     axarr[9-ind].set_ylabel(names[ind])
>>> axarr[-1].legend(loc='lower right')
>>> axarr[-1].set_ylim(-0.1, 0.8)
>>> axarr[-2].set_ylim(-0.1, 0.8)
>>> axarr[-3].set_ylim(-0.1, 0.8)
>>> plt.xlabel('Redshift')
>>> fig.tight_layout()
>>> fig.subplots_adjust(hspace=0)
>>> fig, ax = plt.subplots(10, 1, figsize=(6, 12), sharex=True)
>>> axarr = ax.flatten()
>>> for ind in range(10):
>>>     axarr[9-ind].plot(zarr, base_colours[ind], label='Base model')
>>>     axarr[9-ind].plot(zarr, pky_colours[ind], label='High-EW emission lines')
>>>     axarr[9-ind].plot(zarr, wdy_colours[ind], label='High-blueshift emission lines')
>>>     axarr[9-ind].set_ylabel(names[ind])
>>> axarr[-1].legend(loc='lower right')
>>> axarr[-1].set_ylim(-0.1, 0.8)
>>> axarr[-2].set_ylim(-0.1, 0.8)
>>> axarr[-3].set_ylim(-0.1, 0.8)
>>> plt.xlabel('Redshift')
>>> fig.tight_layout()
>>> fig.subplots_adjust(hspace=0)
```
----------
Exceptions
----------
Quasar_sed may raise one of two exceptions if the input wavelength array is unsuitable:

`Exception('wavlen must be monotonic')`

The input wavelength array is not monotonically increasing. Ensure that the
    input array 'wavlen' is sorted in increasing order.

`Exception('wavlen must cover 4000-5000 A for galaxy normalisation'
           + '\n Redshift is {}'.format(self.z))`

If `gflag=True`, i.e. if the host galaxy component is switched on, then the
    input wavelength array 'wavlen' must cover the 4000-5000A region to ensure
    the correct normalisation of the galaxy component.
    If using `get_colours` or `get_mags`, this exception can occur if the filters
    which have been chosen do not cover 4000-5000A in the rest frame.
    To avoid problems, it is highly recommended to use a dense set of filters
    which cover a wide wavelength range. For example, even if you only want
    LSST ugrizy colours for redshifts 0<z<5, it is preferable to find all
    LSST-Euclid-WISE ugrizyYJHW1W2 colours using `get_colours()`, and then take
    the first 5 columns of the output matrix.

In addition, if you pass a `**kwarg` to `Quasar_sed` which is not recognised, it
will print a warning, and then proceed to ignore the unknown `**kwarg`.
