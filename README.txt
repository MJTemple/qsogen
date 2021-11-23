------------
Introduction
------------

qsogen is a collection of Python code to model quasar colours, magnitudes and SEDs.
It requires numpy, scipy, and astropy; the examples below also use matplotlib. 

The code has been written on a Linux machine running RHEL7 and Python 3.6.7,
with astropy 4.0, matplotlib 3.0.2, numpy 1.16.2 and scipy 1.2.1 from anaconda3.

The code was written by Matthew Temple, with core functionality translated from
earlier FORTRAN code written by Paul Hewett. The structure of the model is
described in Temple, Hewett & Banerji 2021MNRAS.508..737T (arXiv:2109.04472),
which should be read in conjunction with this file.

If you use this code in a scientific publication, please cite
Temple, Hewett & Banerji (2021).

Contact: Matthew.Temple[at]mail.udp.cl

Contents of this README file:
Introduction
File descriptions
A note on measures of luminosity
Example use cases
Exceptions

-----------------
File descriptions
-----------------

The required code files are:
qsosed.py
    Defines a class Quasar_sed, which generates an instance of the model SED.
    To print documentation for Quasar_sed, including a description of input
    parameters, run 'python qsosed.py' from the terminal.
config.py
    Contains a dictionary of parameters which are passed to Quasar_sed.
    All parameters can be overruled by passsing them as **kwargs to Quasar_sed.
model_colours.py
    Defines functions get_colours and get_mags which return arrays of model
    colours and model magnitudes, respectively, for a given set of redshifts.
    Model parameters can be passed as **kwargs to Quasar_sed which overwrite
    settings from config.py.

In addition, for the model to run, it requires some additional input files:
qsosed_emlines.dat
    Emission line templates, with columns
    [wavelength in A, median line template, reference continuum,
    high-EW template, high-blueshift template, narrow line template]
S0_template_norm.sed
    S0 galaxy template from SWIRE, with columns [wavelength in A, f_lambda]
pl_ext_comp_03.sph
    Quasar extinction curve, with columns [wavelength in A, E(lambda-V)/E(B-V)]

If you want to use model_colours.py, you also need filter response files in the
same directory, which must have the form [wavelength in A, filter response].
Those filter response files currently provided are:
DECam_u.filter
DECam_g.filter
DECam_r.filter
DECam_i.filter
DECam_z.filter
DECam_Y.filter
Euclid_H.filter
Euclid_J.filter
Euclid_Y.filter
HSC_g.filter
HSC_r.filter
HSC_r2.filter
HSC_i.filter
HSC_i2.filter
HSC_z.filter
HSC_Y.filter
LSST_g.filter
LSST_i.filter
LSST_r.filter
LSST_u.filter
LSST_y.filter
LSST_z.filter
SDSS_g.filter
SDSS_i.filter
SDSS_r.filter
SDSS_u.filter
SDSS_z.filter
UKIDSS_H.filter
UKIDSS_J.filter
UKIDSS_K.filter
UKIDSS_Y.filter
UKIDSS_Z.filter
VISTA_H.filter
VISTA_J.filter
VISTA_Ks.filter
VISTA_Y.filter
VISTA_Z.filter
WISE_W1.filter
WISE_W2.filter
WISE_W3.filter
WISE_W4.filter

In addition, if you want to compute Vega zero-point magnitudes for any new filters,
you need the Vega spectrum:
vega_2007.lis

--------------------------------
A note on measures of luminosity
--------------------------------

There are two different measures of luminosity which can be input to the model.

LogL3000 controls the output monochromatic 3000A continuum luminosity in erg/s.
    All this does is re-scale the output model SED appropriately.
    This should be used when you care about the absolute value of the flux units
    e.g. when using get_mags to return synthetic magnitudes, or using get_mags
    to fit the model to observed magnitudes and returning the luminosity of the
    observed quasar. [NB. in such cases ebv should also be a free parameter.]
    The default is to normalise the model such that L3000 = 10^46 erg/s.

M_i represents the absolute i-band magnitude at z=2, as defined by Richards+ 
    2006AJ....131.2766R and as adopted by the SDSS quasar catalogues.
    This parameter is used 'under the hood' to control the emission line 
    properties and the relative contribution of the host galaxy component.
    The default settings are to use the average M_i as a function of redshift
    from the sample of 18.6<i_AB<19.1 quasars in SDSS DR16Q. For exploring the
    properties of significantly brighter or fainter populations it is recommended
    to input your own luminosity-redshift relation using the parameter zlum_lumval.
    Alternatively, one can override the default settings and specify M_i exactly
    for individual objects.
    Predicted values of M_i as a function of redshift and apparent magnitude
    can be approximated by 
    M_i(z, i) = -log10(z)*(0.250*(i/20) + 5.050) - (17.40*(20/i) + 6.82)

A very crude approximation for converting between these parameters is
LogLbol = -0.4*MI + 36 = LogL3000 + 0.7

-----------------
Example use cases
-----------------
#
# Create and plot a z=2 quasar model using the default parameters in rest frame
#
>>> from qsosed import Quasar_sed
>>> import matplotlib.pyplot as plt
>>> Quasar2 = Quasar_sed(z=2)
>>> plt.subplots()
>>> plt.plot(Quasar2.wavlen, Quasar2.flux)
>>> plt.xlabel('Rest Wavlength [A]')
>>> plt.ylabel('Flux density per unit wavlength')
#
# Plot the same model in the observed frame
#
>>> plt.subplots()
>>> plt.loglog(Quasar2.wavred, Quasar2.wavred*Quasar2.flux, label='z=2 model')
>>> plt.xlabel('Observed Wavelength [A]')
>>> plt.ylabel('Flux density $\lambda F_\lambda$')
#
# Add a z=4 quasar model to the plot
#
>>> Quasar4 = Quasar_sed(z=4)
>>> plt.loglog(Quasar4.wavred, Quasar4.wavred*Quasar4.flux, label='z=4 model')
>>> plt.legend()
#
# Return intra-band colours for default filters SDSS-UKIDSS-WISE for z=[1,2,3]
#
>>> from model_colours import get_colours
>>> colours = get_colours([1,2,3])
# colours[0] is (u-g, g-r, r-i, i-z, z-Y, Y-J, J-H, H-K, K-W1, W1-W2) for z=1
# colours[1] is the same for z=2
# colours[2] is the same for z=3
#
# get_colours() takes the same **kwargs as Quasar_sed:
>>> redder_colours = get_colours([1,2,3], ebv=0.1)
# redder_colours[0] is an array of colours for z=1 with extinction E(B-V)=0.1
#
# If you want to find colours for a different filter set, use a list of filters
#
>>> LSST_colours = get_colours([1,2,3],
>>>                            filters=['LSST_u_AB', 'LSST_g_AB', 'LSST_r_AB',
>>>                                     'LSST_i_AB', 'LSST_z_AB', 'LSST_y_AB',
>>>                                     'Euclid_Y_AB', 'Euclid_J_AB',
>>>                                     'Euclid_H_AB', 'WISE_W1_Vega',
>>>                                     'WISE_W2_Vega'])
#
# get_mags() has same syntax as get_colours(), but returns synthetic magnitudes
# instead of synthetic colours. See note on measures of luminosity above.
# Example: fit model to observed ugriz photometry to estimate L3000 and E(B-V)
#          for object of known redshift Z<2.3 (i.e. not u-band dropout)
# Note the first arguement of get_mags is an array of redshifts [Z], and it
# returns an array of arrays, so get_mags([2, 3])[0] is an array of magnitudes
# for Z=2 and get_mags([2, 3])[1] is an array of magnitudes for Z=3.
#
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
#
# Explore the effect of dust reddening E(B-V), galaxy contribution, and emission
# lines on colour-redshift tracks
#
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

----------
Exceptions
----------
Quasar_sed may raise an exception if the input wavelength array is unsuitable:

Exception('wavlen must be monotonic')
    The input wavelength array is not monotonically increasing. Ensure that the
    input array 'wavlen' is sorted in increasing order.

Exception('wavlen must cover 4000-5000 A for galaxy normalisation'
           + '\n Redshift is {}'.format(self.z))
    If gflag=True, i.e. if the host galaxy component is switched on, then the
    input wavelength array 'wavlen' must cover the 4000-5000A region to ensure
    the correct normalisation of the galaxy component.
    If using get_colours or get_mags, this exception can occur if the filters
    which have been chosen do not cover 4000-5000A in the rest frame.
    To avoid problems, it is highly recommended to use a dense set of filters
    which cover a wide wavelength range. For example, even if you only want
    LSST ugrizy colours for redshifts 0<z<5, it is preferable to find all
    LSST-Euclid-WISE ugrizyYJHKW1W2 colours using get_colours(), and then take
    the first 5 columns of the output matrix.

In addition, if you pass a **kwarg to Quasar_sed which is not recognised, it
will print a warning, and then proceed to ignore the unknown kwarg.

----------
End README
----------
