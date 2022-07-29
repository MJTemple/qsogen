----------------------
Filter Response Curves
----------------------

Filter response curves are stored in a subfolder here, but need to be in the working directory
when running `model_colours.py` to generate synthetic photometry (magnitudes or colours).

Such filter response files need to have the form [wavelength in A, filter response].

Those filter response files currently provided are:

* DECam_g.filter
* DECam_i.filter
* DECam_r.filter
* DECam_u.filter
* DECam_Y.filter
* DECam_z.filter
* Euclid_H.filter
* Euclid_J.filter
* Euclid_Y.filter
* GALEX_FUV.filter
* GALEX_NUV.filter
* HSC_g.filter
* HSC_i.filter
* HSC_i2.filter
* HSC_r.filter
* HSC_r2.filter
* HSC_Y.filter
* HSC_z.filter
* LSST_g.filter
* LSST_i.filter
* LSST_r.filter
* LSST_u.filter
* LSST_y.filter
* LSST_z.filter
* SDSS_g.filter
* SDSS_i.filter
* SDSS_r.filter
* SDSS_u.filter
* SDSS_z.filter
* UKIDSS_H.filter
* UKIDSS_J.filter
* UKIDSS_K.filter
* UKIDSS_Y.filter
* UKIDSS_Z.filter
* VISTA_H.filter
* VISTA_J.filter
* VISTA_Ks.filter
* VISTA_Y.filter
* VISTA_Z.filter
* WISE_W1.filter
* WISE_W2.filter
* WISE_W3.filter
* WISE_W4.filter

In addition, if you want to compute Vega zero-point magnitudes for any new filters,
you need the Vega spectrum, available in the main code directory:
* ../vega_2007.lis
