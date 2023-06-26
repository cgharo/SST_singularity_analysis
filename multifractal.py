"""
Tools for singularity analysis and computing different multifractal quantities

Functions
---------
singspec        Compute singularity spectrum from a set of singularity expoents
transinv        Compute thje shift necessary to have translational invariance
                of singularity spectra
written by Jordi Isern
"""
import numpy as np

###############################################################################
# singspec: Compute singularity analysis
###############################################################################


def singspec(hin, hmin=False, hmax=False, hbin=0.02, dmax=2):
    """
    h, hist, dh, edh = singspec(hin, hmin=False, hmax=False, hbin=0.02,
    dmax=dim_max)

    Compute the Singularity Spectrum D(h) from an array of singularity
    exponents

    Parameters
    ----------

    hin: array with the singularity exponents

    hmax: maximum singularity exponent
    hmin: minimum singularity exponent
    hbin: bin size
    dmax: gemoetric dimension (1, 2 or 3)

    Return
    ------

    h: output singularity exponents
    hist: histogram of singularity expionents
    dh: singularity spectrum
    edh: error in singularity exponents

    """

    # General properties

    if not hmin:
        hmin = np.min(hin)
    if not hmax:
        hmax = np.max(hin)

    nbins = int((hmax - hmin) / hbin)  # Number of bins
    n = len(hin.flatten())  # Number of points #################### 
    r = 1 / np.sqrt(n)  # Minimum scale

    # Compute the PDF. Compute the historgram and normalize by its integral

    hist, bins = np.histogram(hin.flatten(), range=(hmin, hmax), bins=nbins)
    h = ((bins + np.roll(bins, -1)) / 2)[:-1]

    # Compute the Singularity Spectra

    dh = np.zeros(hist.shape)
    histmax = np.max(hist)
    if histmax == 0:
        histmax = 1
    ind = (hist > 0)
    dh[ind] = dmax - (np.log(hist[ind]) - np.log(histmax)) / np.log(r)

    # Compute the error of dH

    edh = np.zeros(len(hist))
    ind = hist > 0
    edh[ind] = 3 / np.sqrt(hist[ind]) / np.abs(np.log(r))

    return h, hist, dh, edh

###############################################################################
# transinv: Impose translational invarianve
###############################################################################


def transinv(hin, dmax=2,hmin=False, hmax=False, hbin=0.05, return_all=False):
    """

    shift = transinv(hin, dmax=2, hbin=0.05, return_all=False)

    Impose translational invariance

    Parameters
    ----------

    hin:            array with the singularity exponents

    hbin:           bin size
    dmax:           gemoetric dimension (1, 2 or 3)

    return_all:     Returns all information. For testing

    Return
    ------

    shift: shift to be added to singularity exponents

    """

    # Compute the singularity spectrum

    hi, ni, dhi, erri = singspec(hin, hbin=hbin, dmax=dmax,hmin=hmin, hmax=hmax)

    # Get the tangent line

    imax = np.argmax(dhi)
    hp = hi[:imax]
    dhp = dhi[:imax]
    c = dhp - hp

    i1 = np.argmax(c)
    h1 = hp[i1]
    d1 = dhp[i1]
    shift = d1 - dmax - h1

    # Returning

    if return_all:
        return shift, (c[i1], hi, dhi, h1, d1)
    else:
        return shift
