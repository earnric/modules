#!/usr/bin/env py.test
import numpy as np
from numpy.testing import assert_allclose
import extinction


def test_ccm89():
    # NOTE: Test is only to precision of 0.016 because there is a discrepancy
    # of 0.014 for the B band wavelength of unknown origin (and up to 0.002 in
    # other bands).
    #
    # Note that a and b can be obtained with:
    # b = ccm89(wave, 0.)
    # a = ccm89(wave, 1.) - b
    # 
    # These differ from the values tablulated in the original paper.
    # Could be due to floating point errors in the original paper?
    #
    # U, B, V, R, I, J, H, K band effective wavelengths from CCM '89 table 3

    x_inv_microns = np.array([2.78, 2.27, 1.82, 1.43, 1.11, 0.80, 0.63, 0.46])
    wave = 1.e4 / x_inv_microns

    # A(lambda)/A(V) for R_V = 3.1 from Table 3 of CCM '89
    ref_values = np.array([1.569, 1.337, 1.000, 0.751, 0.479, 0.282, 0.190,
                           0.114])

    assert_allclose(extinction.ccm89(wave, 1.0, 3.1), ref_values,
                    rtol=0.016, atol=0.)


def test_odonnell94():
    # NOTE: The tabulated values go to 0.001, but the test is only for matching
    # at the 0.005 level, because there is currently a discrepancy up to 0.0047
    # of unknown origin.
    #
    # Tests od94() at Rv = 3.1 against the widely used values tabulated in
    # Schlegel, Finkbeiner and Davis (1998)
    # http://adsabs.harvard.edu/abs/1998ApJ...500..525S
    #
    # This is tested by evaluating the extinction curve at a (given)
    # effective wavelength, since these effective wavelengths:
    # "... represent(s) that wavelength on the extinction curve 
    # with the same extinction as the full passband."
    #
    # The test does not include UKIRT L' (which, at 3.8 microns) is 
    # beyond the range of wavelengths allowed by the function
    # or the APM b_J filter which is defined in a non-standard way. 
    #
    # The SFD98 tabulated values go to 1e-3, so we should be able to match at
    # that level.

    wave = np.array([3372., 4404., 5428., 6509., 8090.,
                     3683., 4393., 5519., 6602., 8046.,
                     12660., 16732., 22152.,
                     5244., 6707., 7985., 9055.,
                     6993.,
                     3502., 4676., 4127.,
                     4861., 5479.,
                     3546., 4925., 6335., 7799., 9294.,
                     3047., 4711., 5498.,
                     6042., 7068., 8066.,
                     4814., 6571., 8183.])

    ref_values = np.array([1.664, 1.321, 1.015, 0.819, 0.594,
                           1.521, 1.324, 0.992, 0.807, 0.601,
                           0.276, 0.176, 0.112,
                           1.065, 0.793, 0.610, 0.472,
                           0.755,
                           1.602, 1.240, 1.394,
                           1.182, 1.004,
                           1.579, 1.161, 0.843, 0.639, 0.453,
                           1.791, 1.229, 0.996,
                           0.885, 0.746, 0.597,
                           1.197, 0.811, 0.580])

    assert_allclose(extinction.odonnell94(wave, 1.0, 3.1), ref_values,
                    rtol=0.0051, atol=0.)


def test_fitzpatrick99_knots():
    """Test that knots match values in Fitzpatrick (1999) Table 3 for
    fitzpatrick99 function (with R_V = 3.1)"""

    wave = np.array([np.inf, 26500., 12200., 6000., 5470., 4670., 4110.,
                     2700., 2600.])
    x = np.array([0.0, 0.377, 0.820, 1.667, 1.828, 2.141, 2.433, 3.704,
                  3.846])
    
    # A(lambda) values for E(B-V) = 1 or A_V = 3.1
    ref_values = np.array([0.0, 0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265,
                           6.591])

    assert_allclose(extinction.fitzpatrick99(wave, 3.1, unit='aa'), ref_values,
                    rtol=0., atol=0.001)

    # atol = 0.002 because the input values are less precise (e.g., 0.377
    # rather than 1.e4 / 26500.)
    assert_allclose(extinction.fitzpatrick99(x, 3.1, unit='invum'), ref_values,
                    rtol=0., atol=0.002)


def test_fm07():
    wave = np.arange(3000, 9000, 1000)
    ref_values = [ 1.84202329,  1.42645161,  1.13844058,  0.88840962,  0.69220634, 0.54703201]
    assert_allclose(extinction.fm07(wave, 1.), ref_values)

 
def test_calzetti00():
    """Test calzetti against another translation of the same base code"""
    
    wave = np.array([2000., 4000., 8000.])
    flux = np.ones(3)

    new_flux = extinction.apply(extinction.calzetti00(wave, -1., 3.1), flux)
    
    # derived using Julia version of IDL calz_unred
    ref_values = np.array([10.5288, 3.88153, 1.61769])

    assert_allclose(new_flux, ref_values, atol=0.0001)
