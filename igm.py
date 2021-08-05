import numpy as np
import astropy
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from scipy.interpolate import interp1d
import extinction

# Setup the cosmology used in the runs
# This provides comoving distance, lumin distance, etc...
myCosmo = FlatLambdaCDM(H0=71.0, Om0=0.267, Ob0=0.0449, name="myCosmo")


## Remember, reddening happens in the galaxy's rest frame
## The Calzetti model only works down to a rest frame waveln of 1200 Ang
## Since my filter fluxes are all in the observer's frame, need to divide
## by (1+z) to get the absorption (in mags) in the center of each filter.

## hubbFiltCenters = {
filtCenters = {
    "F098M_WFC3": 9800.0,
    "F125W_WFC3": 12500.0,
    "F606W_ACS": 6060.0,
    "F775W_ACS": 7750.0,
    "F435W_ACS": 4350.0,
    "F225W_WFC3": 2250.0,
    "F105W_WFC3": 10500.0,
    "F850LP_ACS": 8500.0,
    "F160W_WFC3": 16000.0,
    "F275W_WFC3": 2750.0,
    "F336W_WFC3": 3360.0,
    ##     }
    ## jwstFiltCenters = {
    "f200w": 20000.0,
    "f182m": 18200.0,
    "f300m": 30000.0,
    "f410m": 41000.0,
    "f335m": 33500.0,
    "f115w": 11500.0,
    "f430m": 43000.0,
    "f460m": 46000.0,
    "f356w": 35600.0,
    "f150w": 15000.0,
    "f162m": 16200.0,
    "f360m": 36000.0,
    "f277w": 27700.0,
    "f140m": 14000.0,
    "f070w": 7000.0,
    "f480m": 48000.0,
    "f444w": 44400.0,
    "f090w": 9000.0,
    "f210m": 21000.0,
    "f250m": 25000.0,
    ##     }
    ## jhkFiltCenters  = {
    "H": 1650.0,
    "J": 1250.0,
    "Ks": 2150.0,
}
## def dustAbsFlux():
##     """ Computes the change in flux at each wavelen for the specified
##     values of A_v[=1.0] and R_v[=4.05]. Returns an interpolating function
##     over the input range, 0.0 outside.
##     """
##     magAbsorpValues = extinction.calzetti00(wavelens,A_v,R_v)
##     absorpFunc = interp1d(np.arange(1200.,25000.),10**(-0.4*(magAbsorpValues-48.6))
##     return absorpFunc

## def dustAbsMag(wavelens,A_v=1.0,R_v=4.05):
##     """ Computes the change in mag at each wavelen for the specified
##     values of A_v[=1.0] and R_v[=4.05]. Returns an interpolating function
##     over the input range, 0.0 outside.
##     """
##     magAbsorpValues = extinction.calzetti00(wavelens,A_v,R_v)
##     absorpFunc = interp1d(np.arange(1200.,25000.),magAbsorpValues)
##     return absorpFunc


def filterFlux(filterFunc, waveln, LperA, z):
    """Computes the flux in a filter function 'filterFunc' at a given redshift.
    Inputs are a filterfuncion that takes anstroms and returns a fraction {0,1},
    a set of wavelens (units assigned to raw input will be ang)
    a set of corresponding luminosities per ang (units will be assigned erg/s/A)
    Returns the flux in units of erg/s/cm^2/Hz
    """
    waveLnU = waveln * u.angstrom  # Assuming waveln is in Angstroms
    LperAU = (
        LperA * u.erg / u.second / u.angstrom
    )  # Assuming lumonisity is in erg/s/A[/M]
    freqU = (waveLnU).to(u.Hz, equivalencies=u.spectral())[
        ::-1
    ]  # Reverse order for freq...
    LperHzU = (LperAU * waveLnU ** 2 / astropy.constants.c).to(u.erg / u.s / u.Hz)[
        ::-1
    ]  # L * \lambda^2/c -> [erg/s/Hz]

    # Redshift the SED and compute the Lyman forest absorption
    rsWaveln, rsLperA = rsSEDwavelen(waveLnU, LperAU, z)
    rsFreq, rsLperHz = rsSEDfreq(freqU, LperHzU, z)
    igmAbsorp = lyTau(
        z
    )  # A function that return transmission % as a function of wavelen.
    lyForFluxHz = rsLperHz * igmAbsorp(rsWaveln[::-1])

    # Compute the flux in the filter:
    # \int {R(\nu) f(\nu)/\nu} / \int R(\nu) /\nu}
    # Note that "/ \nu" is really "/ (h \nu)" so we are dividing out the energy per Hz in the filter
    # Since h is constant we can pull it out and ignore (divides away).
    normFlux = np.trapz(filterFunc(rsWaveln[::-1]) / rsFreq, rsFreq)
    normFlux[
        normFlux <= 0.0
    ] = -1.0  # Incase we have redshifted the filter out of the freq range
    filterFlux = (
        np.trapz(filterFunc(rsWaveln[::-1]) * lyForFluxHz / rsFreq, rsFreq) / normFlux
    )
    return filterFlux.value


def filterFluxC(filterFunc, waveln, LperA, z):
    """Computes the flux in a filter function 'filterFunc' at a given redshift.
    Inputs are a filterfuncion that takes anstroms and returns a fraction {0,1},
    a set of wavelens (units assigned to raw input will be ang)
    a set of corresponding luminosities per ang (units will be assigned erg/s/A)
    Returns the flux in units of erg/s/cm^2/Hz
    """
    waveLnU = waveln * u.angstrom  # Assuming waveln is in Angstroms
    LperAU = (
        LperA * u.erg / u.second / u.angstrom
    )  # Assuming lumonisity is in erg/s/A[/M]
    freqU = (waveLnU).to(u.Hz, equivalencies=u.spectral())[
        ::-1
    ]  # Reverse order for freq...
    LperHzU = (LperAU * waveLnU ** 2 / astropy.constants.c).to(u.erg / u.s / u.Hz)[
        ::-1
    ]  # L * \lambda^2/c -> [erg/s/Hz]

    # Redshift the SED and compute the Lyman forest absorption
    rsWaveln, rsLperA = rsSEDwavelen(waveLnU, LperAU, z)
    rsFreq, rsLperHz = rsSEDfreq(freqU, LperHzU, z)
    igmAbsorp = lyTauC(
        z
    )  # A function that return transmission % as a function of wavelen.
    lyForFluxHz = rsLperHz * igmAbsorp(rsWaveln[::-1])

    # Compute the flux in the filter:
    # \int {R(\nu) f(\nu)/\nu} / \int R(\nu) /\nu}
    # Note that "/ \nu" is really "/ (h \nu)" so we are dividing out the energy per Hz in the filter
    # Since h is constant we can pull it out and ignore (divides away).
    normFlux = np.trapz(filterFunc(rsWaveln[::-1]) / rsFreq, rsFreq)
    normFlux[
        normFlux <= 0.0
    ] = -1.0  # Incase we have redshifted the filter out of the freq range
    filterFlux = (
        np.trapz(filterFunc(rsWaveln[::-1]) * lyForFluxHz / rsFreq, rsFreq) / normFlux
    )
    return filterFlux.value


def lyTau(z):
    """ Generates the optical depth due to Lyman line and continuum absorption
    over a range of wavelengths out to a specified redshift & returns
    an interpolating function that provides the fraction of flux transmitted per A
    at the specified wavelengths.
    Parameter is the redshift (z).
    The function returns a scipy interp1d function which takes an array of
    wavelengths (ang) returning an array of fractions at each lambda.
    NOTE that the interpolating function return 0's outside the range [89.1, 5011872] ang
    """
    lamRange = np.logspace(1.0, 6.7, 15000)
    lineStren = np.array(
        [
            3.6e-3,
            1.7e-3,
            1.2e-3,
            9.3e-4,
            8.0e-4,
            7.0e-4,
            6.2e-4,
            5.7e-4,
            5.2e-4,
            4.8e-4,
            4.5e-4,
            4.2e-4,
            3.9e-4,
            3.7e-4,
            3.5e-4,
            3.3e-4,
            3.16e-4,
        ]
    )
    upE = np.arange(
        2.0, len(lineStren) + 2.0
    )  # Need + 2 since we want to include the end in the sequence
    taus = np.zeros(len(lamRange), dtype=float)
    lamLi = 911.7535 / (1.0 - upE ** -2)  # Wavelength of line absorbers, rest frame

    rs = 1.0 + z
    # Compute line & sum line optical depth
    for lli, ls in zip(lamLi, lineStren):
        for indx, il in enumerate(lamRange):
            if il < rs * lli:
                tauli = ls * (il / lli) ** 3.46
            else:
                tauli = 0.0
            taus[indx] = taus[indx] + tauli

    # Compute continuum optical depth
    for indx, il in enumerate(lamRange):
        xc = il / 911.7535
        if (
            il < 911.7535 * rs
        ):  # If the wavelen < Ly limit in the rest frame... do absorption
            tauc = (
                (xc ** 3.0) / 4.0 * (rs ** 0.46 - xc ** 0.46)
                + 9.4 * xc ** 1.5 * (rs ** 0.18 - xc ** 0.18)
                - 0.7 * xc ** 3 * (xc ** -1.32 - rs ** -1.32)
                - 0.023 * (rs ** 1.68 - xc ** 1.68)
            )
        else:
            tauc = 0.0
        taus[indx] = taus[indx] + tauc

    return interp1d(
        lamRange,
        np.exp(-1.0 * taus),
        kind="slinear",
        bounds_error=False,
        fill_value=0.0,
    )


def lyTauC(z):
    """ Generates the optical depth due to Lyman line and continuum absorption
    over a range of wavelengths out to a specified redshift & returns
    an interpolating function that provides the fraction of flux transmitted per A
    at the specified wavelengths. This version ensures that all flux shortward of the
    lyman limit (in the rest frame) is zero'd out.
    Parameter is the redshift (z).
    The function returns a scipy interp1d function which takes an array of
    wavelengths (ang) returning an array of fractions at each lambda.
    NOTE that the interpolating function return 0's outside the range [89.1, 5011872] ang
    """
    lamRange = np.logspace(1.0, 6.7, 15000)
    lineStren = np.array(
        [
            3.6e-3,
            1.7e-3,
            1.2e-3,
            9.3e-4,
            8.0e-4,
            7.0e-4,
            6.2e-4,
            5.7e-4,
            5.2e-4,
            4.8e-4,
            4.5e-4,
            4.2e-4,
            3.9e-4,
            3.7e-4,
            3.5e-4,
            3.3e-4,
            3.16e-4,
        ]
    )
    upE = np.arange(
        2.0, len(lineStren) + 2.0
    )  # Need + 2 since we want to include the end in the sequence
    taus = np.zeros(len(lamRange), dtype=float)
    lamLi = 911.7535 / (1.0 - upE ** -2)  # Wavelength of line absorbers

    rs = 1.0 + z
    # Compute line & sum line optical depth
    for lli, ls in zip(lamLi, lineStren):
        for indx, il in enumerate(lamRange):
            if il < rs * lli:
                tauli = ls * (il / lli) ** 3.46
            else:
                tauli = 0.0
            taus[indx] = taus[indx] + tauli

    # Compute continuum optical depth
    for indx, il in enumerate(lamRange):
        xc = il / 911.7535
        if il < 911.7535 * rs:
            tauc = 100.0
        else:
            tauc = 0.0
        taus[indx] = taus[indx] + tauc

    xmissFunc = interp1d(
        lamRange, np.exp(-1 * taus), kind="slinear", bounds_error=False, fill_value=0.0
    )

    return xmissFunc


def rsSEDwavelen(inLams, inFlux, z=0):
    """ Redshifts a SED.
    Parameters are wavelength array (ang) and the L/A (luminosity per angstrom). The
    final parameter is the redshift, z. Default is 0.
    The routine returns two arrays: the redshifted wavelength range and the
    corresponding fluxes. The units on the flux are [erg/s/A/cm^2] (assuming input L
    was in [erg/s/A]
    """
    if z < 0:
        print("You must specify z >= 0")
        return
    lamNew = inLams * (1.0 + z)
    fluxNew = (
        inFlux
        / (4.0 * np.pi * myCosmo.luminosity_distance(z).to(u.cm) ** 2.0)
        / (1.0 + z)
    )
    return lamNew, fluxNew.value


def rsSEDfreq(inFreq, inFlux, z=0):
    """ Redshifts a SED.
    Parameters are freq array (Hz) and the L/Hz (luminosity per Hertz). The
    final parameter is the redshift, z. Default is 0.
    The routine returns two arrays: the redshifted freq range and the
    corresponding fluxes. The units on the flux are [erg/s/Hz/cm^2] (assuming input L
    was in [erg/s/Hz]
    """
    if z < 0:
        print("You must specify z >= 0")
        return
    newFreq = inFreq / (1.0 + z)
    fluxNew = (
        inFlux
        / (4.0 * np.pi * myCosmo.luminosity_distance(z).to(u.cm) ** 2.0)
        * (1.0 + z)
    )
    return newFreq, fluxNew.value


def AtoHz(waveLn, LperA):
    """ Compute the equiv freq and LperHz from wavelen and LperA
    Input should not have astropy units attached, but will be assumed
    to be Ang and erg/s/Ang
    Return two arrays: freq (Hz) and LperHz (erg/s/Hz)
    """
    waveLn = waveLn * u.Angstrom
    LperA = LperA * u.erg / u.s / u.Angstrom
    freq = waveLn.to(u.Hz, equivalencies=u.spectral())[::-1]  # Reverse the order
    LperHz = (LperA * (waveLn ** 2 / astropy.constants.c)).to(u.erg / u.s / u.Hz)[::-1]  # Reverse order...
    return freq.value, LperHz.value
