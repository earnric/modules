"""
Routines for computing magnitudes.
Depends on astropy.
"""
import 
def mag(flux):
    """Compute m_ab from flux"""
    return -2.5 * np.log10(flux) - ABref

def compMags(z):
    """Computes M_ab from flux arrays."""
    # Distance modulus ...
    convertToM = -5.0*np.log10((cosmo.luminosity_distance(z)/(10 * u.pc)))
    # Or, -5 log (d/10 * (1+z)) ... Need the (1+z) to get to luminosity distance
    print("Convert to Abs Mag. DM = {:.2f} @ z={:.1f}".format(convertToM,z))
    absMag   = mag(fluxes[z]['1500A'])+convertToM
    pop3Mag  = mag(fluxes[z]['1500A_P3'])+convertToM
    nmAbsMag = mag(fluxes[z]['1500A_NM'])+convertToM
    return absMag, pop3Mag,nmAbsMag

def getHaloMasses(z):
    return fluxes[z]['MstarMsun']

def getHaloP3Masses(z):
    return fluxes[z]['M3StarMsun']

def getHaloRadii(z):
    return fluxes[z]['r_v']
