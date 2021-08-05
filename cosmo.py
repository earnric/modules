import numpy as np
import math as ma
import scipy.integrate as integrate

""" Use -- from astropy.cosmology import FlatLambdaCDM
with the proper parameters:
# Setup the cosmology used in the runs
# This provides comoving distance, lumin distance, etc... 
cosmo = FlatLambdaCDM(H0=71.0, Om0=0.267, Ob0=0.0449,name='myCosmo')
"""

speedOfLight = 3e5 # In km/s ... This will work with H specified in km/s/Mpc
cmPerMpc     = 3.08e24 # Not precise, but this is what RAMSES uses...

def ageAtz(H0=71.0, z=0.0):
    """ageAtz - returns the age, in Myr, at a specified H and redshift
    First parameter is the Hubble const, second is z
    Return Myr
    """
    H = H0 * 1.02201e-6 # in Myr^-1 ... ??? What? Must be inverse seconds
    # Analytic solution to the integral (1/(Om (1+z)^5 + Ol(1+z)^2)^.5)
    # Here the cosmology is Ol=0.733, Om=0.267
    return 1.0/H * ((0.4023577624*(1. + z)**(5/2) * np.sqrt(1. + 2.7453183520/(1. + z)**3)* \
          np.arcsinh(1.65690022393/(1. + z)**(3/2))))/ \
          np.sqrt((1. + z)**2 * (0.733 + 0.267 * (1. + z)**3))


def arcAtz(H0=71.0,z=0.0,om=0.3,ol=0.7):
    """arcAtz - returns the size (proper distance - physical size in rest frame) of an arcsec
    at the specified reshift in kpc.
    Parameters are Hubble const (in km/s/Mpc), redshift, optional Om and Ol
    Return kpc/arcsec
    """
    sol = integrate.quad(lambda z: 1/np.sqrt(om*(1+z)**3+ol), 0,z)
    # Since c/H gives us Mpc, need to multiply final result by 1000
    # to get to kpc
    result = speedOfLight/H0 * sol[0]/(1.+z) * ma.pi/648000.0 # last term is radians per arcsec
    return result * 1000.0
    
def dp(H0=71.0,z=0.0,om=0.3,ol=0.7):
    """ Computer the proper distance to z, in Mpc
    Parameters are Hubble const (in km/s/Mpc), redshift, optional Om and Ol
    Returns Mpc
    """
    sol = integrate.quad(lambda z: 1.0/np.sqrt(om*(1.+z)**3+ol), 0, z)
    return speedOfLight/H0 * sol[0]

    
def dl(H0=71.0,z=0.0,om=0.3,ol=0.7):
    """ Computer the luminosity distance to z, in Mpc
    Parameters are Hubble const (in km/s/Mpc), redshift, optional Om and Ol
    Returns Mpc
    """
    sol = dp(H0, z, om, ol) * (1+z)
    return sol

    
def da(H0=71.0,z=0.0,om=0.3,ol=0.7):
    """Copmute the angular size distance at z in arcsec
    Parameters are Hubble const (in km/s/Mpc), redshift, optional Om and Ol
    Returns Mpc
    """
    sol = dl(H0, z, om, ol)/(1+z)**2
    return sol
    
