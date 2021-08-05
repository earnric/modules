##############################################################
# Rick Sarmento - AST 598, Fall 2014
# Numerically integrating the age fo the universe
# from the FRW eqn. Using the benchmark cosmology
# from class.
##############################################################

from pylab import *
from scipy import integrate

# Conversion factors
kmToCm = 100000.0
MpcToCm = 3.086E24
sInGyr = 3.15e16

# Constants
c  = 29979245800.0 # cm/s
H0 = 71.0 * kmToCm/MpcToCm # H0 in "per sec"
omegaR = 8.4e-5
omegaM = 0.267
omegaL = 0.733
omega0 = omegaR + omegaM + omegaL

##############################################################
# Computes the (comoving) proper distance for a specified
# cosmology (above).
# Parameterized in a, the scale factor.
# All distances returned relative to a = 1 
##############################################################
def dp(a):
    dpe = lambda a: 1.0/(H0 * sqrt(omegaR + omegaM * a + omegaL * a**4.0 + (1.0 - omega0) * a**2.0))
    dp1 = integrate.quadrature(dpe, a, 1.0, maxiter=150)[0]
    return (c * dp1)/MpcToCm
##############################################################

##############################################################
# Utility routine
##############################################################
def zToa(z):
    return (1.0 / (1.0 + z))
##############################################################


def main():

    print "Numerically integrating FRW..."
    print ("H0 = %.3e 1/s" %H0)
    frw = lambda a: 1.0/(H0 * sqrt(omegaR/a**2.0 + omegaM/a + omegaL * a**2.0 + 1.0-omega0))
    age   = integrate.quadrature(frw, 0.0, 1.0,maxiter=150)
    zp5   = integrate.quadrature(frw, 0.0, 1.0/(1.0+0.5),maxiter=150)   # a = 1/(1+z)
    zp75  = integrate.quadrature(frw, 0.0, 1.0/(1.0+0.75),maxiter=150)
    z1    = integrate.quadrature(frw, 0.0, 1.0/(1.0+1.0),maxiter=150)
    z2    = integrate.quadrature(frw, 0.0, 1.0/(1.0+2.0),maxiter=150)
    
    print ("Age of the benchmark universe in Gyr %5.2f" % (age[0]/sInGyr))
    print ("Age of the benchmark universe @ z=0.5 %5.2f" % (zp5[0]/sInGyr))
    print ("Age of the benchmark universe @ z=0.75 %5.2f" % (zp75[0]/sInGyr))
    print ("Age of the benchmark universe @ z=1.0 %5.2f" % (z1[0]/sInGyr))
    print ("Age of the benchmark universe @ z=2.0 %5.2f" % (z2[0]/sInGyr))

    print ("Proper distance to")
    print ("z = 0.5  is %8.2f" % dp(zToa(0.5))) 
    print ("z = 0.75 is %8.2f" % dp(zToa(0.75))) 
    print ("z = 1.0  is %8.2f" % dp(zToa(1.0))) 
    print ("z = 2.0 is %8.2f" % dp(zToa(2.0))) 

    print ("Luminosity distance to") # Comoving distance * (1+z)
    print ("z = 0.5  is %8.2f" % (dp(zToa(0.5)) * (1.0+0.5)))
    print ("z = 0.75 is %8.2f" % (dp(zToa(0.75)) * (1.0+0.75)))
    print ("z = 1.0  is %8.2f" % (dp(zToa(1.0)) * (1.0+1.0)))
    print ("z = 2.0 is %8.2f" % (dp(zToa(2.0)) * (1.0+2.0)))


##############################################################
# Call the main program... 
##############################################################
if __name__ == "__main__":
    main()
