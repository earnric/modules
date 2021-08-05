# This is the HI Photoionization Cross Section

def sigHI(E):
  Eth    = 1.360E+01
  Emax   = 5.000E+04
  E0     = 4.298E-01

  sigma0 = 5.475E+04
  ya     = 3.288E+01
  P      = 2.963E+00
  yw     = 0
  y0     = 0
  y1     = 0
  x      = E/E0-y0
  y      = np.sqrt(x*x+y1*y1)
  F      = ((x-1)*(x-1)+yw*yw)*np.power(y,0.5*P-5.5)*np.power((1+np.sqrt(y/ya)),-P)
  if (E<Eth): 
    F = 0
  if (E>Emax):
    F = 0
  sigma  = sigma0*F*1.0E-18
  return sigma

# This is the HeI Photoionization Cross Section
def sigHeI(E):
  Eth    = 2.459E+01
  Emax   = 5.000E+04
  E0     = 1.361E+01
  sigma0 = 9.492E+02
  ya     = 1.469E+00
  P      = 3.188E+00
  yw     = 2.039E+00
  y0     = 4.434E-01
  y1     = 2.136E+00
  x      = E/E0-y0
  y      = np.sqrt(x*x+y1*y1)
  F      = ((x-1)*(x-1)+yw*yw)*np.power(y,0.5*P-5.5)*np.power((1+np.sqrt(y/ya)),-P)
  if (E<Eth):  
    F = 0
  if (E>Emax):
    F = 0
  sigma  = sigma0*F*1.0E-18
  return sigma

# This is the CVI Photoionization Cross Section
def sigCVI(E):
  Eth    = 4.900E+02
  Emax   = 5.000E+04
  E0     = 1.548E+01
  sigma0 = 1.521E+03
  ya     = 3.288E+01
  P      = 2.963E+00
  yw     = 0
  y0     = 0
  y1     = 0
  x      = E/E0-y0
  y      = np.sqrt(x*x+y1*y1)
  F      = ((x-1)*(x-1)+yw*yw)*np.power(y,0.5*P-5.5)*np.power((1+np.sqrt(y/ya)),-P)
  if (E<Eth):
    F = 0
  if (E>Emax): 
    F = 0
  sigma  = sigma0*F*1.0E-18
  return sigma
