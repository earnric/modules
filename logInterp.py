import numpy as np
from scipy.interpolate import interp1d

def log_interp1d(xx, yy, kind='linear'):
    """ Generates an interpolating function, linear by default,
    over the range of log x & log y.
    Return the interpolating function
    """
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = interp1d(logx, logy, kind=kind,fill_value='extrapolate')
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp
