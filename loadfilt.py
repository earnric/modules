"""This modules provides routines to load JWST filters and convolve (TBD) them with fluxes
in a specified array of wavelengths.
It provides the following routines:
    loadJWSTFilters()
"""
import numpy as np
from scipy.interpolate import interp1d
import glob

def loadJWSTFilters(path="/Users/earnric/Google Drive/ASU/Research-Observability/nircam/",suppress=False):
    """This routine loads the JWST NIRCAM filter curves from the default (or specified) location
    and creates a dictionary of linear functions fitted to the filter curves. The splines
    are returned as a dictionary array (indexed by the filter name) and can be used to determine
    the filters' transmission of light of a given wavelength.
    e.g. - filterResponses = frofilt.loadJWSTFilters()
           filterResponses['f070w']([6000,6005,6010])
    would return an array the transmission fraction for the f070w filter at the list of specified wavelengths
    Parameters: path - the path to the NIRCAM filter files; suppress - flag indicating that the
    routine should (False) or should not (True) print the dictionary keys for the returned set of
    spline fits.
    """
    # Load the filter files
    filePattern = path+'f*.trans'
    filterFiles = glob.glob(filePattern)
    filters     = np.array([np.loadtxt(file,skiprows=1) for file in filterFiles])

    # Create a cubic spline fit for each table of filter data
    # Entry 0 is line #, 1 is lambda in Ang, 2 is response in %
    filterFits = np.array([interp1d(filter[:,1],filter[:,2],kind='linear',bounds_error=False,fill_value=0.0)
                           for filter in filters])

    # Create a dictionary of filter response fits... 
    filterResponse = dict(zip([filename[-11:-6] for filename in filterFiles],filterFits))
    if not suppress:
        print('dict keys: ')
        print('\n'.join(filterResponse.keys()))

    return filterResponse

def loadHubbleFilters(path="/Users/earnric/Google Drive/ASU/Research-Observability/hubble/",suppress=False):
    """This routine loads the Hubble filter curves from the default (or specified) location
    and creates a dictionary of linear functions fitted to the filter curves. The splines
    are returned as a dictionary array (indexed by the filter name) and can be used to determine
    the filters' transmission of light of a given wavelength.
    e.g. - filterResponses = frofilt.loadHubbleFilters()
           filterResponses['f070w']([6000,6005,6010])
    would return an array the transmission fraction for the f070w filter at the list of specified wavelengths
    Parameters: path - the path to the NIRCAM filter files; suppress - flag indicating that the
    routine should (False) or should not (True) print the dictionary keys for the returned set of
    spline fits.
    """
    # Load the filter files
    filePattern = path+'F*.*'
    filterFiles = glob.glob(filePattern)
    filters     = np.array([np.loadtxt(file,skiprows=1) for file in filterFiles])

    # Create a cubic spline fit for each table of filter data
    # Entry 0 is line #, 1 is lambda in Ang, 2 is response in %
    filterFits = np.array([interp1d(filter[:,0],filter[:,1],kind='linear',bounds_error=False,fill_value=0.0)
                           for filter in filters])

    # Create a dictionary of filter response fits...
    filterResponse = dict(zip([filename[filename.rfind('/')+1:filename.rfind('.')] # extract the filter name
                               for filename in filterFiles],filterFits))
    if not suppress:
        print('dict keys: ')
        print('\n'.join(filterResponse.keys()))

    return filterResponse


def loadJHKFilters(path="/Users/earnric/Google Drive/ASU/Codes/PythonCode/modules/flamingo/",suppress=False):
    """This routine loads the J, H & K filter curves from the default (or specified) location
    and creates a dictionary of linear functions fitted to the filter curves. The splines
    are returned as a dictionary array (indexed by the filter name) and can be used to determine
    the filters' transmission of light of a given wavelength.
    e.g. - filterResponses = frofilt.loadJHKFilters()
           filterResponses['J']([6000,6005,6010])
    would return an array with the transmission fraction for the J filter at the list of specified wavelengths
    Parameters: path - the path to the NIRCAM filter files; suppress - flag indicating that the
    routine should (False) or should not (True) print the dictionary keys for the returned set of
    spline fits.
    """
    # Load the filter files
    filePattern = path+'FLAMINGOS.BARR*.*'
    filterFiles = glob.glob(filePattern)
    filters     = np.array([np.loadtxt(file) for file in filterFiles])

    # There are negative numbers in some of these filter specs! Fix 'em
    for filt in filters:
        cond = (filt[:,1] < 0.50)
        filt[:,1][cond] = 0.0 # NOTE THE ORDER is important! Slice then condition...

    # Create a cubic spline fit for each table of filter data
    # Entry 0 is line #, 1 is lambda in Ang, 2 is response in %
    filterFits = np.array([interp1d(filter[:,0],filter[:,1]/100.0,kind='linear',bounds_error=False,fill_value=0.0)
                           for filter in filters])

    # Create a dictionary of filter response fits...
    filterResponse = dict(zip([filename[filename.rfind('RR.')+3:filename.rfind('.ColdWitness')] # extract the filter name
                               for filename in filterFiles],filterFits))
    if not suppress:
        print('dict keys: ')
        print('\n'.join(filterResponse.keys()))

    return filterResponse

