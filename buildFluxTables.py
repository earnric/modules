import numpy as np
from scipy.io import readsav

import re

#from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import astropy 

from matplotlib.pylab import *
import matplotlib.pyplot as plt 
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import NullFormatter

rcParams['figure.figsize'] = (10,8)
rcParams['font.size'] = 22

import sys
sys.path.append('/Users/earnric/Google Drive/ASU/Codes/PythonCode/modules')
import loadfilt as lf
import igm as lyA

import itertools
import os
import subprocess
import glob
import gc
import linecache


#
# Build the file from the flux data and redshift
#
def buildOutFile(z,allAges,ages,wavelns,LperA):
    prevAge = 0.0
    outfile = []
    print("working {} at z={:.1f}".format(fname,z))
    absorb = lyA.lyTauC(z) # Create a lyman forest absorption function for z... ***** NEW truncated absorption
    for age in ages: # Once for each age group in the list
        if age - prevAge < 0.025:
            continue # We don't super-tight spacing in age... 
        print("log age={:.4f}".format(age))
        prevAge = age
        ageCond = (allAges == age) # select records based on log age field... 
        # Convert to freq & Lumin/s/Hz
        freq = (wavelns[ageCond][::-1]).to(u.Hz, equivalencies=u.spectral()) # Reverse array order
        LperHz = (LperA[ageCond] * wavelns[ageCond]**2/astropy.constants.c).to(u.erg/u.s/u.Hz)[::-1] # Reverse array order...

        rsWaveln, rsLperA  = lyA.rsSEDwavelen(wavelns[ageCond], LperA[ageCond], z) # Redshift SED in wavelen
        rsFreq,   rsLperHz = lyA.rsSEDfreq(freq, LperHz, z)      # Redshift SED in freq
        lyForFluxHz = (rsLperHz * absorb(rsWaveln[::-1]))

        # Need to extract only the data for one sed (at a single age)
        # to compute the flux in the filters
        rsWavelnOneAge    = rsWaveln
        rsFreqOneAge      = rsFreq
        lyForFluxHzOneAge = lyForFluxHz
        # Create first part [log age, z, ...]
        # Flux is in erg/s/Hz/cm^2
        jwstNormFlux = np.array([np.trapz(jwstFilters[aFilt](rsWavelnOneAge[::-1])/rsFreqOneAge,rsFreqOneAge).value for aFilt in jwstFilters])
        jwstNormFlux[jwstNormFlux <= 0.0] = -1.0
        jwstFlux = np.array([(np.trapz(jwstFilters[aFilt](rsWavelnOneAge[::-1])*lyForFluxHzOneAge/rsFreqOneAge,rsFreqOneAge)).value
                    for aFilt in jwstFilters])
        jwstFlux[jwstFlux < 1e-90] = 0.0
        jwstFlux   = np.divide(jwstFlux,jwstNormFlux).tolist()

        hubbNormFlux = np.array([np.trapz(hubbleFilters[aFilt](rsWavelnOneAge[::-1])/rsFreqOneAge,rsFreqOneAge).value for aFilt in hubbleFilters])
        hubbNormFlux[hubbNormFlux <= 0.0] = -1.0
        hubbFlux = np.array([(np.trapz(hubbleFilters[aFilt](rsWavelnOneAge[::-1])*lyForFluxHzOneAge/rsFreqOneAge,rsFreqOneAge)).value
                    for aFilt in hubbleFilters])   
        hubbFlux[hubbFlux < 1e-90] = 0.0
        hubbFlux = np.divide(hubbFlux,hubbNormFlux).tolist()

        jhkNormFlux = np.array([np.trapz(jhkFilters[aFilt](rsWavelnOneAge[::-1])/rsFreqOneAge,rsFreqOneAge).value for aFilt in jhkFilters])
        jhkNormFlux[jhkNormFlux <= 0.0] = -1.0
        jhkFlux = np.array([(np.trapz(jhkFilters[aFilt](rsWavelnOneAge[::-1])*lyForFluxHzOneAge/rsFreqOneAge,rsFreqOneAge)).value
                    for aFilt in jhkFilters])
        jhkFlux[jhkFlux < 1e-90] = 0.0
        jhkFlux = np.divide(jhkFlux,jhkNormFlux).tolist()

        aLine = [age] + [z] + jwstFlux + hubbFlux + jhkFlux
        outfile.append(aLine)

    return outfile
#
# Read in Schaerer and SB99 luminosity files and generate flux in filter.
# These tables only consider Ly-forest absorption -- and redshift.
# User needs to apply reddening.
#

def buildFilterFluxFiles():
    """Creates flux-in-filter files from the SED.
    This is untested but a copy of the python notebook that works.
    """
    jwstFilters   = lf.loadJWSTFilters(suppress=True)
    hubbleFilters = lf.loadHubbleFilters(suppress=True)
    
    lamRange      = np.logspace(1.95,5.7,5500)

    schaererPath = '/Users/earnric/Research/Research-Observability/Software-Models/Schaerer/'
    schaererDirs = ['pop3_TA/','pop3_TE/','e-70_mar08/','e-50_mar08/']
    Zs           = [0.0, 0.0, 1.0e-7, 1.0e-5]
    schaererPopFilePattern  = 'pop3_ge0_log?_500_001_is5.[0-9]*' # is5 files have ages in step of 1 Myr
    schaererLowZFilePattern = 'e-?0_sal_100_001_is2.[0-9]*'      # is2 files have ages in step of 0.05 dex

    # Load the schaerer files... 
    # Note that due to spacing (in age), there are sometimes two files that 
    # are SEDS for the same time-stamp! Skip the second one (and third!)
    lastAge = 0.0
    for i, (Z, schaererDir) in enumerate(zip(Zs,schaererDirs)):
        if schaererDir.startswith('pop3'):
            schaererFilePattern = schaererPath + schaererDir + schaererPopFilePattern  # Pop III files, 1 Myr spacing
        else:
            schaererFilePattern = schaererPath + schaererDir + schaererLowZFilePattern # Low Z files, 0.05 dex spacing

        schaererFiles   = glob.glob(schaererFilePattern)  # All the files in the dir... 
        schaererFiles   = [a for a in schaererFiles if not re.search('\.[1-2][0-9][b]*$', a)] # remove .1? and .2? files
        schaererAges    = [linecache.getline(file,13) for file in schaererFiles]    # Get the line with the (log) age... 
        schaererAges    = np.array([float(sa[30:]) for sa in schaererAges],dtype=float)         # Log age starts at position 30

        schaererData    = np.array([np.loadtxt(file,skiprows=16) for file in schaererFiles])
        ageSortIndxes   = schaererAges.argsort()          # Array of indices to sort things by age...

        schaererData    = schaererData[ageSortIndxes]
        schaererAges    = schaererAges[ageSortIndxes]
        print(len(schaererData))
        print(len(schaererAges),schaererAges)

        # Ignore data files with the same age! This occurs in the popIII dirs
        # because the timestep is smaller than the age-resolution printed in the file
        # Hence we get 2 files with different data but the same time stamp
        lastAge = 0.0
        schaererDataGood = []
        schaererAgesGood = []
        for ii,(sd,age) in enumerate(zip(schaererData,schaererAges)):
            if age == lastAge:
                # Remove it
                continue
            lastAge = age
            schaererDataGood.append(sd)
            schaererAgesGood.append(float(age))

        # The following builds an array of arrays (one for each age) with each array's entries:
        # log age, Z, waveln, lum/A
        allSchaererData = [np.insert(sed[:,[0,2]],[0],[[anAge] for ii in range(0,len(sed))], axis=1) 
            for anAge,sed in zip(schaererAgesGood, schaererDataGood)]
        allSchaererData = np.array(allSchaererData).reshape(len(allSchaererData)*len(allSchaererData[0]),3)
        if i == 0:
            pop3TA = allSchaererData # may need a np.copy(...) here... ??
        elif i == 1:
            pop3TE = allSchaererData
        elif i == 2:
            Zem7 = allSchaererData
        elif i == 3:
            Zem5 = allSchaererData
        # We now have:
        # [[log age, waveln, flux], [], ...]

    # Generate the filter-flux data... 
    
    redshifts     = [2.0,3.0,4.0,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12.0,13.0]
    schaererList  = [pop3TA, pop3TE, Zem7, Zem5]
    schaererNames = ["pop3TA", "pop3TE", "Zem7", "Zem5"]
    hdr = 'LogAge, redshift, '
    hdr += ', '.join([jFilt for jFilt in jwstFilters])
    hdr += ', '
    hdr += ', '.join([hFilt for hFilt in hubbleFilters])
    hdr += ', '
    hdr += ', '.join([jFilt for jFilt in jhkFilters])
    fmtStr = '%.4f, %.1f, '
    fmtStr += ', '.join(['%.4e' for ii in arange(len(jwstFilters)+len(hubbleFilters)+len(jhkFilters))]) 
    for schaererData, fname in zip(schaererList,schaererNames):
        ages = np.unique(schaererData[:,0]) # Get the list of ages
        ages = ages[ages <= 9.01] # We don't need flux for stars older than 1.02 Gyr
        # Get the log age, wavelength & Luminosity/s/Ang
        allAges = schaererData[:,0]
        wavelns = schaererData[:,1] * u.Angstrom
        LperA   = schaererData[:,2] * u.erg / u.second / u.angstrom
        for z in redshifts:
            buildOutFile(z,allAges,ages,wavelns,LperA)
            filename = fname + "_" + str(z) + ".gz"
            print('writing {}'.format(filename))
            np.savetxt(filename, outfile, fmt=fmtStr, delimiter=', ', header=hdr)

    print("Done with Schaerer...")
    # Process SB99 files...
    

    # SB99 format:     TIME [YR]    WAVELENGTH [A]   LOG TOTAL  LOG STELLAR  LOG NEBULAR  [ERG/SEC/A]
    # REMEMBER, SB99 data is for a population of 1e6 M_sun
    SB99Path = '/Users/earnric/OneDrive/STARBURST99-runs/' # Home computer dir... 
    SB99Dirs = ['padova0004-op/','padova004-op/','padova008-op/','padova02-op/']
    Zs       = [0.0004, 0.004, 0.008, 0.02]
    SB99FilePat = 'padova*.spectrum1'

    # SB99 format:     TIME [YR]    WAVELENGTH [A]   LOG TOTAL  LOG STELLAR  LOG NEBULAR  [ERG/SEC/A]
    for i, (Z, SB99Dir) in enumerate(zip(Zs,SB99Dirs)):
        SB99FilePattern = SB99Path + SB99Dir + SB99FilePat 
        SB99Files   = glob.glob(SB99FilePattern)  # All the files in the dir... should be one!
        if len(SB99Files) != 1:
            print('Error: too many files in an SB99 dir! - ',SB99Path + SB99Dir)
            sys.exit()
        SB99Data    = np.loadtxt(SB99Files[0],skiprows=6)
        if i == 0:
            SB990004 = np.dstack((np.log10(SB99Data[:,0]),SB99Data[:,1],10**(SB99Data[:,2]-6.0))).reshape(len(SB99Data[:,1]),3)
        elif i == 1:
            SB99004 = np.dstack((np.log10(SB99Data[:,0]),SB99Data[:,1],10**(SB99Data[:,2]-6.0))).reshape(len(SB99Data[:,1]),3)
        elif i == 2:
            SB99008 = np.dstack((np.log10(SB99Data[:,0]),SB99Data[:,1],10**(SB99Data[:,2]-6.0))).reshape(len(SB99Data[:,1]),3)
        elif i == 3:
            SB9902 = np.dstack((np.log10(SB99Data[:,0]),SB99Data[:,1],10**(SB99Data[:,2]-6.0))).reshape(len(SB99Data[:,1]),3)
        # We now have:
        # [[log age, waveln, flux], [], ...]


        # Generate the filter-flux data...
        
    redshifts    = [2.0,3.0,4.0,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12.0,13.0]
    # redshifts    = [14.0,15.0,16.0]
    sb99List  = [SB990004, SB99004, SB99008, SB9902]
    sb99Names = ["SB990004", "SB99004", "SB99008", "SB9902"]
    # sb99List  = [SB9902]
    # sb99Names = ["SB9902"]
    print(hdr)
    for sb99Data, fname in zip(sb99List,sb99Names):
        ages = np.unique(sb99Data[:,0]) # Get the list of ages
        # Get the log age, wavelength & Luminosity/s/Ang
        allAges = sb99Data[:,0]
        wavelns = sb99Data[:,1] * u.Angstrom
        LperA   = sb99Data[:,2] * u.erg / u.second / u.angstrom
        for z in redshifts:
            buildOutFile(z,allAges,ages,wavelns,LperA)
            filename = fname + "_" + str(z) + ".gz"
            print('writing {}'.format(filename))
            np.savetxt(filename, outfile, fmt=fmtStr, delimiter=', ', header=hdr)
        
    print("Done with SB99...")
    return 
