"""
Loads star formation rate data from the file with format:
    z time totSPM totPIIIM totPIIM totCPIIIM totZPM
mapping to
    redshift, proper time, 

Data is accessed by sfrd, sfrdP3, sfrdCP3 fields
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import astropy


class nSFR:
    """
    This class reads in the data generated by the genStarData-RT 
    python program. The constructor expects the size of the box,
    the input file path & nbame, as well as the hubble parameter --
    which defaults to 0.71.

        size : float
            Boxsize in Mpc/h
        prefix : str ("./")
            path to the star particle data summary file
        hubble : float (0.71)
            H0/(100 km/s/Mpc)
    """

    def __init__(self, size, prefix="./", hubble=0.71,verbose=False):
        """
        This constructor sets up the object's cosmology as
        well as the root path to the file loaded with the
        load method.
        """
        self.prefix = prefix
        self.hub = hubble
        self.size = size
        self.cosmo = FlatLambdaCDM(
            H0=hubble * 100.0, Om0=0.267, Ob0=0.0449, name="myCosmo")
        self.cm_vol = (self.size / self.hub) ** 3
        if verbose:
            print("Object h={:.2f} with path prefix {}".format(self.hub, self.prefix))
        self.verbose = verbose
        return

    def load(self, file="starData.txt"):
        """
        Loads data expecting the following columns:
            z time toSPM totPIIIM totPIIM totCPIIIM totZPM

        User specifies the filename to read.
            file : (starData.txt) the file to load

        Once the data is loaded, the user can access the fields below:
            redshift : an array of redshifts at which the SFRD value is computed
            time     : the proper time associated with z, accounting for 
                       my standard cosmology
            sfrdP3   : the SFRD for Pop III stars
            sfrdCP3  : the SFRD for Pop III stars that does not consider subgrid
                        mixing

        """
        self.file = file
        if self.verbose:
            print("Loading new-style stardata file")
            print("Loading file @ {}".format(self.prefix + self.file))
        self.sfrFile = np.loadtxt(self.prefix + self.file, skiprows=1)

        # Determine delta_t between each entry (each entry is for a single time/redshift)
        # Use the latest time as the time associated with the SFRD up to that time
        # Compute delta_mass
        # Compute average SFRD over the preceeding interval: delta_mass/delta_t * nmlze
        # where nmlze is the comoving volume of the simulation

        # z time totSPM totPIIIM totPIIM totCPIIIM totZPM
        # 0   1    2     3         4        5       6
        self.time = self.sfrFile[:, 1]  # Get time in Myr
        self.delta_t = np.diff(self.time) * 1e6  # Time between entries, yrs
        # End time point. Still calling it "mid_time"
        # Ignore the first dt since it isn't correct
        self.mid_time = self.time[1:] # In Myr
        if self.verbose:
            print("Times ", self.mid_time)
        # How much solar mass was created over the previous interval?
        self.delta_m = np.diff(self.sfrFile[:, 2])
        # If we have little or no star formation, but SN, dM could be negative!!
        self.delta_m[self.delta_m < 0.0] = np.nan
        self.sfrd = self.delta_m / self.delta_t / self.cm_vol # M_sun/yr/Mpc^3
        if self.verbose:
            print("SFRD values ", self.sfrd)
        # self.sfrd[np.isnan(self.sfrd)] = 1e-10     # No star formation during delta t, set to small #

        # Compute the redshift associated with the time...
        self.midz = []
        for t in self.mid_time:
            # Using astropy is ok, although it doesn't precisely match the orig file's redshifts
            self.midz.append(astropy.cosmology.z_at_value(
                self.cosmo.age, t * u.Myr))

        # Associate the delta sfrd with the later z
        self.redshift = self.sfrFile[1:, 0]

        self.delta_mP3 = np.diff(self.sfrFile[:, 3])
        self.delta_mP3[self.delta_mP3 < 0.0] = np.nan
        self.sfrdP3 = self.delta_mP3 / self.delta_t / self.cm_vol
        # self.sfrdP3[np.isnan(self.sfrdP3)] = 1e-10
        self.delta_mCP3 = np.diff(self.sfrFile[:, 5])
        self.delta_mCP3[self.delta_mCP3 < 0.0] = np.nan
        self.sfrdCP3 = self.delta_mCP3 / self.delta_t / self.cm_vol
        # self.sfrdCP3[np.isnan(self.sfrdCP3)] = 1e-10

    # Experimental ... use fitting to interpolate the values and compute the
    # SFRD at arbitrary redshifts...
    # def set_baseline(self, base):
    # self.baseline={}
    ##     self.baseline['tot']   = interp1d(base.zData,base.totStarMass,kind='linear')
    ##     self.baseline['pop3']  = interp1d(base.zData,base.popIIIMass, kind='linear')
    ##     self.baseline['cpop3'] = interp1d(base.zData,base.subCritMass,kind='linear')

    # def get_diff_funcs(self,newVal,base=None):
    # if not 'self.baseline' in locals and base==None:
    ##         print("Specify baseline for differencing... ")
    # return
    # if not 'self.baseline' in locals:
    # self.initBaseline(base)

    # self.other={}
    ##     self.other['tot']   = interp1d(newVal.zData,newVal.totStarMass,kind='linear')
    ##     self.other['pop3']  = interp1d(newVal.zData,newVal.popIIIMass, kind='linear')
    ##     self.other['cpop3'] = interp1d(newVal.zData,newVal.subCritMass,kind='linear')
    ##     zmax = min(newVal.zData.max(),self.baseline.zData.max())
    ##     zmin = max(newVal.zData.min(),self.baseline.zData.min())
    ##     theRange = np.linspace(zmin,zmax,10)
    # return theRange,interp1d(theRange, -self.baseline['tot'](theRange) + self.other['tot'](theRange),kind='linear'), \
    ##             interp1d(theRange, -self.baseline['pop3'](theRange) + self.other['pop3'](theRange),kind='linear'), \
    ##             interp1d(theRange, -self.baseline['cpop3'](theRange) + self.other['cpop3'](theRange), kind='linear')

    def loadOld(self, file):
        """
        Loads data expecting the following columns:
            z	tStart	 tEnd	 totStarMass	 totPop3StarMass	 
            totPollStarMass	 totPrimordStarMass	 totGasMass	 \
            totPristGasMass	 totSubcritStarMass	 totNonPrimordStarMass

        User specified the filename to read.

        Once the data is loaded, the user can access the fields below:
            redshift : an array of redshifts at which the SFRD value is computed
            sfrdP3   : the SFRD for Pop III stars
            sfrdCP3  : the SFRD for Pop III stars that does not consider subgrid
                        mixing

        """
        self.file = file
        if self.verbose:
            print("Loading old-style stardata file")
            print("Loading file @ {}".format(self.prefix + self.file))
        self.rawSFRdata = np.loadtxt(self.prefix + self.file, skiprows=1)

        # Determine delta_t between each entry (each entry is for a single time/redshift)
        # Use the time col as the time associated with the SFRD up to that time
        # Compute delta_mass
        # Compute average SFRD over the preceeding interval: delta_mass/delta_t * nmlze
        # where nmlze is the comoving volume of the simulation

        self.time = self.rawSFRdata[:, 2]  # Get all end times, units are years
        self.delta_t = np.diff(self.time)  # Time between entries
        # End time point. Still calling it "mid_time"
        self.mid_time = self.time[1:] / 1e6 # Times in Myr
        if self.verbose:
            print("Times ", self.mid_time)
        # How much solar mass was created over the previous interval?
        self.delta_m = np.diff(self.rawSFRdata[:, 3])
        # If we have little or no star formation, but SN, dM could be negative!!
        self.delta_m[self.delta_m < 0.0] = np.nan
        self.sfrd = self.delta_m / (self.delta_t) / self.cm_vol # M_sun/yr/Mpc^3
        if self.verbose:
            print("SFRD values ", self.sfrd)
        # self.sfrd[np.isnan(self.sfrd)] = 1e-10     # No star formation during delta t, set to small #

        # Compute the redshift associated with the time...
        self.midz = []
        for t in self.mid_time:
            # Using astropy is ok, although it doesn't precisely match the orig file's redshifts
            self.midz.append(
                astropy.cosmology.z_at_value(self.cosmo.age, t * u.Myr)
            )  # New time is Myr

        # Associate the delta sfrd with the later z
        self.redshift = self.rawSFRdata[1:, 0]

        self.delta_mP3 = np.diff(self.rawSFRdata[:, 4])
        self.delta_mP3[self.delta_mP3 < 0.0] = np.nan
        self.sfrdP3 = self.delta_mP3 / (self.delta_t) / self.cm_vol
        # self.sfrdP3[np.isnan(self.sfrdP3)] = 1e-10
        self.delta_mCP3 = np.diff(self.rawSFRdata[:, 9])
        self.delta_mCP3[self.delta_mCP3 < 0.0] = np.nan
        self.sfrdCP3 = self.delta_mCP3 / (self.delta_t) / self.cm_vol
        # self.sfrdCP3[np.isnan(self.sfrdCP3)] = 1e-10
