# Load ramses unit information
# R. Sarmento - 2 Feb 2021
"""
Loads unit information from a given snapshot
returning conversion factors for 
    density
    length
    time
    mass
"""

import numpy as np
import os.path
import sys


class RamUnits:
    """
    Loads ramses unit information from the info_xxxxx.txt file
    and makes the information available to users.
    """

    def __init__(self, num, prefix="./"):
        """
        Loads the ramses info file and stores information in local
        vars
        This ctor updates the conversion factors for use in astro
        applications. 
        """
        self.gPerMsun = 1.988e33   # g/M_sun
        self.cmPerKpc = 3.0857e21  # cm/kpc
        fname = prefix+"output_{:05d}/info_{:05d}.txt".format(num, num)
        if os.path.isfile(fname):
            with open(fname) as f:
                lines = f.readlines()
                # print(lines[15:18])
                # Unit conversion factors go from ramses
                # -- len units to cm
                # -- density units to g/cc
                # -- time units to sec
                self.lenUnit = float(lines[15].split(
                    '=')[1]) / self.cmPerKpc  # cm to kpc
                self.rhoUnit = float(lines[16].split(
                    '=')[1]) / self.gPerMsun * self.cmPerKpc**3  # g/cc to M_sun/kpc^3
                self.timeUnit = float(lines[17].split('=')[1])  # to sec
                self.massUnit = self.rhoUnit * self.lenUnit**3  # in M_sun
        else:
            print("File {} not found".format(fname))
            sys.exit(1)
        return

    def getMass(self, unit="Msun"):
        """
        Returns the conversion factor from RAMSES
        mass units (e.g. - hop Ramses) to solar masses
        Multiply the value returned by ramses by this 
        number to get solar masses
        """
        if unit == "Msun":
            return self.massUnit
        elif unit == "g":
            return self.massUnit * self.gPerMsun
        else:
            print("RamUnits:getLen --UNKNOWN unit")
            return 0

    def getLen(self, unit="kpc"):
        """
        Returns the conversion factor from RAMSES
        length units to kpc, by default. Centimeters
        are also supported (use unit="cm").
        Multiply the values returned by ramses by this 
        number to get kpc or cm.
        """
        if unit == "kpc":
            return self.lenUnit
        elif unit == "cm":
            return self.lenUnit * self.cmPerKpc
        else:
            print("RamUnits:getLen --UNKNOWN unit")
            return 0

    def getTime(self):
        """
        Returns the conversion factor from RAMSES
        time units to sec.
        Multiply the values returned by ramses by this 
        number to get sec.
        """
        return self.timeUnit

    def getDen(self, unit="kpc"):
        """
        Returns the conversion factor from RAMSES
        density units to M_sun/kpc^3. Use unit="g"
        to return g/cc.
        Multiply the values returned by ramses by this 
        number to get either M_sun/kpc^3 or g/cc.
        """
        if unit == "kpc":
            return self.rhoUnit
        elif unit == "g":
            return self.rhoUnit * self.gPerMsun / self.cmPerKpc**3
        else:
            print("RamUnits:getLen --UNKNOWN unit")
            return 0
