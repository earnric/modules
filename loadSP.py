"""
Routines for loading star particle characteristics and halo
data.
"""
import numpy as np
from os import path
import gc

# Module constant
cmInPc = 3.086e18


def loadSP(z, prefix="./", normedLocs=False, gz=False):
    """loadSP - loads the locs, mass, birthTimes, Z, PZ, PPF files from the local dir
    This routine assumes the sp*.txt files are located in the current dir. However,
    you can override that by specifying the path as the second parameter.
    It will also handle gziped files: just override gz=False (default) in the
    parameter list.
    If the version of the sp locations files contained normalized locations
    in the interval (0,1), use normedLocs=True.
    Returns locs, mass, birthTimes, Z, PZ, PPF as in the file:
    If normedLocs: Locations are set to  (-0.5,0.5)
    else:  Corrected locations in kpc, physical centered on 0,0
    Solar masses
    Birth Times are in Myr
    Z SOLAR units in the file...
    PZ SOLAR units in the file...
    PPF Pure fraction...
    """
    if gz:
        gzStr = ".gz"
    else:
        gzStr = ""

    redString = "{:05.2f}".format(z)
    if not path.isfile(prefix + "spLoc_" + redString.format(z) + "_Normed.txt" + ".gz"):
        redString = "{:.2f}".format(z)

    if normedLocs:
        locs = np.loadtxt(
            prefix + "spLoc_" + redString + "_Normed.txt" + gzStr,
            skiprows=1,
        )  # Locations are (0,1)
        locs -= [0.5, 0.5, 0.5]  # Locations are now (-0.5,0.5)
    else:
        locs = np.loadtxt(
            prefix + "spLoc_" + redString + ".txt" + gzStr, skiprows=1
        )  # Corrected locations in kpc, physical
    mass = np.loadtxt(
        prefix + "spMass_" + redString + ".txt" + gzStr, skiprows=1
    )  # Solar masses
    birthTimes = np.loadtxt(
        prefix + "spBT_" + redString + ".txt" + gzStr, skiprows=1
    )  # BirthTimes are in Myr
    Z = np.loadtxt(
        prefix + "spZ_" + redString + ".txt" + gzStr, skiprows=1
    )  # SOLAR units in the file...
    PZ = np.loadtxt(
        prefix + "spPZ_" + redString + ".txt" + gzStr, skiprows=1
    )  # SOLAR units in the file...
    PPF = np.loadtxt(
        prefix + "spPPF_" + redString + ".txt" + gzStr, skiprows=1
    )  # Pure fraction...
    return locs, mass, birthTimes, Z, PZ, PPF


#
# Routines below deal with HOP halo files...
#
def loadHaloGrps(num, prefix="./hop/"):
    """
    loadHaloGrps - loads the halo locations from the group file.

    The routine assumes halos are located in dir './hop/'. This can be
    overridden using prefix.
    num is the outputfile number associated with the group file.
    The data in the HOP file is in the interval (0,1), so we normalize
    to (-0.5, 0.5)
    The information in the file is in the following form:
    #   npart  mass  cont.frac   xc    yc    zc    uc    vc    wc
    Mass units are in terms of the RAMSES mass unit: unit_l^3 * unit_d
    so mass * unit_l^3 * unit_d will give you the actual mass of the halo
    Note however, that we can also use the DM particle mass, easy to compute
    from rho_crit, the sim size and the initial level (levelmin). Hence I
    return the number of particles in the halo.
    Return
       Position in x,y,z in the range -.5,.5
       Number of DM particles in the halo
    """
    # Key to pos file is #   npart,mass,cont.frac,xc,yc,zc,uc,vc,wc
    grpfile = prefix + "grp%05d.pos" % num
    halosRaw = np.loadtxt(grpfile, skiprows=1)
    # Recenter so range is (-.5,.5) vs (0,1) Center is now (0,0)
    halosRaw[:, 4:7] -= 0.5
    # Return position information, normalized to (-0.5,0.5) along with npart
    return (halosRaw[:, 4:7], halosRaw[:, 1])


def loadHaloSizes(num, prefix="./hop/"):
    """loadHaloGrps - loads the halo sizes. This is number of particles
    The routine assumes halos are located in './hop/'. This can be
    overridden by the second parameter
    """
    # Key to pos file is #   npart,mass,cont.frac,xc,yc,zc,uc,vc,wc
    grpfile = prefix + "grp%05d.size" % num
    sizeList = np.loadtxt(
        grpfile,
        dtype={"names": ("index", "count"), "formats": ("i4", "i8")},
        skiprows=3,
    )
    return sizeList  # Return sizes..


#
# Routines below deal with my halo files...
#
def loadHaloXtra(z, prefix="./"):
    """loadHaloXtra - loads the halo locations from the txt file.
    This file only contains a list of halo positions.
    The routine assumes halos are located in './'. This can be
    overridden by the second parameter.
    The interval of this data is (-0.5, 0.5) ...
    """
    # Key to pos file is #   npart,mass,cont.frac,xc,yc,zc,uc,vc,wc
    grpfile = prefix + "extraHalos_{:.1f}.txt".format(z)
    halosRawPos = np.loadtxt(grpfile, skiprows=1)
    return halosRawPos  # Just return position information, normalized to (-0.5,0.5)


def loadAllHalos(z, prefix="./"):
    """loadAllHalos - loads the halo locations from the txt file.
    This file contains a list of halo positions, stellar masses and radii
    The data is expected to be in the form:
    x       y       z       mass      radius
    normalized to the range [-.5,.5] ... with Mass in solar units and
    radius also normalized [0,1]
    The routine assumes halos are located in './'. This can be
    overridden by the second parameter.
    Return positions, masses and radius of halo
    """
    # Key to pos file is x,y,z, mass, rad ... all normalized to -0.5->0.5 for position
    # and 0-1 for radius.
    halofile = prefix + "All_halos_{:.1f}.txt".format(z)
    halosRaw = np.loadtxt(halofile, skiprows=1)
    # Returns position information, normalized to (-0.5,0.5)
    return halosRaw[:, 0:3], halosRaw[:, 3], halosRaw[:, 4]


def loadHMHalos(num, prefix="./"):
    """loadHMHalos - loads the halo locations, masses and r_v from the txt file.
    This file contains a list of halo positions, stellar masses and radii
    The data is expected to be in the form:
    mass,   x,       y,       z,   radius
    which is generated by the extractMass_Rad_Pos.sh script
    Position is normalized to the range [-.5,.5]
    Radius is in range [0, 1]
    Mass is DM mass and converted to solar masses (file contains mass in M_sun/1e11 units)
    The routine assumes halos are located in './'. This can be
    overridden by the second parameter.
    Return positions, masses and radius of halo
    """
    halofile = prefix + "processed_Halos_{:03d}.txt".format(num)
    # The above processed_Halos file is created by extractMass_Rad_Pos.sh
    halosRaw = np.loadtxt(halofile, skiprows=1)
    # Returns position information, normalized to (-0.5,0.5)
    return halosRaw[:, 1:4] - (0.5, 0.5, 0.5), halosRaw[:, 0] * 1e11, halosRaw[:, 4]


def loadHMSubHalos(num, prefix="./"):
    """loadHMHalos - loads the halo locations, masses and r_v from the txt file.
    This file contains a list of halo positions, stellar masses and radii
    The data is expected to be in the form:
    mass,  parentHaloID, x,       y,       z,   radius
    which is generated by the extractMass_Rad_Pos.sh script
    Position is normalized to the range [-.5,.5]
    Radius is in range [0, 1]
    Mass is DM mass and converted to solar masses (file contains mass in M_sun/1e11 units)
    The routine assumes halos are located in './'. This can be
    overridden by the second parameter.
    Return positions, masses and radius of halo
    """
    halofile = prefix + "processed_SubHalos_{:03d}.txt".format(num)
    # The above processed_Halos file is created by extractMass_Rad_Pos.sh
    halosRaw = np.loadtxt(halofile, skiprows=1)
    # Returns position information, normalized to (-0.5,0.5)
    return halosRaw[:, 2:5] - (0.5, 0.5, 0.5), halosRaw[:, 0] * 1e11, halosRaw[:, 5]


def getBoxSize(num, prefix="./"):
    """getBoxSize - returns the size of the RAMSES output file boxsize in cm
    The routine assumes the output files are in the current dir. This can be
    overridden with the second parameter.
    """
    infile = prefix + "output_%05d/info_%05d.txt" % (num, num)
    f = open(infile)
    lines = f.readlines()
    boxsizecm = float(lines[15].split()[2])  # Get the number at the end of the line...
    return boxsizecm


def getAgeSize(filename="snapKeyInfo.txt"):
    """getBoxSize - returns the size of the RAMSES output file boxsize in cm
    The routine assumes the output files are in the current dir. This can be
    overridden with the second parameter.
    """
    theData = np.genfromtxt(filename, delimiter=",", names=True)
    return theData
