"""Routines for plotting and retrieving halo stars
    pltHalo(locs,halosPos,indx,z,size=10)
    getHaloStars(locs, mass, Z, PZ, PPF, halosPos,indx,s=10)
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import math as ma
import numpy as np


def pltView(locs, haloPos, size=10, indx=0, ax=None, Z=None, s=40):
    """pltView - creates an axis object with a plot of the stars visible at a specified
    location. Note that since we are looking at a 2D plane, some of the stars in the field
    may be much nearer or further than the centroid of the halo!!
    Parameters are star locations, halo location, plot size in kpc.
    ** Note that locs and halosPos should be in the same coordinate system! **
    Note that this routine plots the star partciles at their physical location at redshift
    z.
    """
    if not ax:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.gca()
    else:
        fig = plt.gcf()
    # Plot the center location
    # ax.scatter(
    #     haloPos[0], haloPos[1], s=120, facecolors="none", edgecolors="r"
    # )  # Plot halo locs
    if Z is None:
        ax.scatter(
            locs[:, 0],
            locs[:, 1],
            s=s,
            c="b",
            alpha=0.75,
        )  # Plot star particles
    else:
        norm = mpl.colors.LogNorm(vmin=1e-5, vmax=1.0, clip=True)
        mapper = cm.ScalarMappable(norm=norm, cmap=cm.brg)
        ax.scatter(
            locs[:, 0],
            locs[:, 1],
            s=s,
            c=mapper.to_rgba(Z),
            cmap=cm.brg,
            alpha=0.75,
        )  # Plot star particles
        axins1 = inset_axes(
            ax,
            width="50%",  # width = 50% of parent_bbox width
            height="5%",  # height : 5%
            loc="upper right",
        )
        cbar = fig.colorbar(mapper, cax=axins1, orientation="horizontal")
        cbar.ax.set_xlabel("$Z \; [Z_\odot]$", rotation=0, loc="center")
        # cbar.ax.set_ylabel("$Z \; [Z_\odot]$", loc="top")
        # cbar.ax.set_title("$Z \; [Z_\odot]$", fontsize=14)

    ax.grid(b=True, which="major", color="b", linestyle="--", alpha=0.5)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_xlim([haloPos[0] - size / 2, haloPos[0] + size / 2])
    ax.set_ylim([haloPos[1] - size / 2, haloPos[1] + size / 2])

    if indx:
        label = "halo{0}".format(indx)
        ax.annotate(
            label,
            xy=(haloPos[0], haloPos[1]),
            xytext=(-40, 40),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=12,
            bbox=dict(boxstyle="round,pad=0.05", fc="yellow", alpha=0.5),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0"),
        )

    return ax


def getViewStars(locs, mass, ages, Z, PZ, PPF, haloPos, fov=10):
    """getViewStars - Returns the stars that are within the view 2D (!!!) view range
    specified. Note that since we are looking at a 2D plane, some of the stars in the field
    may be much nearer or further than the centroid of the halo!!
    The centroid of the halo is translated to 0,0 ... and all stars within abs(x,y) < (fov/2,fov/2)
    are included in the results.
    Parameters are star particle locations, mass, ages, Z, PZ, PPF, halo location
    and FOV side-length (assumes square). Units are kpc.
    Return locs, mass, ages, Z, PZ, PPF of stars in radius s
    """

    # Find a radius that encompasses a cube with side len 's'
    xydist = fov / 2.0  # Compute how far to look in each direction from center of halo
    centrd = locs - haloPos  # (0,0) is now the center of the halo, z??
    centrd2D = centrd[::, [0, 1]]  # Just keep 2D data...

    cond = (np.abs(centrd2D[::, 0]) <= xydist) & (np.abs(centrd2D[::, 1]) <= xydist)
    haloStars = centrd[
        cond
    ]  # Note that we're getting 0,0 centered coords, but ignoring z!
    halomass = mass[cond]
    haloages = ages[cond]
    haloZ = Z[cond]
    haloPZ = PZ[cond]
    haloPPF = PPF[cond]

    return haloStars, halomass, haloages, haloZ, haloPZ, haloPPF


def getHaloStars_normed(locs, mass, ages, Z, PZ, PPF, haloPos, radius):
    """getHaloStars_normed - returns locations, masses, ages, Z, PZ, PPF of all stars within a specified
    distance of haloPos.
    Parameters are star particle locations, mass, ages, Z, PZ, PPF, halo location
    and radius. Units must be consistent! REMEMBER: r is a radius not a diamter.
    Returned coords are *** relative to haloPos ***. So the center of the view is (0,0)
    Returns locs, mass, ages, Z, PZ, PPF of stars in radius s
    """
    centrd = locs - haloPos  # Recenter the star locations on the halo center
    dists = np.linalg.norm(centrd, axis=1)
    haloStars = centrd[dists <= radius]
    halomass = mass[dists <= radius]
    haloages = ages[dists <= radius]
    haloZ = Z[dists <= radius]
    haloPZ = PZ[dists <= radius]
    haloPPF = PPF[dists <= radius]

    return haloStars, halomass, haloages, haloZ, haloPZ, haloPPF


def getHaloStars_origCoord(locs, mass, ages, Z, PZ, PPF, haloPos, radius):
    """getHaloStars_origCoord - returns locations, masses, ages, Z, PZ, PPF of all stars within a specified
    distance of haloPos.
    Parameters are star particle locations, mass, ages, Z, PZ, PPF, halo location
    and radius. Units must be consistent! REMEMBER: r is a radius not a diamter.
    Returned coords are in the original system.
    Returns locs, mass, ages, Z, PZ, PPF of stars in radius s
    """
    centrd = locs - haloPos  # Recenter the star locations on the halo center
    dists = np.linalg.norm(centrd, axis=1)
    haloStars = centrd[dists <= radius]
    halomass = mass[dists <= radius]
    haloages = ages[dists <= radius]
    haloZ = Z[dists <= radius]
    haloPZ = PZ[dists <= radius]
    haloPPF = PPF[dists <= radius]

    return haloStars + haloPos, halomass, haloages, haloZ, haloPZ, haloPPF


def getHaloStarIndxs(locs, haloPos, radius):
    """getHaloStarIndxs - returns indices into location array for all stars within
    radius.
    Parameters are star particle locations, halo location
    and radius. Units must be consistent! REMEMBER: r is a radius not a diamter.
    Returns List of indices into locs within radius
    """
    centrd = locs - haloPos  # Recenter the star locations on the halo center
    dists = np.linalg.norm(centrd, axis=1) # centrd has same size as locs
    haloStarIndxs = np.where(dists <= radius) # find indices of stars within radius

    return haloStarIndxs[0]

def FindHaloStars(
    locs, haloPos, radius, mass=None, ages=None, Z=None, PZ=None, PPF=None
):
    """FindHaloStars - returns locations all stars within a specified
    distance of haloPos.
    Parameters are star particle locations, halo location and radius.
    Units are must be consistent! REMEMBER: r is a radius not a diameter.
    Returns locs
    """
    rel_coords = locs - haloPos  # Recenter the star locations on the halo center
    dists = np.linalg.norm(rel_coords, axis=1)
    haloStars = rel_coords[dists <= radius]
    if mass:
        mass = mass[dists <= radius]
    if ages:
        ages = ages[dists <= radius]
    if Z:
        Z = Z[dists <= radius]
    if PZ:
        PZ = PZ[dists <= radius]
    if PPF:
        PPF = PPF[dists <= radius]

    # Return stars using orig coords
    return haloStars + haloPos, mass, ages, Z, PZ, PPF
