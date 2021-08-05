"""
Python module to make KEPLER default plots.

(under construction)
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from logged import Logged
import kepdump

# make kepplot base class providing common stuff, like
#
#  - mass coordinates
#
#  - labeling of lines
#
#  - add caption header?
#
# overall design should be like kepler mongo plot, i.e., have a
# container plot with the info on top and having several subplots


# rgbcolors =((1.0D0,1.0D0,1.0D0),
#      1     0.0D0,0.0D0,0.0D0,   0.5D0,0.5D0,0.5D0,
#      2     0.5D0,0.0D0,0.0D0,   1.0D0,0.0D0,0.0D0,
#      3     1.0D0,0.0D0,0.5D0,   0.5D0,0.0D0,0.5D0,
#      4     1.0D0,0.0D0,1.0D0,   1.0D0,0.5D0,1.0D0,
#      5     1.0D0,0.5D0,0.5D0,   1.0D0,0.5D0,0.0D0,
#      6     1.0D0,1.0D0,0.0D0,   1.0D0,1.0D0,0.5D0,
#      7     0.5D0,0.5D0,0.0D0,   0.5D0,1.0D0,0.0D0,
#      8     0.0D0,0.5D0,0.0D0,   0.0D0,1.0D0,0.0D0,
#      9     0.5D0,1.0D0,0.5D0,   0.0D0,1.0D0,0.5D0,
#      o     0.0D0,0.5D0,0.5D0,   0.0D0,1.0D0,1.0D0,
#      1     0.5D0,1.0D0,1.0D0,   0.0D0,0.5D0,1.0D0,
#      2     0.0D0,0.0D0,0.5D0,   0.0D0,0.0D0,1.0D0,
#      3     0.5D0,0.5D0,1.0D0,   0.5D0,0.0D0,1.0D0/
#
#       data ioncolor / 5,22,19,6,4,11,24,16,7,10,
#      1                20,26,14,13,21,2,18,25,9/
#       data isocolor / 22,4,11,16,7,10,20,26,1,18,
#      1                2,24,9,14,13,5,19,25,6,21,
#      2                3,8,15,17,23,
#      3                22,4,11,16,7,10,20,26,1,18,
#      4                2,24,9,14,13,5,19,25,6,21,
#      5                3,8,15,17,23 /
#
#       data ielemcol /
#      &     22,4,19,16,7,10,20,26,18,
#      &     12,24,9,14,13,5,11,25,6,21,
#      &     3,8,15,17,23 /
#
#       data mixcol /16,11,10,9,20,7,5/
#
#       data magcol /13,3/



class compplot(Logged):
    def __init__(self, dump):
        if isinstance(dump, str):
            dump = kepdump.load(dump)
        assert isinstance(dump, kepdump.KepDump)
        self.dump = dump

    def plot(self):
        self.fig = plt.figure(
            figsize = (8,6),
            dpi = 102,
            facecolor = 'white',
            edgecolor = 'white',
            )
        self.ax = self.fig.add_subplot(111)

        d = self.dump
        xrange = d.zm_sun[[0,-2]]
        yrange = np.array([1.e-3, 1])

        ylabel = 'mass fraction'
        xlabel = 'enclosed mass / solar masses'

        xscale = 'linear'
        yscale = 'log'
        self.ax.set_xscale(xscale)
        self.ax.set_yscale(yscale)
        self.ax.set_xlim(*xrange)
        self.ax.set_ylim(*yrange)
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

        abu = d.abu
        m = d.zm_sun
        for n, r in abu.netnum_ranges:
            for i in abu.ions_netnum(n):
                self.ax.plot(
                    m[r[0]:r[1]+1],
                    abu.ion_abu(i)[r[0]:r[1]+1],
                    label = i.mpl,
                    )
        # TODO
        #
        #  - filter duplicates from legend; actually, add only select values to legend
        #
        #  - option top plot thisn on lines as in KEPLER
        #
        #  - get consistent line colors for each ion, take from KEPLER
        #    (similar should be done for Ion/Element
        #
        #  - allow reverse video, similar to KEPLER

        self.ax.legend(loc = 'best', ncol = 2, fontsize = 8)
        self.fig.tight_layout()
        plt.draw()
