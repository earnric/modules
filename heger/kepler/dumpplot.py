import kepdump
import time
import os.path
import numpy as np
from . import plot as kp
from .datainterface import DataInterface

class DumpData(DataInterface):
    def __init__(self, dump = None):
        if not isinstance(dump, kepdump.KepDump):
            dump = kepdump.load(dump)
        self.dump = dump
    @property
    def filename(self):
        return os.path.basename(self.dump.filename)
    @property
    def runpath(self):
        return os.path.dirname(os.path.abspath(self.dump.filename))
    @property
    def datatime(self):
        return time.asctime(time.strptime(self.dump.lastrun,'%H:%M:%S %d%b%y'))
    def __getattr__(self, var):
        o = object()
        val = getattr(self.dump, var, o)
        if val is not o:
            return val
        raise AttributeError(var)
    @property
    def idtcsym(self):
        sym = [x.strip() for x in self.dump.idtcsym]
        return [''] + sym
    @property
    def angd(self):
        return np.zeros((self.qparm.jm+2, 5))
    @property
    def entropies(self):
        """
        return array of entropes

        stored data only contains total entropy
        """
        data = np.zeros((self.qparm.jm+2, 6), dtype = np.float64)
        data[:,0] = self.stot
        return data


def plot(dump = None, plot = 1, **kwargs):
    """
    make kepler plot from dump
    """
    if not isinstance(dump, kepdump.KepDump):
        dump = kepdump.load(dump)
    dd = DumpData(dump)

    # this is to be adjusted for plot types
    p = kp.kepplots[plot](dd, **kwargs)
    return p
