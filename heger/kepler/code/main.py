"""
Kepler main routines

Example:

"""

# STDL import
import numpy as np
import os
import os.path
import sys
import time

from collections import OrderedDict
from math import floor, ceil, log10
from select import select

# alex lib imports (REPLACE)
import keppar
from logged import Logged

# kepler imports
from ..datainterface import DataInterface
from ..plot import kepplots
from ..plot.manager import PlotManager
from ..plot.base import BasePlot
from ..exception import KeplerTerminated

# local dir imports
from .kepbin import _kepler
from .util import b2s, s2b, z2int, z2str
from .interface import Item, CharVar, Parameter

# define local variables
default_plot = []

parm = Parameter(keppar.p, 'p')
qparm = Parameter(keppar.q, 'q')

var = Parameter(OrderedDict([
    ('jm', 0),
    ('ndump', 0),
    ('nsdump', 0),
    ('nedit', 0),
    ('gencom.mode', 4),
    ('vars.jdtc', 0),
    ]))

info = dict()

class KeplerData(DataInterface):
    # TODO - cache values using cycle ID?
    @property
    def filename(self):
        return b2s(_kepler.namecom.nameprob)
    @property
    def nameprob(self):
        return b2s(_kepler.namecom.nameprob)
    @property
    def runpath(self):
        return info['runpath']
    @property
    def qparm(self):
        return qparm
    @property
    def parm(self):
        return parm
    def __getattr__(self, var):
        o = object()
        val = getattr(_kepler.vars, var, o)
        if val is not o:
            return val
        # try loadbuf
        jm = qparm.jm
        buf, label, ierr = _kepler.loadbuf_(var, 1, jm)
        if ierr == 0:
            val = np.zeros(jm+2)
            val[1:-1] = buf[:jm]
            return val
        val = getattr(_kepler.charsave, var, o)
        if val is not o:
            return val
        raise AttributeError(var)
    @property
    def iconv(self):
        jm = qparm.jm
        icon, label, ierr = _kepler.loadbuf_('convect', 1, jm)
        return np.int_(np.round(icon))
    @property
    def angit(self):
        jm = qparm.jm
        it = np.zeros(jm+2)
        it[1:-1] = np.cumsum(self.angi[1:jm+1]*self.xm[1:jm+1])
        return it
    @property
    def jdtc(self):
        return var.jdtc
    @property
    def datatime(self):
        return time.asctime()
    @property
    def idtcsym(self):
        sym = [x.tostring().decode().strip() for x in _kepler.charsave.idtcsym]
        return sym
    @staticmethod # copy from kepdump --> should be in keputils
    def center2face(x, logarithmic=False):
        """
        create zone interface quantity by averaging zone values
        """
        y = x.copy()
        if logarithmic:
            y[:-1] *= x[1:]
            y[1:-2] = np.sqrt(y[1:-2])
        else:
            y[:-1] += x[1:]
            y      *= 0.5
        y[0]    = x[1]
        y[-2]   = x[-2]
        y[-1]   = np.nan
        return y
    @property
    def dnf(self):
        """
        Density at zone interface (flat extrapolation for boundaries) (g/cm**3).
        """
        return self.center2face(self.dn)
    @property
    def entropies(self):
        """
        return array of entropes from call to EOS

        This is done not just using loadbuf as this would reqire
        manymore calls to the equation of state.
        """
        return _kepler.getentropies_(1,qparm.jm)
    # BURN interfce
    @property
    def numib(self):
        return z2int(self.znumib)
    @property
    def netnumb(self):
        """
        here we have to define a function to overwrite that a float
        value would be returned by loadbuf
        """
        return z2int(self.znetnumb)
    @property
    def ionnb(self):
        return z2int(self.zionnb)
    @property
    def ionsb(self):
        return np.array([i.decode('us_ascii') for i in _kepler.charsave.ionsb])
    @property
    def isosym(self):
        return np.array([i.decode('us_ascii') for i in _kepler.charsave.isosym])
    @property
    def ppnb(self):
        x = np.ndarray((qparm.jm, qparm.imaxb), buffer=_kepler.vars.ppnb).transpose()
        x = np.insert(x, [0,qparm.jm],0,axis=1)
        return x
    # APPROX/QSE/NSE interface
    @property
    def ionn(self):
        return z2int(self.zionn)
    @property
    def numi(self):
        return z2int(self.znumi)
    @property
    def ppn(self):
        x = np.ndarray((qparm.jm, qparm.imax), buffer=_kepler.vars.ppn).transpose()
        x = np.insert(x, [0,qparm.jm],0,axis=1)
        return x

kd = KeplerData()
kp = PlotManager(kd)

class Kepler():
    def __init__(self, *args, **kwargs):
        self.start(*args, **kwargs)

    def start(self,
              run = 'xxx',
              gen = 'xxxz',
              parms = '',
              progname = 'python',
              path = '~/kepler/test',
              killburn = False,
              plot = True,
              # here a flag telling whether it runs as server process
              ):

        # set up current path
        self.cwd = os.getcwd()
        if path is None:
            self.path = self.cwd
        else:
            self.path = os.path.expanduser(os.path.expandvars(path))
            os.chdir(self.path)
        self.info['runpath'] = self.path

        commands = np.zeros(21, dtype="S80")

        args = [progname,
                run,
                gen,
                ] + list(parms.split())
        if killburn:
            args += ['k']
        for i, cmd in enumerate(args):
            commands[i] = cmd + ' '*(80-len(cmd))

        S = np.ndarray(80*21, dtype='S1', buffer = commands.data)
        _kepler.start_(len(args)-1, S)
        if plot:
            self.plot()
        self.status = 'started'

    def execute(self, cmdline, remote = False):
        """
        Execute Kepler command as on the terminal.  Pass string
        """
        _kepler.execute_(cmdline + ' '*(132-len(cmdline)), 1)

    def cycle(self, n = 1, remote = False):
        """
        cycle for n steps.
        """
        if self.status is not "started":
            raise KeplerTerminated()
        def getch():
            i, o, e = select([sys.stdin], [], [], 0)
            if i:
                ch = sys.stdin.read(1)
            else:
                ch = ''
            return ch

        _kepler.gencom.mode = 1
        try:
            for i in range(n):
                _kepler.cycle_(False)
                if not remote:
                    x = getch()
                    if x == 's':
                        break
                    elif x != '':
                        x = str(x).strip()
                        self.execute(x)
                if i == n-1 or _kepler.gencom.mode == 0:
                    break
            if len(self.kp.plots) > 0 and self.parm.itvstart < 1:
                # copy from cycle.f
                iseconds = int(time.time())
                if (self.qparm.ilastpl > iseconds):
                    self.qparm.ilastpl = 0
                if ((((self.qparm.ncyc % self.parm.npixedit) == 0) and
                    (iseconds - self.qparm.ilastpl >= self.parm.ipdtmin)) or
                    (i == n - 1)):
                    self.plot(interactive = False)
                    self.qparm.ilastpl = iseconds
        except KeplerTerminated:
            self.status = 'terminated'
            self.closewin()
            raise

        ncyc = qparm.ncyc
        print(' [CYCLE] Stop at cycle {:d}'.format(ncyc))

    #@property # broken autocompletion
    def s(self, n = 1):
        """
        Do a single step.  No need to call as function.
        """
        self.cycle(n)

    #@property # broken autocompletion
    def g(self):
        """
        Do run continuously ('go').
        """
        self.cycle(sys.maxsize)

    def terminate(self, message = '', remote = False):
        self.closeplot()
        os.chdir(self.cwd)
        _kepler.terminate_(message)

    def plot(self, *args, **kwargs):
        self.kp.plot(*args, **kwargs)

    def closewin(self):
        self.kp.closewin()

    def pt(self, mode):
        p = kepplots.get(mode, None)
        if p is not None:
            return p(kd)

    #@property # broken autocompletion
    def run(self):
        """
        Interactive KEPLER mode

        TODO replace with something that allows plots to continue work.
        App?
        """

        KeplerCmdLoop(self, kd).cmdloop()

    @property
    def parm(self):
        return parm

    @property
    def qparm(self):
        return qparm

    @property
    def kd(self):
        return kd

    @property
    def kp(self):
        return kp

    @property
    def info(self):
        return info

    def __getattr__(self, key):
        try:
            return parm[key]
        except:
            pass
        try:
            return qparm[key]
        except:
            pass
        o = object()
        val = getattr(_kepler.vars, key, o)
        if val is not o:
            return val
        # try loadbuf
        jm = self.qparm.jm
        buf, label, ierr = _kepler.loadbuf_(key, 1, jm)
        if ierr == 0:
            val = np.zeros(jm+2)
            val[1:-1] = buf[:jm]
            return val
        raise AttributeError()

    def __setattr__(self, key, value):
        if hasattr(parm, key):
            parm[key] = value
            return
        if hasattr(qparm, key):
            raise AttributeError("read only variable '{}'".format(key))
            qparm[key] = value
            return
        super().__setattr__(key, value)

# def end_kepler():
#     raise Exception('Kepler ended.')

# define and set callbacks
def endkepler(code = None):
    """
    Callback to end KEPLER.
    """
    print(' [endkepler] has been called.')
    raise KeplerTerminated(code)

def plotkepler():
    """
    Here we may need to keep track of kepler plots that were generated
    or need to be updated.
    """
    print(' [kepler] called python plot', parm.ipixtype)
    if parm.itvstart == 0:
        return
    if parm.itvstart == -1:
        parm.itvstart == 0
        kp.closewin()
    if parm.ipixtype == 0:
        return
    p = kepplots.get(parm.ipixtype, None)
    kp.plot(p, interactive = False, duplicates = False)

def getskepler(msg):
    prompt = b2s(_kepler.namecom.nameprob).strip('\x00').strip() + '> '
    s = input(prompt)
    msg[:] = 0
    for i, c in enumerate(s):
        msg[i] = bytes(c, encoding = 'us-ascii')[0]

_kepler.endkepler = endkepler
_kepler.plotkepler = plotkepler
_kepler.getskepler = getskepler
