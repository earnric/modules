import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os.path
import time

import subprocess
import re

from ..datainterface import DataInterface

class BasePlot():
    def __init__(self,
                 data = None,
                 fig = None,
                 ax = None,
                 scale = 1,
                 style = 'kepler',
                 **kwargs):

        assert isinstance(data, DataInterface), 'need to provide data source interface DataInterface'

        self.data = data

        np.seterr(under = 'ignore')
        try:
            plt.style.use(os.path.join(os.path.dirname(__file__), f'{style}.mplstyle'))
        except OSError:
            plt.style.use(style)

        self.scale = scale

        if ax is not None:
            self.ax = ax
            self.fig = ax.figure
        else:
            if fig is not None:
                self.fig = figure
            else:
                self.set_figure()
            self.ax = self.fig.add_subplot(111)

        self.time = 0

    def set_figure(self):
        self.fig = plt.figure(
            figsize = (8,6),
            )
        title = r'KEPLER - {} - {}'.format(self.data.filename,
                                           self.data.runpath)
        self.fig.canvas.set_window_title(title)

        # this part needs to be rewritten based on nucplot.py code
        # if xrandr is not available
        dpi0 = self.fig.get_dpi()
        dpi = 100
        if mpl.get_backend() == 'TkAgg':
            output = subprocess.run(['xrandr'], stdout=subprocess.PIPE).stdout.decode()
            pixels = np.array(re.findall('(\d+)x(\d+) ', output)[0], dtype = np.int)
            mm = np.array(re.findall('(\d+)mm x (\d+)mm', output)[0], dtype = np.int)
            dpi = 25.4 * np.sqrt(0.5 * np.sum(pixels**2) / np.sum(mm**2))
        if dpi > dpi0:
            print(' [Plot] setting dpi to {}'.format(dpi))
            self.fig.set_dpi(dpi * self.scale)

    def update(self, interactive = True):
        # if (not interactive and
        #     time.time() - self.time <  self.data.parm.ipdtmin):
        #     return
        self.time = time.time()
        self.draw()
        self.show()

    def show(self):
        # plt.draw()
        self.fig.canvas.manager.show()
        plt.pause(0.001) # find plot_idle()
        #? self.fig.show(block = False)
        #? delf.fig.canvas.draw()
        #time.sleep(1e-6)

    def draw(self):
        raise NotImplementedError()

    def close(self):
        self.fig.canvas.manager.destroy()

    def clear(self):
        self.ax.lines.clear()
        self.ax.texts.clear()
        self.ax.patches.clear()

        try:
            self.ax2.lines.clear()
            self.ax2.texts.clear()
            self.ax2.patches.clear()
        except:
            pass

        try:
            self.ax.legend_.remove()
        except:
            pass
        try:
            self.ax2.legend_.remove()
        except:
            pass

    def save(self, filename):
        print(self.fig.get_facecolor())
        self.fig.savefig(
            filename,
            facecolor = self.fig.get_facecolor(),
            transparent=True,
            )
