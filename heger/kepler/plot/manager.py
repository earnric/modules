import numpy as np

from . import kepplots
from .base import BasePlot

class PlotManager():
    def __init__(self, data = None, default = None):
        self.data = data
        self.default = default
        self.plots = []

    def plot(self, *args, **kwargs):
        kwargs.setdefault('duplicates', True)
        self.gc()
        for p in args:
            plot = self.add_plot(p, **kwargs)
        if len(self.plots) is 0:
            if self.default is not None:
                plot = self.add_plot(self.default)
            else:
                plot = self.add_plot(self.data.parm.ipixtype, *args, **kwargs)
        if len(args) == 0:
            self.update(**kwargs)

    def update(self, *args, **kwargs):
        self.gc()
        if len(args) == 0:
            args = range(len(self.plots))
        interactive = kwargs.pop('interactive', True)
        for i in args:
            p = self.plots[i]
            if isinstance(p, type):
                self.plots[i] = p(data, **kwargs)
            else:
                p.update(interactive = interactive)

    def gc(self):
        """
        garbage collect closed plots.
        """
        for p in self.plots.copy():
            if isinstance(p, type):
                continue
            if not p.fig.canvas.is_drawable():
                print(f' [kepler.plot] Removing closed plot {p}')
                self.plots.remove(p)

    def add_plot(self, *args, **kwargs):
        if len(args) > 0:
            plot, *args = args
        else:
            plot = kwargs.pop('plot', None)
        plot_ = plot
        duplicates = kwargs.pop('duplicates', False)
        self.gc()
        if plot is None and self.default is not None:
            plot = self.default
        if isinstance(plot, (int, np.int32)):
            plot = kepplots.get(plot, None)
            if plot is None:
                print(f' [add_plot] Plot type {plot_} not supported.')
        if plot is None:
            print(f' [add_plot] No plot added for {plot_}.')
            return
        if not duplicates:
            for p in self.plots:
                if isinstance(p, type):
                    if isinstance(plot, type):
                        if p == plot:
                            return
                    elif isinstance(plot, p):
                        return
                elif isinstance(plot, type):
                    if isinstance(p, plot):
                        return
                elif type(p) == type(plot):
                    return
        if isinstance(plot, type) and issubclass(plot, BasePlot):
            plot = plot(self.data, *args, **kwargs)
        assert isinstance(plot, BasePlot)
        self.plots.append(plot)
        return plot

    def closewin(self):
        for p in self.plots:
            p.close()
        self.plots = []
