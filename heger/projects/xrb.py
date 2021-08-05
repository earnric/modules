"""
Utilities for XRB
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import kepdump
import ionmap
import isotope
import color
import bdat
import datetime
import subprocess
import re
import color

import utils

class FigTemplate(object):
    def get_fig(self, fig = None, ax = None, figsize = (6.4, 4.8), **kwargs):
        if ax is None:
            if fig is None:
                fig = plt.figure(figsize=figsize)

                if mpl.get_backend() == 'TkAgg':
                    try:
                        dpi0 = fig.get_dpi()
                        output = subprocess.run(['xrandr'], stdout=subprocess.PIPE).stdout.decode()
                        pixels = np.array(re.findall('(\d+)x(\d+) ', output)[0], dtype = np.int)
                        mm = np.array(re.findall('(\d+)mm x (\d+)mm', output)[0], dtype = np.int)
                        dpi = 25.4 * np.sqrt(0.5 * np.sum(pixels**2) / np.sum(mm**2))
                        print(f' [Plot] old dpi was {dpi0}')
                        print(f' [Plot] setting dpi to {dpi}')
                        fig.set_dpi(dpi)
                    except:
                        raise
                        print(f' [Plot] could not use xrandr to determine dpi')

            ax = fig.add_subplot(111)
        return fig, ax

    def yscale(self, dump,
               ybase = None,
               ytop = None,
               r = 'y',
               scaley = True,
               **kwargs):

        y = np.append(np.cumsum(dump.xm[-2::-1])[::-1], 0.)

        jbase = np.where(dump.zm > dump.parm.bmasslow)[0][0] + 1
        if ybase is None:
            ybase = y[jbase]
        else:
            jbase = np.where(y > ybase)[0][-1] + 1

        if r == 'radius':
            y = dump.rn[-2] - dump.rn
            y[-1] = y[-2]
            ybase = y[jbase]
            ylabel = r'depth / {}cm'
        elif r == 'logy':
            y = np.log10(np.maximum(dump.y, 1e-99))
            ytop = y[-3] - np.log10(2)
            ybase = y[jbase]
            ylabel = r'$\log($ column depth / $\mathrm{{g}}\,\mathrm{{cm}}^{{-2}}$ $)$'
            scaley = False
        elif r == 'y':
            y = dump.y
            ybase = y[jbase]
            ylabel = r'column depth / {}g$\,$cm$^{{-2}}$'
        else:
            ylabel = 'exterior mass coordinate / {}g'

        if ytop is None:
            ytop = 0

        if scaley:
            # master code for scale formatting
            yl = int(np.floor(np.log10(ybase)))
            scale = np.power(10., -yl)
            ybase *= scale
            ytop *= scale
            y = y * scale
            if yl == 0:
                ylabel = ylabel.format('')
            else:
                ylabel = ylabel.format(fr'$10^{{{yl:d}}}\,$')

        return y, ybase, ytop, ylabel, jbase


class IonMapPlot(FigTemplate):
    def __init__(self,
                 dump,
                 abunlim=1.e-10,
                 cmap = None,
                 vmax = None,
                 mode = 'massfrac',
                 **kwargs,
                 ):
        if isinstance(dump, str):
            dump = kepdump.load(dump)
        self.dump = dump
        if cmap is None:
            cmap = color.ColorBlendBspline(('white',)+tuple([color.ColorScale(color.colormap('viridis_r'), lambda x: (x-0.2)/0.7)]*2) + ('black',), frac=(0,.2,.9,1),k=1)

        b = dump.abub
        if mode == 'massfrac':
            a = ionmap.decay(b, decay = False, isobars = True, molfrac_out = False)
            data = np.log10(np.maximum(a.data, abunlim))
            vmax = 0
            vmin = np.log10(abunlim)
            label = 'mass fraction'
            unit = None

        elif mode in ('ME', 'MElog'):
            B = bdat.BDat()
            Ab = isotope.ufunc_A(b.ions)
            me = B.mass_excess(b.ions) - Ab * B.mass_excess('fe56') / 56
            label = 'mass excess relative to ${{}}^{{56}}$Fe'
            unit = 'MeV/nucleon'
            b.data = b.data * me
            a = ionmap.decay(b, decay = False, isobars = True, molfrac_out = True)
            if mode == 'ME':
                data = a.data
                vmin = 0
                cmap = color.ColorScaleGamma(cmap, -3)
            else:
                data = np.log10(np.maximum(a.data, abunlim))
                vmin = np.min(data)
                label = f'log( {label} )'
            if vmax is None:
                vmax = np.max(data)
            # cmap = color.colormap('viridis_r')
            # cmap = color.ColorScaleGamma(cmap, 1/3)

        if unit is not None:
            label = f'{label} ({unit})'

        A = isotope.ufunc_A(a.ions)

        amax = np.max(A)
        missing, = np.where(~np.in1d(np.arange(amax)+1, A))

        x = np.arange(amax) + 0.5
        y, ybase, ytop, ylabel, jbase = self.yscale(dump, **kwargs)
        data = np.insert(np.transpose(data)[:,:-1], missing, np.zeros(len(y)-1), axis = 0)

        fig, ax = self.get_fig(**kwargs)

        m = ax.pcolormesh(y, x, data, cmap = cmap, vmin=vmin, vmax = vmax)
        cb = fig.colorbar(m, label = label)
        ax.set_xlim(ybase, ytop)
        ax.set_ylim(0, max(x))
        ax.set_ylabel('mass number')
        ax.set_xlabel(ylabel)

        fig.tight_layout()

        self.abu = a
        self.fig = fig
        self.ax = ax


class MEPlot(object):
    def __init__(self,
                 dumps = (4160, 5000, 5340),
                 base = 'run32#{}',
                 labels = None,
                 fig = None,
                 ax = None,
                 mscale = 1e21,
                 norm = 'c12',
                 accmass = None,
                 ):
        """
        TODO - not just plot zone centres but zones as steps
        """

        if ax is None:
            if fig is None:
                fig = plt.figure()
            ax = fig.add_subplot(111)

        B = bdat.BDat()
        norm = isotope.ion(norm)
        menorm = B.mass_excess(norm) / norm.A
        mefe56 = B.mass_excess('fe56')/56 - menorm
        print(f' [MEPlot] Normaslising to {norm} (ME = {menorm} MeV/nucleon)')
        zm0 = None
        ep = None
        for i, dump in enumerate(dumps):
            dump = base.format(dump)
            print(f' [MEPlot] processing {dump}')
            d = kepdump.load(dump)
            b = d.abub
            me = B.mass_excess(b.ions)
            mex = np.sum(b.data * me, axis=1) - menorm
            jblo = np.searchsorted(d.zm, d.bmasslow) + 1
            jj = slice(jblo, -1)
            m = (d.zmm[jj] - d.zmm[jblo]) / mscale
            if labels is not None and len(labels) > i:
                label = labels[i]
            else:
                label = f'{dump} ({str(datetime.timedelta(seconds=int(d.time)))})'
            ex = mex[jj]
            ax.plot(m, ex, label = label)

            if accmass is not None:
                if zm0 is None:
                    zm1 = d.zm[-2]
                    zm0 = zm1 - accmass
                jhi = np.searchsorted(d.zm, zm1) + 1
                ii = slice(jblo, jhi)
                e = np.sum(d.xm[ii] * mex[ii])
                if ep is not None:
                    de = ep - e
                    print(f' [MEPlot] energy difference is {de/accmass} MeV/nucleon')
                ep = e

                jlo = np.searchsorted(d.zm, zm0) + 1
                ii = slice(jlo, jhi)
                e = np.sum(d.xm[ii] * mex[ii])
                print(f' [MEPlot] total is {e/accmass + menorm - mefe56} MeV/nucleon.')


        ax.axhline(mefe56, color = 'gray', ls = '--', label = r'${}^{56}\mathrm{Fe}$')

        surf = d.compsurfb
        mesurf = np.sum(surf * me / isotope.ufunc_A(b.ions)) - menorm
        ax.axhline(mesurf, color = 'gray', ls = ':', label = r'accretion')

        m_surf = np.max(m)
        xlim = (-0.02 * m_surf, 1.02 * m_surf)
        ax.set_xlim(xlim)

        ax.axvspan(xlim[0], 0, color = '#cccccc', label = 'substrate')

        if accmass is not None:
            ax.axvspan(*((np.array([zm0, zm1]) - d.zm[jblo]) / mscale),
                       color = '#ffffbb',
                       zorder = -2,
                       ls = '-',
                       lw = 0.5,
                       label = 'accretion column')

        xlabel = 'Mass above substrate / {}g'
        if mscale is None or mscale == 1:
            xlscale = ''
        else:
            xlscale = int(np.log10(mscale))
            xlscale = fr'$10^{{{xlscale:d}}}\,$'

        ax.set_xlabel(f'Mass above substrate / {xlscale}g')
        if abs(menorm) > 1e-5:
            offset = f'relative to {isotope.ion(norm).mpl()} '
        else:
            offset = ''
        ax.set_ylabel(f'Mass excess per nucleon {offset}/ MeV')

        ax.legend(loc = 'best')
        fig.tight_layout()
        fig.show()

        self.fig = fig
        self.ax = ax


import abuset
import mass_table

class ME(object):
    def __init__(self, mode = 'bdat', silent = False, norm = 'c12'):
        if mode == 'bdat':
            M = bdat.BDat()
            ions = M.ext_ions
            me = M.mass_excess(ions)
        elif mode in ('Ame12', 'Ame03', 'Ame16', ):
            M = mass_table.MassTable(version = mode)
            ions = M.ions
            me = M.mass_excess()

        norm = isotope.ion(norm)

        print(f' [ME] Computing ... ', end = '', flush = True)
        EF = np.linspace(0,25,25001)
        A = ions.A
        Z = ions.Z
        ionsidx = ions.idx
        mena = me[np.argwhere(norm == ions.ions())[0][0]] / norm.A
        x = (me[:,np.newaxis] + Z[:,np.newaxis] * EF[np.newaxis,:]) / A[:,np.newaxis] - mena - norm.Z / norm.A * EF[np.newaxis,:]
        ix = np.array([ionsidx[i] for i in np.argmin(x, axis=0)])
        ii = np.where(ix[1:] != ix[:-1])[0] + 1
        ii = np.insert(ii, 0, 0)
        ij = utils.index1d(ix[ii], ionsidx)
        imax = [(ions[xij], EF[xii], x[xij, xii]) for xij, xii in zip(ij, ii)]
        print(' [ME] done.')

        if not silent:
            print(f'{"Ion":5s}: {"EF/MeV":<6s} {"ME/A":>8s}')
            for i in imax:
                print(f'{i[0]!s:5s}: {float(i[1]):6.3f} {float(i[2]):8.5f}')

        self.imax = imax
        self.M = M
        self.ions = ions
        self.me = me
        self.x = x
        self.ij = ij
        self.ii = ii
        self.EF = EF
        self.mode = mode
        self.norm = norm
        self.mena = mena


    def plot(self, *, ax = None, fig = None):
        if ax is None:
            if fig is None:
                fig = plt.figure(figsize=(6.4,4.8+0.5))
            ax = fig.add_subplot(111)

        ii = np.append(self.ii, len(self.EF)-1)
        for ij in self.ij:
            ax.plot(self.EF, self.x[ij,:], alpha = 0.25)
        ax.set_prop_cycle(None)
        for i,ij in enumerate(self.ij):
            jj = slice(ii[i], ii[i+1])
            ax.plot(self.EF[jj], self.x[ij,jj], lw = 3, label = self.ions[ij].mpl())
        vals = [i[2] for i in self.imax]
        vals.append(self.x[self.ij[-1],-1])
        ymax = np.max(vals)
        ymin = np.min(vals)
        d = ymax - ymin
        ax.set_ylim(ymin - 0.1 * d, ymax + 0.1 * d )
        ax.legend(loc = 'best', title = self.mode)
        ax.set_xlabel('electron Fermi energy / MeV')
        ax.set_ylabel(f'mass excess per nucleon relative to {self.norm.LaTeX()} / MeV')
        fig.tight_layout()
        fig.show()
        self.fig = fig
        self.ax = ax

from abuset import IonList

class Heat(FigTemplate):
    def __init__(self, b = None, **kwargs):
        if not isinstance(b, bdat.BDat):
            b = bdat.BDat()
        ions = IonList(sorted(b.ext_ions, key = lambda x: 1024 * x.A - x.Z))
        A = ions.A
        ai = np.argwhere(A[1:] > A[:-1])[:, 0] + 1
        ai = np.insert(ai, (0, len(ai)), (0, len(A)))
        amax = np.max(ions.A) + 1
        nz = ai[1:] - ai[:-1]
        zmax = np.max(nz)
        me = b.mass_excess(ions)
        a = np.tile(None, (amax, zmax))
        e = np.tile(np.nan, (amax, zmax))
        de = np.tile(0., (amax, zmax))
        zi = np.tile(0, amax)
        zm = np.tile(0, amax)
        for iz, j0, j1 in zip(nz, ai[:-1], ai[1:]):
            Ai = A[j0]
            a[Ai, :iz] = ions[j0:j1]
            e[Ai, :iz] = me[j0:j1]
            zi[Ai] = iz
            de[Ai, :iz-1] = e[Ai, 1:iz] - e[Ai, :iz-1]

            # find first minimum
            ii = np.argwhere(de[Ai,:iz-1] > 0)
            if len(ii) > 0:
                im = ii[0,0]
            else:
                im = 0
            zm[Ai] = im

        mode = kwargs.get('mode', 'single')
        nuloss = kwargs.get('nuloss', 0.35)

        self.A = A
        self.a = a
        self.amax = amax
        self.de = de

        # compute heating for all isotopes
        ecmax = np.max(de)
        heat = np.tile(np.nan, (len(ions), zmax))
        ec = np.tile(np.nan, (len(ions), zmax))
        iion = 0
        for ia in range(amax):
            for iz in range(zi[ia]-1, -1, -1):
                i = iz
                i1 = i + 1
                if de[ia, i] < 0:
                    heat[iion + i, 0] = -de[ia, i]
                    ec[iion + i, 0] = 0.
                    if i1 <= zi[ia]-1:
                        heat[iion + i, 1:] = heat[iion + i1, :-1]
                        ec[iion + i, 1:] = ec[iion + i1, :-1]
                    continue
                if i < zi[ia]-2:
                    dq = 0.
                    if de[ia, i1] > de[ia, i]:
                        eci = de[ia, i]
                    elif mode == 'single':
                        eci = de[ia, i]
                        while i1 < zi[ia]-1 and eci > de[ia, i1]:
                            dq += eci - de[ia, i1]
                            i1 += 1
                    elif mode == 'double':
                        eci = 0.5 * (de[ia, i] + de[ia, i1])
                        i1 += 1
                        while i1 < zi[ia]-1:
                            if de[ia, i1] < eci:
                                dq += eci - de[ia, i1]
                                i1 += 1
                            else:
                                break
                    else:
                        raise Exception('uknown mode')
                    heat[iion + i, 0] = dq
                    ec[iion + i, 0] = eci
                    if i1 <= zi[ia]-1:
                        heat[iion + i, 1:] = heat[iion + i1, :-1]
                        ec[iion + i, 1:] = ec[iion + i1, :-1]
                    continue
                if i == zi[ia]-2:
                    heat[iion + i, 0] = 0.
                    ec[iion + i, 0] = de[ia, i]
                    heat[iion + i, 1:] = heat[iion + i1, :-1]
                    ec[iion + i, 1:] = ec[iion + i1, :-1]
                    continue
                heat[iion + i, 0] = 0.
                ec[iion + i, 0] = ecmax
            iion = iion + zi[ia]

        heat[np.isnan(heat)] = -1e10
        jj = heat > 0
        ii = np.where(jj)
        vals, kk = np.unique(ec[jj], return_inverse = True)
        hmap = np.zeros((len(ions), len(vals)))
        hmap[ii[0], kk] = heat[ii]
        data = (np.cumsum(hmap, axis = -1) / ions.A[:,np.newaxis]) * (1 - nuloss)

        self.heat = heat
        self.ec = ec
        self.ions = ions
        self.data = data
        self.vals = vals

        # compute stable subset projected to A
        staba = [j0 + zm[A[j0]] for j0 in ai[:-1]]
        dataa = data[staba,:]
        missing, = np.where(~np.in1d(np.arange(amax), A[staba]))
        dataa = np.insert(dataa, missing, np.zeros(len(vals)), axis = 0)
        dataa = dataa - dataa[:,0][:,np.newaxis]
        valsa, ii = np.unique(vals, return_index = True)
        dataa = dataa[:, ii]

        self.dataa = dataa
        self.valsa = valsa

    def plot_map(self, **kwargs):
        vals = self.valsa
        data = self.dataa
        amax = self.amax

        x = np.append(vals, vals[-1] * 1.05)
        y = np.arange(amax) + 0.5
        label = 'cummulative heating in MeV / nucleon'
        xlabel = 'Fermi Energy / MeV'
        ylabel = 'mass number'
        vmax = kwargs.get('vmax', 0.3)

        cmap = color.ColorBlendBspline(('white',)+tuple([color.ColorScale(color.colormap('viridis_r'), lambda x: (x-0.2)/0.7)]*2) + ('black',), frac=(0,.2,.9,1),k=1)

        fig, ax = self.get_fig(**kwargs)
        m = ax.pcolormesh(x, y, data[1:,:], cmap = cmap, vmax = vmax)
        cb = fig.colorbar(m, label = label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()
        fig.show()

        self.fig_map = fig
        self.ax_map = ax
        self.cb_map = cb

    def abuheatmap(self, abu, **kwargs):
        vals = self.vals

        try:
            if self.abu is abu:
                ii = self.ii_abu
            else:
                raise
        except:
            ii = utils.index1d(abu.idx(), self.ions.idx)
        heat = self.data[ii, :] * abu.X()[:,np.newaxis]
        heat, A = utils.project(heat, abu.A(), return_values = True)
        amax = np.max(A)
        missing, = np.where(~np.in1d(np.arange(amax-1)+1, A))
        heat = np.insert(heat, missing, np.zeros(len(vals)), axis = 0)

        x = np.append(vals, vals[-1] * 1.05)
        y = np.arange(amax) + 0.5

        cmap = color.colormap('viridis_wrb')

        vmax = None
        label = 'cummulative heating in MeV / nucleon'
        xlabel = 'Fermi Energy / MeV'
        ylabel = 'mass number'

        fig, ax = self.get_fig(**kwargs)
        m = ax.pcolormesh(x, y, heat[1:,:], cmap = cmap, vmax = vmax)
        cb = fig.colorbar(m, label = label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

        self.abu = abu
        self.ii_abu = ii

        self.fig = fig
        self.ax = ax

    def abuheat(self, abu, **kwargs):
        vals = self.vals

        ii = utils.index1d(abu.idx(), self.ions.idx)
        heat = self.data[ii, :] * abu.X()[:,np.newaxis]
        heat = np.sum(heat, axis = 0)

        xlabel = 'Fermi Energy / MeV'
        ylabel = 'cummulative heating in MeV / nucleon'

        fig, ax = self.get_fig(**kwargs)
        ax.plot(vals, heat)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

        self.fig = fig
        self.ax = ax

    def dumpheat(self, dump, **kwargs):
        if isinstance(dump, str):
            dump = kepdump.load(dump)

        x, xbase, xtop, xlabel, jbase = self.yscale(dump, **kwargs)

        iimap = slice(jbase, dump.jm+1)
        b = dump.ppnb[:,iimap]
        ionsb = abuset.IonList(dump.ionsb)
        try:
            if self.dump == dump:
                ii = self.ii_dump
            else:
                raise
        except:
            ii = utils.index1d(ionsb.idx, self.ions.idx)
        A = ionsb.A
        heat = np.sum((self.data[ii, -1] * A)[:,np.newaxis] * b[:,:], 0)

        ylabel = 'heating in MeV / nucleon'

        fig, ax = self.get_fig(**kwargs)

        ax.plot(x[iimap], heat)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xbase, xtop)

        fig.tight_layout()
        fig.show()

        self.ii_dump = ii
        self.dump = dump

        self.fig = fig
        self.ax = ax

    def dumpheatmap(self, dump, **kwargs):
        if isinstance(dump, str):
            dump = kepdump.load(dump)

        x, xbase, xtop, xlabel, jbase = self.yscale(dump, **kwargs)

        iimap = slice(jbase, dump.jm+1)
        b = dump.ppnb[:,iimap]
        ionsb = abuset.IonList(dump.ionsb)
        try:
            if self.dump == dump:
                ii = self.ii_dump
        except:
            ii = utils.index1d(ionsb.idx, self.ions.idx)
        A = ionsb.A

        mode = kwargs.get('mode', 'E')
        if mode == 'E':
            heat = np.sum((self.data[ii, :] * A[:,np.newaxis])[:,np.newaxis,:] * b[:,:,np.newaxis], 0)
            heat = np.transpose(heat[:,:])
            ylabel = 'electron Fermi energy / MeV'
            vals = self.vals
            y = np.append(vals, vals[-1] * 1.05)
        elif mode == 'A':
            heat = (self.data[ii, -1] * A)[:, np.newaxis] * b[:,:]
            heat, A = utils.project(heat, A, return_values = True)
            amax = np.max(A)
            missing, = np.where(~np.in1d(np.arange(amax-1)+1, A))
            heat = np.insert(heat, missing, 0., axis = 0)
            y = np.arange(amax+1) + 0.5
            ylabel = 'mass number'
        else:
            raise Exception('unkown mode {mode}')

        label = 'cummulative heating in MeV / nucleon'
        cmap = color.colormap('viridis_wrb')
        fig, ax = self.get_fig(**kwargs)
        m = ax.pcolormesh(x[jbase:], y, heat, cmap = cmap)
        cb = fig.colorbar(m, label = label)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xbase, xtop)

        fig.tight_layout()
        fig.show()

        self.ii_dump = ii
        self.dump = dump

        self.fig = fig
        self.ax = ax



import apng
import io
import os
import os.path
import matplotlib.pylab as plt

def movie_test(n = 11, filename = os.path.expanduser('~/test.png')):
    f = plt.figure(figsize = (6.4,4.8), dpi = 100)
    ax = f.add_subplot(111)
    ani = apng.APNG()
    for i in range(2 * n - 2):
        if i < n:
            ax.plot([i/(n-1),1-i/(n-1)])
        else:
            ax.plot([(i-n+1)/(n-1),(2*n-2-i)/(n-1)],[1,0])
        ax.set_ylim(0, 1)
        ax.set_xlim(0, 1)
        x = io.BytesIO()
        f.savefig(x, format='png')
        ax.clear()
        # ratio needs to lessim 90 for chrome
        ani.append_file(x, delay=1, delay_den=int(30+10*np.sin(i/n*2*np.pi)))
    ani.save(filename)

import matplotlib.animation as anim

def moviewriter_test(n = 21, filename = os.path.expanduser('~/test.mkv'), dpi = 100):
    fig = plt.figure(figsize = (6.4,4.8), dpi = dpi)
    ax = fig.add_subplot(111)
    moviewriter = anim.FFMpegWriter(
        fps = 30,
        # codec = 'libx264rgb', # mkv
        # extra_args = '-preset veryslow -crf 0'.split(),
        codec = 'apng', # apng
        extra_args = '-plays 0 -preset veryslow'.split(),
        # codec = 'ffv1', # avi
        # extra_args = '-preset veryslow'.split(),
        )
    with moviewriter.saving(fig, filename, dpi):
        for i in range(2 * n - 2):
            ax.clear()
            if i < n:
                ax.plot([i/(n-1),1-i/(n-1)])
            else:
                ax.plot([(i-n+1)/(n-1),(2*n-2-i)/(n-1)],[1,0])
            ax.set_ylim(0, 1)
            ax.set_xlim(0, 1)
            x = io.BytesIO()
            moviewriter.grab_frame()
