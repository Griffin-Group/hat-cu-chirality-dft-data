#!/usr/bin/env python3

"""
STM class copied from ASE
https://wiki.fysik.dtu.dk/ase/_modules/ase/dft/stm.html#STM.scan
https://wiki.fysik.dtu.dk/ase/ase/dft/stm.html

Adapted to accept CHGCAR/PARCHG format instead.
(This makes slightly less features than the original format which had full wavefunction data.)


# LICENSE

Copyright 2025 ASE developers
Copyright 2025 Bernard Field


stm.py is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

stm.py is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along
with stm.py. If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

from pymatgen.io.vasp.outputs import Chgcar

class STM:
    def __init__(self, parchg: str|Chgcar):
        """Scanning tunneling microscope.

        parchg: PARCHG, either file or Pymatgen object
        """

        if isinstance(parchg, str):
            self.parchg = Chgcar.from_file(parchg)
        else:
            self.parchg = parchg
        self.ldos = self.parchg.data['total']
        self.cell = self.parchg.structure.lattice.matrix


    def get_averaged_current(self, z):
        """Calculate avarage current at height z (in Angstrom).

        Use this to get an idea of what current to use when scanning."""

        nz = self.ldos.shape[2]

        # Find grid point:
        n = z / self.cell[2, 2] * nz
        dn = n - np.floor(n)
        n = int(n) % nz

        # Average and do linear interpolation:
        return ((1 - dn) * self.ldos[:, :, n].mean() +
                dn * self.ldos[:, :, (n + 1) % nz].mean())


    def scan(self, current, z0=None, repeat=(1, 1)):
        """Constant current 2-d scan.

        Returns three 2-d arrays (x, y, z) containing x-coordinates,
        y-coordinates and heights.  These three arrays can be passed to
        matplotlib's contourf() or pcolormesh() function like this:

        >>> import matplotlib.pyplot as plt
        >>> plt.contourf(x, y, z)
        >>> plt.show()

        """

        L = self.cell[2, 2]
        nz = self.ldos.shape[2]
        h = L / nz

        ldos = self.ldos.reshape((-1, nz))

        heights = np.empty(ldos.shape[0])
        for i, a in enumerate(ldos):
            heights[i] = find_height(a, current, h, z0)

        s0 = heights.shape = self.ldos.shape[:2]
        heights = np.tile(heights, repeat)
        s = heights.shape

        ij = np.indices(s, dtype=float).reshape((2, -1)).T
        x, y = np.dot(ij / s0, self.cell[:2, :2]).T.reshape((2,) + s)

        return x, y, heights


    def scan2(self, z, repeat=(1, 1)):
        """Constant height 2-d scan.

        Returns three 2-d arrays (x, y, I) containing x-coordinates,
        y-coordinates and currents.  These three arrays can be passed to
        matplotlib's contourf() or pcolormesh() function like this:

        >>> import matplotlib.pyplot as plt
        >>> plt.contourf(x, y, I)
        >>> plt.show()

        """

        nz = self.ldos.shape[2]
        ldos = self.ldos.reshape((-1, nz))

        current = np.empty(ldos.shape[0])

        zp = z / self.cell[2, 2] * nz
        zp = int(zp) % nz

        for i, a in enumerate(ldos):
            current[i] = self.find_current(a, zp)

        s0 = current.shape = self.ldos.shape[:2]
        current = np.tile(current, repeat)
        s = current.shape

        ij = np.indices(s, dtype=float).reshape((2, -1)).T
        x, y = np.dot(ij / s0, self.cell[:2, :2]).T.reshape((2,) + s)

        # Returing scan with axes in Angstrom.
        return x, y, current


    def linescan(self, current, p1, p2, npoints=50, z0=None):
        """Constant current line scan.

        Example::

            stm = STM(...)
            z = ...  # tip position
            c = stm.get_averaged_current(z)
            stm.linescan(c, (1.2, 0.0), (1.2, 3.0))
        """

        heights = self.scan(current, z0)[2]

        p1 = np.asarray(p1, float)
        p2 = np.asarray(p2, float)
        d = p2 - p1
        s = np.dot(d, d)**0.5

        cell = self.cell[:2, :2]
        shape = np.array(heights.shape, float)
        M = np.linalg.inv(cell)
        line = np.empty(npoints)
        for i in range(npoints):
            p = p1 + i * d / (npoints - 1)
            q = np.dot(p, M) * shape
            line[i] = interpolate(q, heights)
        return np.linspace(0, s, npoints), line


    def pointcurrent(self, x, y, z):
        """Current for a single x, y, z position."""

        nx = self.ldos.shape[0]
        ny = self.ldos.shape[1]
        nz = self.ldos.shape[2]

        # Find grid point:
        xp = x / np.linalg.norm(self.cell[0]) * nx
        dx = xp - np.floor(xp)
        xp = int(xp) % nx

        yp = y / np.linalg.norm(self.cell[1]) * ny
        dy = yp - np.floor(yp)
        yp = int(yp) % ny

        zp = z / np.linalg.norm(self.cell[2]) * nz
        dz = zp - np.floor(zp)
        zp = int(zp) % nz

        # 3D interpolation of the LDOS at point (x,y,z) at given bias.
        xyzldos = (((1 - dx) + (1 - dy) + (1 - dz)) * self.ldos[xp, yp, zp] +
                   dx * self.ldos[(xp + 1) % nx, yp, zp] +
                   dy * self.ldos[xp, (yp + 1) % ny, zp] +
                   dz * self.ldos[xp, yp, (zp + 1) % nz])

        return dos2current(xyzldos)


    def find_current(self, ldos, z):
        """ Finds current for given LDOS at height z."""
        nz = self.ldos.shape[2]

        zp = z / self.cell[2, 2] * nz
        dz = zp - np.floor(zp)
        zp = int(zp) % nz

        ldosz = (1 - dz) * ldos[zp] + dz * ldos[(zp + 1) % nz]

        return dos2current(ldosz)



def dos2current(dos):
    # Borrowed from gpaw/analyse/simple_stm.py:
    # The connection between density n and current I
    # n [e/Angstrom^3] = 0.0002 sqrt(I [nA])
    # as given in Hofer et al., RevModPhys 75 (2003) 1287
    return 5000. * dos**2


def interpolate(q, heights):
    qi = q.astype(int)
    f = q - qi
    g = 1 - f
    qi %= heights.shape
    n0, m0 = qi
    n1, m1 = (qi + 1) % heights.shape
    z = (g[0] * g[1] * heights[n0, m0] +
         f[0] * g[1] * heights[n1, m0] +
         g[0] * f[1] * heights[n0, m1] +
         f[0] * f[1] * heights[n1, m1])
    return z


def find_height(ldos, current, h, z0=None):
    if z0 is None:
        n = len(ldos) - 2
    else:
        n = int(z0 / h)
    while n >= 0:
        if ldos[n] > current:
            break
        n -= 1
    else:
        return 0.0

    c2, c1 = ldos[n:n + 2]
    return (n + 1 - (current - c1) / (c2 - c1)) * h


# Plotting routines, for convenience.
# I wrote these ones - Bernard Field 2025.
def plot_scan_basic(stm: STM, current: float, z0=None, repeat=(1,1),
                    show: bool = True, **kwargs) -> tuple[plt.Figure, plt.Axes]:
    """
    Basic plotting of constant current scan. Mostly as an example.
    
    Give it an STM object.
    Specify constant current for scan (nA).
    z0 is middle of the vacuum. Defaults to the top of the unit cell.
    (So if your system is at z=0, you'll need to adjust this.)
    repeat is how many unit cells to plot.

    kwargs are passed to matplotlib.pyplot.pcolormesh
    Good options are
    - cmap ('Gray_r' or 'plasma' are nice).
    - vmin (Very important if there are holes in your structure, to constrain
        the colour range).
    """
    x, y, z = stm.scan(current, z0=z0, repeat=repeat)
    fig, ax = plt.subplots()
    ax.pcolormesh(x, y, z, **kwargs)
    if show:
        plt.show()
    return fig, ax


def plot_scan_nice(stm: STM, current: float, z0: float|None = None,
                   show: bool = True, save: str|None = None,
                   cell_args: dict = {'edgecolor':'g'},
                   ncells: int = 1,
                   dpi = None,
                   **kwargs) -> tuple[plt.Figure, plt.Axes]:
    """
    Plot constant-current STM, with nice display (no whitespace).

    Inputs:
        - stm: STM object, with data loaded. (If you want to rotate the data,
            edit lattice vectors in stm.cell).
        - current: target current for constant current, in nano Amperes.
        - z0: Height in the vacuum region. Defaults to top of the unit cell
            (So if your system is at z=0, you'll need to adjust this).
        - show: Whether to plt.show()
        - save: Filename to save the plot to.
        - cell_args: keyword arguments for matplotlib.patches.Polygon.
            'closed=True' and 'fill=False' are already taken.
            You should specify your line style.
            If empty or None, no unit cell is plotted.
        - ncells: Number of unit cells to plot (default 1).
        - dpi: passed to fig.savefig.
    Further keyword arguments are passed to plt.pcolormesh.
    Recommended keyword arguments:
        - cmap: Colour map. 'Grays_r' and 'plasma' are good options.
        - vmin: Minimum bound for colour map. Useful if there are holes in
            your structure, meaning the scan punches through to somewhere
            far below the actual structure and messes with the contrast.
    """
    # Supercell to fill a rectangular area.
    repeat = (3*ncells,3*ncells)
    # Lattice vectors
    a = stm.cell[0, 0:2]
    b = stm.cell[1, 0:2]
    # Get data
    x, y, z = stm.scan(current, z0=z0, repeat=repeat)
    if z.min() == 0 and 'vmin' not in kwargs:
        print("Warning! Scan went all the way to z=0. You should set vmin to preserve contrast.")
    # Plot
    fig, ax = plt.subplots()
    ax.pcolormesh(x, y, z, **kwargs)
    ax.set_aspect(1) # The dimensions matter. Don't skew the plot.
    # Draw unit cell
    points = np.array([a+b, 2*a+b, 2*(a+b), a+2*b])
    if cell_args:
        polygon = mpatch.Polygon(points, closed=True,
                                fill=False, **cell_args)
        ax.add_artist(polygon)
    # Adjust plotting range
    xmin, ymin = points.min(axis=0) * ncells
    xmax, ymax = points.max(axis=0) * ncells
    # Add a little padding
    padx = (xmax - xmin) * 0.02
    pady = (ymax - ymin) * 0.02
    ax.set_xlim(xmin - padx, xmax + padx)
    ax.set_ylim(ymin - pady, ymax + pady)
    # Remove extraneous annotation
    ax.axis('off')
    # Output
    if save:
        # Save with no whitespace, just the content.
        fig.savefig(save, bbox_inches='tight', pad_inches=0, dpi=dpi)
    if show:
        plt.show()
    return fig, ax
