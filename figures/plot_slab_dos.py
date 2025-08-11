#!/usr/bin/env python3
"""
Plot the PDOS of the achiral and Ag slab systems in one plot.

Copyright 2025 Bernard Field
MIT License
"""

import os

import matplotlib.pyplot as plt

from sumo.cli.dosplot import dosplot
from sumo.plotting import styled_plot, sumo_base_style, sumo_dos_style

SCRIPTDIR = os.path.split(os.path.abspath(__file__))[0]
ROOTDIR = os.path.split(SCRIPTDIR)[0]
FIGDIR = SCRIPTDIR

if __name__ == "__main__":
    # Wrap the call to plt.subplots with the appropriate sumo style.
    fig, axs = styled_plot(sumo_base_style, sumo_dos_style)(plt.subplots)(
            nrows=3,
            figsize=(4.4,3.2),
            gridspec_kw = dict(top=0.98, hspace=0.05, bottom=0.13),
            )
    params = dict(
            elements=dict(C=tuple(), Cu=tuple(), N=tuple()),
            zero_line=True,
            xmin=-1,
            xmax=2.5,
            colours=dict(N=dict(p="#D93B2B"),
                         C=dict(p="#17479E"),
                         Cu=dict(d="#426600"),
                         ),
            )
    # Plot one
    plt.sca(axs[0])
    fname = os.path.join(ROOTDIR,'HAT-Cu-AgSlab5','2b.dos','vasprun.xml.gz')
    if not os.path.exists(fname):
        # The downsampled DOS files have no loss of quality, so I don't need to ask.
        fname = os.path.join(ROOTDIR,'HAT-Cu-AgSlab5','2b.dos','vasprun_ds.xml.gz')
    dosplot(fname,
            **params,
            plot_total=False,
            yscale=3,
            plt=plt,
        )
    ax = plt.gca()
    ax.set_xlabel(None)
    ax.set_xticklabels([])
    ax.set_ylabel(None)
    ax.text(0.5, 0.98, 'Ag slab (5 layer)', transform=ax.transAxes, ha='center', va='top')
    # Plot two
    plt.sca(axs[1])
    fname = os.path.join(ROOTDIR,'HAT-Cu-AgSlab','3b.dos','vasprun.xml.gz')
    if not os.path.exists(fname):
        fname = os.path.join(ROOTDIR,'HAT-Cu-AgSlab','3b.dos','vasprun_ds.xml.gz')
    dosplot(fname,
            **params,
            plot_total=False,
            yscale=3,
            legend_on=False,
            plt=plt,
        )
    ax = plt.gca()
    ax.set_xlabel(None)
    ax.set_xticklabels([])
    ax.text(0.5, 0.98, 'Ag slab (3 layer)', transform=ax.transAxes, ha='center', va='top')
    # Plot 3
    plt.sca(axs[2])
    fname = os.path.join(ROOTDIR,'HAT-Cu_free','3c.dos','vasprun.xml.gz')
    dosplot(fname,
            **params,
            legend_on=False,
            yscale=2,
            plt=plt,
        )
    ax = plt.gca()
    ax.set_ylabel(None)
    ax.text(0.5, 0.98, 'Free-standing', transform=ax.transAxes, ha='center', va='top')
    # Export
    plt.savefig(os.path.join(FIGDIR,'slab_dos_comparison.pdf'))
