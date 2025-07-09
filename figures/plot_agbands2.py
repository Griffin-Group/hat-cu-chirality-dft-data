#!/usr/bin/env python3
"""
Plot HAT-Cu-AgSlab5 band structure with multiplicative blending of fatbands.

Copyright 2025 Bernard Field
MIT License
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
from sumo.plotting.bs_plotter import SBSPlotter
from sumo.electronic_structure.bandstructure import get_projections
from pymatgen.io.vasp import BSVasprun

SCRIPTDIR = os.path.split(os.path.abspath(__file__))[0]
ROOTDIR = os.path.split(SCRIPTDIR)[0]
FIGDIR = os.path.join(SCRIPTDIR, 'band_figures')

# Make sure the script directory is in the path before loading these modules.
sys.path.append(SCRIPTDIR)
import blend_plots
import sumo_band_projections

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--downsampled',action='store_true',
                        help="Use vasprun_ds.xml.gz (reduced precision) instead.")
    args = parser.parse_args()
    if args.downsampled:
        ds_str = '_ds'
    else:
        ds_str = ''
    fname = f'vasprun{ds_str}.xml.gz'
    print("Loading Vasprun")
    vasp = BSVasprun(os.path.join(ROOTDIR, 'HAT-Cu-AgSlab5', '2c.bands', fname),
                     parse_projected_eigen=True,
                     parse_potcar_file=False)
    print("Compiling BandStructure")
    bs = vasp.get_band_structure()
    plotter = SBSPlotter(bs)
    # To blend the plots together, we require them to be EXACTLY the same size
    # with the axes elements in EXACTLY the same place.
    # This means no auto-layout, we're manually placing eveything.
    # And the defaults are too narrow to fit all elements.
    # All measurements are in inches.
    axwidth = 2.0
    axheight = 2.0
    figwidth = 3.2
    figheight = 2.2
    lpad = 0.6
    bpad = 0.1
    style = {"figure.figsize": (figwidth,figheight), 
             "figure.subplot.left": lpad/figwidth,
             "figure.subplot.right": (lpad+axwidth)/figwidth,
             "figure.subplot.bottom": bpad/figheight,
             "figure.subplot.top": (bpad+axheight)/figheight}
    kwargs = dict(mode='stacked', zero_line=True, zero_energy=0.2056,
                  ymin=-1, ymax=2.5, normalise=None, legend=False,
                  circle_size=10, style=style,
                  interpolate_factor=3)
    selections = [['C'], ['N'], ['Cu']]
    colors = ['#ff4444','#4444ff','#44ff44']
    figs = [None, None, None]
    # Create first figure, for C.
    for i in range(3):
        print("Plotting ", selections[i][0])
        proj = get_projections(plotter.bs, selections[i], normalise=None)
        sumo_band_projections.get_plot_with_projections(plotter, proj, colours=[colors[i]],
                                                        **kwargs)
        ax = plt.gca()
        figs[i] = plt.gcf()
        # Replace legend with fake legend
        legend_colors = [[0]*4] * 3
        legend_colors[i] = colors[i]
        sumo_band_projections.fake_legend(ax, legend_colors, ['C','N','Cu'])
        # Remove background
        figs[i].patch.set_facecolor("none")
        ax.patch.set_facecolor("none")
    # Blend them together
    print("Blending")
    img = blend_plots.blend_figures(*figs, mode='multiply')
    # Save
    blend_plots.save_img(os.path.join(FIGDIR, f'HAT-Cu-AgSlab5_blended{ds_str}.png'),
                         img, dpi=400)
