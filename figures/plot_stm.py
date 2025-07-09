#!/usr/bin/env python3
"""
Plot simulated STM plots.

Copyright 2025 Bernard Field
Distributed under the MIT license.
"""

import sys
import os
import argparse

import matplotlib.pyplot as plt
import numpy as np

FIGDIR = os.path.split(os.path.abspath(__file__))[0]
ROOTDIR = os.path.split(FIGDIR)[0]

# Make sure the script directory is in the path before loading these modules.
sys.path.append(FIGDIR)
from stm import STM, plot_scan_nice

def diagnose_stm(stm:STM, current:float, **kwargs):
    """
    Plot a histogram of all the heights values from an STM scan.

    All keyword arguments passed to plt.hist
    """
    _, _, z = stm.scan(current)
    plt.hist(z.flat, **kwargs)
    plt.show()

def rotate_PW_cell(stm:STM):
    """
    Rotate the PW unit cell by 30 degrees, in-place
    Is NOT general purpose. Very specific use case only.
    """
    cell = stm.cell.copy()
    cell[0] = [np.linalg.norm(cell[0]), 0, 0]
    cell[1] = [cell[1,1], cell[1,0], 0]
    stm.cell = cell

if __name__ == "__main__":
    # Allow specifying what figures we calculate.
    all_choices = ['PW_1','PW_2','PW_2a','PW_2b','PW_2c','PW_3','PW_4',
                   'Y_1','Y_2','Y_2a','Y_2b','Y_2c','Y_3','Y_4']
    parser = argparse.ArgumentParser()
    parser.add_argument('include', nargs='*',
                        metavar=','.join(all_choices),
                        help="Which plots to generate (default all).")
    parser.add_argument('-d','--downsampled', action='store_true',
                        help="Use down-sampled PARCHG's instead.")
    args = parser.parse_args()
    if any(arg not in all_choices for arg in args.include):
        parser.error(f"argument include: invalid choice(s): {args.include}. (Choose from {', '.join(all_choices)})")
    if not args.include:
        args.include = all_choices
    # Suffix for down-sampled data and plots
    if args.downsampled:
        ds_str = '_ds'
    else:
        ds_str = ''
    # Locations of data.
    PWDIR = os.path.join(ROOTDIR, 'HAT-Cu_PW', '2d.parchg')
    YDIR = os.path.join(ROOTDIR, 'HAT-Cu_Y', '2d.parchg')
    STMFIGDIR = os.path.join(FIGDIR, 'stm_figures')
    # Configure plotting parameters.
    # Global parameters
    params = dict(current=0.1, cmap='plasma', show=False)
    # Parameters for specific plots (defaults to none)
    spec_params = {i:dict() for i in all_choices}
    for k in spec_params:
        if 'Y' in k:
            spec_params[k]['cell_args'] = {'edgecolor':'magenta'}
        else: # PW in k
            spec_params[k]['cell_args'] = {'edgecolor':'cyan'}
    spec_params['Y_1']['vmin'] = 6.8
    spec_params['Y_2']['vmin'] = 7.2
    for i in ['Y_2a','Y_2b','Y_2c']:
        spec_params[i]['vmin'] = 7.5
        spec_params[i]['vmax'] = 11.522
    spec_params['Y_3']['vmin'] = 7.5
    spec_params['Y_4']['z0'] = 15
    spec_params['PW_1']['vmin'] = 7.5
    for i in ['PW_2a','PW_2b','PW_2c']:
        spec_params[i]['vmin'] = 7.5
        spec_params[i]['vmax'] = 11.454
    spec_params['PW_4']['z0'] = 15
    # Whether to plot a second plot, labelled "tight", with different params.
    tight_params = {'Y_2':dict(vmin=9.3),
                    'Y_3':dict(vmin=9.3),
                    'PW_2':dict(vmin=9.3),
                    'PW_3':dict(vmin=9.3),
                   }
    for k in tight_params:
        if 'Y' in k:
            tight_params[k]['cell_args'] = {'edgecolor':'magenta'}
        else: # PW in k
            tight_params[k]['cell_args'] = {'edgecolor':'cyan'}
    # Load and plot data
    for i in args.include:
        print(f"Loading #{i}")
        num = i.split('_')[1] # Get the numeric part.
        # Get the right data directory.
        if 'PW' in i:
            mydir = PWDIR
        elif 'Y' in i:
            mydir = YDIR
        else:
            raise ValueError("Unexpected band/system ", i)
        stm = STM(os.path.join(mydir, f"PARCHG_{num}{ds_str}.gz"))
        if 'PW' in i:
            rotate_PW_cell(stm)
        plot_scan_nice(stm, **params, **spec_params[i],
                       save=os.path.join(STMFIGDIR, f'band{i}_constant_current{ds_str}.png'))
        if i in tight_params:
            plot_scan_nice(stm, **params, **tight_params[i],
                           save=os.path.join(STMFIGDIR, f'band{i}_constant_current_tight{ds_str}.png'))
        if i in ['PW_2','Y_2']:
            # Plot without the unit cell.
            params2 = spec_params[i].copy()
            del params2['cell_args']
            plot_scan_nice(stm, **params, **params2,
                           cell_args = None,
                           ncells = 4, dpi=300,
                           save=os.path.join(STMFIGDIR, f'band{i}_constant_current_large{ds_str}.png'))

    # Done.
