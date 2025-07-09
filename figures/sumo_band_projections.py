#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Supplemental functions for Sumo to allow plotting bands with different
projections, like specific atomic sites, spin polarised projections, or
arbitrary projections.

Core is `get_plot_with_projections`, where you pass arbitrary weights
for each eigenvalue and selection (i.e. pre-computed projections).
It also includes a few more customisation options, like manually 
specifying colour or drawing to a pre-existing Axes.
It otherwise functions identically to `SBSPlotter.get_projected_plot`.

`get_atomic_projections` allows acquiring projections onto atoms
specified by index.

`get_magnetized_projections` obtains atomic projections polarised 
along a spin axis in a noncollinear calculation. This is useful for
spin-polarised edge states in SOC calculations, for example.

`get_atomic_projected_plot`, `get_spin_projected_plot`, and
`get_two_spin_plot` are convenience wrappers for the above functions.


Original functions taken from sumo (https://github.com/SMTG-Bham/sumo)
Copyright 2017 Alex Ganose

Modifications by Bernard Field
Copyright 2025 Bernard Field
MIT License
"""

import itertools as it
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rcParams
import matplotlib.lines as mlines
from scipy.interpolate import interp1d


from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.io.vasp import Vasprun

from sumo.electronic_structure.bandstructure import (
    force_branches,
)
from sumo.plotting import (
    colorline,
    pretty_plot,
    pretty_subplot,
    styled_plot,
    sumo_base_style,
    sumo_bs_style,
)

from sumo.plotting.bs_plotter import SBSPlotter


def proj_to_proj_by_branches(bs:BandStructure, projections):
    """Returns projections for each branch in a band structure.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
        projections:
            Precomputed projections.
            [{spin:projection[band,kpoint]},selection2,selection3,...]

    Returns:
        list: A ``list`` of projections for each branch of the band
        structure, with the format::

            [ [ {spin: projections} ], [ {spin: projections} ], ... ]

        Where spin is a :obj:`pymatgen.electronic_structure.core.Spin`
        object and projections is a :obj:`numpy.array` of::

            projections[band_index][kpoint_index]

        If there are no projections in the band structure, then an array of
        zeros is returned for each spin.
    """
    spins = bs.bands.keys()

    branches = []
    for b in bs.branches:
        s = b["start_index"]
        e = b["end_index"] + 1

        branch_proj = deepcopy(projections)
        for spin, i in it.product(spins, range(len(projections))):
            branch_proj[i][spin] = projections[i][spin][:, s:e]

        branches.append(branch_proj)
    return branches

@styled_plot(sumo_base_style, sumo_bs_style)
def get_plot_with_projections(
        self:SBSPlotter,
        proj,
        mode="rgb",
        normalise="all",
        interpolate_factor=4,
        color1="#FF0000",
        color2="#0000FF",
        color3="#00FF00",
        colorspace="lab",
        circle_size=150,
        projection_cutoff=0.001,
        zero_energy=None,
        zero_to_efermi=True,
        zero_line=False,
        ymin=-6.0,
        ymax=6.0,
        width=None,
        height=None,
        vbm_cbm_marker=False,
        ylabel="Energy (eV)",
        dpi=400,
        plt=None,
        dos_plotter=None,
        dos_options=None,
        dos_label=None,
        plot_dos_legend=True,
        dos_aspect=3,
        aspect=None,
        fonts=None,
        style=None,
        no_base_style=False,
        spin=None,
        title=None,
        labels=None,
        alpha=1,
        ax=None,
        interpolate=True,
        colours=["#3952A3", "#FAA41A", "#67BC47", "#6ECCDD", "#ED2025"],
        legend=True,
        ):
    """
    This method is copy-paste of SBSPlotter.get_projected_plot, but instead
    of giving a selection of projections you feed it the pre-processed
    projections.
    
    proj:: [{spin:projection[band,kpoint]},selection2,selection3,...]
    Where spin is a :obj:`pymatgen.electronic_structure.core.Spin` object
    and projection is a :obj: `numpy.array`.
    
    REMEMBER TO force_branches(bs) before making the projection!!
    Or, alternatively, make the projection from self.bs, as it has already
    done so. But never nest force_branches, as that just makes more branches.
    
    New arguments:
        labels: list of labels for the different selections.
        alpha: transparency (good for when we have degenerate lines)
        ax: axes to plot to (necessary for plotting to inset axes)
        interpolate: whether or not to interpolate for extra smoothness
        colours: list of colours for 'stacked' mode.
        legend: whether to draw legend.
    """
    nselections = len(proj)
    if mode == "rgb" and nselections > 3:
        raise ValueError("Too many elements/orbitals specified (max 3)")
    elif mode == "solo" and dos_plotter:
        raise ValueError("Solo mode plotting with DOS not supported")
    
    if ax is None:
        if dos_plotter:
            plt = pretty_subplot(
                1,
                2,
                width,
                height,
                sharex=False,
                dpi=dpi,
                plt=plt,
                gridspec_kw={"width_ratios": [dos_aspect, 1], "wspace": 0},
            )
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

    data = self.bs_plot_data(zero_to_efermi=zero_to_efermi)
    if zero_energy is not None:
        data = self._reset_zero_energy(data, zero_energy=zero_energy)

    nbranches = len(data["distances"])

    # Ensure we do spin up first, then spin down
    spins = sorted(self.bs.bands.keys(), key=lambda s: -s.value)
    if spin is not None and len(spins) == 1:
        raise ValueError(
            "Spin-selection only possible with spin-polarised "
            "calculation results"
        )

    if spin is Spin.up:
        spins = [spins[0]]
    elif spin is Spin.down:
        spins = [spins[1]]

    #proj = get_projections_by_branches(self.bs, selection, normalise=normalise)
    proj = proj_to_proj_by_branches(self.bs, proj)

    # nd is branch index
    for spin, nd in it.product(spins, range(nbranches)):
        # mask data to reduce plotting load
        bands = np.array(data["energy"][str(spin)][nd])
        mask = np.where(
            np.any(bands > ymin - 0.05, axis=1)
            & np.any(bands < ymax + 0.05, axis=1)
        )
        distances = data["distances"][nd]
        bands = bands[mask]
        weights = [proj[nd][i][spin][mask] for i in range(nselections)]

        if interpolate and len(distances) > 2:  # Only interpolate if it makes sense to do so
            # interpolate band structure to improve smoothness
            temp_dists = np.linspace(
                distances[0], distances[-1], len(distances) * interpolate_factor
            )
            bands = interp1d(
                distances,
                bands,
                axis=1,
                bounds_error=False,
                fill_value="extrapolate",
            )(temp_dists)
            weights = interp1d(
                distances,
                weights,
                axis=2,
                bounds_error=False,
                fill_value="extrapolate",
            )(temp_dists)
            distances = temp_dists

        else:  # change from list to array if we skipped the scipy interpolation
            weights = np.array(weights)
            bands = np.array(bands)
            distances = np.array(distances)

        # sometimes VASP produces very small negative weights
        weights[weights < 0] = 0

        if mode == "rgb":
            # colours aren't used now but needed later for legend
            colours = [color1, color2, color3]
            
            # If only one ornital, then we just use red
            if len(weights) == 1:
                weights = np.insert(weights, 1, np.zeros(weights[0].shape), axis=0)
                weights = np.insert(weights, 2, np.zeros(weights[0].shape), axis=0)
                colours = [color1]
            # if only two orbitals then just use red and blue
            if len(weights) == 2:
                weights = np.insert(weights, 2, np.zeros(weights[0].shape), axis=0)
                colours = [color1, color2]

            ls = "-" if spin == Spin.up else "--"
            lc = colorline(
                distances,
                bands,
                weights.transpose((1, 2, 0)),
                color1=color1,
                color2=color2,
                color3=color3,
                colorspace=colorspace,
                linestyles=ls,
                linewidth=(rcParams["lines.linewidth"] * 1.25),
            )
            lc.set_alpha(alpha)
            ax.add_collection(lc)

        elif mode == "stacked":
            # TODO: Handle spin

            # extend the specified custom colours with the default colours
            colour_series = rcParams["axes.prop_cycle"].by_key()["color"]
            colours.extend(colour_series)

            # very small circles look crap
            weights[weights < projection_cutoff] = 0

            distances = list(distances) * len(bands)
            bands = bands.flatten()
            zorders = range(-len(weights), 0)
            for w, c, z in zip(weights, colours, zorders):
                ax.scatter(
                    distances,
                    bands,
                    c=c,
                    s=circle_size * w**2,
                    zorder=z,
                    rasterized=True,
                    alpha=alpha,
                )

    # plot the legend
    if labels is None:
        labels = [str(x) for x in range(nselections)]
    elif len(labels) < nselections:
        for x in range(len(labels),nselections):
            labels[x] = str(x)
    for c, label in zip(colours, labels):
        ax.scatter([-10000], [-10000], c=c, s=50, label=label, edgecolors="none")

    if dos_plotter:
        loc = 1
        anchor_point = (-0.2, 1)
    else:
        loc = 2
        anchor_point = (0.95, 1)
    
    if legend:
        ax.legend(
            bbox_to_anchor=anchor_point,
            loc=loc,
            frameon=False,
            handletextpad=0.1,
            scatterpoints=1,
        )

    # finish and tidy plot
    self._maketicks(ax, ylabel=ylabel)
    self._makeplot(
        ax,
        ax.figure,
        data,
        zero_to_efermi=zero_to_efermi,
        zero_line=zero_line,
        vbm_cbm_marker=vbm_cbm_marker,
        width=width,
        height=height,
        ymin=ymin,
        ymax=ymax,
        dos_plotter=dos_plotter,
        dos_options=dos_options,
        dos_label=dos_label,
        plot_dos_legend=plot_dos_legend,
        aspect=aspect,
        spin=spin,
        #title=title, # I must have an older version of sumo as it doesn't recognise this word.
    )
    return plt

def get_atomic_projections(bs:BandStructure, selection:list, normalise:str|None=None) -> list:
    """
    Returns atomic projections from a band structure.

    Args:
        bs: The band structure.
        selection: a list-like of list-like of integers, atomic indices to project
            onto (0-based).
        normalise: 'all', 'select', or None.
    Returns:
        [{spin:projection[band,kpoint]},selection2,selection3,...]
    """
    spins = bs.bands.keys()
    nbands = bs.nb_bands
    nkpts = len(bs.kpoints)
    
    #bs.projections is [spin][bands, kpoints, orbitals, ions]
    
    # Normalisation factor
    sum_proj = {s: np.zeros((nbands, nkpts)) for s in spins}
    if normalise == "all":
        for s in spins:
            sum_proj[s] += bs.projections[s].sum(axis=(2,3))
    
    spec_proj = []
    for atoms in selection:
        spec_proj.append(
                {s: bs.projections[s][:,:,:,atoms].sum(axis=(2,3)) for s in spins}
                )
        if normalise == 'select':
            for s in spins:
                sum_proj[s] += bs.projections[s][:,:,:,atoms].sum(axis=(2,3))
    
    # Normalise
    if normalise:
        # to prevent warnings/errors relating to divide by zero,
        # catch warnings and surround divide with np.nan_to_num
        with np.errstate(divide="ignore", invalid="ignore"):
            for spin, i in it.product(spins, range(len(spec_proj))):
                spec_proj[i][spin] = np.nan_to_num(spec_proj[i][spin] / sum_proj[spin])
    
    return spec_proj

def get_magnetized_projections(vasprun:Vasprun, selection:list, normalise:str|None=None, spin_axis:int=2) -> tuple[list, list]:
    """
    Like get_atomic_projections, but splits the results into spin-polarised
    sectors for noncollinear spin calculations.
    
    Returns two projections: one where m>0 along the spin_axis,
        one where m<0 along the spin_axis.
    """
    bs = vasprun.get_band_structure()
    mag = vasprun.projected_magnetisation[...,spin_axis] # mz magnetization
    # Magnetisation is [kpoints, bands, ions, orbitals, mag]
    # Projection is [bands, kpoints, orbitals, ions]
    # Account for force_branches behaviour.
    dup_ids = get_force_branches_dup_ids(bs)
    mag = mag[dup_ids]
    bs = force_branches(bs)
    # Rearrange mag to match projection
    mag = np.transpose(mag, (1,0,3,2))
    
    spins = bs.bands.keys()
    nbands = bs.nb_bands
    nkpts = len(bs.kpoints)
    # Normalisation factor
    sum_proj = {s: np.zeros((nbands, nkpts)) for s in spins}
    if normalise == "all":
        for s in spins:
            sum_proj[s] += (mag * (mag > 0)).sum(axis=(2,3))
            sum_proj[s] += (-mag * (mag<0)).sum(axis=(2,3))
    
    # Get the projections onto each atom.
    proj_up = []
    proj_dn = []
    
    for atoms in selection:
        proj_up.append(
            {s: (mag[:,:,:,atoms] * (mag[:,:,:,atoms]>0)).sum(axis=(2,3)) for s in spins}
            )
        proj_dn.append(
            {s: (-mag[:,:,:,atoms] * (mag[:,:,:,atoms]<0)).sum(axis=(2,3)) for s in spins}
            )
        if normalise == 'select':
            for s in spins:
                sum_proj[s] += proj_up[-1][s]
                sum_proj[s] += proj_dn[-1][s]
    # Normalise
    # Each band and k-point is normalised separately.
    if normalise:
        # to prevent warnings/errors relating to divide by zero,
        # catch warnings and surround divide with np.nan_to_num
        with np.errstate(divide="ignore", invalid="ignore"):
            for spin, i in it.product(spins, range(len(proj_up))):
                proj_up[i][spin] = np.nan_to_num(proj_up[i][spin] / sum_proj[spin])
                proj_dn[i][spin] = np.nan_to_num(proj_dn[i][spin] / sum_proj[spin])
    
    return proj_up, proj_dn

def get_force_branches_dup_ids(bandstructure:BandStructure):
    """
    Get the dup_ids from sumo.electronic_structure.bandstructure.force_branches
    Sumo is really hard to work with.
    """
    kpoints = np.array([k.frac_coords for k in bandstructure.kpoints])
    labels_dict = {k: v.frac_coords for k, v in bandstructure.labels_dict.items()}

    # pymatgen band structure objects support branches. These are formed when
    # two kpoints with the same label are next to each other. This bit of code
    # will ensure that the band structure will contain branches, if it doesn't
    # already.
    dup_ids = []
    high_sym_kpoints = tuple(map(tuple, labels_dict.values()))
    for i, k in enumerate(kpoints):
        dup_ids.append(i)
        if (
            tuple(k) in high_sym_kpoints
            and i != 0
            and i != len(kpoints) - 1
            and (
                not np.array_equal(kpoints[i + 1], k)
                or not np.array_equal(kpoints[i - 1], k)
            )
        ):
            dup_ids.append(i)
    return dup_ids

def fake_legend(ax:plt.Axes, colours:list, labels:list):
    """
    Overwrites a plot legend with specified colours and labels in the
    style of a sumo bandplot legend.
    """
    # Generate handles
    handles = [
        mlines.Line2D([],[],color=c, marker='o', markersize=np.sqrt(50),
                      markeredgecolor='none', linestyle='none', label=l)
        for c,l in zip(colours, labels)
        ]
    # Plot legend
    ax.legend(
        handles=handles,
        bbox_to_anchor=(0.95, 1),
        loc=2,
        frameon=False,
        handletextpad=0.1,
        scatterpoints=1,
    )

def get_atomic_projected_plot(vasprun:Vasprun, selection:list, kpoints_filename:str=None, **kwargs):
    """
    From Vasprun object, plot bands projected onto selection(s) of atoms.
    
    kwargs are passed to get_plot_with_projections
    """
    bs = vasprun.get_band_structure(kpoints_filename)
    plotter = SBSPlotter(bs)
    proj = get_atomic_projections(plotter.bs, selection, normalise=kwargs.get('normalise', 'all'))
    # Plot
    return get_plot_with_projections(plotter,proj,**kwargs)

def get_spin_projected_plot(vasprun:Vasprun, selection:list, spin:int=0,
                            normalise:str='all', spin_axis:int=2,
                            kpoints_filename:str=None, **kwargs):
    """
    From Vasprun object (noncollinear magnetism), plot bands projected
    onto selection(s) of atoms polarised along the spin_axis in either the
    positive (spin=0) or negative (spin=1) directions.
    
    kwargs are passed to get_plot_with_projections
    """
    proj = get_magnetized_projections(vasprun, selection, normalise, spin_axis)[spin]
    plotter = SBSPlotter(vasprun.get_band_structure(kpoints_filename))
    return get_plot_with_projections(plotter, proj, **kwargs)

def get_two_spin_plot(vasprun:Vasprun, selection:list, normalise:str='all',
                      spin_axis:int=2, kpoints_filename:str=None, **kwargs):
    """
    From Vasprun object (noncollinear magnetism), plot bands projected
    onto a single selection of atoms polarised along the spin_axis,
    plotting both positive and negative spin.
    
    kwargs are passed to get_plot_with_projections
    """
    if isinstance(selection[0], int):
        selection = [selection]
    assert len(selection) == 1, "You gave nested list as selection, but only one selection allowed."
    proj_up, proj_dn = get_magnetized_projections(vasprun, selection, normalise, spin_axis)
    # Each proj_?? is [{spin:projection[band,kpoint]}]
    proj = [proj_up[0], proj_dn[0]]
    # Plot
    plotter = SBSPlotter(vasprun.get_band_structure(kpoints_filename))
    return get_plot_with_projections(plotter, proj, **kwargs)