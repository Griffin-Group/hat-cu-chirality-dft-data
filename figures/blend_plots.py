#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Commutative blending of figures and images.
Useful for plotting projected bands where bands with different projections
are degenerate. Conventional plotting has one overlap the other, giving
misleading results dependent on plotting order. But if you use additive
or multiplicative blending instead of alpha or overlay blending, then
the results are commutative.

See test_plots.py for an example of usage, because it's somewhat involved.

Copyright Bernard Field, 2025
"""

import matplotlib.pyplot as plt

import numpy as np

from PIL import Image


def blend_figures(*figs:plt.Figure, mode:str='additive') -> np.array:
    """
    Takes several figures of the same size, rasterizes them, then blends the images
    
    Based on https://stackoverflow.com/a/26712790/
    
    Note that for many blending modes you must set the background to be 
    transparent (i.e. ax.patch.set_facecolor('none'),
    ax.patch.set_edgecolor('none')).
    
    Parameters
    ----------
    *figs : Multiple Figures
    mode : str, optional
        Blending mode. The default is 'additive'.

    Returns
    -------
    Numpy array representing the figure.

    """
    buffers = []
    for fig in figs:
        # Force the figure to draw
        fig.canvas.draw()
        # Need dimensions because the buffer is just 1D.
        w, h = fig.canvas.get_width_height()
        # Convert the figure to an RGBA array.
        img = np.frombuffer(fig.canvas.buffer_rgba(), np.uint8).reshape(h,w,-1).copy()
        buffers.append(img)
        assert img.shape == buffers[0].shape, "Mismatched image sizes detected!"
    return blend_arrays(*buffers, mode=mode)

def blend_arrays(*buffers:np.ndarray, mode:str='additive') -> np.ndarray:
    """
    Takes 2 or more np.arrays(dtype=np.uint8) that are RGBA.
    Blends them together.
    The arrays must be the same size.
    
    The arrays get mutated.

    Parameters
    ----------
    *buffers : np.array
        2 or more np.arrays, shape (h,w,4), RGBA images, dtype uint8.
    mode : str, optional
        Blending mode. The default is 'additive'.

    Returns
    -------
    img : np.array representing the figure.

    """
    # Do blending
    if mode == 'additive':
        # Set transparent pixels to transparent black.
        # Premultiply by alpha so each pixel is weighted by its opacity.
        for img in buffers:
            tmp = np.asarray(img, dtype=np.int64)
            tmp[:,:,0:3] = tmp[:,:,0:3] * tmp[:,:,[-1]] / 255
            img[:] = tmp
        # Additive adds together all the buffers and maxes at 255.
        img = np.zeros(buffers[0].shape, np.uint16)
        for b in buffers:
            img += b
        img = np.minimum(img, 255).astype(np.uint8)
    elif mode == "multiply":
        # Set transparent pixels to solid white, weighted by alpha.
        for i, img in enumerate(buffers):
            img = np.asarray(img, dtype=np.uint64)
            img = (img * img[:,:,[-1]] + np.full(img.shape, 255, np.uint64) * (255-img[:,:,[-1]])) / 255
            img[:,:,-1] = 255
            buffers[i][:] = img
        img = buffers[0].copy().astype(np.uint64)
        for b in buffers[1:]:
            img = np.asarray(img * b / 255, np.uint64)
        img = np.asarray(img, np.uint8)
    else:
        raise ValueError(f"Don't recognise blending mode {mode}.")
    return img

def raw_imshow(img:np.array, show:bool=True):
    """
    Plot an image array as just the image, no frills.
    Calls plt.show() if show=True.
    """
    figsize = (6.4, 6.4*img.shape[0]/img.shape[1])
    fig, ax = plt.subplots(figsize=figsize)
    plt.imshow(img)
    plt.subplots_adjust(0,0,1,1)
    plt.axis("off")
    if show:
        plt.show()

def save_img(filename:str, img:np.array, format:str=None, dpi:int=None):
    """
    Save an array as an image.
    """
    if isinstance(dpi, int) or isinstance(dpi, float):
        dpi = [dpi]*2 # DPI should be array of two numbers
    Image.fromarray(img).save(filename, format=format, dpi=dpi)