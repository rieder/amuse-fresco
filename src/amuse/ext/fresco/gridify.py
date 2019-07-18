# -*- coding: utf-8 -*-
"""
Maps particle positions to a grid, calculates the luminosity in each grid cell.
"""
import numpy as np
from amuse.units import units


def map_to_grid(
        *args,
        **kwargs
):
    if len(args) == 2:
        return map_to_2d_grid(*args, **kwargs)
    return -1


def map_to_2d_grid(
        x,
        y,
        weights=1,
        image_size=(2048, 2048),
        image_width=10 | units.parsec,
        periodic=False,
        mode="simple",
):
    """
    Returns a grid
    """
    if mode == "simple":
        return map_to_2d_grid_simple(
            x,
            y,
            weights=weights,
            image_size=image_size,
            image_width=image_width,
        )
    return -1


def map_to_2d_grid_simple(
        x,
        y,
        weights=1,
        image_size=(2048, 2048),
        image_width=10 | units.parsec,
):
    try:
        x_size = image_size[0]
        y_size = image_size[1]
    except TypeError:
        x_size = image_size
        y_size = image_size
    try:
        x_width = image_width
        y_width = image_width * (y_size / x_size)
    except TypeError:
        x_width = image_width
        y_width = image_width

    # Calculate the X and Y position in pixels
    # width is centred on 0, so shift by 0.5 image size!
    x_px = x_size * ((x / x_width) + 0.5)
    y_px = y_size * ((y / y_width) + 0.5)

    x_0 = np.floor(x_px)
    y_0 = np.floor(y_px)
    x_1 = x_0 + 1
    y_1 = y_0 + 1
    weight_x0 = (x_1 - x_px)
    weight_y0 = (y_1 - y_px)
    weight_x1 = (x_px - x_0)
    weight_y1 = (y_px - y_0)
    pos_weights = [0, 0, 0, 0]
    pos_weights[0] = weight_x0 * weight_y0
    pos_weights[1] = weight_x0 * weight_y1
    pos_weights[2] = weight_x1 * weight_y0
    pos_weights[3] = weight_x1 * weight_y1

    xbins = np.arange(1 + x_size)
    ybins = np.arange(1 + y_size)

    result = np.zeros(x_size * y_size).reshape(x_size, y_size)
    for i in range(2):
        for j in range(2):
            result += np.histogram2d(
                x_0 + i,
                y_0 + j,
                bins=[xbins, ybins],
                weights=weights * pos_weights[2 * i + j],
            )[0]

    return result
