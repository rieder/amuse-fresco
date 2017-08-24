import numpy as np
from amuse.units import units


def map_to_grid(
        *args,
        **kwargs
        ):
    if len(args) == 2:
        return map_to_2d_grid(*args, **kwargs)
    else:
        return -1


def map_to_2d_grid(
        x,
        y,
        weights=1,
        image_size=(2048, 2048),
        image_width=(10 | units.parsec, 10 | units.parsec),
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
    else:
        return -1


def map_to_2d_grid_simple(
        x,
        y,
        weights=1,
        image_size=(2048, 2048),
        image_width=(10 | units.parsec, 10 | units.parsec),
        ):
    try:
        x_size = image_size[0]
        y_size = image_size[1]
    except:
        x_size = image_size
        y_size = image_size
    try:
        x_width = image_width[0]
        y_width = image_width[1]
    except:
        x_width = image_width
        y_width = image_width
    length_unit = x_width.unit

    xbins = (
            np.linspace(-0.5, 0.5, num=1+x_size)
            * x_width.value_in(length_unit)
            )
    ybins = (
            np.linspace(-0.5, 0.5, num=1+y_size)
            * y_width.value_in(length_unit)
            )
    result = np.histogram2d(
            x.value_in(length_unit),
            y.value_in(length_unit),
            bins=[xbins, ybins],
            weights=weights,
            )

    return result[0]
