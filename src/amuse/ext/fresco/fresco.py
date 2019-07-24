# -*- coding: utf-8 -*-
"""
Fresco creates a "simulated observation" of a set of particles.
Particles can be "stars" (point sources emitting light) or "gas" (emitting,
reflecting and/or obscuring light). Gas may also be displayed with contour
lines.
"""

from __future__ import (
    print_function,
    division,
    absolute_import,
)

import numpy as np

from scipy.ndimage import gaussian_filter

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.datamodel.rotation import rotate

import matplotlib.pyplot as plt

from amuse.ext.fresco.ubvi import rgb_frame
from amuse.ext.fresco.fieldstars import new_field_stars


def evolve_to_age(stars, age, stellar_evolution="SeBa"):
    "Evolve stars to specified age with specified code"
    if stellar_evolution == "SeBa":
        from amuse.community.seba.interface import SeBa
        stellar_evolution = SeBa()
    elif stellar_evolution == "SSE":
        from amuse.community.sse.interface import SSE
        stellar_evolution = SSE()
        # SSE can result in nan values for luminosity/radius
    else:
        raise "No such stellar evolution code %s or no code specified" % (
            stellar_evolution
        )
    stellar_evolution.particles.add_particles(stars)
    if age > 0 | units.yr:
        stellar_evolution.evolve_model(age)
    stars.luminosity = np.nan_to_num(
        stellar_evolution.particles.luminosity.value_in(units.LSun)
    ) | units.LSun

    stars.radius = stellar_evolution.particles.radius
    # prevent zero/nan radius.
    x = np.where(
        np.nan_to_num(
            stars.radius.value_in(units.RSun)
        ) == 0.
    )
    stars[x].radius = 0.01 | units.RSun

    stellar_evolution.stop()
    return


def calculate_effective_temperature(luminosity, radius):
    temp = np.nan_to_num(
        (
            (
                luminosity
                / (
                    constants.four_pi_stefan_boltzmann
                    * radius**2
                )
            )**.25
        ).value_in(units.K)
    ) | units.K
    return temp


def make_image(
        stars=None,
        gas=None,
        converter=None,
        image_width=[
            10. | units.parsec,
            10. | units.parsec,
        ],
        image_size=[1024, 1024],
        percentile=0.9995,
        age=0. | units.Myr,
        sourcebands="ubvri",
        vmax=None,
        calc_temperature=True,
        mapper_code=None,  # "FiMap"
        zoom_factor=1.0,
        psf_type="hubble",
        psf_sigma=1.0,
        extinction=False,
        return_vmax=False,
):
    """
    Makes image from gas and stars
    """
    mode=[]
    if gas is not None:
        mode.append("gas")
    if stars is not None:
        mode.append("stars")
    if mode == []:
        return

    if extinction:
        # Extinction can currently only be handled with FiMap
        mapper_code = "FiMap"

    if mapper_code == "FiMap":
        def mapper():
            from amuse.community.fi.interface import FiMap
            mapper = FiMap(converter, mode="openmp")

            # mapper.parameters.minimum_distance = 1. | units.AU
            mapper.parameters.image_size = image_size
            # mapper.parameters.image_target = image_target

            mapper.parameters.image_width = image_width
            # mapper.parameters.projection_direction = (
            #         (image_target-viewpoint)
            #         / (image_target-viewpoint).length()
            #         )
            # mapper.parameters.projection_mode = projection
            # mapper.parameters.image_angle = horizontal_angle
            # mapper.parameters.viewpoint = viewpoint
            mapper.parameters.extinction_flag = extinction
            return mapper
    else:
        # Gridify as default
        mapper = None
        mapper_code = "gridify"

    if "stars" not in mode:
        image = column_density_map(
            gas,
            image_width=image_width,
            image_size=image_size,
            mapper_factory=mapper,
            mapper_code=mapper_code,
            zoom_factor=zoom_factor,
            psf_type=psf_type,
            psf_sigma=psf_sigma,
            return_vmax=return_vmax,
        )
    else:
        image = image_from_stars(
            stars,
            image_width=image_width,
            image_size=image_size,
            percentile=percentile,
            calc_temperature=calc_temperature,
            age=age,
            sourcebands=sourcebands,
            gas=gas,
            vmax=vmax,
            mapper_factory=mapper,
            mapper_code=mapper_code,
            zoom_factor=zoom_factor,
            psf_type=psf_type,
            psf_sigma=psf_sigma,
            return_vmax=return_vmax,
        )
    return image


def column_density_map(
        gas,
        image_width=10. | units.parsec,
        image_size=[1024, 1024],
        mapper_factory=None,
        mapper_code=None,
        zoom_factor=1.0,
        psf_type="gaussian",
        psf_sigma=10.0,
        return_vmax=False,
):
    if mapper_code == "FiMap":
        if callable(mapper_factory):
            mapper = mapper_factory()

        p = mapper.particles.add_particles(gas)
        p.weight = gas.mass.value_in(units.amu)
        projected = mapper.image.pixel_value
        mapper.stop()
        im = gaussian_filter(
            projected,
            sigma=psf_sigma * zoom_factor,
            order=0,
        )
    else:
        from amuse.ext.fresco.gridify import map_to_grid
        gas_in_mapper = gas.copy()
        gas_in_mapper.weight = gas_in_mapper.mass.value_in(units.amu)
        raw_image = map_to_grid(
            gas_in_mapper.x,
            gas_in_mapper.y,
            weights=gas_in_mapper.weight,
            image_size=image_size,
            image_width=image_width,
        )
        im = gaussian_filter(
            raw_image,
            sigma=psf_sigma * zoom_factor,
            order=0,
        ).T
    if return_vmax:
        return (im, -1)
    return im


def image_from_stars(
        stars,
        image_width=10. | units.parsec,
        image_size=[1024, 1024],
        percentile=0.9995,
        calc_temperature=True,
        age=0. | units.Myr,
        sourcebands="ubvri",
        gas=None,
        vmax=None,
        mapper_factory=None,
        mapper_code=None,
        zoom_factor=1.0,
        psf_type="hubble",
        psf_sigma=1.0,
        return_vmax=False,
):
    if calc_temperature:
        # calculates the temperature of the stars from their total luminosity
        # and radius, calculates those first if needed
        stars.temperature = calculate_effective_temperature(
            stars.luminosity,
            stars.radius,
        )

    vmax, rgb = rgb_frame(
        stars,
        dryrun=False,
        image_width=image_width,
        vmax=vmax,
        multi_psf=False,  # True,
        image_size=image_size,
        percentile=percentile,
        sourcebands=sourcebands,
        mapper_factory=mapper_factory,
        gas=gas,
        mapper_code=mapper_code,
        zoom_factor=zoom_factor,
        psf_type=psf_type,
        psf_sigma=psf_sigma,
    )
    if return_vmax:
        return rgb['pixels'], vmax
    return rgb['pixels']


def initialise_image(
        fig=None,
        dpi=150,
        image_size=[2048, 2048],
        length_unit=units.parsec,
        image_width=5 | units.parsec,
        plot_axes=True,
        subplot=0,
        x_offset=0 | units.parsec,
        y_offset=0 | units.parsec,
        z_offset=0 | units.parsec,
):
    if fig is None:
        if plot_axes:
            left = 0.15
            bottom = 0.15
        else:
            left = 0.
            bottom = 0.
        right = 1.0
        top = 1.0
        figwidth = image_size[0] / dpi / (right - left)
        figheight = image_size[1] / dpi / (top - bottom)
        figsize = (figwidth, figheight)

        xmin = x_offset.value_in(length_unit) - 0.5 * image_width.value_in(length_unit)
        xmax = x_offset.value_in(length_unit) + 0.5 * image_width.value_in(length_unit)
        ymin = y_offset.value_in(length_unit) - 0.5 * image_width.value_in(length_unit)
        ymax = y_offset.value_in(length_unit) + 0.5 * image_width.value_in(length_unit)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, dpi=dpi)
        fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
    else:
        # Simply clear and re-use the old figure
        ax = fig.get_axes()[subplot]
        ax.cla()
    ax.set_xlabel("X (%s)" % (length_unit))
    ax.set_ylabel("Y (%s)" % (length_unit))
    ax.set_aspect(1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_facecolor('black')
    return fig
