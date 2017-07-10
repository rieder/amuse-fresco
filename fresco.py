# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        division
        )

from amuse.units import units, constants, nbody_system
from amuse.lab import Particles
from amuse.io import read_set_from_file
from amuse.community.fi.interface import FiMap

from ubvinew import rgb_frame

import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse


def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-s',
            dest='starsfilename',
            default='stars.hdf5',
            help='file containing the stars [stars.hdf5]',
            )
    parser.add_argument(
            '-g',
            dest='gasfilename',
            default='',
            help='file containing the gas (optional) []',
            )
    parser.add_argument(
            '-o',
            dest='imagefilename',
            default='test.png',
            help='write image to this file [test.png]',
            )
    parser.add_argument(
            '-b',
            dest='sourcebands',
            default='ubvri',
            help='colour bands to use [ubvri]',
            )
    parser.add_argument(
            '-x',
            dest='plot_axes',
            action='store_true',
            default=False,
            help='plot axes [False]',
            )
    return parser.parse_args()


def calculate_effective_temperature(luminosity, radius):
    temp = (
            (
                luminosity / (constants.four_pi_stefan_boltzmann*radius**2)
                )**.25
            ).in_(units.K)
    return temp


def make_image(
        mode,
        converter,
        gas,
        stars,
        image_width=10. | units.parsec,
        image_size=[1024, 1024],
        percentile=0.9995,
        age=0. | units.Myr,
        sourcebands="ubvri",
        vmax=None,
        calc_temperature=True,
        ):
    def mapper():
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
        mapper.parameters.extinction_flag =\
            True if mode == "stars+gas" else False
        return mapper

    if mode == "gas":
        image = column_density_map(mapper, gas)
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
                mapper_factory=mapper
                )
    return image


def column_density_map(mapper, part):
    if callable(mapper):
        mapper = mapper()

    p = mapper.particles.add_particles(part)
    p.weight = part.mass.value_in(units.amu)
    projected = mapper.image.pixel_value
    # print part.weight.sum()
    # print projected.min(),projected.max(),projected.mean(),projected.sum()
    mapper.stop()
    im = np.transpose(projected)
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
        mapper_factory=None,
        ):
    if calc_temperature:
        # calculates the temperature of the stars from their total luminosity
        # and radius, calculates those first if needed
        try:
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
        except:
            print(
                    "Calculating luminosity/temperature for %s old stars..."
                    % (age)
                    )
            from amuse.community.sse.interface import SSE
            se = SSE()
            se.particles.add_particles(stars)
            if age > 0 | units.Myr:
                se.evolve_model(age)
            stars.luminosity = se.particles.luminosity
            stars.radius = se.particles.radius
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
            se.stop()

    vmax, rgb = rgb_frame(
            stars,
            dryrun=False,
            image_width=image_width,
            multi_psf=False,  # True,
            image_size=image_size,
            percentile=percentile,
            sourcebands=sourcebands,
            mapper_factory=mapper_factory,
            gas=gas,
            )
    return rgb['pixels']


if __name__ == "__main__":
    args = new_argument_parser()
    starsfilename = args.starsfilename
    gasfilename = args.gasfilename
    imagefilename = args.imagefilename

    plot_axes = args.plot_axes
    sourcebands = args.sourcebands

    length_unit = units.parsec
    dpi = 600
    image_width_arcsec = 160

    # FIXME this should depend on the distance!
    # size = fixed at nr of arcmin
    image_width = 10. | units.parsec

    # FIXME the psf is fixed pixel size, so the pixels in the image here
    # reflects how much pixels will be spread!  In principle, this should be a
    # fixed number, equal to the number of pixels of the ccd for which the psf
    # was made.  Also, the image should reflect the distance to the "observed"
    # cluster.
    image_size = [2048, 2048]
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

    age = 500. | units.Myr
    percentile = 0.9995  # for determining vmax

    xmin = -0.5 * image_width.value_in(length_unit)
    xmax = 0.5 * image_width.value_in(length_unit)
    ymin = -0.5 * image_width.value_in(length_unit)
    ymax = 0.5 * image_width.value_in(length_unit)

    stars = read_set_from_file(
            starsfilename,
            "amuse",
            close_file=True,
            )
    if gasfilename:
        gas = read_set_from_file(
                gasfilename,
                "amuse",
                close_file=True,
                )
    else:
        gas = Particles()
    gas.h_smooth = 0.02 | units.parsec

    # FIXME: add these features
    # - Rotate so that xy = observed x/y axes of figure
    # - Scale positions to desired ra/dec (script Alison)
    # - calculate vmax based on nr of photons/exposure time

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, dpi=dpi)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel("[%s]" % (length_unit))
    ax.set_ylabel("[%s]" % (length_unit))
    ax.set_aspect(1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    converter = nbody_system.nbody_to_si(
            stars.total_mass(),
            image_width,
            )

    image = make_image(
            "gas_and_stars",
            converter,
            gas,
            stars,
            image_width=image_width,
            image_size=image_size,
            percentile=percentile,
            calc_temperature=True,
            age=age,
            sourcebands=sourcebands,
            )

    plt.imshow(
            image,
            origin='lower',
            extent=[
                xmin,
                xmax,
                ymin,
                ymax,
                ],
            )

    plt.savefig(
            imagefilename,
            dpi=dpi,
            )
