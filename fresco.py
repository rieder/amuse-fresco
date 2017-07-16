# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        division
        )

from amuse.units import units, constants, nbody_system
from amuse.lab import Particles
from amuse.io import read_set_from_file
from amuse.community.fi.interface import FiMap

from fresco.ubvinew import rgb_frame
from fresco.fieldstars import new_field_stars

import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse


def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--filetype',
            dest='filetype',
            default='amuse',
            help='filetype [amuse], valid are amuse,starlab,txt,...',
            )
    parser.add_argument(
            '-s',
            dest='starsfilename',
            default='',
            help='file containing stars (optional) []',
            )
    parser.add_argument(
            '-g',
            dest='gasfilename',
            default='',
            help='file containing gas (optional) []',
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
            '-a',
            dest='age',
            default=100.,
            type=float,
            help='age of the stars in Myr [100]',
            )
    parser.add_argument(
            '-w',
            dest='width',
            default=5.,
            type=float,
            help='image width in parsec [5]',
            )
    parser.add_argument(
            '-x',
            dest='plot_axes',
            action='store_true',
            default=False,
            help='plot axes [False]',
            )
    parser.add_argument(
            '--seed',
            dest='seed',
            default=1701,
            type=int,
            help='random seed',
            )
    parser.add_argument(
            '--vmax',
            dest='vmax',
            default=0,
            type=float,
            help='vmax value',
            )
    parser.add_argument(
            '--field',
            dest='n_fieldstars',
            default=0,
            type=int,
            help='add N field stars (optional) [0]',
            )
    return parser.parse_args()


def evolve_to_age(stars, age, se="SeBa"):
    if se == "SeBa":
        from amuse.community.seba.interface import SeBa
        se = SeBa()
    elif se == "SSE":
        from amuse.community.sse.interface import SSE
        se = SSE()
        # SSE can result in nan values for luminosity/radius
    se.particles.add_particles(stars)
    if age > 0 | units.yr:
        se.evolve_model(age)
    stars.luminosity = np.nan_to_num(
            se.particles.luminosity.value_in(units.LSun)
            ) | units.LSun
    # Temp fix: add one meter to radius of stars, to prevent zero/nan radius.
    # TODO: Should fix this a better way, but it's ok for now.
    stars.radius = (1 | units.m) + (np.nan_to_num(
            se.particles.radius.value_in(units.RSun)
            ) | units.RSun)
    se.stop()
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
                vmax=vmax,
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
        vmax=None,
        mapper_factory=None,
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
            )
    return rgb['pixels']


if __name__ == "__main__":
    args = new_argument_parser()
    starsfilename = args.starsfilename
    gasfilename = args.gasfilename
    imagefilename = args.imagefilename
    stellar_evolution = True
    vmax = args.vmax if args.vmax > 0 else None
    n_fieldstars = args.n_fieldstars
    filetype = args.filetype
    np.random.seed(args.seed)
    se_code = "SSE"

    plot_axes = args.plot_axes
    sourcebands = args.sourcebands
    age = args.age | units.Myr
    # FIXME this should depend on the distance!
    # size = fixed at nr of arcmin
    image_width = args.width | units.parsec

    length_unit = units.parsec
    dpi = 600
    image_width_arcsec = 160

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

    percentile = 0.9995  # for determining vmax

    xmin = -0.5 * image_width.value_in(length_unit)
    xmax = 0.5 * image_width.value_in(length_unit)
    ymin = -0.5 * image_width.value_in(length_unit)
    ymax = 0.5 * image_width.value_in(length_unit)

    if starsfilename:
        stars = read_set_from_file(
                starsfilename,
                filetype,
                close_file=True,
                )
        if stellar_evolution:
            print(
                    "Calculating luminosity/temperature for %s old stars..."
                    % (age)
                    )
            evolve_to_age(stars, age, se=se_code)
        com = stars.center_of_mass()
        stars.position -= com
    else:
        stars = Particles()
        com = np.array([0, 0, 0]) | units.parsec

    if n_fieldstars:
        minage = 400 | units.Myr
        maxage = 12 | units.Gyr
        fieldstars = new_field_stars(
                n_fieldstars,
                width=image_width,
                height=image_width,
                )
        fieldstars.age = (
                minage
                + (
                    np.random.sample(n_fieldstars)
                    * (maxage - minage)
                    )
                )
        evolve_to_age(fieldstars, 0 | units.yr, se=se_code)
        # TODO: add distance modulus
        stars.add_particles(fieldstars)

    if gasfilename:
        gas = read_set_from_file(
                gasfilename,
                filetype,
                close_file=True,
                )
        gas.position -= com
    else:
        gas = Particles()
    gas.h_smooth = 0.05 | units.parsec

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
            "stars+gas",
            converter,
            gas,
            stars,
            image_width=image_width,
            image_size=image_size,
            percentile=percentile,
            calc_temperature=True,
            age=age,
            vmax=vmax,
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
