#!/usr/bin/env python
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

import argparse

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.datamodel.rotation import rotate

from amuse.ext.fresco.fieldstars import new_field_stars
from amuse.ext.fresco.fresco import (
    evolve_to_age,
    make_image,
    column_density_map,
    initialise_image,
)


def new_argument_parser():
    "Parse command line arguments"
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
        '-f',
        dest='followfilename',
        default=None,
        help=\
        'file containing star keys to center on (optional) []\n'
        '  implies --com',
    )
    parser.add_argument(
        '-o',
        dest='imagefilename',
        default=None,
        help='write image to this file [test]',
    )
    parser.add_argument(
        '--imagetype',
        dest='imagetype',
        default='png',
        help='image file type [png]',
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
        '--ext',
        dest='calculate_extinction',
        action='store_true',
        default=False,
        help='include extinction by dust [False]',
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
    parser.add_argument(
        '--ax',
        dest='angle_x',
        default=0,
        type=float,
        help='Rotation step around x-axis in deg [0]',
    )
    parser.add_argument(
        '--ay',
        dest='angle_y',
        default=0,
        type=float,
        help='Rotation step around y-axis in deg [0]',
    )
    parser.add_argument(
        '--az',
        dest='angle_z',
        default=0,
        type=float,
        help='Rotation step around z-axis in deg [0]',
    )
    parser.add_argument(
        '--frames',
        dest='frames',
        default=1,
        type=int,
        help='Number of frames (>1: rotate around x,y,z) [1]',
    )
    parser.add_argument(
        '--px',
        dest='pixels',
        default=2048,
        type=int,
        help='Number of pixels along each axis [2048]',
    )
    parser.add_argument(
        '--psf',
        dest='psf_type',
        default='hubble',
        help='PSF type (valid: [hubble], gaussian)',
    )
    parser.add_argument(
        '--sigma',
        dest='psf_sigma',
        default=1.0,
        type=float,
        help='PSF sigma (if PSF type is gaussian)',
    )
    parser.add_argument(
        '--contours',
        dest='contours',
        action='store_true',
        default=False,
        help='Plot gas contour lines [False]',
    )
    parser.add_argument(
        '--com',
        dest='use_com',
        action='store_true',
        default=False,
        help='Center on center of mass [False]',
    )
    parser.add_argument(
        '--xo',
        dest='x_offset',
        default=0.0,
        type=float,
        help='X offset (in parsec)',
    )
    parser.add_argument(
        '--yo',
        dest='y_offset',
        default=0.0,
        type=float,
        help='Y offset (in parsec)',
    )
    parser.add_argument(
        '--zo',
        dest='z_offset',
        default=0.0,
        type=float,
        help='Z offset (in parsec)',
    )
    return parser.parse_args()


def main():
    # Fixed settings
    stellar_evolution = True
    se_code = "SeBa"
    length_unit = units.parsec
    dpi = 600
    percentile = 0.9995  # for determining vmax

    # Parse arguments
    args = new_argument_parser()
    starsfilename = args.starsfilename
    gasfilename = args.gasfilename
    followfilename = args.followfilename
    imagefilename = args.imagefilename
    imagetype = args.imagetype
    vmax = args.vmax if args.vmax > 0 else None
    n_fieldstars = args.n_fieldstars
    filetype = args.filetype
    contours = args.contours
    np.random.seed(args.seed)
    plot_axes = args.plot_axes
    angle_x = args.angle_x | units.deg
    angle_y = args.angle_y | units.deg
    angle_z = args.angle_z | units.deg
    sourcebands = args.sourcebands
    psf_type = args.psf_type.lower()
    psf_sigma = args.psf_sigma
    age = args.age | units.Myr
    image_width = args.width | units.parsec
    pixels = args.pixels
    frames = args.frames
    if followfilename is not None:
        use_com = True
    else:
        use_com = args.use_com
    x_offset = args.x_offset | units.parsec
    y_offset = args.y_offset | units.parsec
    z_offset = args.z_offset | units.parsec
    extinction = args.calculate_extinction

    # Derived settings
    if psf_type not in ["hubble", "gaussian"]:
        print(("Invalid PSF type: %s" % psf_type))
        exit()
    image_size = [pixels, pixels]
    # If the nr of pixels is changed, zoom the PSF accordingly.
    zoom_factor = pixels / 2048.

    if starsfilename:
        stars = read_set_from_file(
            starsfilename,
            filetype,
            close_file=True,
        )
        if stellar_evolution and (age > 0 | units.Myr):
            print((
                "Calculating luminosity/temperature for %s old stars..."
                % (age)
            ))
            evolve_to_age(stars, age, stellar_evolution=se_code)
        if use_com:
            if followfilename is not None:
                followstars = read_set_from_file(
                    followfilename, filetype, close_file=True,
                )
                center_on_these_stars = followstars.get_intersecting_subset_in(
                    stars,
                )
            else:
                center_on_these_stars = stars
            com = center_on_these_stars.center_of_mass()
            x_offset, y_offset, z_offset = com
        stars.x -= x_offset
        stars.y -= y_offset
        stars.z -= z_offset
    else:
        stars = Particles()

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
        evolve_to_age(fieldstars, 0 | units.yr, stellar_evolution=se_code)
        stars.add_particles(fieldstars)

    if gasfilename:
        gas = read_set_from_file(
            gasfilename,
            filetype,
            close_file=True,
        )
        if use_com:
            if stars.is_empty():
                com = gas.center_of_mass()
                x_offset, y_offset, z_offset = com
        gas.x -= x_offset
        gas.y -= y_offset
        gas.z -= z_offset
        # Gadget and Fi disagree on the definition of h_smooth.
        # For gadget, need to divide by 2 to get the Fi value (??)
        gas.h_smooth *= 0.5
        gas.radius = gas.h_smooth

        # Select only the relevant gas particles (plus a margin)
        minx = (1.1 * -image_width/2)
        maxx = (1.1 * image_width/2)
        miny = (1.1 * -image_width/2)
        maxy = (1.1 * image_width/2)
        gas_ = gas.select(
            lambda x, y:
            x > minx
            and x < maxx
            and y > miny
            and y < maxy,
            ["x", "y"]
        )
        gas = gas_
    else:
        gas = Particles()
    # gas.h_smooth = 0.05 | units.parsec

    converter = nbody_system.nbody_to_si(
        stars.total_mass() if not stars.is_empty() else gas.total_mass(),
        image_width,
    )

    # Initialise figure and axes
    fig = initialise_image(
        dpi=dpi,
        image_size=image_size,
        length_unit=length_unit,
        image_width=image_width,
        plot_axes=plot_axes,
        x_offset=x_offset,
        y_offset=y_offset,
        z_offset=z_offset,
    )
    ax = fig.get_axes()[0]
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    for frame in range(frames):
        fig = initialise_image(fig)

        if (frame != 0) or (frames == 1):
            if not stars.is_empty():
                rotate(stars, angle_x, angle_y, angle_z)
            if not gas.is_empty():
                rotate(gas, angle_x, angle_y, angle_z)

        image, vmax = make_image(
            stars=stars if not stars.is_empty() else None,
            gas=gas if not gas.is_empty() else None,
            converter=converter,
            image_width=image_width,
            image_size=image_size,
            percentile=percentile,
            calc_temperature=True,
            age=age,
            vmax=vmax,
            sourcebands=sourcebands,
            zoom_factor=zoom_factor,
            psf_type=psf_type,
            psf_sigma=psf_sigma,
            return_vmax=True,
            extinction=extinction,
        )

        if not stars.is_empty():
            ax.imshow(
                image,
                origin='lower',
                extent=[
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                ],
            )
            if contours and not gas.is_empty():
                gascontours = column_density_map(
                    gas,
                    zoom_factor=zoom_factor,
                    image_width=image_width,
                    image_size=image_size,
                )
                gascontours[np.isnan(gascontours)] = 0.0
                vmax = np.max(gascontours) / 2
                # vmin = np.min(image[np.where(image > 0.0)])
                vmin = vmax / 100
                levels = 10**(
                    np.linspace(
                        np.log10(vmin),
                        np.log10(vmax),
                        num=5,
                    )
                )[1:]
                # print(vmin, vmax)
                # print(levels)
                ax.contour(
                    origin='lower',
                    levels=levels,
                    colors="white",
                    linewidths=0.1,
                    extent=[
                        xmin,
                        xmax,
                        ymin,
                        ymax,
                    ],
                )
        else:
            image = column_density_map(
                gas,
                image_width=image_width,
                image_size=image_size,
            )

            ax.imshow(
                image,
                origin='lower',
                extent=[
                    xmin,
                    xmax,
                    ymin,
                    ymax,
                ],
                cmap="gray",
            )

        if frames > 1:
            savefilename = "%s-%06i.%s" % (
                imagefilename if imagefilename is not None else "test",
                frame,
                imagetype,
            )
        else:
            savefilename = "%s.%s" % (
                imagefilename if imagefilename is not None else "test",
                imagetype,
            )
        plt.savefig(
            savefilename,
            dpi=dpi,
        )


if __name__ == "__main__":
    main()
    print(
        "------------------\n"
        "\n"
        "Fresco is built on the tools below."
        " If you use Fresco for your publication, please cite the following"
        " references:\n"
        "\n"
        "Fresco:"
        " Steven Rieder & Inti Pelupessy."
        " rieder/Fresco (2019)."
        " doi:10.5281/zenodo.3362342\n"
        "AMUSE:"
        " Simon Portegies Zwart & Steve McMillan."
        " Astrophysical Recipes; The art of AMUSE (2019)."
        " ADS bibcode: 2018araa.book.....P."
        " doi:10.1088/978-0-7503-1320-9\n"
        "Matplotlib:"
        " John D. Hunter."
        " Matplotlib: A 2D Graphics Environment."
        " doi:10.1109/MCSE.2007.55\n"
        "Numpy:"
        " Stéfan van der Walt, S. Chris Colbert and Gaël Varoquaux."
        " The NumPy Array: A Structure for Efficient Numerical Computation."
        " doi:10.1109/MCSE.2011.37\n"
        "Python:"
        " Guido van Rossum."
        " Extending and Embedding the Python Interpreter."
        " May 1995. CWI Report CS-R9527.\n"
    )
