# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        )
import os

import numpy

# from numpy.fft import fft2, ifft2
# from numpy import log

from amuse.units import units, nbody_system
from amuse.community.fi.interface import FiMap

# import logging
# logging.basicConfig(level=logging.DEBUG)

from filters3 import filter_band_flux, filter_data
from xyz import xyz_data
from blackbody import B_lambda
from color_converter import (
        ColorConverter,
        XYZ_to_sRGB_linear, sRGB_linear_to_sRGB
        )

# import time

from astropy.convolution import (
        # convolve,
        convolve_fft,
        )

import astropy.io.fits as pyfits


def Convolve(
        image,
        kernel,
        ):
    result = convolve_fft(image, kernel, boundary='fill')
    return result


psf = dict()

for band in "ubvri":
    for i in range(4):
        instrument = "WFPC_II_WFC3"
        this_dir = os.path.split(__file__)[0]
        if this_dir == "":
            this_dir = "."
        datadir = this_dir + "/data/" + instrument + "/"
        # FIXME use some relative dir
        f = pyfits.open(datadir+band+"%2.2i.fits" % i)
        psf[band+str(i)] = numpy.array(f[0].data)


def rgb_frame(
        stars,
        dryrun=False,
        vmax=None,
        percentile=0.9995,
        multi_psf=False,
        sourcebands="ubvri",
        image_width=12. | units.parsec,
        image_size=[1024, 1024],
        ):

    print("luminosities..")

    for band in sourcebands:
        setattr(
                stars,
                band+"_band",
                4 * numpy.pi * stars.radius**2 *
                filter_band_flux(
                    "bess-" + band + ".pass",
                    lambda x: B_lambda(x, stars.temperature),
                    ),
                )

    print("..raw images..")

    conv = nbody_system.nbody_to_si(stars.total_mass(), image_width)
    mapper = FiMap(conv, mode="openmp", redirection="none")

    mapper.parameters.image_width = image_width
    mapper.parameters.image_size = image_size

    mapper.particles.add_particles(stars)
    # mapper.parameters.viewpoint_x   = 0|units.parsec
    # mapper.parameters.viewpoint_y   = 0|units.parsec
    # mapper.parameters.viewpoint_z   = 10|units.parsec
    mapper.parameters.projection_direction = [0, 0, 1]
    mapper.parameters.upvector = [0, -1, 0]

    raw_images = dict()
    for band in sourcebands:
        mapper.particles.weight = getattr(
                stars,
                band+"_band"
                ).value_in(units.LSun)
        im = mapper.image.pixel_value
        raw_images[band] = im

    mapper.stop()

    convolved_images = dict()

    print("..convolving..")

    if multi_psf:
        a = numpy.arange(image_size[0])/float(image_size[0]-1)
        b = numpy.arange(image_size[1])/float(image_size[1]-1)
        w1 = numpy.outer(a, b)
        w2 = numpy.outer(1.-a, b)
        w3 = numpy.outer(a, 1.-b)
        w4 = numpy.outer(1.-a, 1.-b)
        for key, val in raw_images.items():
            im1 = Convolve(val, psf[key+'0'])
            im2 = Convolve(val, psf[key+'1'])
            im3 = Convolve(val, psf[key+'2'])
            im4 = Convolve(val, psf[key+'3'])
            convolved_images[key] = w1*im1+w2*im2+w3*im3+w4*im4
    else:
        for key, val in raw_images.items():
            im1 = Convolve(val, psf[key+'0'])
            convolved_images[key] = im1

    print("..conversion to rgb")

    source = [filter_data['bess-'+x+'.pass'] for x in sourcebands]

    target = [xyz_data['x'], xyz_data['y'], xyz_data['z']]

    conv = ColorConverter(source, target)

    ubv = numpy.array([convolved_images[x] for x in sourcebands])

    xyz = numpy.tensordot(conv.conversion_matrix, ubv, axes=(1, 0))

    conv_xyz_to_lin = XYZ_to_sRGB_linear()

    srgb_l = numpy.tensordot(
            conv_xyz_to_lin.conversion_matrix,
            xyz,
            axes=(1, 0),
            )

    if dryrun or vmax is None:
        flat_sorted = numpy.sort(srgb_l.flatten())
        n = len(flat_sorted)
        vmax = flat_sorted[int(1.-3*(1.-percentile)*n)]
        print("vmax:", vmax)
    if dryrun:
        return vmax

    conv_lin_to_sRGB = sRGB_linear_to_sRGB()

    srgb = conv_lin_to_sRGB.convert(srgb_l/vmax)

    # r = srgb[0,:,:].transpose()
    # g = srgb[1,:,:].transpose()
    # b = srgb[2,:,:].transpose()
    r = numpy.fliplr(srgb[0, :, :])
    g = numpy.fliplr(srgb[1, :, :])
    b = numpy.fliplr(srgb[2, :, :])

    rgb = numpy.zeros(image_size[0] * image_size[1] * 3)
    rgb = rgb.reshape(image_size[0], image_size[1], 3)
    rgb[:, :, 0] = r
    rgb[:, :, 1] = g
    rgb[:, :, 2] = b

    image = dict(
        pixels=rgb,
        size=rgb.size,
        )

    return vmax, image
