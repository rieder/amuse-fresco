# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        )
import os

import numpy

# from numpy.fft import fft2, ifft2
# from numpy import log

from amuse.lab import (
        Particles,
        units,
        nbody_system,
        )
from amuse.community.fi.interface import FiMap

# import logging
# logging.basicConfig(level=logging.DEBUG)

from filters3 import filter_band_flux, filter_data, filter_band_lambda
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


def assign_weights_and_opacities(
        band,
        mapper_stars,
        mapper_gas,
        stars,
        gas,
        dust_to_gas=0.01,
        dust_extinction_cross_section_v_band=4.9e-22 | units.cm**2,
        dust_albedo_v_band=.01,
        Nstar=50,
        ):

    mapper_stars.weight = getattr(stars, band+"_band").value_in(units.LSun)

    if len(gas) == 0:
        return

    lambda_eff = filter_band_lambda("bess-"+band+".pass")
    lambda_v = filter_band_lambda("bess-v.pass")

    f_H = 0.7
    dust_cross_section = (
            dust_extinction_cross_section_v_band
            * (lambda_eff/lambda_v)**-1.5
            )
    dust_albedo = dust_albedo_v_band*(lambda_eff/lambda_v)**0.

    mapper_gas.opacity_area = (
            f_H
            * gas.mass.value_in(units.amu)
            * dust_cross_section
            )

    weight = numpy.zeros(len(mapper_gas))
    stars_ordered_by_lum = stars.sorted_by_attribute(band+"_band")
    for star in stars_ordered_by_lum[-Nstar:]:
        d2 = (
                (gas.x-star.x)**2
                + (gas.y-star.y)**2
                + (gas.z-star.z)**2
                + 0.25 * gas.h_smooth**2
                )
        flux = getattr(star, band+"_band") / (4*numpy.pi*d2)
        # flux = getattr(star, "v_band") / (4*numpy.pi*d2)
        weight += (
                flux
                * mapper_gas.opacity_area
                * dust_albedo
                ).value_in(units.LSun)
    mapper_gas.weight = weight
    return


def rgb_frame(
        stars,
        dryrun=False,
        vmax=None,
        percentile=0.9995,
        multi_psf=False,
        sourcebands="ubvri",
        image_width=12. | units.parsec,
        image_size=[1024, 1024],
        mapper_factory=None,
        gas=None,
        ):
    if gas is None:
        gas = Particles()

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

    if mapper_factory:
        mapper = mapper_factory()
    else:
        conv = nbody_system.nbody_to_si(stars.total_mass(), image_width)
        mapper = FiMap(conv, mode="openmp", redirection="none")
        mapper.parameters.image_width = image_width
        mapper.parameters.image_size = image_size

    stars_in_mapper = mapper.particles.add_particles(stars)
    gas_in_mapper = mapper.particles.add_particles(gas)

    mapper.parameters.projection_direction = [0, 0, 1]
    mapper.parameters.upvector = [0, -1, 0]

    raw_images = dict()
    for band in sourcebands:
        assign_weights_and_opacities(
                band,
                stars_in_mapper,
                gas_in_mapper,
                stars,
                gas,
                Nstar=50,
                )
        # mapper.particles.weight = getattr(
        #         stars,
        #         band+"_band"
        #         ).value_in(units.LSun)
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
