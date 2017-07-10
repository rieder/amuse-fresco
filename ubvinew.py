import os

import numpy

from numpy.fft import fft2, ifft2
from numpy import log

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
        MinPad=True,
        pad=True,
        ):
    result = convolve_fft(image, kernel, boundary='fill')
    # result = convolve(image,kernel)
    return result


def _Convolve(
        image1,
        image2,
        MinPad=True,
        pad=True,
        ):
    """ Not so simple convolution """
    # Just for comfort:
    FFt = fft2
    iFFt = ifft2

    # The size of the images:
    r1, c1 = image1.shape
    r2, c2 = image2.shape

    # MinPad results simpler padding,smaller images:
    if MinPad:
        r = r1+r2
        c = c1+c2
    else:
        # if the Numerical Recipies says so:
        r = 2 * max(r1, r2)
        c = 2 * max(c1, c2)

    # For nice FFT, we need the power of 2:
    if pad:
        pr2 = int(log(r)/log(2.0) + 1.0)
        pc2 = int(log(c)/log(2.0) + 1.0)
        rOrig = r
        cOrig = c
        r = 2**pr2
        c = 2**pc2

    # numpy fft has the padding built in, which can save us some steps here.
    # The thing is the s(hape) parameter:
    # fftimage = FFt(image1,s=(r,c)) * FFt(image2,s=(r,c))
    fftimage = (
            FFt(
                image1.astype(dtype=numpy.float32),
                shape=(r, c),
                )
            * FFt(
                image2[::-1, ::-1].astype(dtype=numpy.float32),
                shape=(r, c),
                )
            )

    if pad:
        return (iFFt(fftimage))[:rOrig, :cOrig].real
    else:
        return (iFFt(fftimage)).real


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

    print "luminosities.."

    allstars = stars
    # stars = stars.select(lambda x: x > 0.10|units.LSun, ["luminosity"])
    print "Using %i stars out of %i" % (len(stars), len(allstars))

    for band in sourcebands:
        # setattr(
        #        stars,
        #        band+"_band",
        #        4 * numpy.pi * stars.radius**2 *
        #        __filter_band_flux__(
        #            "bess-" + band + ".pass",
        #            lambda x: B_lambda(x, stars.temperature),
        #            ),
        #        )
        for star in stars:
            setattr(
                    star,
                    band+"_band",
                    4 * numpy.pi * star.radius**2 *
                    filter_band_flux(
                        "bess-" + band + ".pass",
                        lambda x: B_lambda(x, star.temperature),
                        ),
                    )

    print "..raw images.."

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
    # print mapper.parameters

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

    print "..convolving.."

    if multi_psf:
        a = numpy.arange(image_size[0])/float(image_size[0]-1)
        b = numpy.arange(image_size[1])/float(image_size[1]-1)
        w1 = numpy.outer(a, b)
        w2 = numpy.outer(1.-a, b)
        w3 = numpy.outer(a, 1.-b)
        w4 = numpy.outer(1.-a, 1.-b)
        for key, val in raw_images.items():
            # xpad,ypad   = psf[key+'0'].shape
            im1 = Convolve(val, psf[key+'0'])  # [xpad/2:-xpad/2,ypad/2:-ypad/2]
            im2 = Convolve(val, psf[key+'1'])  # [xpad/2:-xpad/2,ypad/2:-ypad/2]
            im3 = Convolve(val, psf[key+'2'])  # [xpad/2:-xpad/2,ypad/2:-ypad/2]
            im4 = Convolve(val, psf[key+'3'])  # [xpad/2:-xpad/2,ypad/2:-ypad/2]
            convolved_images[key] = w1*im1+w2*im2+w3*im3+w4*im4
    else:
        for key, val in raw_images.items():
            # xpad,ypad   = psf[key+'0'].shape
            im1 = Convolve(val, psf[key+'0'])  # [xpad/2:-xpad/2,ypad/2:-ypad/2]
            convolved_images[key] = im1

    print "..conversion to rgb"

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
        print "vmax:", vmax
    if dryrun:
        return vmax

    conv_lin_to_sRGB = sRGB_linear_to_sRGB()

    srgb = conv_lin_to_sRGB.convert(srgb_l/vmax)
    # print srgb.shape()

    # r = srgb[0,:,:].transpose()
    # g = srgb[1,:,:].transpose()
    # b = srgb[2,:,:].transpose()
    r = numpy.fliplr(srgb[0, :, :])
    g = numpy.fliplr(srgb[1, :, :])
    b = numpy.fliplr(srgb[2, :, :])

    rgb = numpy.zeros(image_size[0] * image_size[1] * 3)
    rgb = rgb.reshape(image_size[0], image_size[1], 3)
    # print image_size[0], image_size[1], rgb.shape, r.shape, g.shape, b.shape
    rgb[:, :, 0] = r
    rgb[:, :, 1] = g
    rgb[:, :, 2] = b

    image = dict(
        pixels=rgb,
        size=rgb.size,
        )

    return vmax, image
