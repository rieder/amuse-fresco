"""
Convolve raw images with a point-spread function
"""
import os
import numpy as np
from scipy.ndimage import zoom, gaussian_filter

from astropy.io import fits
from astropy.convolution import (
    # convolve,
    convolve_fft,
)


def convolve(
    image,
    kernel,
):
    result = convolve_fft(image, kernel, boundary="fill")
    return result


def get_psf(
    instrument="WFC3_TINYTIM",
    zoom_factor=1.0,
):
    psf = {}
    this_dir = os.path.split(__file__)[0]
    datadir = os.path.join(this_dir, "data", instrument)

    for band in "ubvri":
        for i in range(1):
            psf_filename = os.path.join(datadir, band + f"{i:02d}.fits")
            f = fits.open(psf_filename)
            if zoom_factor != 1.0:
                psf[band + str(i)] = zoom(
                    np.array(f[0].data),
                    zoom_factor,
                    order=3,
                    mode="constant",
                    cval=0.0,
                    prefilter=True,
                )
            else:
                psf[band + str(i)] = np.array(f[0].data)
    return psf


def convolve_image_with_psf(
    raw_image,
    colorband,
    psf_type="hubble",
    psf_file=None,
    psf_sigma=None,
    zoom_factor=None,
    multi_psf=False,
    image_size=None,
    verbose=False,
):
    if psf_type == "file":
        print(f"PSF: {psf_file}")
        with fits.open(psf_file) as psf:
            im1 = convolve(raw_image, psf[0].data)
            convolved_image = im1
    elif psf_type == "file":
        print(f"PSF: {psf_file}")
        with fits.open(psf_file) as psf:
            im1 = convolve(raw_image, psf[0].data)
            convolved_image = im1
    elif psf_type == "hubble":
        psf = get_psf(zoom_factor=zoom_factor)
        if multi_psf:
            a = np.arange(image_size[0]) / float(image_size[0] - 1)
            b = np.arange(image_size[1]) / float(image_size[1] - 1)
            w1 = np.outer(a, b)
            w2 = np.outer(1.0 - a, b)
            w3 = np.outer(a, 1.0 - b)
            w4 = np.outer(1.0 - a, 1.0 - b)
            im1 = convolve(raw_image, psf[colorband + "0"])
            im2 = convolve(raw_image, psf[colorband + "1"])
            im3 = convolve(raw_image, psf[colorband + "2"])
            im4 = convolve(raw_image, psf[colorband + "3"])
            convolved_image = w1 * im1 + w2 * im2 + w3 * im3 + w4 * im4
        else:
            im1 = convolve(raw_image, psf[colorband + "0"])
            convolved_image = im1
    elif psf_type == "gaussian":
        im1 = gaussian_filter(raw_image, sigma=psf_sigma, order=0)
        convolved_image = im1
    else:
        raise ValueError(f"Unknown PSF type: {psf_type}")

    return convolved_image
