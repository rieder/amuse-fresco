"""
Contains helper functions to create images from particle positions and
luminosities.
"""

import os
import gzip

import numpy as np
from numpy import pi
from scipy.ndimage import zoom, gaussian_filter

from amuse.datamodel import Particles
from amuse.units import units
from amuse.units.quantities import VectorQuantity
from amuse.units.trigo import cos, sin, arccos, arctan2

# import logging
# logging.basicConfig(level=logging.DEBUG)

from .filters import (
    filter_band_flux,
    get_filter_data,
    filter_band_lambda,
)
from .xyz import xyz_data
from .blackbody import B_lambda
from .color_converter import ColorConverter, XYZ_to_sRGB_linear, sRGB_linear_to_sRGB
from .convolve_psf import convolve_image_with_psf
from .monochrome_images_to_color import monochrome_images_to_color


def assign_weights_and_opacities(
    band,
    mapper_stars,
    mapper_gas,
    stars,
    gas,
    dust_to_gas=0.01,
    dust_extinction_cross_section_v_band=4.9e-22 | units.cm**2,
    dust_albedo_v_band=0.01,
    Nstar=50,  # number of stars to use for reflected light on gas/dust
):
    mapper_stars.weight = getattr(stars, band + "_band").value_in(units.LSun)

    if gas.is_empty():
        return

    lambda_eff = filter_band_lambda("bess-" + band + ".pass")
    lambda_v = filter_band_lambda("bess-v.pass")

    f_H = 0.7
    dust_cross_section = (
        dust_extinction_cross_section_v_band * (lambda_eff / lambda_v) ** -1.5
    )
    dust_albedo = dust_albedo_v_band * (lambda_eff / lambda_v) ** 0.0

    mapper_gas.opacity_area = f_H * gas.mass.value_in(units.amu) * dust_cross_section

    weight = np.zeros(len(mapper_gas))
    stars_ordered_by_lum = stars.sorted_by_attribute(band + "_band")
    Nstar = min(Nstar, len(stars))
    for star in stars_ordered_by_lum[-Nstar:]:
        d2 = (
            (gas.x - star.x) ** 2
            + (gas.y - star.y) ** 2
            + (gas.z - star.z) ** 2
            + 0.25 * gas.radius**2
        )
        flux = getattr(star, band + "_band") / (4 * pi * d2)
        weight += (flux * mapper_gas.opacity_area * dust_albedo).value_in(units.LSun)
    mapper_gas.weight = weight
    return


def rgb_frame(
    stars,
    gas=None,
    dryrun=False,
    vmax=None,
    percentile=0.9995,
    multi_psf=False,
    sourcebands=None,
    image_size=None,
    mapper_factory=None,
    psf_type="hubble",
    psf_file=None,
    psf_sigma=1.0,
    visualisation_mode="visual",
    zoom_factor=1.0,
    distance=None,
    flightpath=None,
    verbose=False,
):
    if sourcebands is None:
        sourcebands = ["u", "b", "v", "r", "i"]
    if image_size is None:
        image_size = [1024, 1024]

    if gas is None:
        gas = Particles()

    if verbose:
        print("luminosities..")

    for band in sourcebands:
        if not hasattr(stars, band + "_band"):
            setattr(
                stars,
                band + "_band",
                4
                * pi
                * stars.radius**2
                * filter_band_flux(
                    "bess-" + band + ".pass",
                    lambda x: B_lambda(x, stars.temperature),
                ),
            )
    if verbose:
        print("..raw images..")

    # Use mapper to make raw (pre-convolved) images
    mapper = mapper_factory()
    print(f"Adding {len(stars)} stars and {len(gas)} gas particles")
    stars_in_mapper = mapper.particles.add_particles(stars)
    gas_in_mapper = mapper.particles.add_particles(gas)
    # mapper.parameters.projection_direction = [0, 0, 1]
    # mapper.parameters.upvector = [0, -1, 0]

    raw_images = {}
    for band in sourcebands:
        raw_images[band] = []

    for i, viewpoint in enumerate(flightpath):
        image_target = mapper.parameters.image_target
        mapper.parameters.viewpoint = viewpoint
        mapper.parameters.projection_direction = (
            (image_target-viewpoint)
            / (image_target-viewpoint).length()
        )
        print(
            f"Angle {i}: rel.viewpoint {viewpoint - image_target} and "
            f"projection direction {mapper.parameters.projection_direction}"
        )

        for band in sourcebands:
            assign_weights_and_opacities(
                band,
                stars_in_mapper,
                gas_in_mapper,
                stars,
                gas,
                Nstar=1000,
            )
            # mapper.particles.weight = getattr(
            #         stars,
            #         band+"_band"
            #         ).value_in(units.LSun)

            print(f"Image target: {image_target}")

            im = mapper.image.pixel_value
            # with gzip.GzipFile(f"flight5_raw_{band}_{i+100:03}.npy.gz", "wb") as f:
            #     np.save(file=f, arr=im)

            # raw_images is a dictionary of lists, and looks similar to this:
            # raw_images = {
            #     "u": [u_raw1, u_raw2, u_raw3],
            #     "b": [b_raw1, b_raw2, b_raw3],
            #     ...etc
            # }
            raw_images[band].append(np.fliplr(im))

    if verbose:
        print("..convolving..")

    convolved_images = {}
    for color_band in sourcebands:
        convolved_images[color_band] = []
        for raw_image in raw_images[color_band]:
            convolved_images[color_band].append(
                convolve_image_with_psf(
                    raw_image,
                    color_band,
                    psf_type=psf_type,
                    psf_file=psf_file,
                    psf_sigma=psf_sigma,
                    zoom_factor=zoom_factor,
                    multi_psf=multi_psf,
                    image_size=image_size,
                    verbose=verbose,
                )
            )

    images = []
    for i in range(len(flightpath)):
        convolved_images_i = {}
        for color_band in sourcebands:
            convolved_images_i[color_band] = convolved_images[color_band][i]

        images.append(
            monochrome_images_to_color(
                convolved_images_i,
                sourcebands=sourcebands,
                xyz_data=xyz_data,
                image_size=image_size,
                visualisation_mode=visualisation_mode,
                percentile=percentile,
                vmax=vmax,
                dryrun=dryrun,
                verbose=verbose,
            )
        )

    mapper.stop()
    return vmax, images
