"""
Contains helper functions to create images from particle positions and
luminosities.
"""

import os
import gzip

import numpy as np
from numpy import pi
from scipy.ndimage import zoom, gaussian_filter

from astropy.convolution import (
    # convolve,
    convolve_fft,
)

from astropy.io import fits

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


def ease_out_sine(x: float) -> float:
    return sin((x * pi) / 2)


def ease_in_sine(x: float) -> float:
    return 1 - cos((x * pi) / 2)


def ease_in_out_sine(x: float) -> float:
    return -(cos(pi * x) - 1) / 2


def get_viewpoint(
    i: int,
    waypoints: list,  # t -> int, x -> float, y -> float, z -> float
    transitions: list,  # easing functions, length is one less than waypoints
) -> VectorQuantity:
    nsteps = len(waypoints)
    if i == 0:
        y, z, x = waypoints[0][1:]
        vec = [
            y.value_in(units.parsec),
            z.value_in(units.parsec),
            x.value_in(units.parsec),
        ] | units.parsec
        return vec
    # Find last waypoint before step i
    for j in range(nsteps - 1):
        if waypoints[j+1][0] > i:
            break
    t = (i - waypoints[j][0]) / (waypoints[j+1][0] - waypoints[j][0])

    # Linear interpolation
    x = waypoints[j][1] + (waypoints[j+1][1] - waypoints[j][1]) * transitions[j](t)
    y = waypoints[j][2] + (waypoints[j+1][2] - waypoints[j][2]) * transitions[j](t)
    z = waypoints[j][3] + (waypoints[j+1][3] - waypoints[j][3]) * transitions[j](t)

    vec = [
        y.value_in(units.parsec),
        x.value_in(units.parsec),
        z.value_in(units.parsec),
    ] | units.parsec
    return vec


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
    if this_dir == "":
        this_dir = "."
    # FIXME use some relative dir
    datadir = this_dir + "/data/" + instrument + "/"

    for band in "ubvri":
        for i in range(4):
            psf_filename = datadir + band + "%2.2i.fits" % i
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
        flux = getattr(star, band + "_band") / (4 * np.pi * d2)
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
    image_width=12.0 | units.parsec,
    image_size=None,
    mapper_factory=None,
    psf_type="hubble",
    psf_file=None,
    psf_sigma=1.0,
    visualisation_mode="visual",
    zoom_factor=1.0,
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
                * numpy.pi
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
        viewpoint_origin = mapper.parameters.viewpoint
        image_target = mapper.parameters.image_target
        upvector = mapper.parameters.upvector

        print(f"Image target: {image_target}")

        # rotate around upvector in 1 deg steps
        distance = -100 | units.pc
        # initial_viewpoint = (
        #     [
        #         0, 0, 5*distance.value_in(units.pc),
        #     ] | units.pc
        # )
        initial_viewpoint = (
            [
                0,
                0,
                distance.value_in(units.pc),
            ] | units.pc
        )
        mid_viewpoint = (
            [
                0.4 * distance.value_in(units.pc),
                0,
                0,
            ] | units.pc
        )
        final_viewpoint = (
            [
                0,
                0,
                -distance.value_in(units.pc),
            ] | units.pc
        )
        nsteps = 1
        waypoints = [
            (0, initial_viewpoint[0], initial_viewpoint[1], initial_viewpoint[2]),
            (nsteps//2, mid_viewpoint[0], mid_viewpoint[1], mid_viewpoint[2]),
            (nsteps, final_viewpoint[0], final_viewpoint[1], final_viewpoint[2]),
        ]
        transitions = [
            ease_in_sine,
            ease_out_sine
        ]
        # for i in range(0, nsteps, 1):
        for i in range(nsteps):
            relative_viewpoint = get_viewpoint(i, waypoints, transitions)

            viewpoint = image_target + relative_viewpoint
            mapper.parameters.viewpoint = viewpoint
            mapper.parameters.projection_direction = (
                (image_target-viewpoint)
                / (image_target-viewpoint).length()
            )
            print(f"Angle {i}: rel.viewpoint {relative_viewpoint} and projection direction {mapper.parameters.projection_direction}")
            print(f"{image_target-viewpoint}, {relative_viewpoint}, {(image_target-viewpoint).length()}")
            im = mapper.image.pixel_value
            # with gzip.GzipFile(f"flight5_raw_{band}_{i+100:03}.npy.gz", "wb") as f:
            #     np.save(file=f, arr=im)
            if i == 0:
                raw_images[band] = np.fliplr(im)
    mapper.stop()

    convolved_images = {}

    if verbose:
        print("..convolving..")

    if psf_type == "file":
        print(f"PSF: {psf_file}")
        with fits.open(psf_file) as psf:
            for key, val in list(raw_images.items()):
                im1 = convolve(val, psf[0].data)
                convolved_images[key] = im1
    elif psf_type == "hubble":
        psf = get_psf(zoom_factor=zoom_factor)
        if multi_psf:
            a = np.arange(image_size[0]) / float(image_size[0] - 1)
            b = np.arange(image_size[1]) / float(image_size[1] - 1)
            w1 = np.outer(a, b)
            w2 = np.outer(1.0 - a, b)
            w3 = np.outer(a, 1.0 - b)
            w4 = np.outer(1.0 - a, 1.0 - b)
            for key, val in list(raw_images.items()):
                im1 = convolve(val, psf[key + "0"])
                im2 = convolve(val, psf[key + "1"])
                im3 = convolve(val, psf[key + "2"])
                im4 = convolve(val, psf[key + "3"])
                convolved_images[key] = w1 * im1 + w2 * im2 + w3 * im3 + w4 * im4
        else:
            for key, val in list(raw_images.items()):
                im1 = convolve(val, psf[key + "0"])
                convolved_images[key] = im1
    elif psf_type == "gaussian":
        for key, val in list(raw_images.items()):
            im1 = gaussian_filter(val, sigma=psf_sigma, order=0)
            convolved_images[key] = im1

    if visualisation_mode == "visual":
        if verbose:
            print("..conversion to rgb")
        filter_data = get_filter_data()
        source = [filter_data["bess-" + x + ".pass"] for x in sourcebands]

        target = [xyz_data["x"], xyz_data["y"], xyz_data["z"]]

        conv = ColorConverter(source, target)

        ubv = np.array([convolved_images[x] for x in sourcebands])

        xyz = np.tensordot(conv.conversion_matrix, ubv, axes=(1, 0))
        conv_xyz_to_lin = XYZ_to_sRGB_linear()
        srgb_l = np.tensordot(
            conv_xyz_to_lin.conversion_matrix,
            xyz,
            axes=(1, 0),
        )
    elif visualisation_mode == "rgb":
        if len(sourcebands) != 3:
            print("Error: rgb visualisation requires exactly 3 bands")
        srgb_l = []
        for x in sourcebands:
            band = convolved_images[x]
            srgb_l.append(band / band.max())
        srgb_l = np.array(srgb_l)

    if dryrun or vmax is None:
        flat_sorted = np.sort(srgb_l.flatten())
        n = len(flat_sorted)
        vmax = flat_sorted[int(1.0 - 3 * (1.0 - percentile) * n)]
        if verbose:
            print("vmax: %f" % vmax)
    if dryrun:
        return vmax

    conv_lin_to_sRGB = sRGB_linear_to_sRGB()

    srgb = conv_lin_to_sRGB.convert(srgb_l / vmax)

    # r = np.fliplr(srgb[0, :, :])
    # g = np.fliplr(srgb[1, :, :])
    # b = np.fliplr(srgb[2, :, :])

    rgb = np.zeros(image_size[0] * image_size[1] * 4)
    rgb = rgb.reshape(image_size[0], image_size[1], 4)
    rgb[:, :, 0] = srgb[0, :, :].T
    rgb[:, :, 1] = srgb[1, :, :].T
    rgb[:, :, 2] = srgb[2, :, :].T
    rgb[:, :, 3] = 1  # alpha channel

    image = {
        "pixels": rgb,
        "size": rgb.size,
    }

    # convert from double to uint8
    # image["pixels"] = np.clip(image["pixels"] * 256, 0, 255).astype("uint8")
    # convert from double to uint16
    image["pixels"] = np.clip(image["pixels"] * 65536, 0, 65535).astype("uint16")

    return vmax, image
