"""
Fresco creates a "simulated observation" of a set of particles.
Particles can be "stars" (point sources emitting light) or "gas" (emitting,
reflecting and/or obscuring light). Gas may also be displayed with contour
lines.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from amuse.units import units, constants
from amuse.community.fi import FiMap
from .ubvi import rgb_frame
from .gridify import map_to_grid


def evolve_to_age(stars, age, stellar_evolution="seba"):
    "Evolve stars to specified age with specified code"
    if stellar_evolution.lower() == "seba":
        from amuse.community.seba import Seba

        stellar_evolution = Seba()
    elif stellar_evolution.lower() == "sse":
        from amuse.community.sse import Sse

        stellar_evolution = Sse()
        # SSE can result in nan values for luminosity/radius
    else:
        raise ValueError(
            f"No such stellar evolution code {stellar_evolution} or no code "
            f"specified"
        )

    stellar_evolution.particles.add_particles(stars)
    if age > 0 | units.yr:
        stellar_evolution.evolve_model(age)
    stars.luminosity = (
        np.nan_to_num(stellar_evolution.particles.luminosity.value_in(units.LSun))
        | units.LSun
    )

    stars.radius = stellar_evolution.particles.radius
    # prevent zero/nan radius.
    location_of_nan_radius = np.where(
        np.nan_to_num(stars.radius.value_in(units.RSun)) == 0.0
    )
    stars[location_of_nan_radius].radius = 0.01 | units.RSun

    stellar_evolution.stop()
    return


def calculate_effective_temperature(luminosity, radius):
    temp = (
        np.nan_to_num(
            (
                (luminosity / (constants.four_pi_stefan_boltzmann * radius**2))**0.25
            ).value_in(units.K)
        )
        | units.K
    )
    return temp


def make_image(
    maps=None,
    stars=None,
    gas=None,
    converter=None,
    image_width=10.0 | units.parsec,
    image_height=None,
    distance=100 | units.pc,
    padding_fraction=0.1,
    image_size=[1024, 1024],
    percentile=0.9995,
    age=0.0 | units.Myr,
    sourcebands=["u", "b", "v", "r", "i"],
    vmax=None,
    calc_temperature=True,
    zoom_factor=1.0,
    psf_type="hubble",
    psf_file=None,
    psf_sigma=1.0,
    extinction=False,
    return_vmax=False,
    origin=None,
    visualisation_mode="visual",
    viewpoint=None
):
    """
    Makes image from gas and stars
    """
    if origin is None:
        origin = [0, 0, 0] | units.pc
    image_target = origin
    image_target[0] += image_width / 2
    image_target[1] += image_width / 2
    viewpoint = image_target.copy()
    viewpoint[2] -= distance
    if image_height is None:
        image_height = image_width
    mode = []
    if gas is not None:
        mode.append("gas")
        # Only keep the relevant gas
        gas_all = gas
        x_lower = origin[0] - padding_fraction * image_width
        x_upper = origin[0] + (padding_fraction + 1) * image_width
        y_lower = origin[1] - padding_fraction * image_height
        y_upper = origin[1] + (padding_fraction + 1) * image_height
        gas_relevant_indices = (
            (gas_all.x >= x_lower)
            & (gas_all.x <= x_upper)
            & (gas_all.y >= y_lower)
            & (gas_all.y <= y_upper)
        )
    if stars is not None:
        mode.append("stars")
    if not mode:
        raise ValueError(
            "No gas or stars specified, or stars/gas out of bounds"
        )

    perspective = True

    def mapper():
        print("Setting up FiMap")
        mapper = FiMap(converter, mode="openmp", )

        # mapper.parameters.minimum_distance = 1. | units.au
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
        if perspective:
            # NOTE: in perspective mode, the image may rotate if the axes are not exactly clear
            mapper.parameters.projection_mode = "perspective"
            # image_target = [0, 0, 0] | units.pc  # coordinate center?
            mapper.parameters.image_target = image_target
            # Warning: when projection direction and upvector align, need
            # to do a "meridian flip"!
            mapper.parameters.upvector = [0, -1, 0]
            mapper.parameters.viewpoint = viewpoint
            mapper.parameters.projection_direction = (
                (image_target-viewpoint)
                / (image_target-viewpoint).length()
            )
        mapper.parameters.extinction_flag = extinction
        with open(f"mapper_width_{image_width}a.txt", "w") as f:
            f.write(str(mapper.parameters))
        print(mapper.parameters)
        print(f"viewpoint: {viewpoint}")
        print(f"image_target: {image_target}")
        print(f"length: {(image_target-viewpoint).length()}")
        return mapper

    if "stars" not in mode:
        image = column_density_map(
            gas,
            image_width=image_width,
            image_size=image_size,
            mapper_factory=mapper,
            zoom_factor=zoom_factor,
            psf_type=psf_type,
            psf_sigma=psf_sigma,
            return_vmax=return_vmax,
        )
    else:
        print("vmax: ", vmax)
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
            zoom_factor=zoom_factor,
            psf_type=psf_type,
            psf_file=psf_file,
            psf_sigma=psf_sigma,
            return_vmax=return_vmax,
            visualisation_mode=visualisation_mode,
        )
    return image


def column_density_map(
    gas,
    image_width=10.0 | units.parsec,
    image_size=[1024, 1024],
    mapper_factory=None,
    zoom_factor=1.0,
    psf_type="gaussian",
    psf_sigma=10.0,
    return_vmax=False,
):
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
    if return_vmax:
        return (im, -1)
    return im


def image_from_stars(
    stars,
    image_width=10.0 | units.parsec,
    image_size=[1024, 1024],
    percentile=0.9995,
    calc_temperature=False,
    age=0.0 | units.Myr,
    sourcebands=["u", "b", "v", "r", "i"],
    gas=None,
    vmax=None,
    mapper_factory=None,
    zoom_factor=1.0,
    psf_type="hubble",
    psf_file=None,
    psf_sigma=1.0,
    return_vmax=False,
    visualisation_mode="visual",
):
    print(f"Return_vmax: {return_vmax}")
    if (
        not hasattr(stars, "temperature")
        and hasattr(stars, "radius")
        and hasattr(stars, "luminosity")
    ):
        calc_temperature = True
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
        zoom_factor=zoom_factor,
        psf_type=psf_type,
        psf_file=psf_file,
        psf_sigma=psf_sigma,
        visualisation_mode=visualisation_mode,
    )
    if return_vmax:
        return rgb["pixels"], vmax
    return rgb["pixels"]


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
    wcs=None,
):
    if wcs is not None:
        fig, ax = plt.subplots(subplot_kw=dict(projection=wcs))
        return fig
    if fig is None:
        if plot_axes:
            left = 0.2
            bottom = 0.2
        else:
            left = 0.0
            bottom = 0.0
        right = 1.0
        top = 1.0
        figwidth = image_size[0] / dpi / (right - left)
        figheight = image_size[1] / dpi / (top - bottom)
        figsize = (figwidth, figheight)

        xmin = x_offset.value_in(length_unit) - 0.5 * image_width.value_in(length_unit)
        xmax = x_offset.value_in(length_unit) + 0.5 * image_width.value_in(length_unit)
        ymin = y_offset.value_in(length_unit) - 0.5 * image_width.value_in(length_unit)
        ymax = y_offset.value_in(length_unit) + 0.5 * image_width.value_in(length_unit)

        fig = plt.figure(figsize=figsize, dpi=dpi)
        print(f"using wcs projection {wcs}")
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
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
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_facecolor("black")
    return fig
