"""
Convert monochrome images of different color bands to a color image
"""
import numpy as np
from .filters import get_filter_data

from .color_converter import (
    ColorConverter, XYZ_to_sRGB_linear, sRGB_linear_to_sRGB,
)


def monochrome_images_to_color(
    convolved_images,
    sourcebands,
    xyz_data,
    image_size,
    visualisation_mode="visual",
    percentile=0.99,
    vmax=None,
    dryrun=False,
    verbose=False,
):
    if visualisation_mode == "visual":
        # Convert color bands to a "natural" visual representation
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
        # Convert color bands to a false color representation
        if len(sourcebands) != 3:
            print("Error: rgb visualisation requires exactly 3 bands")
        srgb_l = []
        for x in sourcebands:
            band = convolved_images[x]
            srgb_l.append(band)
        srgb_l = np.array(srgb_l)

    if dryrun or vmax is None:
        flat_sorted = np.sort(srgb_l.flatten())
        n = len(flat_sorted)
        vmax = flat_sorted[int(1.0 - 3 * (1.0 - percentile) * n)]
        if verbose:
            print(f"vmax: {vmax}")
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

    # convert from double to uint8 or uint16
    bit_depth = 16
    image["pixels"] = np.clip(
        image["pixels"] * 2**bit_depth,
        0,
        2**bit_depth - 1
    ).astype(f"uint{bit_depth}")

    return image
