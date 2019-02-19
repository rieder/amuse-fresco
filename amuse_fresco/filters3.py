import os

import numpy

from amuse.units import units, quantities

# from matplotlib import pyplot

# from blackbody import (
#         energy_flux2,
#         # total_bolometric_flux,
#         B_lambda,
#         )


def get_filter_data(
        instrument="WFPC_II_WFC3",
):
    this_dir, this_filename = os.path.split(__file__)
    if this_dir == "":
        this_dir = "."
    data_dir = this_dir + "/data/" + instrument + "/"

    filter_data = dict()

    bessellfilters = ["bess-u.pass",
                      "bess-b.pass",
                      "bess-v.pass",
                      "bess-r.pass",
                      "bess-i.pass"]  # Bessell 1990 filters

    for bessellfilter in bessellfilters:

        f = open(data_dir + bessellfilter, "r")

        lines = f.readlines()

        wavelength = [] | units.angstrom
        throughput = []

        for l in lines:
            line = l.split()
            wavelength.append(float(line[0]) | units.angstrom)
            throughput.append(float(line[1]))
        throughput = numpy.array(throughput)

        filter_data[bessellfilter] = dict(
            wavelength=wavelength,
            throughput=throughput,
        )
    return filter_data


def filter_band_flux(
        fdata,
        source,
        wavelength=None,
):
    if fdata.__class__.__name__ == "str":
        fdata = get_filter_data()[fdata]
    if wavelength is None:
        wavelength = fdata['wavelength']
        throughput = fdata['throughput']
    else:
        xp = fdata['wavelength']
        fp = fdata['throughput']
        throughput = numpy.interp(
            wavelength.value_in(units.angstrom),
            xp=xp.value_in(units.angstrom),
            fp=fp,
            left=0.,
            right=0.,
        )
    src = throughput * source(wavelength)
    src = quantities.to_quantity(src)
    return numpy.trapz(
        src.number,
        x=wavelength.number,
    ) | (src.unit*wavelength.unit)


def filter_band_lambda(fdata):
    return (
        filter_band_flux(fdata, lambda x: x)
        / filter_band_flux(fdata, lambda x: 1.)
    ).in_(units.angstrom)


# def plot_filters():
#
#     n = 1000
#     for bessellfilter in filters:
#         wavelength = (
#                 3000. + (9200.-3000.)
#                 * numpy.array(range(n+1))
#                 / n
#                 ) | units.angstrom
#         xp = filter_data[bessellfilter]['wavelength']
#         fp = filter_data[bessellfilter]['throughput']
#         f = numpy.interp(
#                 wavelength.value_in(units.nano(units.m)),
#                 xp=xp.value_in(units.nano(units.m)),
#                 fp=fp,
#                 left=0.,
#                 right=0.,
#                 )
#
#         pyplot.plot(lam.value_in(units.nano(units.m)), f)
#
#     pyplot.show()


# if __name__ == "__main__":
#
#     plot_filters()
#
#     T = 5000. | units.K
#
#     fb = filter_band_flux(
#             filter_data["bess-u.pass"],
#             lambda x: B_lambda(x, T),
#             ).in_(units.W/units.m**2)
#
#     print fb
#     print (
#             fb
#             * (1. | units.RSun)**2
#             / (1. | units.AU)**2
#             ).in_(units.W/units.m**2)
#     print (
#             energy_flux2(5778. | units.K)
#             * (1. | units.RSun)**2
#             / (1. | units.AU)**2
#             ).in_(units.W/units.m**2)
