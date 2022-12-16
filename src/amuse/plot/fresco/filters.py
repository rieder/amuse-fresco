import os

import numpy

from amuse.units import units, quantities

import matplotlib.pyplot as plt

# from blackbody import (
#         energy_flux2,
#         # total_bolometric_flux,
#         B_lambda,
#         )

bessellfilters = {
    "u": "bess-u.pass",
    "b": "bess-b.pass",
    "v": "bess-v.pass",
    "r": "bess-r.pass",
    "i": "bess-i.pass",
}  # Bessell 1990 filters


def get_filter_data(
    instrument="WFPC_II_WFC3",
    filters=bessellfilters
):
    this_dir, this_filename = os.path.split(__file__)
    if this_dir == "":
        this_dir = "."
    data_dir = this_dir + "/data/" + instrument + "/"

    filter_data = dict()

    for band in filters:
        wavelength = [] | units.angstrom
        throughput = []
        with open(data_dir + filters[band], "r") as f:
            lines = f.readlines()
            for line in lines:
                line = line.split()
                wavelength.append(float(line[0]) | units.angstrom)
                throughput.append(float(line[1]))
            throughput = numpy.array(throughput)

            filter_data[filters[band]] = dict(
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


def plot_filters(
    filters=bessellfilters,
    min_wavelength=3000 | units.angstrom,
    max_wavelength=9200 | units.angstrom,
    wavelength_unit=units.nano(units.m),
):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    filter_data = get_filter_data(filters=filters)

    n = 1000
    for band in filters:
        wavelength = (
            min_wavelength + (max_wavelength-min_wavelength)
            * (
                numpy.array(range(n+1))
                / n
            )
        )
        xp = filter_data[filters[band]]['wavelength']
        fp = filter_data[filters[band]]['throughput']
        f = numpy.interp(
            wavelength.value_in(wavelength_unit),
            xp=xp.value_in(wavelength_unit),
            fp=fp,
            left=0.,
            right=0.,
        )

        ax.plot(wavelength.value_in(wavelength_unit), f)
        ax.set_xlabel(f'wavelength ({wavelength_unit})')
        ax.set_ylabel('transmission')
        

    plt.show()


if __name__ == "__main__":

    plot_filters(filters=bessellfilters)

    T = 5000. | units.K

    fb = filter_band_flux(
            filter_data["bess-u.pass"],
            lambda x: B_lambda(x, T),
            ).in_(units.W/units.m**2)

    print(fb)
    print(
        (
            fb * (1. | units.RSun)**2
            / (1. | units.AU)**2
        ).in_(units.W/units.m**2)
    )
    print(
        (
            energy_flux2(5778. | units.K)
            * (1. | units.RSun)**2
            / (1. | units.AU)**2
        ).in_(units.W/units.m**2)
    )
