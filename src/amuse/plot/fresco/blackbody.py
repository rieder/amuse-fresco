# -*- coding: utf-8 -*-
"""Functions for calculating a black body spectrum.

"""

from __future__ import (
        print_function,
        division,
        )
import numpy
from numpy import (
    pi, e,
)

from amuse.units import units
from amuse.units.constants import (
    kB, h, c,
)
from amuse.units.constants import Rydberg_constant as Ry
from amuse.units.constants import Stefan_hyphen_Boltzmann_constant as sigma


def B_nu(nu, t):
    return 2*h*nu**3/c**2 * 1. / (e**(h*nu/kB/t)-1)


def B_lambda(wavelength, temperature):
    x = (
        2 * h * c**2
        / wavelength**5
        / (
            e**(
                h * c
                / (
                    (numpy.outer(
                        temperature.value_in(units.K),
                        wavelength.value_in(units.m)
                        ) | units.K * units.m)
                    * kB
                    )
                ) - 1
            )
        )
    return x


def energy(nu):
    return h*nu


def freq(e):
    return e/h


def freq_from_wavenumber(k):
    return c*k


def wavelength(nu):
    return c/nu


def freq_from_wavelength(wavelength):
    return c/wavelength


def wiens_lambda_max(T):
    b = 2897768.6 | units.nano(units.m) * units.K
    return b/T


def wiens_T_from_lambda_max(wavelength):
    b = 2897768.6 | units.nano(units.m) * units.K
    return b/wavelength


def energy_flux(T, lowfreq=0. | units.s**-1, N=1000):
    nu = (numpy.arange(N+1)+1.)/N*(kB*T)/h*25. + lowfreq
    b = pi*B_nu(nu, T)
    return (b[1:]+b[:-1]).sum()/2*(nu[1]-nu[0])


def energy_flux2(T, lambdas=None, throughput=1., N=1000):
    if lambdas is None:
        lmax = wiens_lambda_max(T)
        lambdas = lmax * 10**(-2. + 4. * numpy.arange(N+1)/float(N))
    b = pi*throughput*B_lambda(lambdas, T)
    return numpy.trapz(b.number, x=lambdas.number) | (b.unit*lambdas.unit)


def photon_flux(T, lowfreq=0. | units.s**-1, N=100000):
    nu = (numpy.arange(N+1)+1.)/N*(kB*T)/h*25. + lowfreq
    n = pi*B_nu(nu, T)/energy(nu)
    return (n[1:]+n[:-1]).sum()/2*(nu[1]-nu[0])


def total_bolometric_flux(T):
    return sigma*T**4


if __name__ == "__main__":
    T = 10000. | units.K

    print((wiens_lambda_max(T)))

    print((energy_flux(T).in_(units.W * units.m**-2)))
    print((energy_flux2(T).in_(units.W * units.m**-2)))
    print(((sigma*T**4).in_(units.W * units.m**-2)))

    print()
    nf = photon_flux(T, lowfreq=freq_from_wavenumber(Ry))
    print((numpy.log10(nf.value_in(units.cm**-2 * units.s**-1))))
    print()

    a = photon_flux(T)
    print((numpy.log10(a.value_in(units.cm**-2 * units.s**-1))))
    b = sigma*T**4 / (kB*T)/2.7
    print((numpy.log10(b.value_in(units.cm**-2 * units.s**-1))))

    print((b/a))
    print((nf/b))

    print((wiens_T_from_lambda_max(300. | units.nano(units.m))))
    print((wiens_T_from_lambda_max(610. | units.nano(units.m))))
    print((wiens_T_from_lambda_max(920. | units.nano(units.m))))
