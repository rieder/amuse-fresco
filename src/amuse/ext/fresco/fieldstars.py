# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        division,
        )
import numpy
from amuse.datamodel import Particles
from amuse.units import units
from amuse.ic.salpeter import new_salpeter_mass_distribution


def new_field_stars(
    N,
    width=10 | units.parsec,
    height=10 | units.parsec,
    depth=100 | units.parsec,
    massdistribution="salpeter",
    agespread=3 | units.Gyr,
    seed=1701,
):
    np.random.seed(seed)
    stars = Particles(N)
    stars.x = (numpy.random.random(N)-0.5) * width
    stars.y = (numpy.random.random(N)-0.5) * height
    stars.z = (numpy.random.random(N)-0.02) * depth
    if massdistribution == "salpeter":
        stars.mass = new_salpeter_mass_distribution(N)

    return stars
