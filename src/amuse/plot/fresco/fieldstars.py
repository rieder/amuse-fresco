import numpy as np
from amuse.datamodel import Particles
from amuse.units import units
from amuse.ic.salpeter import new_salpeter_mass_distribution


def new_field_stars(
    number_of_stars,
    width=10 | units.parsec,
    height=10 | units.parsec,
    depth=100 | units.parsec,
    massdistribution="salpeter",
    agespread=3 | units.Gyr,
    seed=1701,
):
    np.random.seed(seed)
    stars = Particles(number_of_stars)
    stars.x = (np.random.random(number_of_stars) - 0.5) * width
    stars.y = (np.random.random(number_of_stars) - 0.5) * height
    stars.z = (np.random.random(number_of_stars) - 0.02) * depth
    if massdistribution == "salpeter":
        stars.mass = new_salpeter_mass_distribution(number_of_stars)

    return stars
