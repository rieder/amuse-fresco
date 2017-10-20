from amuse.units import units, nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_powerlaw_mass_distribution
from amuse.io import write_set_to_file

Nstars = 1000
Ngas = 1000

mass = new_powerlaw_mass_distribution(
        Nstars,
        mass_max=20. | units.MSun,
        )
conv = nbody_system.nbody_to_si(
        mass.sum(),
        1. | units.parsec,
        )
p = new_plummer_model(
        Nstars,
        conv,
        )
p.mass = mass
write_set_to_file(p, "stars.hdf5", "amuse", append_to_file=False)

p = new_plummer_model(
        Ngas,
        conv,
        )
p.h_smooth = 0.1 | units.parsec
write_set_to_file(p, "gas.hdf5", "amuse", append_to_file=False)
