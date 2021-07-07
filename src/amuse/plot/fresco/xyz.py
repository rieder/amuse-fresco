import os

import numpy

from amuse.units import units

this_dir, this_filename = os.path.split(__file__)
if this_dir == "":
    this_dir = "."
instrument = "WFPC_II_WFC3"
data_dir = this_dir + "/data/" + instrument + "/"

f       = open(data_dir + "ciexyz31.csv","r")

lines   = f.readlines()

cielam  = [] | units.nano(units.m)
x       = []
y       = []
z       = []

for l in lines:
    line    = l.split(",")
    cielam.append( float(line[0]) | units.nano(units.m) )
    x.append( float(line[1]) )
    y.append( float(line[2]) )
    z.append( float(line[3]) )
x   = numpy.array(x)
y   = numpy.array(y)
z   = numpy.array(z)

xyz_data        = dict()
xyz_data['x']   = dict(wavelength=cielam, throughput=x)
xyz_data['y']   = dict(wavelength=cielam, throughput=y)
xyz_data['z']   = dict(wavelength=cielam, throughput=z)
