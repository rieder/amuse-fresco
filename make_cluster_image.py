import os,sys

import numpy as np

from amuse.units import units,constants,nbody_system
from amuse.io import read_set_from_file
#from amuse.ext.job_server import JobServer
#from PIL import Image
#from ubvivariant import rgb_frame
from ubvinew import rgb_frame

import matplotlib
import matplotlib.pyplot as plt

def calculate_effective_temperature(luminosity, radius):
    return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)

def image_from_stars(
        stars,
        image_width         = 10. | units.parsec,
        image_size          = [1024,1024],
        percentile          = 0.9995,
        calc_temperature    = True,
        age                 = 0.|units.Myr,
        ):
    if calc_temperature:
        try:
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
        except:
            print "Calculating luminosity/temperature for %s old stars..."%(age)
            from amuse.community.sse.interface import SSE
            se = SSE()
            se.particles.add_particles(stars)
            if age > 0|units.Myr:
                print "Evolve start"
                se.evolve_model(age)
                print "Evolve done"
            stars.luminosity = se.particles.luminosity
            stars.radius = se.particles.radius
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
            se.stop()

    print "Start making RGB image"
    vmax, rgb = rgb_frame(
            stars,
            dryrun      = False,
            image_width = image_width,
            multi_psf   = True,
            image_size  = image_size,
            percentile  = percentile,
            sourcebands = "ubvri",
            )
    return rgb['pixels']
    

if __name__=="__main__":
    filename    = sys.argv[1]
    
    length_unit = units.parsec
    dpi         = 1200
    image_width_arcsec = 160

    image_width = 10. | units.parsec #FIXME this should depend on the distance!
                                     # size = fixed at nr of arcmin
    image_size  = [2048,2048] #FIXME the psf is fixed pixel size, so the pixels
    # in the image here reflects how much pixels will be spread!
    # In principle, this should be a fixed number, equal to the number of
    # pixels of the ccd for which the psf was made.
    # Also, the image should reflect the distance to the "observed" cluster.

    age         = 500.|units.Myr
    percentile  = 0.9995 ## for determining vmax

    xmin        = -0.5 * image_width.value_in(length_unit)
    xmax        =  0.5 * image_width.value_in(length_unit)
    ymin        = -0.5 * image_width.value_in(length_unit)
    ymax        =  0.5 * image_width.value_in(length_unit)

    stars       = read_set_from_file(
            filename,
            "amuse",
            close_file  = True,
            )

    # FIXME: add these features
    # - Rotate so that xy = observed x/y axes of figure
    # - Scale positions to desired ra/dec (script Alison)
    # - calculate vmax based on nr of photons/exposure time

    imagefilename = sys.argv[1]+".png"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_xlabel("[%s]"%(length_unit))
    ax.set_ylabel("[%s]"%(length_unit))

    image = image_from_stars(
            stars,
            image_width         = image_width,
            image_size          = image_size,
            percentile          = percentile,
            calc_temperature    = True,
            age                 = age,
            )

    plt.imshow(
            image,
            origin  = 'lower',
            extent  = [
                xmin,
                xmax,
                ymin,
                ymax,
                ],
            )
    
    plt.savefig(
            imagefilename,
            dpi = dpi,
            )
