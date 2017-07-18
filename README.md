# FRESCO
Fresco aims to simulate observations of particle-based simulations, such as
those of a star cluster. It creates a Hubble WFC3-like image from an AMUSE-type
hdf5 file of stars or gas particles (other filetypes, like the Starlab format,
are also supported). 

For stars, the temperature and radius are calculated using the SeBa (default)
or SSE stellar evolution code, if it is not available in the dataset.

Gas/dust particles may also be read, and will cause extinction and reflection
of light from background- and nearby stars, respectively.

A random field of back/foreground stars may be added to the image, as a way to
make the image more natural looking or to provide a background that may be
obscured by the gas/dust particles.

![Example image](test.png)

## Requirements

- Python 2.7 (3.x experimental)
- Numpy
- Matplotlib
- AMUSE (https://github.com/amusecode/amuse)
  - FIMap
  - SSE (optional)
- Astropy (or Pyfits)

## Authors

Fresco is written by Inti Pelupessy and Steven Rieder
