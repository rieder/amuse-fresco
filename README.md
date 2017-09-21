# Fresco
Fresco aims to simulate observations of particle-based simulations, such as
those of a star cluster. It creates an observation-like image from a list of
stars and/or gas particles. Supported filetypes include AMUSE-type hdf5 files,
Starlab files and plaintext files.

For stars, the temperature and radius are calculated using a stellar evolution
code, if these are not already present in the dataset.

Gas particles may also be read. In combination with stars, these will cause
reflection from nearby stars and optionally extinction of light from
background. Without stars, Fresco will make a density plot of the gas.
Optionally, the gas may also be indicated with contour lines.

A random field of background and foreground stars may be added to the image, as
a way to make the image more natural looking and/or to provide a background
that may be obscured by the gas/dust particles.

![Example image](test.png)

## Requirements

- Python 2.7 (3.x experimental)
- Numpy
- Matplotlib
- AMUSE (https://github.com/amusecode/amuse)
  - FIMap (optional, for extinction)
  - SSE or SeBa (optional, for calculating stellar luminosities and radii)
- Astropy (or Pyfits)

## Authors

Fresco is developed by Inti Pelupessy and Steven Rieder
