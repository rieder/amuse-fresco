# FRESCO
Fabricate Reality-Emulating Star Cluster Observations

Fresco creates a Hubble WFC3-like image from an AMUSE-type hdf5 file of stars
(other filetypes, like the Starlab format, are also supported). The positions
and masses of the stars are always used, if information on the temperature and
radius is not available this is calculated using the SSE stellar evolution
code.
Gas particles may also be read, and will cause extinction and reflection of
light from background- and nearby stars, respectively.

![Example image](test.png)

## Requirements

- Python 2.7 (3.x expiremental)
- Numpy
- Matplotlib
- AMUSE (https://github.com/amusecode/amuse)
  - FIMap
  - SSE (optional)
- Astropy (or Pyfits)
