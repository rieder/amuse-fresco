# FRESCO
Fabricate Reality-Emulating Star Cluster Observations

Fresco requires AMUSE (http://amusecode.org; https://github.com/amusecode/amuse) and Matplotlib to run. It creates a Hubble WFC3-like image from an AMUSE-type hdf5 file of stars (other filetypes, like the Starlab format, are also supported). The positions and masses of the stars are always used, if information on the temperature and radius is not available this is calculated using the SSE stellar evolution code.
Gas particles may also be read, and will cause extinction and reflection of light from background- and nearby stars, respectively.

The current version is an initial but working version.

![Example image](test.png)
