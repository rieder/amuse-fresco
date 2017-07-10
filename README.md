# FRESCO
Fabricate Reality-Emulating Star Cluster Observations

Fresco (currently: make_cluster_image.py) requires AMUSE (http://amusecode.org; https://github.com/amusecode/amuse) and Matplotlib to run. It creates a Hubble WFC3-like image from an AMUSE-type hdf5 file of stars. The positions and masses of the stars are always used, if information on the temperature and radius is not available this is calculated using the SSE stellar evolution code.

The current version is an initial but working version.

![Example image](test.png)
