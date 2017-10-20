#!/bin/sh
python test_set.py
python ../fresco.py -s stars.hdf5 -o test_stars
python ../fresco.py -g gas.hdf5 -o test_gas
python ../fresco.py -s stars.hdf5 -g gas.hdf5 -o test_starsgas
