#!/bin/sh
python test_set.py
python ../fresco.py -s stars.amuse -o test_stars
python ../fresco.py -g gas.amuse -o test_gas
