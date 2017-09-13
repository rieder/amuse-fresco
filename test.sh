#!/bin/sh
python test_set.py
python fresco.py -s stars.amuse -o test_stars.png
python fresco.py -g gas.amuse -o test_gas.png
