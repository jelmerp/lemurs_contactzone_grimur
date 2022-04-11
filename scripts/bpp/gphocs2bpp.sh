#!/bin/bash

INPUT_GPHOCS=results/gphocs/input/hz.mur3gri2c.gphocsInput.txt
INPUT_BPP=results/bpp/input/mur3gri2c.txt

sed 's/.*fa//' $INPUT_GPHOCS > $INPUT_BPP
