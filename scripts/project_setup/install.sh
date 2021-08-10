#!/bin/bash

conda env create --force --file envs/admixture.yml
conda env create --force --file envs/plink.yml
conda env create --force --file envs/vcftools.yml
