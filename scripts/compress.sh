#!/bin/bash
#SBATCH -p CN_BIOT
#SBATCH --ntasks=1
cmpress Rfam/covariance-models/concatenated/Rfam.cm > compress.log 2>&1
