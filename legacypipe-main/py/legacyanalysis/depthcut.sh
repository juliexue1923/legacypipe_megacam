#!/bin/bash

export PYTHONPATH="${PYTHONPATH}:/central/home/zxue/unwise_psf-master/py"
export PYTHONPATH="${PYTHONPATH}:$HOME/astrometry_install/lib/python"
export LEGACY_SURVEY_DIR=$HOME/jet/megacam-dir
export GAIA_CAT_DIR=$HOME/jet/gaia_hp
export DELVE_CAT_DIR=$HOME/jet/delve_hp
export DUST_DIR=$HOME/sfdmaps
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/lib

#export LEGACY_SURVEY_DIR=/global/projecta/projectdirs/cosmo/work/legacysurvey/dr8/DECaLS/
export LEGACY_SURVEY_DIR=~/jet/megacam-dir

outdir=~/jet/megacam-dir

brick="$1"

bri=$(echo $brick | head -c 3)
echo $brick
mkdir -p $outdir/logs/$bri

# Shifter
# cd /src/legacypipe/py

python legacyanalysis/depth-cut.py --outdir $outdir --margin 1 $brick > $outdir/logs/$bri/$brick.log 2>&1
