#!/usr/bin/env python
__doc__ = """Prepare qsub command to align against given Polio genome.

"""

import sys
import os
import glob

def main():
    limsid = sys.argv[1]
    outfile = sys.argv[2]
    refstrain = sys.argv[3]
    outdir = sys.argv[4]

    # use limsID -> as called is a fasta file
    bash5 = limsid

    cmdTemp = """#!/bin/bash
export LD_LIBRARY_PATH_OLD=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="" 
export SEYMOUR_HOME=/mnt/secondary/Smrtpipe/builds/smrtpipe_v2.116526;
. $SEYMOUR_HOME/etc/setup.sh;
export WORKFLOWS=/mnt/secondary/Smrtpipe/martin/8081/data/workflows/;

cd %s

compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=8  --noXML --h5mode=w \
--h5fn=%s \
-x -bestn 1 \
--tmpDir=/scratch --debug \
%s \
%s

errFromCmph5.py %s | sort -n -k 2 > %s.error
"""

    cmd  = cmdTemp % (outdir, outfile, bash5, refstrain, outfile, outfile)

    print cmd

main()
