#!/usr/bin/env python
__doc__ = """clustering consensus for mixed populations

ConsensusClusterSubset.py \
--runDir runinittest_FULL_mix159_mple17_mplf5 \
--fasta runinittest_mix159_mple17_mplf5.fasta \
--subids runinittest_mix159_mple17_mplf5.group1.ids \
--ref quiverResult.consensus.fasta \
--spanThreshold=6400 \
--entropyThreshold=1.0

Work:

- input: fasta file and optional subset list of ids

- subset the reads if necessary

- estimate a quiver consensus on the reads

- cluster the reads using the quiver consensus

- break the read ids into groups based on cluster

"""

import sys
from optparse import OptionParser
import subprocess
import os
import re
import glob
import datetime

################################
def runit( cmd ):
    sys.stderr.write("runit %s\n" % cmd)
    return subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='/bin/bash').communicate()

################################
def ConsensusClusterSubset(*argv, **options):

    options["runDir"] = os.path.abspath(options["runDir"])
    options["fasta"] = os.path.abspath(options["fasta"])
    options["ref"] = os.path.abspath(options["ref"])
    options["basfofn"] = os.path.abspath(options["basfofn"])

    sys.stderr.write("LOG %s : running ConsensusClusterSubset %s\n" % (str(datetime.datetime.now()), str(options.items())))

    if not os.path.exists(options["runDir"]):
        cmd = "mkdir %s" % options["runDir"]
        print cmd
        runit(cmd)

    ################################
    # possibly subset reads
    inputfasta= options["fasta"]
    if not options["subids"] == None:
        options["subids"] = os.path.abspath(options["subids"])

        # does the subset id file exist?
        if not os.path.exists(options["subids"]):
            sys.stderr.write("ERROR: the given subset id file does not exist. Exiting. %s\n" % options["subids"])
            sys.exit(1)

        # generate the index for subsetting
        indexfile = options["fasta"]+".index"
        if not os.path.exists(indexfile):
            cmd = "export SEYMOUR_HOME=/mnt/secondary/Smrtpipe/builds/smrtpipe_v2.116526; . $SEYMOUR_HOME/etc/setup.sh; fastaindex --fasta %s --index %s" % (options["fasta"],indexfile)
            runit(cmd)

        inputfasta = "%s/%s.fasta" % (options["runDir"], os.path.basename(options["subids"]))
        if not os.path.exists(inputfasta):
            cmd = "export SEYMOUR_HOME=/mnt/secondary/Smrtpipe/builds/smrtpipe_v2.116526; . $SEYMOUR_HOME/etc/setup.sh; fastafetch --fasta %s --index %s -F -q %s > %s" % (options["fasta"],indexfile,options["subids"],inputfasta)
            runit(cmd)
    else:
        sys.stderr.write("LOG: using full fasta file with no subset\n")

    ################################
    # estimate quiver consensus
    # this puts computation onto cluster, so just run locally
    if not os.path.exists("%s/quiver.done" % options["runDir"]):
        cmd = "runQuiverFastaBas.py --runDir %s --fasta %s --ref %s --basfofn %s" % (options["runDir"],inputfasta,options["ref"],options["basfofn"])
        dat = runit(cmd)
        sys.stderr.write(dat[0])
        sys.stderr.write("\n")
        sys.stderr.write(dat[1])
        sys.stderr.write("\n")

    ################################
    # cluster the reads and break into groups
    # this puts computation onto cluster, so just run locally
    if not os.path.exists("%s/aac.hclust.png" % options["runDir"]):
        cmd = "alignAndClusterMaxIns.py --runDir %s --limsID %s --ref %s/quiverResult.consensus.fasta --spanThreshold %s --entropyThreshold %s" % (options["runDir"],inputfasta,options["runDir"],options["spanThreshold"],options["entropyThreshold"])
        dat = runit(cmd)
        sys.stderr.write(dat[0])
        sys.stderr.write("\n")
        sys.stderr.write(dat[1])
        sys.stderr.write("\n")

################################

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. hivqwf_2450154-0017/")
    parser.add_option("--fasta", type="string", dest="fasta", help="the fasta files containing HIV genome reads")
    parser.add_option("--subids", type="string", dest="subids", help="a file of newline separated read ids to use")
    parser.add_option("--ref", type="string", dest="ref", help="the generic reference. HIVemory.fasta")
    parser.add_option("--spanThreshold", type="string", dest="spanThreshold", help="how much of the reference must be spanned in order to keep read =6400")
    parser.add_option("--entropyThreshold", type="string", dest="entropyThreshold", help="for clustering the minimum entropy needed in a column to be kept =1.0")
    parser.add_option("--basfofn", type="string", dest="basfofn", help="the bas.h5 fofn. HIV.bash5.fofn")

    (options, args) = parser.parse_args()

    if not options.runDir:
        print "ConsensusClusterSubset.py:"
        parser.print_help()
        sys.exit(1)

    ConsensusClusterSubset(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
