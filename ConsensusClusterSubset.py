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
    if not (options["subids"] == None): options["subids"] = os.path.abspath(options["subids"])

    if not os.path.exists(options["runDir"]):
        cmd = "mkdir %s" % options["runDir"]
        print cmd
        runit(cmd)

    sys.stderr.write("LOG %s : running ConsensusClusterSubset %s\n" % (str(datetime.datetime.now()), str(options.items())))
    fp = open("%s/options.items" % options["runDir"],"w")
    for (k,v) in options.items():
        fp.write("%s\t%s\n" % (k,v))
    fp.close()

    ################################
    # possibly subset reads
    inputfasta= options["fasta"]
    if not (options["subids"] == None):
        # does the subset id file exist?
        if not os.path.exists(options["subids"]):
            sys.stderr.write("ERROR: the given subset id file does not exist. Exiting. %s\n" % options["subids"])
            sys.exit(1)

        inputfasta = "%s/%s.fasta" % (options["runDir"], os.path.basename(options["subids"]))
        if not os.path.exists(inputfasta):
            # cmd = "export SEYMOUR_HOME=%s; . $SEYMOUR_HOME/etc/setup.sh; fastafetch --fasta %s --index %s -F -q %s > %s" % (os.environ['SEYMOUR_HOME'],options["fasta"],indexfile,options["subids"],inputfasta)
            cmd = "fastasub.py %s %s > %s" % (options["fasta"],options["subids"],inputfasta)
            runit(cmd)
    else:
        if options["CCS"]=="0":
            sys.stderr.write("LOG: using full fasta file with no subset\n")
        else:
            # pbalign won't parse CCS fasta ids from reads_of_insert! ERROR. Could not parse title m131210_065552_sherri_c100581732550000001823093904021452_s1_p0/10/ccs
            inputfasta = "%s/all.fasta" % (options["runDir"])
            # NOTE: the "0_1" is important, otherwise pbalign will refuse to align the fasta!!!!!
            cmd = "cat %s | sed 's/ccs$/0_1/' > %s" % (options["fasta"],inputfasta)
            runit(cmd)

    ################################
    # estimate quiver consensus
    if not os.path.exists("%s/quiver.done" % options["runDir"]):
        if options["CCS"]=="0":
            cmd = "runQuiverFastaBas.py --runDir %s --fasta %s --ref %s --basfofn %s --nproc %s" % (options["runDir"],inputfasta,options["ref"],options["basfofn"],options["nproc"])
            dat = runit(cmd)
            sys.stderr.write(dat[0])
            sys.stderr.write("\n")
            sys.stderr.write(dat[1])
            sys.stderr.write("\n")
        else:
            # for CCS use reference unaltered without estimation
            cmd = "cp %s %s/quiverResult.consensus.fasta; touch %s/quiver.done" % (options["ref"], options["runDir"], options["runDir"])
            runit(cmd)

    ################################
    # cluster the reads and break into groups
    # this puts computation onto cluster, so just run locally
    if not os.path.exists("%s/aac.hclust.png" % options["runDir"]):
        # check to see of spanThreshold is a percentage %
        if options["spanThreshold"][-1]== "%":
            dat = runit("fastalength %s/quiverResult.consensus.fasta" % options["runDir"])
            reflength = int(dat[0].split(" ")[0])
            frac = float(options["spanThreshold"][:-1])/100.0
            newSpanThreshold = int(frac*reflength+0.5)
            sys.stderr.write("newSpanThreshold= %d = %d*%f\n" % (newSpanThreshold,reflength,frac))
            options["spanThreshold"] = str(newSpanThreshold)

        cmd = "null"
        if options["CCS"]=="0":
            cmd = "alignAndClusterMaxIns.py --runDir %s --limsID %s --ref %s/quiverResult.consensus.fasta --spanThreshold %s --entropyThreshold %s --nproc %s --doOverlap %s --CCS %s" % (options["runDir"],options["basfofn"],options["runDir"],options["spanThreshold"],options["entropyThreshold"],options["nproc"],options["doOverlap"],options["CCS"])
        else:
            cmd = "alignAndClusterMaxIns.py --runDir %s --limsID %s --ref %s/quiverResult.consensus.fasta --spanThreshold %s --entropyThreshold %s --nproc %s --doOverlap %s --CCS %s" % (options["runDir"],inputfasta,options["runDir"],options["spanThreshold"],options["entropyThreshold"],options["nproc"],options["doOverlap"],options["CCS"])

        dat = runit(cmd)
        sys.stderr.write(dat[0])
        sys.stderr.write("\n")
        sys.stderr.write(dat[1])
        sys.stderr.write("\n")

    ################################
    # all done file
    runit("touch %s/all.done" % options["runDir"])

################################

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. hivqwf_2450154-0017/")
    parser.add_option("--fasta", type="string", dest="fasta", help="the fasta files containing HIV genome reads")
    parser.add_option("--subids", type="string", dest="subids", help="a file of newline separated read ids to use")
    parser.add_option("--ref", type="string", dest="ref", help="the generic reference. HIVemory.fasta")
    parser.add_option("--spanThreshold", type="string", dest="spanThreshold", help="how much of the reference must be spanned in order to keep read =6400 or percentage of ref=99.2%")
    parser.add_option("--entropyThreshold", type="string", dest="entropyThreshold", help="for clustering the minimum entropy needed in a column to be kept =1.0")
    parser.add_option("--basfofn", type="string", dest="basfofn", help="the bas.h5 fofn. HIV.bash5.fofn")
    parser.add_option("--nproc", type="string", dest="nproc", help="the number of processors to use when computing alignments. 1")
    parser.add_option("--doOverlap", type="string", dest="doOverlap", help="compute distances only on overlapping interval 1=yes 0=no. 0")
    parser.add_option("--CCS", type="string", dest="CCS", help="Use CCS fasta reads rather than raw from bas.h5. 1=yes 0=no. 0")

    (options, args) = parser.parse_args()

    if not options.runDir:
        print "ConsensusClusterSubset.py:"
        parser.print_help()
        sys.exit(1)

    if not options.nproc:
        options.nproc = "1"

    if not options.doOverlap:
        options.doOverlap = "0"

    if not options.CCS:
        options.CCS = "0"

    ConsensusClusterSubset(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
