#!/usr/bin/env python
import sys
from optparse import OptionParser
import subprocess
import os
import time
import re
import glob

################################
def runit( cmd ):
    sys.stderr.write("runit %s\n" % cmd)
    return subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True, executable='/bin/bash').communicate()[0]
    # NOTE: shell inherits environment from the shell running python!

################################
# fakes running qsub so it's not needed
def qsubWait( runDir, cmdFile ):

    cmd = "cd %s; chmod 777 %s" % (runDir, cmdFile)
    runit(cmd)
    cmd = "source %s/%s" % (runDir,cmdFile)
    print cmd
    runit(cmd)

################################

def alignAndCluster(*argv, **options):

    options["runDir"] = os.path.abspath(options["runDir"])
    options["ref"] = os.path.abspath(options["ref"])

    # align and cluster 

    if not os.path.exists(options["runDir"]):
        cmd = "mkdir %s" % options["runDir"]
        print cmd
        runit(cmd)

    ## submit the alignment job
    if not os.path.exists("%s/clusalignments.cmd" % options["runDir"]):

        # use pbalign rather than compareSequences.py. I can use the
        # subsetting done when generating the sample Quiver reference
        # here without recomputing it

        template= """#!/bin/bash
cd %s;
export SEYMOUR_HOME=%s;
source $SEYMOUR_HOME/etc/setup.sh;
export PATH=%s
"""
        cmd = template % (options["runDir"],os.environ['SEYMOUR_HOME'],os.environ['PATH'])

        # run the alignment
        templateList = ["pbalign.py %s \\" % options["inputseq"],
                        "%s \\" % options["ref"],
                        "%s \\" % "alignments.cmp.h5",
                        "--regionTable=infastaWhitelist.filter.fofn \\",
                        "--hitPolicy=allbest \\",
                        "--nproc=%s \\" % options["nproc"],
                        "--seed=1234 \\",
                        "--minLength 2048 \\",
                        "--minAccuracy 10 \\",
                        "--maxDivergence 90 \\",
                        "--noSplitSubreads \\",
                        "--algorithmOptions \"-useQuality \" \\",
                        "--maxHits 1",
                        "\n",
                        "errFromCmph5.py %s | sort -n -k 2 > %s.error" % ("alignments.cmp.h5", "alignments.cmp.h5"),
                        "\n"]
        cmd = cmd + "\n".join(templateList)
        if options["useQuality"]=="0":
            cmd = cmd.replace("--algorithmOptions \"-useQuality \" \\\n","")

        # if CCS then run only the fasta CCS basecalls. TODO: this could be the .ccs.h5 to get quality
        if ".fasta" in options["inputseq"]:
            cmd = cmd.replace("--regionTable=infastaWhitelist.filter.fofn \\\n","")
            cmd = cmd.replace("--minLength 2048 \\\n","")

        fp = open("%s/clusalignments.cmd" % options["runDir"],"w")
        fp.write("%s\n" % cmd)
        fp.close()

        qsubWait( options["runDir"], "clusalignments.cmd" )

    sys.stderr.write("got clusalignments.cmd\n")
    sys.stderr.write("got .cmp.h5.error\n")

    if not os.path.exists("%s/alignments.filterFull" % options["runDir"]):
        cmd = "filterFull.py %s/alignments.cmp.h5.error %s > %s/alignments.filterFull" % (options["runDir"], options["spanThreshold"], options["runDir"])
        dat = runit(cmd)
    sys.stderr.write("got .filterFull\n")        

    #### cmph5ToMSA with filtered reads
    if not os.path.exists("%s/aac.msa" % options["runDir"]):
        reffile = options["ref"]

        cmd = "cd %s; export SEYMOUR_HOME=%s; source $SEYMOUR_HOME/etc/setup.sh; export PATH=%s ; cmph5ToMSAMaxInserts.py --cmph5=alignments.cmp.h5 --reffile=%s --msafile=aac.msa --idfile=aac.id --colfile=aac.col --rangeLow=0 --rangeHigh=999999 --seqLow=0 --seqHigh=999999 --filter=alignments.filterFull --maxInsert=4 >msajob.stout 2>msajob.stderr" % (options["runDir"],os.environ['SEYMOUR_HOME'],os.environ['PATH'],reffile)
        cmdFile = "%s/msajob.sh" % (options["runDir"])
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "msajob.sh" )
    sys.stderr.write("got aac.msa\n")

if __name__ == "__main__":
    parser = OptionParser("alignAndCluster")
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. ac_mix3_sabin1_24500100-0045/")
    parser.add_option("--inputseq", type="string", dest="inputseq", help="input sequences file. basfofn or fasta")
    parser.add_option("--ref", type="string", dest="ref", help="The fasta reference. sabin1")
    parser.add_option("--spanThreshold", type="string", dest="spanThreshold", help="How much of the gnome each read must span to be kept.")
    parser.add_option("--nproc", type="string", dest="nproc", help="the number of processors to use when computing alignments. 1")
    parser.add_option("--useQuality", type="string", dest="useQuality", help="use Quality in alignment of reads against quiver. 1=yes 0=no. 0")

    (options, args) = parser.parse_args()

    if not options.runDir:
        parser.print_help()
        sys.exit(1)

    if not options.nproc:
        options.nproc = "1"

    if not options.useQuality:
        options.useQuality = "0"

    alignAndCluster(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
