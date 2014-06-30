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
    # NOTE: shell inherits environment!

################################

def runQuiverFastaBas(*argv, **options):

    options["runDir"] = os.path.abspath(options["runDir"])
    options["fasta"] = os.path.abspath(options["fasta"])
    options["ref"] = os.path.abspath(options["ref"])

    if not os.path.exists(options["runDir"]):
        cmd = "mkdir %s" % options["runDir"]
        print cmd
        runit(cmd)

    inputfasta = options["fasta"]

    ################################
    # align stratified reads against generic HIV reference
    # inputfasta = the reads
    # outh5 = output cmp.h5

    outh5 = "%s/hivalign.cmp.h5" % (options["runDir"])
        
    if not os.path.exists("%s/align.done" % options["runDir"]):

        # use pbalign rather than compareSequences.py which is
        # inconsistent.  to do this, I have to get a whitelist of
        # fasta ids to use as pbalign must point directly to the
        # bas.h5 for quiver

        template= """cd %s;
export SEYMOUR_HOME=%s;
source $SEYMOUR_HOME/etc/setup.sh;
export PATH=%s
"""
        cmd = template % (options["runDir"],os.environ['SEYMOUR_HOME'],os.environ['PATH'])

        # get list of fasta ids into infastaWhitelist.txt
        dat = open(options["fasta"]).read().splitlines()
        ofp = open("%s/infastaWhitelist.txt" % options["runDir"], "w")
        for ll in dat:
            if ll[0]==">":
                # output fasta id with run/zmw ... and not /subread
                ff = ll[1:].split("/")
                ofp.write("%s\n" % "/".join(ff[:2]))
        ofp.close()
          
        # get region table for input into pbalign
        template= """filter_plsh5.py %s \
--outputDir=%s \
--outputFofn=%s \
--outputSummary=filter_plsh5.outputSummary \
--debug \
--logFile=filter_plsh5.log \
--trim=False \
--filter="ReadWhitelist=%s"\n""" 
        cmd = cmd + template % (options["basfofn"], options["runDir"], "infastaWhitelist.filter.fofn", "infastaWhitelist.txt")

        # run the alignment
        templateList = ["pbalign.py %s \\" % options["basfofn"],
                        "%s \\" % options["ref"],
                        "%s \\" % outh5,
                        "--regionTable=infastaWhitelist.filter.fofn \\",
                        "--forQuiver \\",
                        "--hitPolicy=allbest \\",
                        "--nproc=%s \\" % options["nproc"],
                        "--seed=1234 \\",
                        "--minLength 2048 \\",
                        "--minAccuracy 10 \\",
                        "--maxDivergence 90 \\",
                        "--noSplitSubreads \\",
                        "--maxHits 1",
                        "\n"]
        cmd = cmd + "\n".join(templateList)

        fp = open("%s/alignments.cmd" % options["runDir"],"w")
        fp.write("%s\n" % cmd)
        fp.close()

        # run locally rather than qsub
        cmd = "cd %s; chmod 777 alignments.cmd" % (options["runDir"])
        print cmd
        runit(cmd)
        cmd = "source %s/alignments.cmd" % (options["runDir"])
        print cmd
        print runit(cmd)

        cmd = "touch %s/align.done" % options["runDir"]
        runit(cmd)

        sys.stderr.write("got align.done\n")        

    ################################
    # Compute quiver consensus
    if not os.path.exists("%s/quiver.done" % options["runDir"]):

        outvar = "%s/quiverResult" % (options["runDir"])

        template= """cd %s;
export SEYMOUR_HOME=%s;
source $SEYMOUR_HOME/etc/setup.sh;
export PATH=%s
"""
        cmd = template % (options["runDir"],os.environ['SEYMOUR_HOME'],os.environ['PATH'])


        template = """variantCaller.py \
-vv  \
-j8 --algorithm=quiver \
%s \
-r %s \
-o %s.gff -o %s.consensus.fasta -o %s.consensus.fastq \
> %s.vc.out 2> %s.vc.err
samtools faidx %s.consensus.fasta
"""

        # --parameter=AllQVsModel.C2

        cmd = cmd+template % (outh5, options["ref"], outvar, outvar, outvar, outvar, outvar, outvar)

        fp = open("%s/var.cmd" % options["runDir"],"w")
        fp.write("%s\n" % cmd)
        fp.close()

        # run locally rather than qsub
        cmd = "cd %s; chmod 777 var.cmd" % (options["runDir"])
        print cmd
        runit(cmd)
        cmd = "source %s/var.cmd" % (options["runDir"])
        print cmd
        runit(cmd)

        cmd = "touch %s/quiver.done" % options["runDir"]
        runit(cmd)

        sys.stderr.write("got quiver.done\n")        

################################

if __name__ == "__main__":
    parser = OptionParser("hivBarcodeWorkflow")
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. hivqwf_2450154-0017/")
    parser.add_option("--fasta", type="string", dest="fasta", help="the fasta files containing HIV genome reads")
    parser.add_option("--ref", type="string", dest="ref", help="the generic reference. HIVemory.fasta")
    parser.add_option("--basfofn", type="string", dest="basfofn", help="the bas.h5 list. HIV.bash5.fofn")
    parser.add_option("--nproc", type="string", dest="nproc", help="the number of processors to use when computing alignments. 1")

    (options, args) = parser.parse_args()

    if not options.runDir:
        print "HIV quiver workflow:"
        parser.print_help()
        sys.exit(1)

    if not options.nproc:
        options.nproc = "1"

    runQuiverFastaBas(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
