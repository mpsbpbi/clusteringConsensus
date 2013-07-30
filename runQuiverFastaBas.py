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
    return subprocess.Popen( cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]

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

        template= """cd %s;
export SEYMOUR_HOME=%s;
. $SEYMOUR_HOME/etc/setup.sh;
"""
        cmd = template % (options["runDir"],os.environ['SEYMOUR_HOME'])

        template= """compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=%s  --noXML --h5mode=w \
--h5fn=%s \
-x -bestn 1 \
--debug \
%s \
%s
sleep 2
"""
        toadd = template % (options["nproc"],outh5,inputfasta,options["ref"])
        cmd = cmd+toadd
        
        template="""loadPulses %s \
%s \
-metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread
sleep 2
"""
        toadd = template % (options["basfofn"],outh5)
        cmd = cmd+toadd

        # -d debug option removed!!
        template= """cmph5tools.py sort --deep --inPlace %s
"""
        toadd = template % (outh5)
        cmd = cmd+toadd

        fp = open("%s/alignments.cmd" % options["runDir"],"w")
        fp.write("%s\n" % cmd)
        fp.close()

        # run locally rather than qsub
        cmd = "cd %s; chmod 777 alignments.cmd" % (options["runDir"])
        print cmd
        runit(cmd)
        cmd = "bash %s/alignments.cmd" % (options["runDir"])
        print cmd
        runit(cmd)

        cmd = "touch %s/align.done" % options["runDir"]
        runit(cmd)

        sys.stderr.write("got align.done\n")        

    ################################
    # Compute quiver consensus
    if not os.path.exists("%s/quiver.done" % options["runDir"]):

        outvar = "%s/quiverResult" % (options["runDir"])

        template = """variantCaller.py \
-vv  \
-j8 --algorithm=quiver \
%s \
-r %s \
-o %s.gff -o %s.consensus.fasta -o %s.consensus.fastq \
> %s.vc.out 2> %s.vc.err
"""

        # --parameter=AllQVsModel.C2

        cmd = template % (outh5, options["ref"], outvar, outvar, outvar, outvar, outvar)

        fp = open("%s/var.cmd" % options["runDir"],"w")
        fp.write("%s\n" % cmd)
        fp.close()

        # run locally rather than qsub
        cmd = "cd %s; chmod 777 var.cmd" % (options["runDir"])
        print cmd
        runit(cmd)
        cmd = "bash %s/var.cmd" % (options["runDir"])
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
