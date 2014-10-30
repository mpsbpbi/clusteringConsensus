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
# fakes running qsub so it's not needed
def qsubWait( runDir, cmdFile ):

    cmd = "cd %s; chmod 777 %s" % (runDir, cmdFile)
    runit(cmd)
    cmd = "source %s/%s" % (runDir,cmdFile)
    print cmd
    runit(cmd)

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

    sys.stderr.write("got quiver.done\n")

    ################################
    # align reads to quiver consensus, take only those spanning, produce MSA
    if not os.path.exists("%s/aac.msa" % options["runDir"]):
        # check to see of spanThreshold is a percentage %
        if options["spanThreshold"][-1]== "%":
            dat = runit("fastalength %s/quiverResult.consensus.fasta" % options["runDir"])
            reflength = int(dat[0].split(" ")[0])
            frac = float(options["spanThreshold"][:-1])/100.0
            newSpanThreshold = int(frac*reflength+0.5)
            sys.stderr.write("newSpanThreshold= %d = %d*%f\n" % (newSpanThreshold,reflength,frac))
            options["spanThreshold"] = str(newSpanThreshold)
            # TODO: write updated options.items

        if options["CCS"]=="0":
            inputseq = options["basfofn"]
        else:
            inputseq = inputfasta

        cmd = "msaAlignSpanning.py --runDir %s --inputseq %s --ref %s/quiverResult.consensus.fasta --spanThreshold %s --nproc %s --useQuality %s" % (options["runDir"],inputseq,options["runDir"],options["spanThreshold"],options["nproc"],options["useQuality"])

        dat = runit(cmd)
        sys.stderr.write(dat[0])
        sys.stderr.write("\n")
        sys.stderr.write(dat[1])
        sys.stderr.write("\n")

    sys.stderr.write("got aac.msa\n")

    ################################
    # identify variant positions using either entropy or basis chisq
    if not os.path.exists("%s/distjob.usecols" % options["runDir"]):

        ########
        if options["chisqThreshold"] != None:
            # compute using chisq
            
            cmd = "echo basisVar\n"
            ## get the data into HP-region count form
            if not os.path.exists("%s/msaobs-hp-set.counts" % options["runDir"]):
                cmd = cmd + "cd %s; msaobs-hp-set.py aac.msa all > msaobs-hp-set.counts\n" % (options["runDir"])

            # use basis to find minor variants across all
            if not os.path.exists("%s/distjob.usecols" % options["runDir"]):
                cmd = cmd + "cd %s; source /home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/virtscipy/bin/activate; cat msaobs-hp-set.counts | minorMsaObs.py %s %s %s\n" % (options["runDir"], options["chisqThreshold"], options["fisherlThreshold"], options["surround"])

            fp = open("%s/basisVar.cmd" % options["runDir"],"w")
            fp.write("%s\n" % cmd)
            fp.close()
            qsubWait( options["runDir"], "basisVar.cmd" )

        ########
        if options["entropyThreshold"] != None:
            # compute using entropy

            # get the size of the msa from the aac.msa.info file
            dat = runit("cat %s/aac.msa.info" % options["runDir"])[0].strip()
            ff = dat.split(" ")
            mycol = ff[14]
            myrow = ff[16]
            myrowhalf = max(32,int(myrow)/100) # Require 1/100 to be non-empty with minimum of 32

            if options["doOverlap"]=="1":
                cmd = "cd %s; entropyVariants aac.msa %s %s %d %s %d overlap > entropyVariants.stdout 2>entropyVariants.stderr" % (options["runDir"], myrow, mycol, myrowhalf, options["entropyThreshold"], 4) # 4 is the maxInsert size to identify match columns
            else:
                cmd = "cd %s; entropyVariants aac.msa %s %s %d %s %d > entropyVariants.stdout 2>entropyVariants.stderr" % (options["runDir"], myrow, mycol, myrowhalf, options["entropyThreshold"], 4) # 4 is the maxInsert size to identify match columns

            cmdFile = "%s/distjob.sh" % options["runDir"]
            fp = open(cmdFile,"w")
            fp.write("%s\n" % cmd)
            fp.close()
            qsubWait( options["runDir"], "distjob.sh" )

    sys.stderr.write("got distjob.usecols\n")

    ################################

    # cluster / phase the reads based on the identified variant
    # positions using either complete linkage single base
    # agreeFraction, HP-region basis clustering, or HP-region basis
    # phasing. Yields either sets of stratified IDs or possibly phased
    # genomes if clustering not clear
    # TODO:
    if not os.path.exists("%s/cluster.done" % options["runDir"]):

        ################################
        if options["clusterMethod"]=="agreeFracCluster":
            cmd = "agreeFracCluster.py --runDir %s --entropyThreshold %s --doOverlap %s" % (options["runDir"],options["entropyThreshold"],options["doOverlap"])

            dat = runit(cmd)
            sys.stderr.write(dat[0])
            sys.stderr.write("\n")
            sys.stderr.write(dat[1])
            sys.stderr.write("\n")

        ################################
        if options["clusterMethod"]=="basisPhase":

            cmd = "echo basisPhase\n"
            # first get compound variant counts
            if not os.path.exists("%s/basisPhase-msaobs-hp-set.counts" % options["runDir"]):
                cmd = cmd+ "cd %s; msaobs-hp-set.py aac.msa distjob.usecols collapse sanitize %s %s > basisPhase-msaobs-hp-set.counts\n" % (options["runDir"], options["surround"], options["addHPContext"])

            # do the phasing
            if not os.path.exists("%s/basisPhase.output" % options["runDir"]):
                cmd = cmd+"cd %s; source /home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/virtscipy/bin/activate; basisPhase.py basisPhase-msaobs-hp-set.counts %s > basisPhase.output\n" % (options["runDir"], options["basisPhaseThreshold"])

            fp = open("%s/basisPhase.cmd" % options["runDir"],"w")
            fp.write("%s\n" % cmd)
            fp.close()
            qsubWait( options["runDir"], "basisPhase.cmd" )

        cmd = "touch %s/cluster.done" % (options["ref"])

    sys.stderr.write("got cluster.done\n")

    ################################
    # summary and all done file
    runit("cd %s; resultsSummary.py > resultsSummary.html" % options["runDir"])
    runit("touch %s/all.done" % options["runDir"])

################################

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. hivqwf_2450154-0017/")
    parser.add_option("--fasta", type="string", dest="fasta", help="the fasta files containing HIV genome reads")
    parser.add_option("--subids", type="string", dest="subids", help="a file of newline separated read ids to use")
    parser.add_option("--ref", type="string", dest="ref", help="the generic reference. HIVemory.fasta")
    parser.add_option("--spanThreshold", type="string", dest="spanThreshold", help="how much of the reference must be spanned in order to keep read =6400 or percentage of ref=99.2%")
    parser.add_option("--entropyThreshold", type="string", dest="entropyThreshold", help="exclusive from chisqThreshold. for clustering the minimum entropy needed in a column to be kept =1.0")
    parser.add_option("--chisqThreshold", type="string", dest="chisqThreshold", help="exclusive from entropyThreshold. for variant positions the maximum p-value to be kept for chi stat=1.0e-100")
    parser.add_option("--fisherlThreshold", type="string", dest="fisherlThreshold", help="exclusive from entropyThreshold. for variant positions the maximum p-value to be kept for fisherl stat=1.0e-100")
    parser.add_option("--addHPContext", type="string", dest="addHPContext", help="for selected HP context add the HP context downstream. 'addHPConext'=yes ''=no. ''")
    parser.add_option("--clusterMethod", type="string", dest="clusterMethod", help="for clustering/phasing (agreeFracCluster, basisCluster, basisPhase) =basisCluster")
    parser.add_option("--basisPhaseThreshold", type="string", dest="basisPhaseThreshold", help="for basisPhase threshold to add new haplotype = 1.0E-10")
    parser.add_option("--basfofn", type="string", dest="basfofn", help="the bas.h5 fofn. HIV.bash5.fofn")
    parser.add_option("--nproc", type="string", dest="nproc", help="the number of processors to use when computing alignments. 1")
    parser.add_option("--doOverlap", type="string", dest="doOverlap", help="compute distances only on overlapping interval 1=yes 0=no. 0")
    parser.add_option("--CCS", type="string", dest="CCS", help="Use CCS fasta reads rather than raw from bas.h5. 1=yes 0=no. 0")
    parser.add_option("--useQuality", type="string", dest="useQuality", help="use Quality in alignment of reads against quiver. 1=yes 0=no. 0")
    parser.add_option("--surround", type="string", dest="surround", help="use surround'ing bases in minorMsaObs to id variant positions. 'surround'=yes ''=no. ''")

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

    if not options.useQuality:
        options.useQuality = "0"

    if not options.surround:
        options.surround=""
    else:
        options.surround="surround"

    if not options.chisqThreshold:
        options.useQuality = "1.0e-100"

    if not options.fisherlThreshold:
        options.useQuality = "1.0e-100"

    if not options.addHPContext:
        options.addHPContext=""
    else:
        options.addHPContext="addHPContext"

    ConsensusClusterSubset(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
