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

################################
# fakes running qsub so it's not needed
def qsubWait( runDir, cmdFile ):

    cmd = "cd %s; chmod 777 %s" % (runDir, cmdFile)
    runit(cmd)
    cmd = "bash %s/%s" % (runDir,cmdFile)
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

        cmdTemp = """#!/bin/bash
export SEYMOUR_HOME=%s;
. $SEYMOUR_HOME/etc/setup.sh;

cd %s

compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=%s  --noXML --h5mode=w \
--h5fn=%s \
-x -bestn 1 \
--tmpDir=/scratch --debug \
%s \
%s

errFromCmph5.py %s | sort -n -k 2 > %s.error
"""

        cmd  = cmdTemp % (os.environ['SEYMOUR_HOME'], options["runDir"], options["nproc"], "alignments.cmp.h5", options["limsID"], options["ref"], "alignments.cmp.h5", "alignments.cmp.h5")

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

    #sys.exit(0)

    #### cmph5ToMSA with filtered reads
    if not os.path.exists("%s/aac.msa" % options["runDir"]):
        # get ref file from ref repository
        # OLD: reffile = glob.glob("%s/*.fa" % options["ref"])[0]
        reffile = options["ref"]

        cmd = "cd %s; export SEYMOUR_HOME=%s; source $SEYMOUR_HOME/etc/setup.sh; cmph5ToMSAMaxInserts.py --cmph5=alignments.cmp.h5 --reffile=%s --msafile=aac.msa --idfile=aac.id --colfile=aac.col --rangeLow=0 --rangeHigh=999999 --seqLow=0 --seqHigh=999999 --filter=alignments.filterFull --maxInsert=4 >msajob.stout 2>msajob.stderr" % (options["runDir"],os.environ['SEYMOUR_HOME'],reffile)
        cmdFile = "%s/msajob.sh" % (options["runDir"])
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "msajob.sh" )
    sys.stderr.write("got aac.msa\n")

    #sys.exit(0)

    #### compute distances with 1/2 read coverage and entropy
    if not os.path.exists("%s/aac.msaToDist" % options["runDir"]):
        # get the size of the msa from the aac.msa.info file
        dat = runit("cat %s/aac.msa.info" % options["runDir"]).strip()
        ff = dat.split(" ")
        mycol = ff[14]
        myrow = ff[16]
        myrowhalf = int(myrow)/2 # TODO: this was /2 or /8 for smaller
        cmd = "cd %s; pairwiseAlignDist aac.msa %s %s %d %s %d > aac.msaToDist 2>distjob.err" % (options["runDir"], myrow, mycol, myrowhalf, options["entropyThreshold"], 4) # 4 is the maxInsert size to identify match columns

        cmdFile = "%s/distjob.sh" % options["runDir"]
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "distjob.sh" )

    sys.stderr.write("got aac.msaToDist\n")

    #### Compute the hclust and plot it in R
    if not os.path.exists("%s/aac.hclust.png" % options["runDir"]) or not os.path.exists("%s/.RData" % options["runDir"]):
        Rscript= """
inputs = scan("Rscript.input",what="character")
# num objects
# msaToDist file
# plot file

num = as.numeric(inputs[1])
mydist = matrix(nrow=num, ncol=num)
dd = read.table(inputs[2], head=F,sep="\t")
for (ii in 1:dim(dd)[1]){
  myx = dd$V1[ii]+1
  myy = dd$V2[ii]+1
  myv = dd$V3[ii]
  mydist[myx,myy]=myv
  mydist[myy,myx]=myv
}
for (ii in 1:num){
  mydist[ii,ii]=0.0
}
myhc = hclust(as.dist(mydist))
png(inputs[3],width=800,height=800,type="cairo")
plot(myhc,main="")
#abline(h=0.6,col="red")
dev.off()

# break the reads ids into cluster groups
# this is not the right threshold. TODO: optimize cut points
myct = cutree(myhc, h=0.75)

myid = read.table("aac.id", head=F,sep="\\t")

# the first sequence is always the reference the way I'm running
# it. And get rid of the subread information

jj = cbind(as.character(gsub("\\\\/[0-9_]+$","",myid$V1[2:length(myct)])), myct[2:length(myct)])

writegroup = function(group,name){
  file = paste("clustergroup_",name,".ids",sep="")
  cat(jj[ jj[,2]==group, 1], file=file, sep="\\n")
}

bysize = names(sort(table(jj[,2]),decr=T))
for (ii in 1:length(bysize)){
  writegroup(bysize[ii],ii)
}

"""
        fp = open("%s/Rscript.R" % options["runDir"],"w")
        fp.write(Rscript)
        fp.close()

        dat = runit("cat %s/aac.msa.info" % options["runDir"]).strip()
        ff = dat.split(" ")
        myrow = ff[16]

        fp = open("%s/Rscript.input" % options["runDir"],"w")
        fp.write("%s\n" % myrow)
        fp.write("aac.msaToDist\n")
        fp.write("aac.hclust.png\n")
        fp.close()
        
        cmd = "cd %s; R CMD BATCH Rscript.R" % (options["runDir"])
        cmdFile = "%s/Rjob.sh" % options["runDir"]
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "Rjob.sh" )

    sys.stderr.write("got aac.hclust.png\n")

if __name__ == "__main__":
    parser = OptionParser("alignAndCluster")
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. ac_mix3_sabin1_24500100-0045/")
    parser.add_option("--ref", type="string", dest="ref", help="The fasta reference. sabin1")
    parser.add_option("--limsID", type="string", dest="limsID", help="The LIMS ID to run. 24500100-0045")
    parser.add_option("--spanThreshold", type="string", dest="spanThreshold", help="How much of the gnome each read must span to be kept. 5500")
    parser.add_option("--entropyThreshold", type="string", dest="entropyThreshold", help="Minimum entropy a MSA column must have to be included in distance. 1.0")
    parser.add_option("--nproc", type="string", dest="nproc", help="the number of processors to use when computing alignments. 1")

    (options, args) = parser.parse_args()

    if not options.runDir:
        print "Align Run Against Reference and Cluster:"
        parser.print_help()
        sys.exit(1)

    if not options.nproc:
        options.nproc = "1"

    alignAndCluster(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
