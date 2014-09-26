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

    # align and cluster 

    if not os.path.exists(options["runDir"]):
        cmd = "mkdir %s" % options["runDir"]
        print cmd
        runit(cmd)

    #### compute distances with distjob.usecols identified positions
    if not os.path.exists("%s/aac.msaToDist" % options["runDir"]):
        # get the size of the msa from the aac.msa.info file
        dat = runit("cat %s/aac.msa.info" % options["runDir"]).strip()
        ff = dat.split(" ")
        mycol = ff[14]
        myrow = ff[16]
        myrowhalf = max(32,int(myrow)/100) 
        if options["doOverlap"]=="1":
            cmd = "cd %s; agreeFracDist aac.msa %s %s %d %s %d overlap > aac.msaToDist 2>agreeFrac.err" % (options["runDir"], myrow, mycol, myrowhalf, options["entropyThreshold"], 4) # 4 is the maxInsert size to identify match columns
        else:
            cmd = "cd %s; agreeFracDist aac.msa %s %s %d %s %d > aac.msaToDist 2>agreeFrac.err" % (options["runDir"], myrow, mycol, myrowhalf, options["entropyThreshold"], 4) # 4 is the maxInsert size to identify match columns

        cmdFile = "%s/agreeFrac.sh" % options["runDir"]
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "agreeFrac.sh" )

    sys.stderr.write("got aac.msaToDist\n")

    #### Compute the hclust and plot it in R
    if not os.path.exists("%s/Rplots.pdf" % options["runDir"]) or not os.path.exists("%s/.RData" % options["runDir"]):
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
pdf(inputs[3])
#png(inputs[3],width=800,height=800,type="cairo")
plot(myhc,main="")
#abline(h=0.6,col="red")
dev.off()

# break the reads ids into cluster groups
# compute threshold based on inputs
############
# goal: compute clustering cutoff bound based on simple binomial
# get numPos from number of rows in distjob.usecols TODO: compute first before clustering
# get totalPos from aac.msa.info columns= [-3] or wc
# get numObjects from aac.msa.info rows= [-1] or wc
#

obsErr = function(totalPos){
  # given we take most variable numPos positions with mean=0.26 (raw pairwise) error, what err to we get?
  f = function(x){ (1.0 - exp(totalPos*pnorm(x,mean=0.26,sd=0.08,log=T))) -exp(-1)}
  uniroot( f , interval=c(0,1))$root
}

thresholdGivenPositions=function(numPos){
  useErr = obsErr(totalPos)
  k=seq(0,numPos)
  dist = pbinom(k,size=numPos,prob=useErr)
  ur = uniroot(approxfun(k/numPos,dist-(1-2*0.05/(numObj*numObj))),interval=c(0,1))
  return(ur$root)
}

numPos = as.numeric(scan(pipe("cat distjob.usecols | wc -l")))
totalPos = as.numeric(scan(pipe("cat aac.col | wc -l")))
numObj = as.numeric(scan(pipe("cat aac.id | wc -l")))

mythresh = thresholdGivenPositions(numPos)
write(mythresh, "clusterThreshold.txt")
write.table(data.frame(mythresh,numPos,totalPos,numObj, obsErr(totalPos), 2*0.05/(numObj*numObj)), "clusterThreshold.support",row.names=F)
############

if (mythresh<0.99){

myct = cutree(myhc, h=mythresh)

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
        fp.write("Rplots.pdf\n")
        fp.close()
        
        cmd = "cd %s; R CMD BATCH Rscript.R" % (options["runDir"])
        cmdFile = "%s/Rjob.sh" % options["runDir"]
        fp = open(cmdFile,"w")
        fp.write("%s\n" % cmd)
        fp.close()
        qsubWait( options["runDir"], "Rjob.sh" )

    sys.stderr.write("got Rplots.pdf\n")

if __name__ == "__main__":
    parser = OptionParser("agreeFracCluster")
    parser.add_option("--runDir", type="string", dest="runDir", help="where to store results. ac_mix3_sabin1_24500100-0045/")
    parser.add_option("--entropyThreshold", type="string", dest="entropyThreshold", help="Minimum entropy a MSA column must have to be included in distance. 1.0")
    parser.add_option("--doOverlap", type="string", dest="doOverlap", help="compute distances only on overlapping interval 1=yes 0=no. 0")

    (options, args) = parser.parse_args()

    if not options.runDir:
        parser.print_help()
        sys.exit(1)

    if not options.doOverlap:
        options.doOverlap = "0"

    alignAndCluster(**options.__dict__) # object to dict for kwargs, TODO: i guess this is right
