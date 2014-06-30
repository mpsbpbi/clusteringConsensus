#!/usr/bin/env python
__doc__ = """Generate RESULT.html that summarizes all the results.
python /home/UNIXHOME/mbrown/mbrown/workspace2013Q1/pacbioCode-viral-clusteringConsensus-v1/code/resultsSummary.py > ~/resultsSummary.html
"""

import glob
import os

################################
mapFileToDescriptionDat = """options.items	initial program options	0
infastaWhitelist.txt	list of read ids to include in analysis	0
alignments.cmd	command line to generate initial alignments	0
*.rgn.h5	read whitelisting data files for bas/bax.h5 files	0
filter_plsh5.outputSummary	read whitelisting summary	0
infastaWhitelist.filter.fofn	read whitelisting files	0
filter_plsh5.log	read whitelisting log	0
hivalign.cmp.h5	initial alignments	0
align.done	inital alignment checkpoint file	0
var.cmd	command line to generate consensus	0
quiverResult.vc.out	quiver stdout	0
quiverResult.consensus.fasta	quiver fasta result	1
quiverResult.consensus.fastq	quiver fastq result	1
quiverResult.gff	quiver gff result	0
quiverResult.vc.err	quiver stderr	1
quiverResult.consensus.fasta.fai	index on quiver fasta result	0
quiver.done	quiver checkpoint file	0
clusalignments.cmd	commandline to generate alignments against quiver for clustering	0
alignments.cmp.h5	alignments for clustering	0
alignments.cmp.h5.error	hit summaries for all clustering alignments	1
alignments.filterFull	clustering alignments that passed filter	1
msajob.sh	command line to dump MSA for clustering alignment	0
msajob.stderr	msa stderr	0
aac.msa.info	msa information rows / cols	1
aac.msa	dump of msa clustering alignments with rows and cols as given in aac.msa.info	1
aac.msa.inserts	for every insert greater than 4 bases wrt the reference (rowNum, baseColumnNum Identity)	1
aac.id	summary stats for every row in the msa	1
aac.col	summary stats for every column in the msa	1
msajob.stout	msa stdout	0
distjob.sh	commandline for computing variant positions and pairwise distances	0
distjob.usecols	columns in the msa to use for clustering	1
distjob.err	distance job stderr	0
aac.msaToDist	pairwise distances (coli, colj, dist, fullDiffs, fullTotal, overlapDiffs, overlapTotal)	0
Rscript.R	R code to compute clustering and cutpoints	0
Rscript.input	R code inputs	0
Rjob.sh	command line to run the R code	0
Rplots.pdf	clustering plot from the R code	1
clusterThreshold.txt	threshold used to cut clustering	0
clusterThreshold.support	data supporting cluster threshold	1
clustergroup_*.ids	list of read ids in each cluster. greater than 1 cluster is evidence of subspecies	1
.RData	data store from the R code	0
Rscript.Rout	output from the R code	0
all.done	checkpoint for program done	0
"""

mapFileToDescription = []
for ll in mapFileToDescriptionDat.split("\n"):
  if len(ll)>1:  
      mapFileToDescription.append(ll.split("\t"))

################################

def short( name ):
    if len(name)>40:
        return(name[:15]+"..."+name[-15:])
    else:
        return(name)

def linksize( ff ):
    return( "<a href=\"%s\">%s</a> %d bytes" % (ff, ff, os.path.getsize(ff)))

print "<html>"
print "#resultsSummary-version=1.0"
print "<h1>Results for CluCon in %s </h1>" % os.getcwd()
print "Gives all results for the CluCon run"
print "<pre>"
print "----"
print "- quiver reference length: %s" % open("quiverResult.consensus.fasta.fai").read().split("\t")[1]
print "- quiver reference: %s" % linksize("quiverResult.consensus.fasta")

print "- number of passing alignments to quiver reference: %d" % len(open("alignments.filterFull").read().splitlines())

dd = open("aac.msa.info").read().strip().split(" ")
print "- size of MSA: rows=%s cols=%s" % (dd[-1],dd[-3])

print "- number of variant columns in MSA: %d" % len(open("distjob.usecols").read().splitlines())

print "- clustering plot: <a href=\"Rplots.pdf\">Rplots.pdf</a>"

print "- number of subclusters: %d" % len(glob.glob("clustergroup_*.ids"))

print "----"
print "</pre>"

################################
print "<h2>options.items</h2>"
print "<table border=\"1px\">"
for ll in open("options.items").read().splitlines():
    ff=ll.split("\t")
    print "<tr> <td>%s</td> <td>%s</td> </tr>" % (ff[0],ff[1])
print "</table>"

################################
print "<h2>All Result Files</h2>"
print "Green is a more final result."

print "<table border=\"1px\">"
print "<tr bgcolor=\"#dd99ee\"> <td>num</td> <td>%s</td> <td>%s</td> <td>%s</td> </tr>" % ("filename","description", "size")
for ii in range(len(mapFileToDescription)):
    ff = mapFileToDescription[ii]
    if ff[2]=="1":
        bgcolor="#00dd00"
    else:
        bgcolor="#ffffff"

    if "*" not in ff[0]:
        print "<tr bgcolor=\"%s\"> <td>%d</td> <td><a href=\"%s\">%s</a></td> <td>%s</td> <td>%d</td> </tr>" % (bgcolor, ii,ff[0],short(ff[0]),ff[1],os.path.getsize(ff[0]))
    else:
        for file in glob.glob(ff[0]):
            print "<tr bgcolor=\"%s\"> <td>%d</td> <td><a href=\"%s\">%s</a></td> <td>%s</td> <td>%d</td> </tr>" % (bgcolor,ii,file,short(file),ff[1],os.path.getsize(file))
print "</table>"
