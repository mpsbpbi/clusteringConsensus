#!/usr/bin/env python
__doc__ = """python cluconRecurse.py cluconFrom4k-1-runallgivenref-AC222_D18_assembly.fa-222-18_15Mar~A01_1/ > cluconFrom4k-1-runallgivenref-AC222_D18_assembly.fa-222-18_15Mar~A01_1.cluconRecurse

Convention: foo -> foo-~_D1C1, foo-~_D1C2
"""

import sys
import os
import glob

threshNumber = 16

children = [(sys.argv[1],0)]

tree = { sys.argv[1]: {"name": sys.argv[1], "children": list() } }

template = "export SEYMOUR_HOME=/mnt/secondary/Smrtanalysis/opt/smrtanalysis;  source $SEYMOUR_HOME/etc/setup.sh; export PATH=/home/UNIXHOME/mbrown/mbrown/workspace2013Q1/pacbioCode-viral-clusteringConsensus-v1/code:$PATH; cd %s; time ConsensusClusterSubset.py --nproc=8 --runDir %s --fasta %s --subids %s --ref %s --spanThreshold=%s --entropyThreshold=%s --basfofn %s"
# runDir fasta subids ref spanThreshold entropyThreshold basfofn

while len(children)>0:
    current = children.pop()

    print "================================"
    print "current\t%s" % str(current)

    if not os.path.exists("%s/all.done" % current[0]):
        print "not done yet. wait"
        continue

    clusterfiles = glob.glob("%s/clustergroup_*.ids" % current[0])

    if len(clusterfiles)<2:
        print "DONE. one cluster in %s" % current[0]
        continue

    for cc in clusterfiles:
        print ">>>", cc

        dat = open(cc).read().splitlines()
        if len(dat)<threshNumber:
            print "child %s too small with %d reads" % (cc,len(dat))
            continue

        # look for existing child
        mydepth = current[1]+1
        mynum = os.path.basename(cc)
        mynum = mynum.replace("clustergroup_","")
        mynum = mynum.replace(".ids","")
        mynum = int(mynum)
        mychild = "%s-~_D%dC%d" % (current[0],mydepth,mynum)
        
        if os.path.exists(mychild):
            # child exists, work on it next iteration
            print "push %s" % mychild
            children.append( (mychild, mydepth) )
            newnode = {"name": mychild, "children": list()}
            tree[current[0]]["children"].append(newnode)
            tree[mychild] = newnode

        else:
            # we need to run it
            print "runit %s" % mychild

            # get the info for the run
            options = dict()
            dat = open("%s/options.items" % current[0]).read().splitlines()
            for ll in dat:
                ff=ll.split("\t")
                options[ff[0]]=ff[1]

            # runDir fasta subids ref spanThreshold entropyThreshold basfofn
            cmd = template % (os.getcwd(), mychild, options["fasta"], cc, "%s/quiverResult.consensus.fasta" % current[0], options["spanThreshold"], options["entropyThreshold"], options["basfofn"])
            print cmd
 
print "tree=%s" % tree[ sys.argv[1] ]
