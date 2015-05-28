#!/home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/virtscipy/bin/python2.7
__doc__ = """Compute probabilities that positions contain a minor variant.
rewrite minorDecompData.py to minorMsaObs.py to take output of
msaobs-hp-set.py (on all positions) and compute chi-squared rather
than summed log likelihood.
From /home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/minorMsaObs.py
"""

import sys
sys.path.append("/home/UNIXHOME/mbrown/mbrown/workspace2013Q1/pacbioCode-viral-clusteringConsensus-v1/code")

import math
import numpy as np
import probderive as pd
import scipy
import scipy.stats
import os

def monmoff( truth, ss ):
  ss=ss.upper()
  ss=ss.replace("-","")
  ss=ss.replace(".","")
  tc = ss.count(truth)
  nc = len(ss)-tc
  return( tc, nc )

################################
def scoreCounts( modelname, allobs, sumcounts):
    # remove all counts < 0.005 (1/2% ...3) TODO: threshold
    thresh = 0.005*sumcounts
    for (k,v) in allobs.items():
        if v<thresh: del(allobs[k])

    aos = sorted(allobs, key=allobs.get, reverse=True)
    observed = np.zeros( len(aos) )
    ii=0
    for kk in aos:
        observed[ii]=allobs[kk]
        ii+=1
        # print kk, allobs[kk]
    # C 1267, CT 683,  76, CTT 57, T 40, AC 34, CC 21, CTG 16, CG 15, CTTT 14, TC 12, A 8, GC 6, AAC 5, G 5, CA 5, CCT 4, AAAC 4, CTTTT 4, GCT 3, CCC 2, CGG 2, TGC 2, GGC 2, CTGT 2, CGT 2, AGCAC 1, AGC 1, CAT 1, AAACTTT 1, GCCTT 1, CGGGCAATA 1, CAAAC 1, GAAC 1, AACTT 1, TTGC 1, ACC 1, TA 1, GG 1, GAAGC 1, TTTCG 1, GCC 1, CCCC 1, TGGTCC 1, ACTT 1, TTTT 1, GGGGC 1, CTCC 1, CTGTT 1, TGTGC 1

    ################################
    # score data
    
    ## get expected probabilities
    if modelname not in mypd.models:
        return([1.0,1.0])

    expectedProb = np.zeros( len(aos) )
    ii=0
    for kk in aos:
        # kk, allobs[kk] = C 1267
        expectedProb[ii] = ( math.pow(2.0, -mypd.pderive(modelname,kk) ) )
        # print "# expectedProb", modelname, kk, expectedProb[ii]
        ii += 1
    expected = np.round(sum(observed) * expectedProb / sum(expectedProb) )

    bestfisher = 1.0
    debug = True
    if debug: print "#### #id", modelname, sum(observed)
    for ii in range(len(aos)):
        table = [[int(observed[ii]+0.5), int(sum(observed))],[ int(expected[ii]+0.5), int(sum(observed))]]
        (ftodds, pval) = scipy.stats.fisher_exact(table)
        if debug: print "# fisher:", ii, aos[ii], table, pval
        if observed[ii]>expected[ii]:
            if pval<bestfisher: bestfisher=pval

    if sum(observed)>1:
      result = scipy.stats.chi2_contingency( np.vstack( (observed, expected) ) )
      return( [result[1], bestfisher] )
    else:
      return([1.0,1.0])

################################

def scoreIt(dat):

    rMap = dict(zip("ACGTacgt-","TGCAtgca-"))

    ################################
    # read in all observations of singleHP region, "singleHP.msaobs-hp.txt" with fw and rc separated:
    # 5-15-G.C.1.T+	{ 'C>': 1144,'C<': 1078,'CTTG>': 1 };
    (truth, obs) = dat.split("\t")
    (begin,end,context) = truth.split("-")
    (tleft,tmid,tmlen,tright) = context.split(".")

    allobs = {}
    allobsrc = {}
    modelname = "%s%s" % (tmid,tmlen) # C1
    modelnamerc = "%s%s" % (rMap[tmid],tmlen) # G1
    mysum = 0
    mysumrc = 0

    for (key,counts) in [x.split(": ") for x in obs[2:-3].split(",")]: # obs discard "{ "..." };"

        key = key[1:-1] # discard single quotes
        counts = int(counts)

        if key[-1] == ">":
            # forward
            isFW = True
            key = key[:-1]
        else:
            # RC
            isFW = False
            key = key[:-1]
            key = "".join([rMap[x] for x in key[::-1]]) # reverse-complement

        if isFW:
          allobs[ key ] = allobs.get(key,0) + counts
          mysum += counts
        else:
          allobsrc[ key ] = allobsrc.get(key,0) + counts
          mysumrc += counts

    (chifw, fisherfw) = scoreCounts(modelname, allobs, mysum)
    (chirc, fisherrc) = scoreCounts(modelnamerc, allobsrc, mysumrc)

    return( ( begin, end, context, min(chifw,chirc), min(fisherfw,fisherrc), chifw,fisherfw, chirc,fisherrc) )

if __name__ == "__main__":

    chithreshold = float(sys.argv[1])
    fisherthreshold = float(sys.argv[2])

    input= sys.stdin.read().splitlines()

    # compute and write results

    ################################
    # get the model data
    mypd = pd.probderive()
    mypd.readallmodels()
    #mypd.models["A1"][1,0]

    keep = []
    if not os.path.exists("distjob.ranking"):
        ofp = open("distjob.ranking","w")
        for ii in input:
            res = scoreIt(ii)

            ofp.write("\t".join([ str(x) for x in res]))
            # res[3] for chi, res[4] for max fisher
            if (res[3]< chithreshold) and (res[4] < fisherthreshold):
                keep.append(res)
                ofp.write("\tkeep")
            else:
                ofp.write("\t")
            ofp.write("\n")
        ofp.close()
    else:
        for ll in open("distjob.ranking").read().splitlines():
            # res[3] for chi, res[4] for max fisher
            if (float(ll.split("\t")[3]) < chithreshold) & (float(ll.split("\t")[4]) < fisherthreshold):
                keep.append(ll.split("\t"))

    # write out kept columns
    ofp = open("distjob.usecols","w")
    for ii in range(len(keep)):

        ## output entire range
        # ff = int(keep[ii][0])
        # tt = int(keep[ii][1])
        # rr = "\t".join([str(x) for x in range((ff+1),tt)])
        # ofp.write("%s\n" % rr)

        # the single left-most interior match:
        ofp.write("%d\n" % ( int(keep[ii][0])+5 ))

    ofp.close()
                  
