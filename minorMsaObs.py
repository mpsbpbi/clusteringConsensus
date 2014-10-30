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

def scoreIt(dat):

    ################################
    # read in all observations of singleHP region, "singleHP.msaobs-hp.txt" (not .decomp.)
    # old# 10	25	C.T.2.G	{ 'C....T....T....G': 605,'                ': 196,'     T....T....G': 83,'          T....G': 70,'               G': 21,'C....T....Tg...G': 6,'C....-....T....-': 1 };
    # new# 12415-12440-T.G.4.A+	{ 'T....G....G....G....G....A+': 1647,'T....-....G....G....G....A+': 394,...: 1 };
    allobs = {}
    ff = dat.split("\t")
    gg = ff[0][0:-1].split("-")[2].split(".")
    if not doSurround:
        modelname = "%s%s" % (gg[1],gg[2])
    else:
        modelname = "%s1+%s%s+%s1" % (gg[0],gg[1],gg[2],gg[3]) # try surround context

    mysum = 0
    for (key,vv) in [x.split(": ") for x in ff[1][2:-3].split(",")]: # ff[1] discard "{ "..." };"
        if not doSurround:
            key = key[2:-3] # discard surrounding quotes and first and last base and the tailing "+"
        else:
            key = key[1:-2] # discard surrounding quotes and the tailing "+". try surround context
        if " " in key: continue # don't look at end-incomplete data
        key=key.upper()
        key=key.replace("-","")
        key=key.replace(".","")
        key=key.replace(" ","")
        if len(key)>0:
            allobs[ key ] = allobs.get(key,0) + int(vv)
            mysum += int(vv)

    # remove all counts < 0.005 (1/2% ...3) TODO: threshold
    thresh = 0.005*mysum
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

    # could use "matix" nx1 as in README_NIAID-closehiv-simplephase.html or chi-squared

    ### chisquared
    
    ## get expected probabilities
    if not doSurround:
        if modelname not in mypd.models:
            sys.stderr.write( "modelname %s not in np.models\n" % modelname)
            tmp = ff[0].split("-")
            tmp.append(1.0)
            tmp.append(1.0)
            return(tmp)

    expectedProb = np.zeros( len(aos) )
    ii=0
    for kk in aos:
        # kk, allobs[kk] = C 1267
        expectedProb[ii] = ( math.pow(2.0, -mypd.pderive(modelname,kk) ) )
        ii += 1

    #todo what is sumof expectedProb?
    expected = np.round(sum(observed) * expectedProb / sum(expectedProb) )

    bestfisher = 1.0
    debug = True
    if debug: print "#id", ff[0], sum(observed), "sum(expectedProb)", sum(expectedProb)
    for ii in range(len(aos)):
        table = [[int(observed[ii]+0.5), int(sum(observed))],[ int(expected[ii]+0.5), int(sum(observed))]]
        (ftodds, pval) = scipy.stats.fisher_exact(table)
        if debug: print aos[ii], observed[ii], expected[ii], expectedProb[ii], pval
        if observed[ii]>expected[ii]:
            if pval<bestfisher: bestfisher=pval

    result = scipy.stats.chi2_contingency( np.vstack( (observed, expected) ) )
    tmp = ff[0].split("-")
    tmp.append(result[1])
    tmp.append(bestfisher)
    return( tmp )

if __name__ == "__main__":

    chithreshold = float(sys.argv[1])
    fisherthreshold = float(sys.argv[2])

    doSurround=False
    if len(sys.argv)>3:
        if sys.argv[3]=="surround":
            doSurround=True

    input= sys.stdin.read().splitlines()

    # compute and write results

    # sortof messy here, but not read in each time
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
                  
