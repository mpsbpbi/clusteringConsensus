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

def scoreIt(dat):

    ################################
    # get the model data
    mypd = pd.probderive()
    mypd.readallmodels()
    #mypd.models["A1"][1,0]

    ################################
    # read in all observations of singleHP region, "singleHP.msaobs-hp.txt" (not .decomp.)
    # 10	25	C.T.2.G	{ 'C....T....T....G': 605,'                ': 196,'     T....T....G': 83,'          T....G': 70,'               G': 21,'C....T....Tg...G': 6,'C....-....T....-': 1 };
    allobs = {}
    ff = dat.split("\t")
    gg= ff[2].split(".")
    modelname = "%s%s" % (gg[1],gg[2])

    mysum = 0
    for (key,vv) in [x.split(": ") for x in ff[3][2:-3].split(",")]: # ff[3] discard "{ "..." };"
        key = key[2:-2] # discard surrounding quotes and first and last base
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
        #print kk, allobs[kk]
    # C 1267, CT 683,  76, CTT 57, T 40, AC 34, CC 21, CTG 16, CG 15, CTTT 14, TC 12, A 8, GC 6, AAC 5, G 5, CA 5, CCT 4, AAAC 4, CTTTT 4, GCT 3, CCC 2, CGG 2, TGC 2, GGC 2, CTGT 2, CGT 2, AGCAC 1, AGC 1, CAT 1, AAACTTT 1, GCCTT 1, CGGGCAATA 1, CAAAC 1, GAAC 1, AACTT 1, TTGC 1, ACC 1, TA 1, GG 1, GAAGC 1, TTTCG 1, GCC 1, CCCC 1, TGGTCC 1, ACTT 1, TTTT 1, GGGGC 1, CTCC 1, CTGTT 1, TGTGC 1

    ################################
    # score data

    # could use "matix" nx1 as in README_NIAID-closehiv-simplephase.html or chi-squared

    ### chisquared
    
    ## get expected probabilities
    if modelname not in mypd.models:
        sys.stderr.write( "modelname %s not in np.models\n" % modelname)
        return((ff[0],ff[1],ff[2],1.0))

    expectedProb = np.zeros( len(aos) )
    ii=0
    for kk in aos:
        # kk, allobs[kk] = C 1267
        expectedProb[ii] = ( math.pow(2.0, -mypd.pderive(modelname,kk) ) )
        ii += 1

    expected = np.round(expectedProb * sum(observed))

    if False:
        print "id observed expected expectedProb over", sum(observed)
        for ii in range(len(aos)):
            print aos[ii], observed[ii], expected[ii], expectedProb[ii]

    result = scipy.stats.chi2_contingency( np.vstack( (observed, expected) ) )
    return( (ff[0],ff[1],ff[2],result[1]) )

if __name__ == "__main__":

    threshold = float(sys.argv[1])

    input= sys.stdin.read().splitlines()

    # compute and write results
    keep = []
    ofp = open("distjob.ranking","w")
    for ii in input:
        res = scoreIt(ii)
        ofp.write("\t".join([ str(x) for x in res]))
        if res[3] < threshold:
            keep.append(res)
            ofp.write("\tkeep")
        ofp.write("\n")
    ofp.close()

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
                  
