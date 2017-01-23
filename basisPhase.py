#!/home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/virtscipy/bin/python2.7
"""use basis system matrix to phase multiple SNPs
"""

import sys
sys.path.append("/home/UNIXHOME/mbrown/mbrown/workspace2013Q1/pacbioCode-viral-clusteringConsensus-v1/code")

import probderive as pd
import math
import numpy as np
import scipy
import scipy.stats

debug = False

################################
def monmoff( truth, ss ):
  ss=ss.upper()
  ss=ss.replace("-","")
  ss=ss.replace(".","")
  tc = ss.count(truth)
  nc = len(ss)-tc
  return( tc, nc )

def HPify( ss ):
    out = []
    base = ss[0]
    basec = 1
    ii=1
    while ii<len(ss):
        if ss[ii]!=base:
            out.append("%s%d" % (base,basec))
            base = ss[ii]
            basec = 1
        else:
            basec+=1
        ii+=1
    # last
    out.append("%s%d" % (base,basec))
    return("+".join(out))

################################
def multipleScore( indatfile, threshold ):
    __doc__="Example to score a multiple HP regions separated by +'s"

    mypd = pd.probderive()
    mypd.readallmodels()

    # read in all observations of singleHP region, "msaobs-hp-set.output"
    # 11950-11960-T.G.1.T+12085-12095-G.C.1.A+	{ 'T....G....T+G....C....A+': 294,'T....A....T+G....G....A+': 271, ... ,'T....G....T+Gt...C....-+': 1 };
    
    # get number of regions
    (refhpdat, obs) = open(indatfile).read().split("\t")
    refhp = refhpdat.split("+")
    refhpnum = (len(refhp)-1)
    allobs = [ list() for x in range(refhpnum) ]

    # get number of observations
    numobs = 0 # for all: len(obs[2:-3].split(","))

    # break apart observations.
    observed=[]
    haplotypes=[]
    for (key,vv) in [x.split(": ") for x in obs[2:-3].split(",")]: # obs discard "{ "..." };"

        if int(vv)<5: continue # threshold all counts less than 5
        if " " in key: continue # " " is missing from ends so discard

        numobs+=1
        observed.append( int(vv) )
        haplotypes.append(key)

        if debug: print "keyxxx", key
        key = key[1:-1] # discard single quotes
        keyhp = key.split("+")
        for ii in range(refhpnum):
            key = keyhp[ii]
            # old from full obs including edges
            #key = key[2:-2] # discard surrounding quotes and first and last base
            #key=key.upper()
            #key=key.replace("-","")
            #key=key.replace(".","")
            allobs[ii].append(key)

    # NOW: allobs[hpnum][obsnum] = count

    # for each region, compute system matrix independently. TODO: if
    # two regions are adjacent then concatenate because not
    # independent

    system = np.zeros( (numobs, numobs) )
    for ii in range(numobs):
        for jj in range(numobs):
            prodnlp=1.0
            for hp in range(refhpnum):
                if len(allobs[hp][ii])==0:
                    # null cant derive anything
                    nlp=0.0
                else:
                    nlp = math.pow(2.0,-mypd.pderive( HPify(allobs[hp][ii]), allobs[hp][jj]))
                prodnlp *= nlp # TODO: work in logs

                if (ii<4) and (jj<4):
                    if debug: print "for", ii, jj, hp, allobs[hp][ii], allobs[hp][jj], nlp

            system[jj,ii] = prodnlp
            if (ii<4) and (jj<4):
                if debug: print "system for", ii, jj, haplotypes[ii],haplotypes[jj], prodnlp

    if debug: print "refhpnum", refhpnum, "numobs", numobs, "allobs", allobs, "observed", observed

    if debug:
      print "system= cbind( "
      for ii in range(numobs):
              print "c(", ",".join([ str(x) for x in system[ii,] ]), "),"
      print ")"

    ################################
    # cycle through the top starting from 1 until nothing is added

    converged = False
    toadd = [0]

    print "================================ top 0"
    print haplotypes[0], [observed[0], sum(observed)]

    maxiter = min(30,numobs)

    for topnum in range(1,maxiter):
        print "================================ top %d" % topnum
        hyp = np.zeros(numobs)
        mysum = 0.0
        for ii in toadd:
            mysum = mysum + observed[ii]
        for ii in toadd:
            hyp[ii] = observed[ii]/(mysum) # fractions sum to 1.0 and not counts: sum(observed)*

        added=False
        hypObs = np.dot(system,hyp)
        hypObs = sum(observed)*hypObs/sum(hypObs) # fractions to total

        print "hyp= c(%s)" % ",".join([str(x) for x in hyp])
        print "hypObs= c(%s)" % ",".join([str(x) for x in hypObs])
        for ii in range(max(toadd)+1, maxiter):
            table= [[int(hypObs[ii]+0.5), int(sum(observed))],[ int(observed[ii]+0.5), int(sum(observed))]]
            print haplotypes[ii], table, 
            (ftodds, pval) = scipy.stats.fisher_exact(table)
            print pval
            justadd=False
            if pval<0:
              # TODO: hack to get around numeric crap library!
              print >>sys.stderr, "horrible negative pval! taking abs", pval
              pval = abs(pval)
              if table[0][0]<100 and table[1][0]>100: justadd=True

            if justadd or (pval<threshold and observed[ii]>hypObs[ii]):
                print "should add!", ii, haplotypes[ii]
                toadd.append(ii)
                added=True
                break
        else:
            # didn't break so didn't add anything in this topnum iteration
            # toadd explains all the data
            print "================================"
            print "*** converged! the following haplotypes and fractions explain all the data"
            print "*** referencePositions: %s" % refhpdat
            converged = True
            for ii in toadd:
                print "***", ii, haplotypes[ii], hyp[ii]
            break

    if not converged:
        print "*** NOT converged! Still adding after %d!" % maxiter

################################
if __name__ == "__main__":
    multipleScore( sys.argv[1], float(sys.argv[2]))
