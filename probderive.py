#!/usr/bin/env python
"""compute probability of Observing bases given Truth given homopolmer
basis sets: pderive( "A2+G1+A2", "AAGAA" ).

From /home/UNIXHOME/mbrown/mbrown/workspace2014Q3/decon-system-matrix/probderive.py
"""

import sys
import math
import numpy as np
# import fisher

################################
def choose(n,k):
    f=math.factorial
    return f(n)/f(k)/f(n-k)

################################
# NOTE: error when trying with class method: TypeError: pderive() takes exactly 3 arguments (2 given)
def memoize(f):
    """ Memoization decorator for functions taking one or more arguments. """
    class memodict(dict):
        def __init__(self, f):
            self.f = f
        def __call__(self, *args):
            return self[args]
        def __missing__(self, key):
            ret = self[key] = self.f(*key)
            return ret
    return memodict(f)

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
class probderive:
#    maxrows=14
#    maxcols=18

    maxrows=26
    maxcols=27

    ################################
    def __init__(self):
        self.models = dict()
        self.memo=dict()

    ################################
    # @memoize
    def pderive( self, HPstr, obs):
        __doc__ = "prob(HP->obs) by trying all prefixes"

        kkey = HPstr+"+"+obs
        if kkey in self.memo:
            return(self.memo[kkey])

        ## print "pderive entry ", HPstr, obs

        HP = HPstr.split("+")
        if len(HP)==1:
            # table lookup
            model = HP[0]
            (mon,moff) = monmoff( model[0], obs)
            if model in self.models:
                if mon>self.maxrows: mon=self.maxrows
                if moff>self.maxcols: moff=self.maxcols
                result=self.models[model][mon][moff]
                result = result + (-math.log(pow(3.0,-moff))/math.log(2.0))
                result = result + (-math.log(pow(choose(mon+moff,moff),-1))/math.log(2.0)) # prob, 3^moff off string, 1/choose to place moffs, TODOaccounts for duplicates GG????
            else:
                result = 30 # no model assume 2^-30 prob
            self.memo[kkey]=result
            return(result)
        else:
            logmysum = 99E+99
            bestlogp = 99E+99
            bestsplit= -1
            for split in range( len(obs)+1 ):
                left  = self.pderive( HP[0], obs[0:split] )
                right = self.pderive(  "+".join(HP[1:]), obs[split:]  )
                ## print "pderive split=%d %f + %f = %f" % (split, left,right, left+right)
                logmysum = -math.log(math.pow(2.0,-logmysum) + math.pow(2.0,-(left+right))) / math.log(2.0)
                # if (left+right)<bestlogp:
                #     bestlogp=left+right
                #     bestsplit=split
            ## print "pderive best ", HPstr,obs, bestlogp, bestsplit
            self.memo[kkey]=logmysum
            return(logmysum)

    ###############################
    def constructLogP( self, ll ):
        ff = ll.split("\t")
        rows=int(ff[4])
        cols=int(ff[5])
        mat = [[0.0 for x in range(self.maxcols+1)] for x in range(self.maxrows+1)]

        dd = ff[6].split(",")
        tt=1
        mysum=0.0
        for rr in range(rows+1):
            for cc in range(cols+1):
                mat[rr][cc]=float(dd[tt])
                mysum+=float(dd[tt])
                tt+=1
        mysum += 0.01*self.maxcols*self.maxrows # pseudo-count

        for rr in range(self.maxrows+1):
            for cc in range(self.maxcols+1):
                mat[rr][cc]= -math.log((mat[rr][cc]+0.01)/mysum)/math.log(2)

        return((ff[1],ff[2],np.array(mat))) # A,1,mat

    def readallmodels(self):
        #for hp in range(1,8):
        #    dat = open("/home/UNIXHOME/mbrown/mbrown/workspace2014Q3/basis-variantid/D100.msaobs-hp.decomp.halforder.%d.txt" % hp).read().splitlines()
            for ll in open("/home/UNIXHOME/mbrown/mbrown/workspace2014Q4/clucon-better-model/basis-full-DNA.mmmatrix").read().splitlines():
                model = self.constructLogP(ll)
                self.models["%s%s" % (model[0],model[1])] = model[2]

    def modelForR(self, name):
        out = "modelForR = matrix( c( "
        numbs = []
        for rr in range(self.maxrows+1):
            for cc in range(self.maxcols+1):
                numbs.append( "%f" % self.models[name][rr][cc] )
        out += ", ".join(numbs)
        out += "), nrow=%d, ncol=%d, byrow=T)" % (self.maxrows+1, self.maxcols+1)
        return(out)


    ################################
    def test1():
        probX= self.pderive( "A4", "AAGAA" )
        probY= self.pderive( "A2+G1+A2", "AAGAA" )
        print "probX, probY", probX, probY
        print "probX, probY 5.23696614586 1.18144788777"
        print "without 3^-k old= probX, probY 3.65200364514 1.15015473386"

        print HPify("AAGAA")
        print "A2+G1+A2"

################################
def singleScore():
    __doc__="Example to score a single HP region"
    mypd = probderive()

    mypd.readallmodels()

#     # read in all observations of singleHP region, "singleHP.msaobs-hp.txt" (not .decomp.)
#     allobs = {}
#     dat = open(sys.argv[1]).read().splitlines()
#     ff = dat[0].split("\t")
#     for (key,vv) in [x.split(": ") for x in ff[3][2:-3].split(",")]: # ff[3] discard "{ "..." };"
#         key = key[2:-2] # discard surrounding quotes and first and last base
#         key=key.upper()
#         key=key.replace("-","")
#         key=key.replace(".","")
#         if len(key)>0:
#             allobs[ key ] = allobs.get(key,0) + int(vv)

#     aos = sorted(allobs, key=allobs.get, reverse=True)
#     observed = np.zeros( len(aos) )
#     ii=0
#     for kk in aos:
#         observed[ii]=allobs[kk]
#         ii+=1
#         print kk, allobs[kk]
#     # C 1267, CT 683,  76, CTT 57, T 40, AC 34, CC 21, CTG 16, CG 15, CTTT 14, TC 12, A 8, GC 6, AAC 5, G 5, CA 5, CCT 4, AAAC 4, CTTTT 4, GCT 3, CCC 2, CGG 2, TGC 2, GGC 2, CTGT 2, CGT 2, AGCAC 1, AGC 1, CAT 1, AAACTTT 1, GCCTT 1, CGGGCAATA 1, CAAAC 1, GAAC 1, AACTT 1, TTGC 1, ACC 1, TA 1, GG 1, GAAGC 1, TTTCG 1, GCC 1, CCCC 1, TGGTCC 1, ACTT 1, TTTT 1, GGGGC 1, CTCC 1, CTGTT 1, TGTGC 1

#     # compute system matrix: system[0][1] = prob(C->CT)
#     system = np.zeros( (len(aos),len(aos)) )
#     for ii in range(len(aos)):
#         for jj in range(len(aos)):
#             if len(aos[ii])==0:
#                 # null cant derive anything
#                 nlp=0.0
#             else:
#                 nlp = math.pow(2.0,-mypd.pderive( HPify(aos[ii]), aos[jj] ))
#             system[ii,jj] = nlp


#     estimate = np.linalg.solve(system, observed)
#     #estimate = np.linalg.lstsq(system, observed)
#     print "observed", observed
#     print "estimate", estimate
#     print "errors", (estimate-observed)


#     print "================================ top 1"
#     hyp = np.zeros(len(aos))
#     hyp[0] = sum(observed)
#     hypObs = np.dot(system,hyp)
#     for ii in range(len(aos)):
#         table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
#         print aos[ii], table, 
#         ft=fisher.FishersExactTest(table)
#         pval = ft.two_tail_p()
#         print pval
#         if pval<1.0E-05 and observed[ii]>hypObs[ii]:
#             print "should add!", aos[ii]
# # CT [[90, 2239], [683, 2239]] 4.53302972746e-99
# # should add! CT
# # CTT [[12, 2239], [57, 2239]] 3.29962440949e-08
# # should add! CTT
# # T [[7, 2239], [40, 2239]] 9.78189685198e-07
# # should add! T

#     print "================================ top 2"
#     hyp = np.zeros(len(aos))
#     hyp[0] = sum(observed)*observed[0]/(observed[0]+observed[1])
#     hyp[1] = sum(observed)*observed[1]/(observed[0]+observed[1])
#     hypObs = np.dot(system,hyp)
#     for ii in range(len(aos)):
#         table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
#         print aos[ii], table, 
#         ft=fisher.FishersExactTest(table)
#         pval = ft.two_tail_p()
#         print pval
#         if pval<1.0E-05 and observed[ii]>hypObs[ii]:
#             print "should add!", aos[ii]


#     print "================================ top 3"
#     hyp = np.zeros(len(aos))
#     hyp[0] = sum(observed)*observed[0]/(observed[0]+observed[1]+observed[2])
#     hyp[1] = sum(observed)*observed[1]/(observed[0]+observed[1]+observed[2])
#     hyp[2] = sum(observed)*observed[2]/(observed[0]+observed[1]+observed[2])
#     hypObs = np.dot(system,hyp)
#     for ii in range(len(aos)):
#         table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
#         print aos[ii], table, 
#         ft=fisher.FishersExactTest(table)
#         pval = ft.two_tail_p()
#         print pval
#         if pval<1.0E-05 and observed[ii]>hypObs[ii]:
#             print "should add!", aos[ii]

################################
def multipleScore():
    __doc__="Example to score a multiple HP regions separated by +'s"

    mypd = probderive()
    mypd.readallmodels()

    # # read in all observations of singleHP region, "msaobs-hp-set.output"
    # # 11950-11960-T.G.1.T+12085-12095-G.C.1.A+	{ 'T....G....T+G....C....A+': 294,'T....A....T+G....G....A+': 271, ... ,'T....G....T+Gt...C....-+': 1 };
    
    # # get number of regions
    # (refhpdat, obs) = open(sys.argv[1]).read().split("\t")
    # refhp = refhpdat.split("+")
    # refhpnum = (len(refhp)-1)
    # allobs = [ list() for x in range(refhpnum) ]

    # # get number of observations
    # numobs = 0 # for all: len(obs[2:-3].split(","))

    # # break apart observations.
    # observed=[]
    # haplotypes=[]
    # for (key,vv) in [x.split(": ") for x in obs[2:-3].split(",")]: # obs discard "{ "..." };"

    #     if int(vv)<5: continue # threshold all counts less than 5
        
    #     numobs+=1
    #     observed.append( int(vv) )
    #     haplotypes.append(key)

    #     print "keyxxx", key
    #     key = key[1:-1] # discard single quotes
    #     keyhp = key.split("+")
    #     for ii in range(refhpnum):
    #         key = keyhp[ii]
    #         # old from full obs including edges
    #         #key = key[2:-2] # discard surrounding quotes and first and last base
    #         #key=key.upper()
    #         #key=key.replace("-","")
    #         #key=key.replace(".","")
    #         allobs[ii].append(key)

    # # NOW: allobs[hpnum][obsnum] = count

    # # for each region, compute system matrix independently. TODO: if
    # # two regions are adjacent then concatenate because not
    # # independent

    # system = np.zeros( (numobs, numobs) )
    # for ii in range(numobs):
    #     for jj in range(numobs):
    #         prodnlp=1.0
    #         for hp in range(refhpnum):
    #             if len(allobs[hp][ii])==0:
    #                 # null cant derive anything
    #                 nlp=0.0
    #             else:
    #                 nlp = math.pow(2.0,-mypd.pderive( HPify(allobs[hp][ii]), allobs[hp][jj]))
    #             prodnlp *= nlp # TODO: work in logs

    #             if (ii<4) and (jj<4):
    #                 print "for", ii, jj, hp, allobs[hp][ii], allobs[hp][jj], nlp

    #         system[jj,ii] = prodnlp
    #         if (ii<4) and (jj<4):
    #             print "system for", ii, jj, haplotypes[ii],haplotypes[jj], prodnlp

    # print "refhpnum", refhpnum, "numobs", numobs, "allobs", allobs, "observed", observed
    # # print "system"
    # # for ii in range(numobs):
    # #     for jj in range(numobs):
    # #         print " ", system[ii,jj], 
    # #     print

    # print "================================ top 1"
    # hyp = np.zeros(numobs)
    # hyp[0] = sum(observed)
    # hypObs = np.dot(system,hyp)
    # for ii in range(10):
    #     table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
    #     print haplotypes[ii], table, 
    #     ft=fisher.FishersExactTest(table)
    #     pval = ft.two_tail_p()
    #     print pval
    #     if pval<1.0E-05 and observed[ii]>hypObs[ii]:
    #         print "should add!", haplotypes[ii]

    # print "================================ top 2"
    # hyp = np.zeros(numobs)
    # hyp[0] = sum(observed)*observed[0]/(observed[0]+observed[1])
    # hyp[1] = sum(observed)*observed[1]/(observed[0]+observed[1])
    # hypObs = np.dot(system,hyp)
    # for ii in range(10):
    #     table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
    #     print haplotypes[ii], table, 
    #     ft=fisher.FishersExactTest(table)
    #     pval = ft.two_tail_p()
    #     print pval
    #     if pval<1.0E-05 and observed[ii]>hypObs[ii]:
    #         print "should add!", haplotypes[ii]

    # print "================================ top 3"
    # hyp = np.zeros(numobs)
    # hyp[0] = sum(observed)*observed[0]/(observed[0]+observed[1]+observed[2])
    # hyp[1] = sum(observed)*observed[1]/(observed[0]+observed[1]+observed[2])
    # hyp[2] = sum(observed)*observed[2]/(observed[0]+observed[1]+observed[2])
    # hypObs = np.dot(system,hyp)
    # for ii in range(10):
    #     table= [[int(hypObs[ii]), int(sum(observed))],[ int(observed[ii]), int(sum(observed))]]
    #     print haplotypes[ii], table, 
    #     ft=fisher.FishersExactTest(table)
    #     pval = ft.two_tail_p()
    #     print pval
    #     if pval<1.0E-05 and observed[ii]>hypObs[ii]:
    #         print "should add!", haplotypes[ii]

################################
if __name__ == "__main__":
    #singleScore()
    multipleScore()
