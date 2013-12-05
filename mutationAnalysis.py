#!/usr/bin/env python
__doc__= """python mutationAnalysis.py aac.msa 591 44795

I want to read in the alginemnts "aac.msa" with 44795 columns and 591
rows assuming it is headed by the reference, step through the columns
counting characters, for each of the bases compute a Fisher's exact
test against some null fraction and output summary stats:

column IDtested IDisRef pval frac numer denom nullFrac nullNumer nullDenom

This essentially replaces codonMutAnalysisALLPOS.py because the bulk
of what it does is step through the alignment in memory which is
already dumped here. codonMutAnalysisALLPOS.py might still be useful
for very large alignments where writing explicitly to disk might be
burdensome.
"""

import sys
from scipy.stats import fisher_exact

print "column\tIDtested\tRefBase\tisgreater\tpval\tfrac\tnumer\tdenom\tnullFrac\tnullNumer\tnullDenom"

rows = int(sys.argv[2])
cols = int(sys.argv[3])+1 # +1 FOR NEWLINE

dat = open(sys.argv[1]).read()
if len(dat)!=rows*cols:
    sys.stderr.write("ERROR. len(dat)!=rows*cols %d %d %d %d\n" % (len(dat),rows,cols, rows*cols))
    sys.exit(1)
else:
    sys.stderr.write("LOG. good len(dat)==rows*cols %d %d %d %d\n" % (len(dat),rows,cols, rows*cols))

# step through columns. TODO: step in blocks rather than single
myc = 0
while myc<(cols-1):
    RefBase = dat[myc] # assumes at 0th row

    obs = dict()
    for myr in range(1, rows):
        myind = myr*cols+myc
        obs[dat[myind]] = obs.get(dat[myind],0)+1

    denom = float(rows - obs.get(" ",0))

    # print "LOG myc=%d obs=%s denom=%f" % (myc,str(obs),denom)

    # not enough data
    if denom<64:
        myc+=5
        continue

    for (IDtested, numer) in obs.items():

        if IDtested==" ": continue

        frac = numer/denom

        nullNumer = 15.0/5.0
        nullDenom = 100.0
        nullFrac = nullNumer/nullDenom

        if frac > nullFrac:
            isgreater = 1
        else:
            isgreater = 0

        pval = fisher_exact( [[numer,denom], [nullNumer,nullDenom]])[1]

        print "%d\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (myc, IDtested, RefBase, isgreater, pval, frac, numer, denom, nullFrac, nullNumer, nullDenom)

    myc+=5
