#!/usr/bin/env python
from pbpy.io.cmph5 import  CmpH5AlnHit
from pbpy.io import cmph5
from pbpy.io.cmph5 import CmpH5

import numpy as np
import sys
import math

from optparse import OptionParser

nmap = dict( zip( (0,1,2,3), ('a','c','g','t') ) ) 
rMap = dict(zip("ACGTacgt-","TGCAtgca-"))
UC = dict(zip("ACGTacgt-","ACGTACGT-"))
alphabet = "ACGTacgt-."

################################

def cmph5ToMSA(*argv, **options):

    # read in reference template
    f = open(options["reffile"])
    refSeq = ""
    for l in f:
        l = l.strip()
        if l[0] == ">": 
            continue
        else:
            refSeq += l
    f.close()

    # read in the query_ids to keep
    queryidsToKeep = {}
    f = open(options["filter"])
    dat = f.read().splitlines()
    f.close()
    for ll in range(1,len(dat)):
        ff = dat[ll].split("\t")
        queryidsToKeep[ff[0]] = 1
    sys.stderr.write("len(queryidsToKeep) %d\n" % len(queryidsToKeep))

    # what part of the template do we want the alignment for [low_0,high) ?
    rangeLow = int(options["rangeLow"])
    rangeHigh = int(options["rangeHigh"])
    if rangeHigh>(len(refSeq)+1):
        rangeHigh=(len(refSeq)+1)
    rangeLength = rangeHigh-rangeLow

    # what part of the template do we want the alignment for [low_1,high) ?
    seqLow=int(options["seqLow"])
    seqHigh=int(options["seqHigh"])

    alignment = []
    for ii in range(rangeLength):
        alignment.append({})

    allqueryids = []

    # get alignment data
    cmph5f = cmph5.factory.create(options["cmph5"], 'r')

    # put in reference
    currentHit = 0
    query_id = "reference"
    allqueryids.append(query_id)
    # I guess cmp.h5 is 1-based forcing to start at 1??
    for qp, tp, qb, tb  in zip(range(1,len(refSeq)+1), range(1,len(refSeq)+1), refSeq, refSeq):
        if tp<rangeLow or tp>(rangeHigh-1): continue
        tIdx = tp-rangeLow
        alignment[tIdx][query_id] = UC[qb]

    # step through alignments correctly handling fwd/rc sense
    for h in cmph5f.alnHitIterator():
        # continue if not in queryidsToKeep from filtering
        if not(h.query_id in queryidsToKeep): continue

        currentHit += 1
        # print "working on currentHit %d" % currentHit

        # read in all sequences in region to get correct spacing and
        # then only output selected range below...

        # shortcut on sequence if spacing between blocks not required to be consistent
        # if currentHit>(seqHigh-1) or currentHit<seqLow: continue

        allqueryids.append(h.query_id)

        ts = h.target_start
        te = h.target_end
        qs = h.query_start
        qe = h.query_end
        strand = h.target_strand

        # is the hit overlapping our range? b2<e1 and b1<e2
        if not (ts<rangeHigh and te>rangeLow): continue

        if strand == "+":
            tM =  [ 1 if x != '-' else 0 for x in h.alignedTarget] 
            qM =  [ 1 if x != '-' else 0 for x in h.alignedQuery] 
            qSeq = h.alignedQuery
            tSeq = h.alignedTarget
            tOffset = ts + np.cumsum(tM)
            qOffset = qs + np.cumsum(qM)
        else:
            tM =  [ 1 if x != '-' else 0 for x in h.alignedTarget[::-1]] 
            qM =  [ 1 if x != '-' else 0 for x in h.alignedQuery[::-1]] 
            qSeq = [rMap[c] for c in h.alignedQuery[::-1]]
            tSeq = [rMap[c] for c in h.alignedTarget[::-1]]
            tOffset = ts + np.cumsum(tM)
            qOffset = qs + np.cumsum(qM)

        for qp, tp, qb, tb  in zip(qOffset, tOffset, qSeq, tSeq):
            if tp<rangeLow or tp>(rangeHigh-1): continue
            tIdx = tp-rangeLow
            tmp = alignment[tIdx].get(h.query_id, "")
            # print tp, tIdx, tmp, qb, "%s%s" % (tmp,qb)
            if tmp=="":
                out = UC[qb]
            else:
                out = "%s%s" % (tmp,qb)
            alignment[tIdx][h.query_id] = out

    print "cmp.h5 has %d reads. genome has %d bases." % ( currentHit, len(refSeq))

    # take care of seqHigh more than number in alignment
    if (currentHit+1)<seqHigh:
        seqHigh = currentHit+1
    seqLength = seqHigh-seqLow

    # RESULT: alignment[templatePosition][query_id] = s. len(s)>=1
    # s[0] = "-" (delete) or base (match). s[>0] = base (insert)

    ## I don't have to do this if I'm eliminating inserts.
    # step through and find the maximum length of each aligned column,
    # so I can space out the multiple alignment. This is based on the
    # entire multiple alignment including the reference at [0]

    columnSpace = [0]*rangeLength
    for cc in range(rangeLength):
        columnSpace[cc] = 1+options["maxInsert"]
        ## for inserts it takes a bit more work without fixed size
        # vv = alignment[cc].values()
        # if len(vv)==0:
        #     columnSpace[cc] = 0
        # else:
        #     columnSpace[cc] = max([len(ss) for ss in vv])

    cscs=np.array([0])
    cscs = np.append(cscs,np.cumsum(columnSpace))

    ### construct 2d array of filtered amino acid alignment with
    ### sum(columnSpace) columns and seqLength rows. ("acgt" 4 bases,
    ### "-" Deleted wrt reference " "=not present "."=spaced insert)
    print "MSA %d sequences %d columns for range %d - %d" % (seqLength, sum(columnSpace), rangeLow, rangeHigh)
    
    # check to see if empty in the sequence range
    totallength = 0
    for rr in range(seqLow, seqHigh):
        for cc in range(rangeLength):
            sa = alignment[cc].get(allqueryids[rr],"")
            totallength += len(sa)
    if(totallength==0):
        msaf = open(options["msafile"],"w")
        msaf.write("# cmph5= %s reffile= %s rangeLow= %d rangeHigh= %d seqLow= %d seqHigh= %d columns= %d rows= %d\n" % ( options["cmph5"], options["reffile"], rangeLow, rangeHigh, seqLow, seqHigh, sum(columnSpace), seqLength))
        msaf.write("EMPTY")
        msaf.close()
        sys.exit(0)

    msa = np.empty([seqLength, sum(columnSpace)], np.character)
    rowDist = []
    for ii in range(seqLength):
        rowDist.append({})
    colDist = []
    for ii in range(sum(columnSpace)):
        colDist.append({})

    msafs = open("%s.info" % options["msafile"],"w")
    msafs.write("# cmph5= %s reffile= %s rangeLow= %d rangeHigh= %d seqLow= %d seqHigh= %d columns= %d rows= %d\n" % ( options["cmph5"], options["reffile"], rangeLow, rangeHigh, seqLow, seqHigh, sum(columnSpace), seqLength))
    msafs.close()

    msaf = open(options["msafile"],"w")
    msainsf = open("%s.inserts" % options["msafile"],"w")
    ##msafullseq = open("%s.fullseq" % options["msafile"],"w")

    for rr in range(seqLow, seqHigh):
        ##thisfullseq = []
        for cc in range(rangeLength):
            tofill = columnSpace[cc] # should be constant = options["maxInsert"]
            sa = alignment[cc].get(allqueryids[rr],"")
            ##thisfullseq.append(sa)

            # get rid of inserts that are too long
            mymi = options["maxInsert"]
            if len(sa)>(mymi+1):
                # record the full insert for later
                msainsf.write("%d %d %s\n" % (rr, cc, sa))

                sa = sa[0:(mymi+1)]
                
            # for inserts    
            if (len(sa)==0):
                sa = " "*tofill
            else:
                sa = sa + "."*(tofill-len(sa))

            for tt in range(cscs[cc], cscs[cc+1]):
                myc = sa[tt-cscs[cc]]
                msa[(rr-seqLow), tt] = myc
                colDist[tt][myc] = (colDist[tt].get(myc,0))+1
                rowDist[(rr-seqLow)][myc] = (rowDist[(rr-seqLow)].get(myc,0)) +1

        msaf.write( "".join(msa[(rr-seqLow),]))
        msaf.write( "\n")
        ##fullseq = "".join(thisfullseq)
        ##fullseq = fullseq.replace("-","").upper()
        ##msafullseq.write(">%s\n%s\n" % (allqueryids[rr],fullseq))
    msaf.close()
    msainsf.close()
    ##msafullseq.close()

    # write out the sequence query ids used in the block. append the frequency of observed characters (ACGTacgt-.)
    idf = open(options["idfile"],"w")
    for rr in range(seqLow, seqHigh):
        idf.write("%s" % allqueryids[rr])
        myv = []
        for myc in alphabet:
            idf.write("\t%s\t%d" % (myc, rowDist[rr-seqLow].get(myc,0)))
            myv.append(rowDist[rr-seqLow].get(myc,0)+0.001) # get in correct order
        idf.write("\tmax\t%s" % sorted( zip(tuple(alphabet),myv), key= lambda x: -x[1])[0][0])
        z = sum(myv)
        idf.write("\tz\t%f" % (z))
        e = sum([ -(x/z)*math.log(x/z)/math.log(2.0) for x in myv])
        idf.write("\tentropy\t%f\n" % (e))
    idf.close()

    # write out the column stats in the block. append the frequency of observed characters (ACGTacgt-.)
    idf = open(options["colfile"],"w")
    for cc in range(sum(columnSpace)):
        idf.write("%d" % cc)
        myv = []
        for myc in alphabet:
            idf.write("\t%s\t%d" % (myc, colDist[cc].get(myc,0)))
            myv.append(colDist[cc].get(myc,0)+0.001) # get in correct order
        ss = sorted( zip(tuple(alphabet),myv), key= lambda x: -x[1])
        idf.write("\tmax\t%s" % ss[0][0])
        z = sum(myv)
        idf.write("\tz\t%f" % (z))
        e = sum([ -(x/z)*math.log(x/z)/math.log(2.0) for x in myv])
        idf.write("\tentropy\t%f" % (e))
        idf.write("\tsecond\t%s" % (ss[1][0]))
        idf.write("\tsecondP\t%f\n" % (ss[1][1]/z))
        
        
    idf.close()

    cmph5f.close()

if __name__ == "__main__":
    parser = OptionParser("usage")
    parser.add_option("--cmph5", type="string", dest="cmph5", help="cmp.h5 alignment of ccs reads against reference")
    parser.add_option("--reffile", type="string", dest="reffile", help="file containing reference")
    parser.add_option("--msafile", type="string", dest="msafile", help="the file to store multiple alignments")
    parser.add_option("--idfile", type="string", dest="idfile", help="the file of the ids of the reads")
    parser.add_option("--colfile", type="string", dest="colfile", help="the file of the column statistics")
    parser.add_option("--rangeLow", type="string", dest="rangeLow", help="low end of range of reference to generate MSA for")
    parser.add_option("--rangeHigh", type="string", dest="rangeHigh", help="high end of range of reference to generate MSA for")
    parser.add_option("--seqLow", type="string", dest="seqLow", help="low end of range of seqs")
    parser.add_option("--seqHigh", type="string", dest="seqHigh", help="high end of range of seqs")
    parser.add_option("--filter", type="string", dest="filter", help="file giving query_ids to keep in 1st column (row 0 = header)")
    parser.add_option("--maxInsert", type="int", dest="maxInsert", help="how many columns to space for inserts, 4")

    (options, args) = parser.parse_args()

    if not options.cmph5:
        print "error cmph5 must be specified"
        parser.print_help()
        sys.exit(1)

    if not options.maxInsert:
        print "maxInsert set to 4"
        options.maxInsert = 4

    cmph5ToMSA(**options.__dict__) # object to dict for kwargs
