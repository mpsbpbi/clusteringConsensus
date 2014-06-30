#!/usr/bin/env python
# NOTE: this is depracated in 2.3: TODO: replace with pbcore.io

from pbpy.io.cmph5 import CmpH5 # this was from pbpy.io import CmpH5IO
import sys

def errFromCmph5(cmph5fn):
  print "query_id\terr\tnumerr\trlen\ttarget_id\ttargetStart\ttargetEnd"

  cmph5f = CmpH5.factory.create(cmph5fn,'r')
  for h in cmph5f.alnHitIterator():
    # print dir(h)
    # ['__contains__', '__doc__', '__getitem__', '__init__', '__len__', '__module__', '__str__', '_aln2query', '_aln2target', '_alnDataFromQueryId', '_buildQueryName', '_clearMatchInfo', '_colToValue', '_computeAlnPosTranslations', '_isReadOverhanging', '_parseReadId', '_query2aln', '_target2aln', 'alignedLength', 'alignedQuery', 'alignedTarget', 'alignmentArray', 'alignmentId', 'alignmentIndexData', 'alnToQueryPos', 'alnToTargetPos', 'constructBuffer', 'fullyAligned', 'getNumAllErrors', 'hasAlignment', 'hasAlignments', 'hasFullOverlap', 'hasMatchInfo', 'hasPulseInfo', 'identity', 'indexCols', 'loadFromCmpH5Dataset', 'loadPulseInfo', 'mapQV', 'moleculeId', 'nDel', 'nIns', 'nMatch', 'nMismatch', 'parseAgar', 'pulseInfo', 'queryContained', 'queryExtendsTargetLeft', 'queryExtendsTargetRight', 'queryToAlnPos', 'query_end', 'query_id', 'query_length', 'query_start', 'query_strand', 'readable', 'readableLine', 'score', 'setMatchInfo', 'stripGaps', 'switchTargetQuery', 'targetToAlnPos', 'target_end', 'target_id', 'target_length', 'target_start', 'target_strand', 'validateHit', 'zScore']

    numerr=(h.nMismatch + h.nIns + h.nDel)
    rlen = (h.query_end-h.query_start)
    if numerr>0:
      err = float(numerr)/float(rlen)
    else:
      err = 1.0/float(rlen+1)
    if h.target_strand=="-":
      tid = h.target_id+"rc"
    else:
      tid = h.target_id
    print("%s\t%f\t%d\t%d\t%s\t%d\t%d" % (h.query_id, err,numerr, rlen,tid, h.target_start, h.target_end))

################################
def main():
# argv[1]="CCS.subreads.cmp.h5"
  errFromCmph5(sys.argv[1])

if __name__ == "__main__":
    main()
