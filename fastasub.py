#!/usr/bin/env python
__doc__ = """fastasub.py foo.fasta listofids.txt

note: PBI annoying change!
listofids.txt: m131017_075620_sherri_c100581822550000001823093904021432_s1_p0/21274
foo.fasta:    >m131017_075620_sherri_c100581822550000001823093904021432_s1_p0/21274/3592_23677
"""

import sys
import re

# get ids into hash
toget = dict()
dat = open(sys.argv[2]).read().splitlines()
for ll in dat:
    # get rid of subread info if there
    if ll.count("/")==2:
        ll = re.sub("/[0-9_]*$","",ll)
    toget[ll+"/"] = 1

# now go through fasta
name = ''
realname = ''
seqBuffer = []
with open(sys.argv[1]) as fp:
    for line in fp:
        if len(line)<1: continue
        line = line.rstrip()
        if len(line)==0: continue
        if line[0]==">":
            if len(seqBuffer)>0:
                if name in toget:
                    print (">%s\n%s" % (realname, "".join(seqBuffer)))
            seqBuffer = []
            if len(line)==1:
                name = ''
                realname=''
            else:
                realname = line[1:]
                # handle sub
                if "ccs" in name:
                    name = re.sub("ccs$","",realname)
                else:
                    name = re.sub("[0-9_]+$","",realname)
        else:
            seqBuffer.append( line )
    if len(seqBuffer)>0:
        if name in toget:
            print (">%s\n%s" % (realname, "".join(seqBuffer)))

            

