#!/usr/bin/env python
_doc__ = """Filter reads that span the boundries (with slop) and have less than 50% error

filterFull.py hiv.2450117-0178.cmp.h5.error > hiv.2450117-0178.cmp.h5.error.fullfiltered

The boundries were found by aligning the genes against the reference:

Here are the elements where the primers are placed in our reference:
gag 209:1688
pol 1492:4492
vif 4436:5015
nef 8164:8788

first = dat$targetStart<(209+100) & dat$targetEnd>(4492-100)
second = dat$targetStart<(4436+100) & dat$targetEnd>(8788-100)

Input:
query_id	err	numerr	rlen	target_id	targetStart	targetEnd
m120224_171206_42142_c100283442390000001523008607041295_s1_p0/1425/1812_1881	0.014493	1	69	C.ZM.2003.ZM246F_flD5.ZM246Frc	2179	2247
m120224_171206_42142_c100283442390000001523008607041295_s1_p0/53469/2119_2186	0.014925	1	67	C.ZM.2003.ZM246F_flD5.ZM246F	4135	4201
m120224_171206_42142_c100283442390000001523008607041295_s1_p0/25078/4825_4882	0.017241	0	57	C.ZM.2003.ZM246F_flD5.ZM246Frc	181	238
"""

import sys

infile = sys.argv[1]
spanT = int(sys.argv[2])

slop = 100

mybegin = 209
myend=8788

dat = open(infile).read().splitlines()

ll = 0
numGood = 0
print dat[ll]
ll += 1
while ll<len(dat):
    ff = dat[ll].split("\t")
    err = float(ff[1])
    ts = int(ff[5])
    te = int(ff[6])
    rlen = int(ff[3])

    #if ts<(mybegin+slop) and te>(myend-slop):
    span = te-ts
    if span>spanT:
        print dat[ll]
        numGood += 1

    ll += 1

sys.stderr.write("numGood %d\n" % numGood)
