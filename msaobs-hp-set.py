#!/usr/bin/env python
"""break multiple sequence alignment data into homopolymer region counts
From /home/UNIXHOME/mbrown/mbrown/workspace2014Q3/NIAID-closehiv/msaobs-hp-set.py
"""

import sys

# the input MSA
dat = open(sys.argv[1]).read().splitlines()

# the list of msa positions that need to be accounted for. if "all"
# then all hp regions are tallied independently
msapos = dict()
if sys.argv[2] != "all":
    for ll in open(sys.argv[2]).read().splitlines():
        # print "## looking for msapos", int(ll)
        msapos[int(ll)] = 1

################################
# find x+y^n+z ranges
ref = dat[0]
ranges=[]

# make sure begin and end aren't HPs
begin=5
while ref[begin]==ref[begin+5]:
    begin+=5
end=len(ref)-5
while ref[end]==ref[end-5]:
    end-=5

rangebegin = begin
this=rangebegin+5
next=this+5
while (this<end):
    while (ref[next]==ref[this]) and (this<end):
        this+=5
        next+=5
    rangeend=next
    if sys.argv[2] == "all":
        ranges.append( (rangebegin,rangeend) )
    else:
        # cycle through inefficiently to see if range contains desired position
        for pp in range(rangebegin+1,rangeend):
            if msapos.get(pp,0)==1:
                ranges.append( (rangebegin,rangeend) )
                #print "## hp range desired", rangebegin, rangeend

    rangebegin=this
    this=rangebegin+5
    next=this+5

sys.stderr.write("#ranges=%s\n" %  repr(ranges))
# collapse hp regions if overlapping
newranges= [ranges[0]]
for ii in range(1,len(ranges)):
    if newranges[-1][1]>=ranges[ii][0]:
        newranges[-1]= (newranges[-1][0],ranges[ii][1])
        sys.stderr.write("#collapse %d = %d\n" % ( (ii),len(newranges)-1 ))
    else:
        newranges.append(ranges[ii])
sys.stderr.write("#newranges=%s\n" %  repr(newranges))
ranges=newranges

################################
# collect the data

hpdat=[]

if sys.argv[2]=="all":
    # look at all HP ranges
    for rr in ranges:
        rangelen = rr[1]/5-rr[0]/5+1
        truth = "%s.%s.%d.%s" % (ref[rr[0]], ref[rr[0]+5], rangelen-2, ref[rr[1]])

        store = {}
        res = "%d\t%d\t%s\t" % (rr[0],rr[1],truth)

        for ii in range(1,len(dat)):
            ll = dat[ii] 
            region = ll[rr[0]:(rr[1]+1)]
            store[region] = store.get(region,0)+1

        ss = sorted(store.items(), key=lambda x: x[1], reverse=True)

        out = []
        for sss in ss:
            out.append("'%s': %d" % (sss))

        res = res+"{ %s };" % ",".join(out)
        hpdat.append(res)

else:
    # look only at the positions that are specified
    hphaplo = [""]*len(dat)

    for ii in range(len(dat)):
        ll = dat[ii] 
        for rr in ranges:
            rangelen = rr[1]/5-rr[0]/5+1
            if ii==0:
                # reference
                hphaplo[ii] += "%d-%d-%s.%s.%d.%s+" % (rr[0],rr[1],ref[rr[0]], ref[rr[0]+5], rangelen-2, ref[rr[1]])
            else:
                # old: everything including ends
                #hphaplo[ii]+= "%s+" % ll[rr[0]:(rr[1]+1)]

                # only bases without ends
                key = "%s+" % ll[ (rr[0]+1) : (rr[1]+1-1) ]
                key=key.upper()
                key=key.replace("-","")
                key=key.replace(".","")
                hphaplo[ii]+= key

    # count up unique occur
    store = {}
    for ii in range(1,len(dat)):
        store[hphaplo[ii]] = store.get(hphaplo[ii],0)+1

    ss = sorted(store.items(), key=lambda x: x[1], reverse=True)
    out = []
    for sss in ss:
        out.append("'%s': %d" % (sss))

    hpdat.append("%s\t{ %s };" % (hphaplo[0], ",".join(out)))

# hpdat now contains all the counts for all the ranges.
#5	15	G.C.1.T	{ 'G....C....T': 572,'           ': 287,'          T': 83,'     C....T': 39,'G....C....-': 1 };
#10	25	C.T.2.G	{ 'C....T....T....G': 605,'                ': 196,'     T....T....G': 83,'          T....G': 70,'               G': 21,'C....T....Tg...G': 6,'C....-....T....-': 1 };
# or all on a single line for linked
    # "11950-11960-T.G.1.T+12085-12095-G.C.1.A+\t{ 'T....G....T+G....C....A+': 294,'T....A....T+G....G....A+': 271,'T....A....T+G....-....A+': 37,'T....G....T+-....C....A+': 18,'T....-....T+G....G....A+': 14,'T....G....-+G....C....A+': 13,'T....G....T+Gc...C....A+': 11,'T....Ta...T+G....G....A+': 10,'-....G....T+G....C....A+': 10,'T....G....T+Gg...C....A+': 10,'T....-....T+G....C....A+': 9,'T....G....T+G....Cc...A+': 9,'Ta...A....T+G....G....A+': 9,'T....A....T+G....Gg...A+': 9,'T....A....T+G....-....G+': 9,'T....A....T+G....A....A+': 9,'T....G....T+G....C....-+': 9,'T....G....T+G....-....A+': 8,'T....A....-+G....G....A+': 7,'T....Gg...T+G....C....A+': 7,'T....G....T+G....Ca...A+': 7,'Ta...T....T+G....G....A+': 7,'T....A....T+G....C....A+': 6,'T....Gt...T+G....C....A+': 5,'T....G....T+Gt...C....A+': 5,'T....A....T+G....Cg...A+': 4,'Tg...G....T+G....C....A+': 4,'T....G....T+G....Caa..A+': 3,'T....Aa...T+G....G....A+': 3,'T....A....T+G....-....-+': 3,'A....-....T+G....-....A+': 3,'T....T....T+G....C....A+': 3,'T....A....T+G....T....A+': 3,'A....-....T+G....G....A+': 3,'T....Gc...T+G....C....A+': 3,'-....A....T+G....G....A+': 3,'-....G....T+-....C....A+': 2,'T....G....T+C....C....A+': 2,'T....-....T+G....G....-+': 2,'T....C....T+G....G....A+': 2,'T....-....T+-....C....A+': 2,'T....G....C+G....C....A+': 2,'Ta...A....T+G....-....G+': 2,'Ta...G....T+G....-....A+': 2,'T....Ga...T+G....C....A+': 2,'T....G....T+G....G....A+': 2,'Tc...G....T+G....C....A+': 2,'T....A....T+G....G....G+': 2,'Tt...G....T+G....C....A+': 2,'T....G....T+Ga...A....A+': 2,'T....G....T+G....Ct...A+': 2,'T....A....T+-....-....A+': 2,'T....A....T+G....G....-+': 2,'Ta...A....T+G....A....A+': 2,'T....-....A+G....G....A+': 2,'T....A....T+Ga...G....A+': 2,'Ta...Gc...T+G....G....A+': 1,'T....A....T+Gt...C....-+': 1,'T....Gcg..T+G....Cg...A+': 1,'-....G....-+G....C....A+': 1,'C....G....T+Gg...Cc...A+': 1,'Tg...G....T+G....Ca...A+': 1,'T....-....T+G....T....A+': 1,'Tg...G....T+Gc...C....A+': 1,'Tt...Gg...T+G....C....A+': 1,'Ta...G....T+G....G....A+': 1,'T....A....T+G....Ct...A+': 1,'T....Gg...T+G....-....A+': 1,'T....G....T+Gcc..C....A+': 1,'T....Aa...T+G....Ag...A+': 1,'G....G....T+G....Ca...A+': 1,'Tca..T....T+G....G....A+': 1,'Ta...A....T+G....-....-+': 1,'T....G....T+G....A....A+': 1,'T....Ga...T+G....-....A+': 1,'-....Ga...T+G....G....A+': 1,'T....Gt...T+G....Ca...A+': 1,'Ta...A....T+G....-....A+': 1,'T....A....-+G....Gg...A+': 1,'T....-....A+G....A....A+': 1,'T....G....-+G....-....A+': 1,'-....-....A+G....-....G+': 1,'T....-....T+G....C....C+': 1,'T....Ga...T+G....G....A+': 1,'T....Ta...T+G....Cg...A+': 1,'T....Gt...T+A....C....A+': 1,'T....G....T+Gt...Ct...A+': 1,'C....G....T+Gc...C....A+': 1,'T....Taa..T+G....Ag...A+': 1,'T....Gc...T+G....Ca...A+': 1,'Tca..G....T+G....C....A+': 1,'Tct..G....T+G....C....A+': 1,'T....G....G+Tt...Cc...A+': 1,'T....Gg...T+G....Cc...A+': 1,'T....A....T+G....Ggg..A+': 1,'T....A....T+Gt...A....A+': 1,'Tt...G....T+G....-....A+': 1,'T....-....-+T....Cgtg.A+': 1,'T....Gtca.T+G....C....A+': 1,'T....Aa...T+G....Gg...A+': 1,'Ta...G....T+G....C....A+': 1,'T....G....T+G....Cg...A+': 1,'T....G....T+A....Ca...A+': 1,'T....A....T+GgtgtCcgttA+': 1,'T....G....T+T....C....A+': 1,'T....G....-+G....C....-+': 1,'Tct..G....T+G....-....A+': 1,'T....G....T+G....C....G+': 1,'T....G....T+-....T....A+': 1,'T....Gg...T+Gg...C....A+': 1,'T....A....T+G....Cgg..A+': 1,'-....G....T+G....-....A+': 1,'Tgg..G....T+Gc...Ca...A+': 1,'T....A....C+G....G....A+': 1,'T....-....-+G....G....A+': 1,'T....-....-+G....C....A+': 1,'-....-....-+G....G....A+': 1,'T....G....T+Gc...Ca...A+': 1,'T....Ca...T+G....G....A+': 1,'T....Ca...T+G....C....A+': 1,'T....A....T+GatttG....A+': 1,'T....A....T+G....Tg...A+': 1,'T....A....T+Gg...C....A+': 1,'T....-....T+G....-....A+': 1,'A....A....T+G....G....A+': 1,'T....-....T+G....Gg...A+': 1,'Ttt..Gg...T+G....Ctct.A+': 1,'T....C....T+G....-....A+': 1,'T....A....T+Ggg..Cc...A+': 1,'T....A....T+Ga...C....-+': 1,'A....-....T+G....Cgg..A+': 1,'T....A....T+Gta..C....A+': 1,'T....G....T+-....-....A+': 1,'T....Gta..T+Gat..CtagtA+': 1,'C....G....T+G....C....A+': 1,'-....G....T+Gg...C....-+': 1,'T....-....T+G....A....A+': 1,'Tat..G....T+G....G....A+': 1,'T....Gc...T+T....C....A+': 1,'T....Ta...T+G....-....A+': 1,'G....G....T+G....C....A+': 1,'T....G....T+-....C....-+': 1,'T....A....T+G....Gggg.A+': 1,'T....T....T+G....G....A+': 1,'T....-....T+G....C....-+': 1,'T....A....T+Gg...C....-+': 1,'T....A....T+Gc...Cg...A+': 1,'T....G....T+Ga...G....A+': 1,'T....A....T+Ga...A....A+': 1,'T....G....T+Gggg.C....A+': 1,'A....A....T+Gt...G....A+': 1,'-....-....T+G....C....A+': 1,'T....Gt...T+T....C....A+': 1,'Ta...A....T+G....T....A+': 1,'-....G....T+A....A....A+': 1,'T....-....T+Gc...Cg...A+': 1,'T....G....T+G....C....T+': 1,'T....-....T+-....-....A+': 1,'T....G....T+-....C....C+': 1,'T....G....T+Gt...C....-+': 1 };"
for hh in hpdat:
    print hh
