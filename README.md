clusteringConsensus (CluCon)
============================

Goal
----

Your biological sample contains a number of distinct genomes and you
want to know the number of distinct genomes, their exact identities,
and their relative abundances.

CluCon is analysis code that clusters reads and estimates consensus
from PacBio sequencing data from mixed populations to answer these
questions. It is designed for mixed genomes (or genome regions) where
PacBio reads can cover the entire genome.

Example
-------

Here is an example data analysis run that examines a mixture of
near-full length HIV genomes (9kb long): [README_HIV-three-clones.html]("https://s3.amazonaws.com/files.pacb.com/Users/mbrown/HIV-three-clones/README_HIV-three-clones.html")

PacBio reads can sequence entire HIV genomes from single molecules as
single, continuous 9kb+ reads. Given a simple containing an unknown
number of HIV genomes, one PacBio sequencing chip in three hours, and
the CluCon software, we were able to determine that the sample
contained three full-length HIV species in a (60%/20%/20%) mixture and
got the exact genomic identity of all three species.

The calling sequence is simple taking the sequence data and a generic
reference to seed the process:

    ConsensusClusterSubset.py \
    --nproc=1 \
    --runDir=bound-clucon-dna622 \
    --fasta=raw-niaid-dna622.fasta \
    --ref=hiv_hxb2_whole_genome-covered.fasta \
    --spanThreshold=99.0% \
    --entropyThreshold=1.0 \
    --basfofn=dna622.baxh5.fofn

----