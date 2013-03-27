clusteringConsensus
===================

goal
----

Provide calling code to run clustering consensus on PacBio sequencing
from mixed populations.

methods
-------

This is a simple pull of my calling scripts into github so they can be
shared and tracked.

After installation here is the calling sequence:

    export SEYMOUR_HOME=/opt/smrtanalysis
    source /opt/smrtanalysis/etc/setup.sh
    export PATH=/home/mbrown/Desktop/smrtanalysis/code:$PATH
    time ConsensusClusterSubset.py \
    --runDir testit \
    --fasta 2450417-0003.fasta \
    --ref HIVemory.fasta \
    --spanThreshold=6400 \
    --entropyThreshold=1.0 \
    --basfofn 2450417-0003.bas.fofn \
    > 2450417-0003.workflow.output 2>&1 &
