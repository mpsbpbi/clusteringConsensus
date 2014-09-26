all: pairwiseAlignDist agreeFracDist entropyVariants

pairwiseAlignDist: pairwiseAlignDist.cc
	g++ -std=c++0x -O2 -o pairwiseAlignDist pairwiseAlignDist.cc -lm

agreeFracDist: agreeFracDist.cc
	g++ -std=c++0x -O2 -o agreeFracDist agreeFracDist.cc -lm

entropyVariants: entropyVariants.cc
	g++ -std=c++0x -O2 -o entropyVariants entropyVariants.cc -lm
