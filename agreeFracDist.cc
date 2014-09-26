#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <algorithm>
#include <string.h>

using namespace std;

int rowcolToIndex(int numcol, int row, int col){
  return(row*numcol+col);
}

int main (int argc, const char** argv) {

  /* argv[1] = "monoNS5ACon1.2450129-0010.dccs.512.msa"
     argv[2] = numrow
     argv[3] = numcol
     argv[4] = threshold (Minimum number of characters in column to keep it.)
     argv[5] = entropy threshold
     argv[6] = max insert size in alignment. from cmph5ToMSAMaxInserts.py --maxInsert=4
     argv[7] = "overlap" or "full" (or null). is the distance based on overlap region or full.
   */

  int length;
  char *buffer;

  ifstream is;
  is.open(argv[1], ios::binary );

  // get length of file:
  is.seekg(0, ios::end);
  length = is.tellg();
  is.seekg(0, ios::beg);

  // allocate memory:
  buffer = new char [length];

  // read data as a block: TODO mmap it
  is.read(buffer,length);
  is.close();

  int numrow = atoi(argv[2]);
  int numcol = atoi(argv[3])+1;
  int maxInsert = atoi(argv[6]);

  int doOverlap;
  if (argc<8){
    doOverlap=0;
  } else {
    if (strcmp(argv[7],"overlap")==0){
      doOverlap=1;
    } else {
      doOverlap=0;
    }
  }
  cerr << "doOverlap: " << doOverlap << endl;

  cerr << argv[1] << endl;
  cerr << "length= " << length << " should be " << numrow << "*" << numcol << "=" << (numrow*numcol) << endl;
  if (length != (numrow*numcol)){
    exit(1);
  }

  ////////////////////////////////
  // look for the file ""distjob.usecols". if exists read in the
  // columns to be used, otherwise do the work to find high entropy
  // columns and write the results

  int threshold = atoi(argv[4]);
  float ethreshold = atof(argv[5]);
  int numAboveThreshold = 0;
  int numHP = 0;

  cerr << "at ethreshold=\t" << ethreshold << "\ttreshold=\t" << threshold <<endl;

  vector<int> filteredColumns;

  ifstream usecolsfile("distjob.usecols");
  if (usecolsfile.good()){
    ////////////////////////////////
    // read in list of integers specifying the columns to be used.
    cerr << "LOG: reading in list of columns for clustering from 'distjob.usecols'" <<endl;
    int tmpint;
    while (usecolsfile>>tmpint){
      filteredColumns.push_back(tmpint);
    }
  } else {
    cerr << "ERROR! distjob.usecols is not readable" << endl;
    exit(1);
  }

  ////////////////////////////////
  // for all pairs of sequences, compute distances on all filteredColumns
  // All delete mismatches  count as matches.
  unordered_map<string, double> dist;  // "<from><to>"
  dist["AA"] = 0.0;
  dist["CA"] = 1.0;
  dist["GA"] = 1.0;
  dist["TA"] = 1.0;
  dist["-A"] = 1.0;
  dist[" A"] = 0.1;
  
  dist["AC"] = 1.0;
  dist["CC"] = 0.0;
  dist["GC"] = 1.0;
  dist["TC"] = 1.0;
  dist["-C"] = 1.0;
  dist[" C"] = 0.1;
  
  dist["AG"] = 1.0;
  dist["CG"] = 1.0;
  dist["GG"] = 0.0;
  dist["TG"] = 1.0;
  dist["-G"] = 1.0;
  dist[" G"] = 0.1;
  
  dist["AT"] = 1.0;
  dist["CT"] = 1.0;
  dist["GT"] = 1.0;
  dist["TT"] = 0.0;
  dist["-T"] = 1.0;
  dist[" T"] = 0.1;
  
  dist["aa"] = 0.0;
  dist["ca"] = 1.0;
  dist["ga"] = 1.0;
  dist["ta"] = 1.0;
  dist[".a"] = 1.0;
  dist[" a"] = 0.1;
  
  dist["ac"] = 1.0;
  dist["cc"] = 0.0;
  dist["gc"] = 1.0;
  dist["tc"] = 1.0;
  dist[".c"] = 1.0;
  dist[" c"] = 0.1;
  
  dist["ag"] = 1.0;
  dist["cg"] = 1.0;
  dist["gg"] = 0.0;
  dist["tg"] = 1.0;
  dist[".g"] = 1.0;
  dist[" g"] = 0.1;
  
  dist["at"] = 1.0;
  dist["ct"] = 1.0;
  dist["gt"] = 1.0;
  dist["tt"] = 0.0;
  dist[".t"] = 1.0;
  dist[" t"] = 0.1;
  
  dist["A-"] = 1.0;
  dist["C-"] = 1.0;
  dist["G-"] = 1.0;
  dist["T-"] = 1.0;
  dist["--"] = 0.0;
  dist[" -"] = 0.1;
  
  dist["a."] = 1.0;
  dist["c."] = 1.0;
  dist["g."] = 1.0;
  dist["t."] = 1.0;
  dist[".."] = 0.0;
  dist[" ."] = 0.1;
  
  dist["A "] = 0.1;
  dist["C "] = 0.1;
  dist["G "] = 0.1;
  dist["T "] = 0.1;
  dist["a "] = 0.1;
  dist["c "] = 0.1;
  dist["g "] = 0.1;
  dist["t "] = 0.1;
  dist["- "] = 0.1;
  dist[". "] = 0.1;
  dist["  "] = 0.1;

  int thisc;
  double doverlap,doverlapZ,dglobal,dglobalZ;
  char cc1;
  char cc2;
  string key;
  double ratio;
  double numer, denom;

  for (int ii=0; ii<(numrow-1); ii++){
    for (int jj=(ii+1); jj<numrow; jj++){

      doverlap = 0.0;
      doverlapZ = 0.0;
      dglobal = 0.0;
      dglobalZ = 0.0;

      for (int cc=0; cc<filteredColumns.size(); cc++){
	thisc = filteredColumns[cc];
	key = "";
	cc1 = buffer[rowcolToIndex(numcol,ii,thisc)];
	cc2 = buffer[rowcolToIndex(numcol,jj,thisc)];
	key += cc1;
	key += cc2;
	if ((cc1 != ' ') && (cc2 != ' ')){
	  // information at both
	  doverlap += dist[key];
	  doverlapZ += 1.0;
	}
	dglobal += dist[key];
	dglobalZ += 1.0;
      }

      if (doOverlap != 1){
	// global
	numer = dglobal;
	denom = dglobalZ;
	ratio = dglobal/dglobalZ;
      } else {
	// overlap
	numer=doverlap;
	denom = doverlapZ;
	if (doverlapZ == 0.0){
	  // no information, assume identical
	  ratio=0.0;
	} else {
	  ratio = doverlap/doverlapZ;
	}
      }
      cout << ii << "\t" << jj << "\t" << ratio << "\t" << dglobal << "\t" << dglobalZ << "\t" << doverlap << "\t" << doverlapZ << endl;
    }
  }

  delete[] buffer;
  return 0;
}
