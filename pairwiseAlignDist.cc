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
  //cout.write (buffer,length);

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
    ////////////////////////////////
    if (0){
      // Now go through the alignment counting characters
      int numbases;
      char myc;
      for (int cc=0; cc<(numcol-1); cc++){
	numbases = 0;
	for (int rr=0; rr<numrow; rr++){
	  myc = buffer[rowcolToIndex(numcol,rr,cc)];
	  if ( ((int)myc >= 65) && ((int)myc <= 122) ){
	    numbases ++;
	  }
	}
	cout << cc << "\t" << numbases << endl;
      }
    }

    ////////////////////////////////
    // How many columns contain the minimum number of characters?
    int numbases;
    char myc;
    int numAboveThreshold = 0;

    //// hash
    // unordered_map<int,int> goodColumns;
    //     goodColumns[cc] = 1;
    // for (unordered_map<int,int>::iterator it = goodColumns.begin(); it != goodColumns.end(); it++){
    //   cout << it->first << " " << it->second << endl;
    // }
    // cout << goodColumns[15034] << endl; // 1
    // cout << goodColumns[4] << endl;     // 0

    vector<int> goodColumns;
    double counts[5]; // ACGT-
    double freqs[5]; // ACGT-
    double sumCounts;
    double entropy;

    // output ranking of all positions so I can compute ROC of feature
    // selection
    ofstream featureRanking("distjob.ranking");
    featureRanking << "column\tentropy\tnumbases\tfreq\taccepted" << endl;

    int doMatchOnly = 0; // Discard all insert columns in the alignment and only look at match.

    for (int cc=0; cc<(numcol-1); cc++){
      if ( (doMatchOnly ==1) && ( (cc % (maxInsert+1)) != 0) ){ continue; } // Skip all non-match positions

      numbases = 0;
      for (int ii=0; ii<5; ii++){ counts[ii]=0.01;}

      for (int rr=0; rr<numrow; rr++){
	myc = buffer[rowcolToIndex(numcol,rr,cc)];

	if (myc=='A'||myc=='a'){ counts[0] += 1.0; }
	if (myc=='C'||myc=='c'){ counts[1] += 1.0; }
	if (myc=='G'||myc=='g'){ counts[2] += 1.0; }
	if (myc=='T'||myc=='t'){ counts[3] += 1.0; }
	if (myc=='-'||myc=='.'){ counts[4] += 1.0; }

	if ( ((int)myc >= 65) && ((int)myc <= 122) ){
	  numbases ++;
	}
      }

      // compute entropy
      sumCounts=0.0;
      int tot=5; // 4 for bases and 5 to include delete in entropy

      for (int ii=0; ii<tot; ii++){ sumCounts += counts[ii];}
      for (int ii=0; ii<tot; ii++){ freqs[ii] = counts[ii]/sumCounts;}
      entropy = 0.0;
      for (int ii=0; ii<tot; ii++){ entropy += -freqs[ii]*log(freqs[ii])/log(2.0);}

      // if (cc==3964){
      //   cerr << entropy << "::" << numbases << "::" << counts[0] << "," << counts[1] << "," << counts[2] << "," << counts[3] << "," << counts[4] << endl;
      // }

      // sort to get minor frequency 
      int elements = sizeof(freqs)/sizeof(freqs[0]);
      sort(freqs, freqs+elements);
      
      // output to ranking
      featureRanking << cc << "\t" << entropy << "\t" << numbases << "\t" << freqs[4];

      if ( (numbases>threshold) && (entropy>ethreshold) ){
	// if (freqs[3] > 0.08){ // 2nd minor frequency > XX%, rather than entropy
	numAboveThreshold++;
	cerr << "cc " << cc << " numbases " << numbases << " entropy " << entropy << " freqs[4] " << freqs[4] << endl;
	goodColumns.push_back(cc);
	featureRanking << "\t1" << endl;
      } else {
	featureRanking << "\t0" << endl;
      }
    }
    featureRanking.close();

    ////////
    // remove all columns that are homopolymer positions in the
    // reference (row==0), so cut down on HP errors affecting the
    // distances.

    unordered_map<int,int> HPCol;
    //     goodColumns[cc] = 1;
    // for (unordered_map<int,int>::iterator it = goodColumns.begin(); it != goodColumns.end(); it++){
    //   cout << it->first << " " << it->second << endl;
    // }
    // goodColumns.count(x)

    char prev,curr,next;
    for (int pp = 1; pp < (numcol-1); pp++){
      prev = buffer[rowcolToIndex(numcol,0,pp-1)];
      curr = buffer[rowcolToIndex(numcol,0,pp)];
      next = buffer[rowcolToIndex(numcol,0,pp+1)];

      if ((curr == prev) || (curr==next)){
	HPCol[pp] = 1;
	numHP += 1;
      }
    }

    for (int i =0; i < goodColumns.size(); i++){
      //    if (HPCol.count(goodColumns[i]) ==0){
      // Take all columns not just non-HP columns now that deletes count less
      if (1){
	filteredColumns.push_back(goodColumns[i]);
      }
    }

    // write the results to "distjob.usecols"
    ofstream outusecolsfile("distjob.usecols");
    for (int i =0; i < filteredColumns.size(); i++){
      outusecolsfile << filteredColumns[i] << endl;
    }
    outusecolsfile.close();

    cerr << "numAboveThreshold=\t" << numAboveThreshold << "\tnumHP=\t" << numHP << endl;
  }

  // got the columns either by reading the file or computing it
  cerr << "filteredColumns.size()=\t"  << filteredColumns.size() << endl;
  for (int i =0; i < filteredColumns.size(); i++){
    cerr << filteredColumns[i] << " ";
  }
  cerr << endl;

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
