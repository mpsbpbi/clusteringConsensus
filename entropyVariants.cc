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

  int threshold = atoi(argv[4]);
  float ethreshold = atof(argv[5]);
  int numAboveThreshold = 0;
  int numHP = 0;

  cerr << "at ethreshold=\t" << ethreshold << "\ttreshold=\t" << threshold <<endl;

  vector<int> filteredColumns;

  if (1==1){
    ////////////////////////////////
    // How many columns contain the minimum number of characters?
    int numbases;
    char myc;
    int numAboveThreshold = 0;

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

    for (int i =0; i < goodColumns.size(); i++){
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

  delete[] buffer;
  return 0;
}
