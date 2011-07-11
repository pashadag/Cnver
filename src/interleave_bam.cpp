#include<cassert>
#include<cmath>
#include<cstring>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<cstdarg>
#include<algorithm>


#include "include/defs.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "bamtools/BamWriter.h"

using namespace std;
using namespace BamTools;

BamReader reader1;
BamReader reader2;
BamWriter writer;

void usage(int argc, char * argv[]) {
	cout << "Usage: " << argv[0] << " bamfile1 bamfile2 outbamfile\n";
	cout << "\tInterleaves two bamfiles.  So first it takes all the mappings for the first read in bamfile1, then all the mappings for the first read in bamfile2, ";
	cout << "then all the mappings for the second read in bamfile1, and so on.  The bamfiles must have the same amount of reads. The name of each mapping is appended with a suffix of _1 or _2, depending on which file it comes from.\n";
	exit(1);
}

bool process_block(BamReader & reader, BamAlignment & leftover, string nameSuffix) {
	BamAlignment cur = leftover;

	//read the first alignment if there are no leftovers from last block
	if (cur.Name == "2paclives") {
		if (!reader.GetNextAlignment(cur)) {
			cerr << "Error: empty BAM file.\n";
			exit(1);
		}
	}
	string firstName = cur.Name;
	bool stillgoing = true;
	do {
		cur.Name = cur.Name + nameSuffix;
		writer.SaveAlignment(cur);
		stillgoing = reader.GetNextAlignment(cur);
	} while (stillgoing && cur.Name == firstName);
	
	leftover = cur;
	return stillgoing;
}

int main(int argc, char * argv[]) {

	if (argc != 4) usage(argc, argv);
	reader1.Open(argv[1]);
	reader2.Open(argv[2]);
	writer.Open(argv[3], reader1.GetHeaderText(), reader1.GetReferenceData());

	BamAlignment leftover1;
	BamAlignment leftover2;
	leftover1.Name = "2paclives";
	leftover2.Name = "2paclives";
	bool stillgoing = true;

	while (stillgoing) {
		stillgoing = process_block(reader1, leftover1, "_1");
		process_block(reader2, leftover2, "_2");
	}
	reader1.Close();
	reader2.Close();
	writer.Close();
	return 0;
}

