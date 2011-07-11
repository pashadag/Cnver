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
#include<queue>
#include<cstdarg>
#include<algorithm>

using namespace std;
ostringstream oMsg;
string sbuf;
string filename, baseFilename;


void usage(int argc, char * argv[]) {
	cout << "Usage: " << argv[0] << " <genome fasta file> <start (zero-based)> <end(inclusive)> <readlen>" << endl;
	exit(1);
}

int main(int argc, char * argv[]) {

	if (argc != 5) usage(argc, argv);
	filename=argv[1];
	int startPos = atoi(argv[2]);
	int endPos = atoi(argv[3]);
	int readLen = atoi(argv[4]);

	ifstream inFile;
	inFile.open(filename.c_str());

	char buf[256];
	char c;
	int curPos=0;
	string genome;

	inFile.get(buf, 256);

	while ((c = inFile.rdbuf()->sbumpc()) != EOF) {
		if (c != '\n') {
			if (curPos >= startPos) {
				genome.push_back(c);
			}
			curPos++;
			if (curPos > (endPos + readLen - 1)) break;
		}
	}
	inFile.close();

	int i;
	for (i = 0; i < genome.size() - readLen + 1; i++) {
		cout << ">SIM_LOC" << startPos + i << endl;
		//cerr << "genome = " << genome << endl;
		cout << genome.substr(i, readLen) << endl;
		//this might be aproblem in the end of the genome
	}


}

