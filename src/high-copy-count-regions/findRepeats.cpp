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

// takes as stdin a file with every line a used location, using 0-based. Not necessarily sorted
// outputs all intervals which do not contain any used locations
// intervals are 1-based and closed;
// the output is also sorted 

void usage(int argc, char * argv[]) {
	cout << "Usage: " << argv[0] << " <numLocations> " << endl;
	exit(1);
}

int main(int argc, char * argv[]) {
	if (argc != 2) usage(argc, argv);
	int numLocs = atoi(argv[1]);
	vector<bool> used;
	used.resize(numLocs, false);
	while (getline(cin, sbuf)) {
		int num = atoi(sbuf.c_str());
		if ((num < 0) || (num >= numLocs)) {
			cerr << "Read line: " << sbuf << endl;
			cerr << "The extracted number " << num << " is invalid." << endl;
			cout << "numLocs = " << numLocs << endl;
			cerr << "This line will be ignored." << endl;
		} else {
			used[num] = true;
		}
	}
	bool lastUsed = true;
	int curStart = 0;
	for (int i = 0; i < used.size(); i++) {
		if (used[i]) {
			if (lastUsed) {
				continue;
			} else {
				cout << curStart << "\t" << i - 1 << endl;
			}
		} else {
			if (!lastUsed) {
				continue;
			} else {
				curStart = i;
			}
		}
		lastUsed = used[i];
	}
	if (!lastUsed) {
		cout << curStart << "\t" << numLocs - 1 << endl;
	}
}

