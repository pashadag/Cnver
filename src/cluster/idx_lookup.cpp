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

#include "../include/dbtypes.h"

using namespace std;

bool operator==(const Map_t & x, const Map_t & y) {
	return x.read_id  == y.read_id;
}

bool operator<(const Map_t & x, const Map_t & y) {
	return x.read_id  < y.read_id;
}

#include "index.h"
int main(int argc, char * argv[]) {
	if (argc < 4) {
		cerr << "Usage: idx_lookup <normal|double|rmap|bamnames> [indexFile] <dataFile> <entry>.\n";
		exit(1);
	}

	string type = argv[1];

	if (type == "double") {
		cout << binary_search_lookup(argv[2], argv[3], atoi(argv[4]));
	} else if (type == "normal") {
		cout << lookup(argv[2], argv[3], atoi(argv[4]));
	} else if (type == "rmap") {
		Map_t val;
		val.read_id = atol(argv[3]);
		cout << binary_search_lookup_in_binary_file(argv[2], val);
	} else if (type == "bamnames") {
		char readname[MAX_BAM_READNAME];
		ifstream dataFile;
		open_file_binary(dataFile, argv[2]);
		dataFile.seekg(MAX_BAM_READNAME * atol(argv[3]));
		dataFile.read((char *) &readname, MAX_BAM_READNAME);
		cout << readname << endl;
	} else {
		cerr << "idx_lookup: unknown parameter: " << type << endl;
		exit(1);
	}

	return 0;

}

