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
		cerr << "Usage: idx_lookup <normal|double|rmap> [indexFile] <dataFile> <entry>.\n";
		exit(1);
	}

	if (strcmp(argv[1], "double") == 0) {
		cout << binary_search_lookup(argv[2], argv[3], atoi(argv[4]));
	} else if (strcmp(argv[1], "normal") == 0) {
		cout << lookup(argv[2], argv[3], atoi(argv[4]));
	} else if (strcmp(argv[1], "rmap") == 0) {
		Map_t val;
		val.read_id = atol(argv[3]);
		cout << binary_search_lookup_in_binary_file(argv[2], val);
	} else {
		cerr << "idx_lookup: unknown parameter: " << argv[1] << endl;
		exit(1);
	}

	return 0;

}

