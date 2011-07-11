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

using namespace std;
#include "include/defs.h"


void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " fasta_file ref_name ";
	exit(1);
}

int main(int argc, char * argv[]) {
	if (argc != 3) usage(argc, argv);
	string genome = read_genome(argv[1]);
	string ref_name = argv[2];
	bool inside = false;
	int start = 0;


	//the interval file should be 1-based
	for (long i = 0; i < genome.length(); i++) {
		long one_based_i = i + 1;
		
		if (genome[i] == 'n' || genome[i] == 'N') {
			if (inside) {
				continue;
			} else {
				inside = true;
				start = one_based_i;
			}
		} else {
			if (!inside) {
				continue;
			} else {
				cout << ref_name << "\t" << start << "\t" << one_based_i - 1 << endl;
				inside = false;
			}
		}
	}







}

