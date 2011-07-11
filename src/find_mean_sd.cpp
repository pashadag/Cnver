#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


#include "include/dbtypes.h"
#include "include/defs.h"
#include "include/InputReader.h"
#include "bamtools/BamReader.h"
#include "bamtools/BamAux.h"

using namespace BamTools;

string sbuf;
MeanSD calculator;
int lowd, highd, read_len;


void usage( int argc, char ** argv) {
	cerr << "Usage: " << argv[0] << " mapping file input_format(bam|txt) lower_limit upper_limit read_len.\n";
	cerr << "Estimates the mean and sd from a single mapping file.\n";
	cerr << "Matepairs mapping uniquely within the given distance limits are used for the estimate.\n";
	printf("Usage: %s mappings\n", argv[0]);

	cout << endl;
	exit(1);
}

void judge(Map_t & map1, Map_t & map2) {
	int dist;
	if (strcmp(map1.ref_name, map2.ref_name) == 0 && map1.orientation != map2.orientation) {
		dist = map2.ref_pos + read_len - map1.ref_pos;
		if (dist < highd && dist > lowd) calculator.addVal(dist);
		//cout << "[ADDED]\t" << map1 << "\t" << map2 << endl;
	} else {
		//cout << "[DROPPED_1]\t" << map1 << "\t" << map2 << endl;
	}
	return;
}

int main ( int argc, char ** argv) {
	if (argc != 6) usage(argc, argv);
	string map_file = argv[1];
	string input_format = argv[2];
	lowd = atoi(argv[3]);
	highd = atoi(argv[4]);
	read_len = atoi(argv[5]);

	//read in the read mappings, grouping them into blocks with the same read_id.  
	//Each block is then handled in process_block
	uint64_t lastReadId = -1;
	Map_t mapt;
	deque<Map_t> block;
	string lastReadName = "2paclives";
	InputReader reader(input_format, map_file);
	Map_t map1, map2;
	bool foundLeft = false;
	while (reader.getNext(mapt)) {
		if (lastReadId != mapt.read_id) {
			if (lastReadId != -1) { //process  block
				if (lastReadId % 2 == 0) { //block of first mate done
					if (block.size() == 1 && 0 != strcmp(block[0].ref_name, "-1")) { //potentialy good unique matepair
						map1 = block[0];
						foundLeft = true;
					} else {
						foundLeft = false;
					}
				} else { //block of second mate done
					if (foundLeft && lastReadId == map1.read_id + 1 && block.size() == 1 && 0 != strcmp(block[0].ref_name, "-1")) { //potentialy good unique matepair
						map2 = block[0];
						judge(map1, map2);
					}
					foundLeft = false;
				}
			}

			block.clear();
		} 
		lastReadId = mapt.read_id;
		block.push_back(mapt);
	}

	cout << "Mean:\t" << calculator.getMean() << endl;
	cout << "Standard error:\t" << calculator.getSD() << endl;
	cout << "Number of samples:\t" <<  calculator.numpoints << endl;

	return 0;	
}
