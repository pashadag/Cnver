#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "../include/defs.h"
#include "../include/interval.h"
#include "../include/block.h"


string sbuf;
vector<string> row;
vector<Block> blocks;
vector<IntervalFromBlock> intervals;
vector<double> blockFlow;
int ploidy;

void process_interval(int start, int end) {
	for (int i = start; i <= end; i++){
		int blockIdx = intervals[i].blockIdx;
		tokenize(intervals[i].intv.label, row);
		double cov = atof(row[0].c_str());
		double exp = atof(row[1].c_str());
		int    len = atoi(row[2].c_str());
		cout << intervals[i].intv.start << "\t" << intervals[i].intv.end << "\t";
		cout << blockFlow[blockIdx] << "\t";
		cout << blocks[blockIdx].intv.size() * ploidy << "\t";
		cout << intervals[i].intv.end - intervals[i].intv.start + 1 << "\t";
		cout << len << "\t";
		cout << cov / exp << "\t";
		cout << intervals[i].intv.start - intervals[start].intv.start << "\t";
		if (blocks[blockIdx].inv[intervals[i].intIdx]) cout << "-";
		cout << blockIdx << endl;
	}
}



int main(int argc, char ** argv) {

	if (argc != 6 && argc != 7) {
		fprintf(stderr, "usage: %s <block_file> <solution_file> ploidy [interval]\n", argv[0]);
		fprintf(stderr, "       interval can be without an end\n", argv[0]);
		return -1;
	}

	ploidy = atoi(argv[3]);
	//read in blocks
	ifstream inf;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	inf.close();
	getIntervalsFromBlocks(blocks, intervals);


	//load flow of blocks from solution file
	blockFlow.resize(blocks.size(), 0);
	open_file(inf, argv[2]);
	while (get_row_whitespace(inf, row, sbuf)) {
		if (row[0] == "block") {
			double flow  = atof(row[1].c_str());
			int blockIdx = atoi(row[2].c_str());
			blockFlow.at(blockIdx) = int(flow);
		}
	}
	inf.close();

	int start, end, dummy;
	find_pos(intervals, atoi(argv[5]), start, dummy);
	if (argc == 7) {
		find_pos(intervals, atoi(argv[6]), end,   dummy);
	} else {
		end = intervals.back().intv.end;	
	}
	process_interval(start, end);

	return 0;
}
