#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/defs.h"
#include "include/interval.h"
#include "include/block.h"

string sbuf;
int ploidy;
vector<int> blockFlow;
vector<int> tempFlow;
vector<Block> blocks;
vector<IntervalFromBlock> intervals;
vector<int> walkCounts;

//removes as much flow along a walk as possible
//assumes that at least one unit can be removed
int removeMaxFlow(int start, int len, int call) {
	vector<int> walkBlocks;
	for (int i = start; i < start + len; i++) {
		int blockIdx = intervals[i].blockIdx;
		walkBlocks.push_back(blockIdx);
	}

	//fill in the relevant fields of walkCounts...the point is that we don't just
	//want to clear this array because it might be huge.
	for (int i = 0; i < walkBlocks.size(); i++)  walkCounts[walkBlocks[i]] = 0;
	for (int i = 0; i < walkBlocks.size(); i++)  walkCounts[walkBlocks[i]]++;
	
	//remove duplicates
	sort(walkBlocks.begin(), walkBlocks.end());
	walkBlocks.erase(unique(walkBlocks.begin(), walkBlocks.end()), walkBlocks.end());

	bool stillgood = true;
	int unitsRemoved = 0;
	while (stillgood) {
		for (int i = 0; i < walkBlocks.size(); i++) {
			//start; i < start + len; i++) {
			//int blockIdx = intervals[i].blockIdx;
			int blockIdx = walkBlocks[i];
			int refCount = blocks[blockIdx].intv.size() * ploidy;

			if (call > 0) {
				assert (blockFlow.at(blockIdx) >= refCount);
			} else if (call < 0) {
				assert (blockFlow.at(blockIdx) <= refCount); 
			}

			int deb = blockFlow.at(blockIdx);
			blockFlow.at(blockIdx) -= call * walkCounts[blockIdx];

			if (call > 0) {
				assert (blockFlow.at(blockIdx) >= refCount);
			} else if (call < 0) {
				assert (blockFlow.at(blockIdx) <= refCount); 
			}

			if (call > 0) {
				if (blockFlow.at(blockIdx) - walkCounts[blockIdx] < refCount) stillgood = false;
			} else if (call < 0) {
				if (blockFlow.at(blockIdx) + walkCounts[blockIdx] > refCount) stillgood = false;
			}

		}
		unitsRemoved++;
	}

	return unitsRemoved;
}


void restoreUnitFlow(int start, int len, int change) {
	for (int i = start; i < start + len; i++) {
		int blockIdx = intervals[i].blockIdx;
		blockFlow[blockIdx] += change;
	}
}

//finds the longest walk, leaving the flow's unchanged
void findLongestPath(int & maxStart, int & maxBlockLen, int & maxBpLen, int & maxCall) {
	maxBpLen = 0;
	int start;
	int bpLen = 0;
	int lastCall = 0;
	for (int i = 0; i < intervals.size(); i++) {
		int blockIdx = intervals[i].blockIdx;
		int curFlow = blockFlow.at(blockIdx);
		int refCount = blocks[blockIdx].intv.size() * ploidy;
		int curCall;
		if (curFlow == refCount)     
			curCall = 0;
		else if (curFlow > refCount) 
			curCall = 1;
		else if (curFlow < refCount) 
			curCall = -1;

		if (curCall != lastCall && lastCall != 0) { //end path 
			restoreUnitFlow(start, i - start, lastCall); //restore original flow
			if (bpLen > maxBpLen) {
				maxBpLen = bpLen;
				maxBlockLen = i - start;
				maxStart = start;
				maxCall = lastCall;
			}
		} 

		if (curCall != 0) { //either start a new path or continue an old one
			if (curCall != lastCall) { //start a new path
				start = i;
				bpLen = 0;
			} 
			blockFlow.at(blockIdx) -= curCall;
			bpLen += intervals[i].intv.end - intervals[i].intv.start + 1;
		}
		lastCall = curCall;
	}

	//fix up last path
	if (lastCall != 0) {
		restoreUnitFlow(start, intervals.size() - start, lastCall); 
	}


}

int main(int argc, char ** argv) {

	if (argc != 5) {
		fprintf(stderr, "usage: %s <block_file> <solution_file> <ploidy> <min_path_len>\n", argv[0]);
		return -1;
	}

	ploidy = atoi(argv[3]);
	int minPathLen = max(1, atoi(argv[4]));

	//read in blocks
	ifstream inf;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	inf.close();
	getIntervalsFromBlocks(blocks, intervals);


	//load flow of blocks from solution file
	blockFlow.resize(blocks.size(), 0); 
	walkCounts.resize(blocks.size(), 0); 
	tempFlow.resize(blocks.size(), 0); 
	open_file(inf, argv[2]);
	vector<string> row;
	while (get_row_whitespace(inf, row, sbuf)) {
		if (row[0] == "block") {
			double flow  = atof(row[1].c_str());
			int blockIdx = atoi(row[2].c_str());
			if (int(flow) != flow) { //get rid of half-integral flow by moving it closer to the ploidy
				if (flow > ploidy) {
					flow -= 0.5;
				} else {
					flow += 0.5;
				}
			}
			blockFlow.at(blockIdx) = int(flow);
		}
	}

	//find long walks and report them
	int start, blockLen, bpLen, call;
	findLongestPath(start, blockLen, bpLen, call);
	while (bpLen >= minPathLen) {

		//decrease flow along path and find the size of the call;
		int amountToCall = removeMaxFlow(start, blockLen, call);

		//output path
		cout << intervals[start].intv.chr << " " << intervals[start].intv.start << "\t";
		cout << intervals[start + blockLen - 1].intv.end << " " << call * amountToCall << endl;

		//and do it all again
		findLongestPath(start, blockLen, bpLen, call);
	}

	return 0;
}
