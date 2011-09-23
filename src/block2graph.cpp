#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <functional>
#include "include/defs.h"
#include "include/interval.h"
#include "include/block.h"


int block2node(int blockIndex, int end) { 
	if (end == 0) 
		return blockIndex * 2 + 1; 
	else 
		return blockIndex * 2 + 2;
}



//an interval along with a pointer into its location in the blocks
class IntervalFromBlock {
public:
	Interval intv;
	int blockIdx;
	int intIdx;
	IntervalFromBlock(Interval _intv, int _blockIdx, int _intIdx) : intv(_intv), blockIdx(_blockIdx), intIdx(_intIdx) {}
	bool operator< (const IntervalFromBlock & other) const {
		return intv < other.intv;
	}
};


vector<Block> blocks;
vector<IntervalFromBlock> intervals;

bool comp_lt(const IntervalFromBlock & ifb, const int & pt ) {
	return (ifb.intv.end < pt);
}

void find_pos(int pos, int & block, int & end) {
	vector<IntervalFromBlock>::iterator it = lower_bound(intervals.begin(), intervals.end(), pos, comp_lt);
	assert(it != intervals.end());
	if (it->intv.start == pos) {
		end = 0;
	} else if (it->intv.end == pos) {
		end = 1;
	} else {
		assert(0);
	}
	block = it->blockIdx;
	bool inv = blocks.at(it->blockIdx).inv.at(it->intIdx);
	if (inv) end = abs(1 - end);
}

int main(int argc, char ** argv) {

	if (argc != 4) {
		fprintf(stderr, "usage: %s <block_file> <link_file> <problem_file>\n", argv[0]);
		return -1;
	}

	//read in blocks
	ifstream inf;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	inf.close();

	//read in links
	vector<Link> links;
	open_file(inf, argv[2]);
	read_links(inf, links);
	inf.close();


	//output nodes
	ofstream outf;
	open_file(outf, argv[3]);
	int numNodes = 2 * blocks.size();
	outf << numNodes << endl;
	for (int i = 0; i < blocks.size(); i++) {
		outf << "NODE\t" << block2node(i,0) << "\t0" << endl;
		outf << "NODE\t" << block2node(i,1) << "\t0" << endl;
	}


	//output sequence edges
	for (int i = 0; i < blocks.size(); i++) {
		Block b = blocks[i];
		int totLen = 0;
		double totExp = 0;
		double totCov = 0;
		for (int j = 0; j < b.intv.size(); j++) {
			istringstream is(b.intv[j].label);
			int curLen;
			double curExp, curCov;
			is >> curCov >> curExp >> curLen;
			totLen += curLen;
			totExp += curExp;
			totCov += curCov;
		}
		outf << "ARC\t" << block2node(i,0) << "\t" << block2node(i,1) << "\t" << totCov << "\t" << totLen << "\t" << totExp << endl; //not sure about the last one
	}

	//output reference adjacency edges
	for (int i = 0; i < blocks.size(); i++) {
		for (int j = 0; j < blocks[i].intv.size(); j++) {
			intervals.push_back(IntervalFromBlock(blocks[i].intv[j], i, j));
		}
	}
	sort(intervals.begin(), intervals.end());
	for (int i = 0; i < intervals.size(); i++) {
		int curBlock  = intervals[i].blockIdx;
		int curInv    = blocks.at(curBlock).inv.at(intervals[i].intIdx) ? 1 : 0;
		int nextBlock  = intervals[(i+1) % intervals.size()].blockIdx;
		int nextInv    = blocks.at(nextBlock).inv.at(intervals[(i+1) % intervals.size()].intIdx) ? 1 : 0;
		if (i == intervals.size() - 1) { //the last arc is cirulcation arc
			outf << "ARC\t" << block2node(curBlock, abs(curInv - 1) ) << "\t" << block2node(nextBlock, nextInv) << "\t0\t-1\t0" << endl;
		} else {
			outf << "ARC\t" << block2node(curBlock, abs(curInv - 1) ) << "\t" << block2node(nextBlock, nextInv) << "\t0\t0\t0" << endl;
		}
	}


	//output donor adjacency edges (links)
	for (int i = 0; i < links.size(); i++) {
		Link link = links[i];
		int block1, end1, block2, end2;
		find_pos(link.from, block1, end1);
		find_pos(link.to, block2, end2);
		block1 *= end1 * 2 - 1;
		block2 *= -1 * (end2 * 2 - 1);
		outf << "ARC\t" << block1 << "\t" << block2 << "\t0\t0\t0" << endl;
	}

	return 0;
}
