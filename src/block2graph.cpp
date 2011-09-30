#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/defs.h"
#include "include/interval.h"
#include "include/block.h"


#define MAX_ADJ_ARC_CAPACITY 8888
#define EPSILON_FLOW_COST 1

int block2node(int blockIndex, int end) { 
	if (end == 0) 
		return blockIndex * 2 + 1; 
	else 
		return blockIndex * 2 + 2;
}


vector<Block> blocks;
vector<IntervalFromBlock> intervals;

int main(int argc, char ** argv) {

	if (argc != 5) {
		fprintf(stderr, "usage: %s <block_file> <link_file> <problem_file> <base_flow>\n", argv[0]);
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
	ofstream placedLinksOut;
	open_file(placedLinksOut, make_string(argv[2]) + ".placed");

	int baseFlow = atoi(argv[4]);

	//output nodes
	ofstream outf;
	open_file(outf, argv[3]);
	for (int i = 0; i < blocks.size(); i++) {
		outf << "n " << block2node(i,0) << " 0" << endl;
		outf << "n " << block2node(i,1) << " 0" << endl;
	}


	//output sequence edges
	for (int i = 0; i < blocks.size(); i++) {
		Block b = blocks[i];
		double totLen = 0;
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
		double avgLen = totLen / b.intv.size();
		double expDocPerPos;
		if (totLen > 0) {
			expDocPerPos = totExp / (baseFlow * totLen);
		} else {
			expDocPerPos = 0;
		}
		outf << "convex_arc " << block2node(i,0) << " " << block2node(i,1) << " ";
		outf << totCov << " " << avgLen << " " << expDocPerPos << endl; 
	}

	//output reference adjacency edges
	getIntervalsFromBlocks(blocks, intervals);
	for (int i = 0; i < intervals.size(); i++) {
		int curBlock  = intervals[i].blockIdx;
		int curInv    = blocks.at(curBlock).inv.at(intervals[i].intIdx) ? 1 : 0;
		int nextBlock = intervals[(i+1) % intervals.size()].blockIdx;
		int nextInv   = blocks.at(nextBlock).inv.at(intervals[(i+1) % intervals.size()].intIdx) ? 1 : 0;
		int fromNode  = block2node(curBlock,  abs(curInv - 1)) * (curInv  * -2 + 1); 
		int toNode    = block2node(nextBlock, abs(nextInv))    * (nextInv * -2 + 1);
		outf << "a " << fromNode << " " << toNode << " ";
		if (i == intervals.size() - 1) { //the last arc is cirulcation arc
			outf << baseFlow << " " << baseFlow << " " << EPSILON_FLOW_COST  << endl;
		} else {
			outf << "0 " << MAX_ADJ_ARC_CAPACITY << " " << EPSILON_FLOW_COST << endl;
		}
	}


	//output donor adjacency edges (links)
	for (int i = 0; i < links.size(); i++) {
		Link link = links[i];
		int block1, end1, block2, end2;
		find_pos(blocks, intervals, link.from, block1, end1);
		find_pos(blocks, intervals, link.to, block2, end2);
		int fromNode  = block2node(block1, end1);
		int toNode    = block2node(block2, end2);
		if (end1 == 0) fromNode *= -1;
		if (end2 == 1) toNode   *= -1;
		outf << "a " << fromNode << " " << toNode << " ";
		outf << "0 " << MAX_ADJ_ARC_CAPACITY << " " << EPSILON_FLOW_COST << endl;
		placedLinksOut << links[i] << "\t" << block1 << "\t" << block2 << "\t" << end1 << "\t" << end2 << endl;
	}

	outf.close();
	placedLinksOut.close();
	return 0;
}
