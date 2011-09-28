#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/defs.h"
#include "include/interval.h"

class Block {
public:
	vector<Interval> intv;
	vector<bool>     inv;
};
		

void getBlockBps(vector<Block> & blocks, vector<int> & bps) {
	for (int i = 0; i < blocks.size(); i++) {
		Block b = blocks[i];
		for (int idx = 0; idx < b.intv.size(); idx++) {
			bps.push_back(b.intv[idx].start); 
			bps.push_back(b.intv[idx].end + 1); 
		}
	}
}


void split_block(Block b, vector<int> & bps, vector<Block> & result, int main_idx) {
	if (main_idx >= b.intv.size()) {
		result.push_back(b);
		return;
	}

	int bp_start, bp_end;
	Interval range = b.intv[main_idx];
	searchRange(bps, range, bp_start, bp_end);
	for (int i = bp_start; i < bp_end; i++) {
		int bp = bps[i];
		//assume alignments have no gaps
		//break off the beginning
		if (bp == b.intv[main_idx].start) continue;
		Block newblock = b;
		for (int sec_idx = 0; sec_idx < b.intv.size(); sec_idx++) {
			if (sec_idx == main_idx) continue;
			b.intv[sec_idx].start += bp - b.intv[main_idx].start;
			newblock.intv[sec_idx].end = b.intv[sec_idx].start - 1;
		}
		b.intv[main_idx].start = bp;
		newblock.intv[main_idx].end = b.intv[main_idx].start - 1;
		result.push_back(newblock);
	}
	result.push_back(b);
}

void split_blocks(vector<Block> & blocks, vector<int> &bps) {
	// Assume that there are no gaps (this case would be difficult to handle properly and would require more information in the alignment)
	vector<Block> newblocks;
	sort(bps.begin(), bps.end());
	int maxBlockSize = 0;
	for (int i = 0; i < blocks.size(); i++) maxBlockSize = max(maxBlockSize, (int) blocks[i].intv.size());

	for (int main_idx = 0; main_idx < maxBlockSize; main_idx++) {
		for (int i = 0; i < blocks.size(); i++) {
			split_block(blocks[i], bps, newblocks, main_idx);
		}
		blocks = newblocks;
		newblocks.clear();
	}

}


void read_blocks(ifstream & inf, vector<Block> & blocks) {
	string lastIndex = "-1";
	vector<string> row;
	while (get_row_whitespace(inf, row)) {
		ostringstream o;
		bool inv = (row[3] == "0" ? false : true);
		string curSize = row[4]; //not really used
		string curIndex = row[5];
		if (row.size() > 6) { //trailing label -- lets store it in the interval label
			o << row[6];
			for (int i = 7; i < row.size(); i++) o << "\t" << row[i];
		}
		Interval intv(row[0], atoi(row[1].c_str()), atoi(row[2].c_str()), o.str());
		if (lastIndex != curIndex) { // start of new block
			Block b;
			b.intv.push_back(intv);
			b.inv.push_back(inv);
			blocks.push_back(b);
		} else {
			blocks.back().intv.push_back(intv);
			blocks.back().inv.push_back(inv);
		}
		lastIndex = curIndex;
	}
	return;
}

void print_block(ostream & out, Block & b, int index) {
	for (int i = 0; i < b.intv.size(); i++) {
		out << b.intv[i].chr << "\t" << b.intv[i].start << "\t" << b.intv[i].end  << "\t" << b.inv[i] << "\t" << b.intv.size() << "\t" << index;
		if (b.intv[i].label != "") out << "\t" << b.intv[i].label;
		out << endl;
	}
}

//an Interval along with a pointer into its location in a Block vector
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

void getIntervalsFromBlocks(vector<Block> & blocks, vector<IntervalFromBlock> & intervals) {
	for (int i = 0; i < blocks.size(); i++) {
		for (int j = 0; j < blocks[i].intv.size(); j++) {
			intervals.push_back(IntervalFromBlock(blocks[i].intv[j], i, j));
		}
	}
	sort(intervals.begin(), intervals.end());
}


/*
   void split_intervals(vector<Interval> & intervals, vector<int> & bps) {
   vector<Interval> newIntervals;
   sort(bps.begin(), bps.end());
   for (int i = 0; i < intervals.size(); i++) {
   int bp_start, bp_end;
   Interval intv = intervals[i];
   searchRange(bps, intv, bp_start, bp_end);
   for (int i = bp_start; i < bp_end; i++) {
   int bp = bps[i];
   if (bp == intv.start) continue;
   Interval newIntv = intv;
   intv.start = bp;
   newIntv.end = intv.start - 1;
   newIntervals.push_back(newIntv);
   }
   newIntervals.push_back(intv);
   }
   intervals = newIntervals;
   return;
   }
 */


