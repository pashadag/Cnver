#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <set>
#include "defs.h"
#include "interval.h"

class Block {
public:
	vector<Interval> intv;
	vector<bool>     inv;
};
		
void getBlockBps(vector<Block> & blocks, set<int> & bps) {
	for (int i = 0; i < blocks.size(); i++) {
		Block b = blocks[i];
		for (int idx = 0; idx < b.intv.size(); idx++) {
			bps.insert(b.intv[idx].start); 
			bps.insert(b.intv[idx].end + 1); 
		}
	}
}
/*
void getBlockBps(vector<Block> & blocks, vector<int> & bps) {
	for (int i = 0; i < blocks.size(); i++) {
		Block b = blocks[i];
		for (int idx = 0; idx < b.intv.size(); idx++) {
			bps.push_back(b.intv[idx].start); 
			bps.push_back(b.intv[idx].end + 1); 
		}
	}
}
void split_block(int blockIdx, vector<int> & bps, vector<int> & newbps, vector<Block> & blocks, int main_idx) {
	assert(blockIdx < blocks.size());
	if (main_idx >= blocks[blockIdx].intv.size()) {
		return;
	}

	int bp_start, bp_end;
	Interval range = blocks[blockIdx].intv[main_idx];
	searchRange(bps, range, bp_start, bp_end);
	for (int i = bp_start; i < bp_end; i++) {
		int bp = bps[i];
		Block * b = &blocks[blockIdx];
		//assume alignments have no gaps
		//break off the beginning
		if (bp == b->intv[main_idx].start) continue;
		Block newblock = *b;
		for (int sec_idx = 0; sec_idx < b->intv.size(); sec_idx++) {
			if (sec_idx == main_idx) continue;
			b->intv[sec_idx].start += bp - b->intv[main_idx].start;
			newblock.intv[sec_idx].end = b->intv[sec_idx].start - 1;
			newbps.push_back(b->intv[sec_idx].start);
		}
		b->intv[main_idx].start = bp;
		newblock.intv[main_idx].end = b->intv[main_idx].start - 1;
		blocks.push_back(newblock);
	}
}


void split_blocks(vector<Block> & blocks, vector<int> &bps, vector<int> &newbps) {
	// Assume that there are no gaps (this case would be difficult to handle properly and would require more information in the alignment)
	sort(bps.begin(), bps.end());
	int maxBlockSize = 0;
	for (int i = 0; i < blocks.size(); i++) maxBlockSize = max(maxBlockSize, (int) blocks[i].intv.size());

	for (int main_idx = 0; main_idx < maxBlockSize; main_idx++) {
		int numBlocks = blocks.size();
		for (int i = 0; i < numBlocks; i++) {
			split_block(i, bps, newbps, blocks, main_idx);
		}
	}

}

*/

void split_block(int blockIdx, set<int> & bps, set<int> & newbps, vector<Block> & blocks, int main_idx) {
	assert(blockIdx < blocks.size());
	if (main_idx >= blocks[blockIdx].intv.size()) {
		return;
	}

	int bp_start, bp_end;
	Interval range = blocks[blockIdx].intv[main_idx];
	//searchRange(bps, range, bp_start, bp_end);
	for ( set<int>::iterator bp = bps.lower_bound(range.start); bp != bps.end() && *bp <= range.end; bp++) {
		Block * b = &blocks[blockIdx];
		//assume alignments have no gaps
		//break off the beginning
		if (*bp == b->intv[main_idx].start) continue;
		Block newblock = *b;
		for (int sec_idx = 0; sec_idx < b->intv.size(); sec_idx++) {
			if (sec_idx == main_idx) continue;
			b->intv[sec_idx].start += *bp - b->intv[main_idx].start;
			newblock.intv[sec_idx].end = b->intv[sec_idx].start - 1;
			bps.insert(b->intv[sec_idx].start);
			newbps.insert(b->intv[sec_idx].start);
		}
		b->intv[main_idx].start = *bp;
		newblock.intv[main_idx].end = b->intv[main_idx].start - 1;
		blocks.push_back(newblock);
	}
}


void split_blocks(vector<Block> & blocks, set<int> &bps, set<int> & newbps) {
	// Assume that there are no gaps (this case would be difficult to handle properly and would require more information in the alignment)
	int maxBlockSize = 0;
	for (int i = 0; i < blocks.size(); i++) maxBlockSize = max(maxBlockSize, (int) blocks[i].intv.size());

	for (int main_idx = 0; main_idx < maxBlockSize; main_idx++) {
		int numBlocks = blocks.size();
		for (int i = 0; i < numBlocks; i++) {
			split_block(i, bps, newbps, blocks, main_idx);
		}
	}

}


void split_blocks(vector<Block> & blocks, set<int> &bps) {
	set<int> newbps;
	split_blocks(blocks, bps, newbps);
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

bool comp_lt(const IntervalFromBlock & ifb, const int & pt ) {
	return (ifb.intv.end < pt);
}

void find_pos(vector<IntervalFromBlock> & intervals, int pos, int & interval, int & end) {
	vector<IntervalFromBlock>::iterator it = lower_bound(intervals.begin(), intervals.end(), pos, comp_lt);
	interval = it - intervals.begin();
	assert(it != intervals.end());
	if (it->intv.start == pos) {
		end = 0;
	} else if (it->intv.end == pos) {
		end = 1;
	} else {
		end = -1;
	}
}


void find_pos(vector<Block> & blocks, vector<IntervalFromBlock> & intervals, int pos, int & blockIdx, int & end) {
	int interval;
	find_pos(intervals, pos, interval, end);
	blockIdx = intervals[interval].blockIdx;
	int intIdx   = intervals[interval].intIdx;
	bool inv = blocks.at(blockIdx).inv.at(intIdx);
	if (inv && end != -1) end = abs(1 - end);
}
/*
void selectBlocks(vector<Block> & blocks, vector<int> & bps, vector<int> & blockIdxs) {
	UnionFindClass uf(blocks.size());
	vector<IntervalFromBlock> intervals;
	getIntervalsFromBlocks(blocks, intervals);

	for (int i = 0; i < blocks.size(); i++) {
		for (int intIdx = 0; intIdx < blocks[i].intv.size(); intIdx++) {
			int bp_start, bp_end;
			Interval range = blocks[i].intv[intIdx];
			searchRange(bps, range, bp_start, bp_end);
			if (bp_start != bp_end) {
				blockIdx.push_back(i);
			}
		}
	}
}



void split_blocks2(vector<Block> & blocks, vector<int> &bps, vector<int> &newbps) {
	// Assume that there are no gaps (this case would be difficult to handle properly and would require more information in t he alignment)
	sort(bps.begin(), bps.end());
	int maxBlockSize = 0;
	for (int i = 0; i < blocks.size(); i++) maxBlockSize = max(maxBlockSize, (int) blocks[i].intv.size());

	for (int main_idx = 0; main_idx < maxBlockSize; main_idx++) {
		int numBlocks = blocks.size();
		for (int i = 0; i < numBlocks; i++) {
			split_block(i, bps, newbps, blocks, main_idx);
		}
	}

}

*/

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


