#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <functional>
#include "include/graphtypes.h"
#include "include/defs.h"
#include "include/interval.h"
#include "include/union.h"
#include "include/block.h"

/* Glue class and operators */
class Glue {
	public:
		Interval intv[2];
		bool inv;
};

ostream & operator << (ostream & out, const Glue & g) {
	out << g.intv[0] << "\t" << g.intv[1] << "\t" << g.inv;
	return out;
}


bool comp_le1(const Glue & g1, const Glue &g2) {
	return comp_le(g1.intv[0], g2.intv[0]);
}

bool comp_le2(const Glue & g1, const Glue &g2) {
	return comp_le(g1.intv[1], g2.intv[1]);
}


ostream & operator << (ostream & out, const chr_pos_t & x) {
	out << x.contigname << "\t" << x.contigstart << "\t" << x.contigend << "\t" << x.strand;
	return out;
}

istream & operator >> (istream &in, chr_pos_t x) {
	in >> x.contigname >> x.contigstart >> x.contigend >> x.strand;
	return in;
}

ostream & operator << (ostream & out, const pos_entry_t& x) {
	out << x.pos << "\t" << x.edgeindex;
	return out;
}

istream & operator << (istream & in, pos_entry_t x) {
	in >> x.pos >> x.edgeindex;
	return in;
}

/* End Glue class definitions */




bool operator== (const Interval & i1, const Interval & i2) {
	return (i1.chr == i2.chr && i1.start == i2.start && i1.end == i2.end && i1.label == i2.label);
}

/* Takes, as input, a glue file. Produces, as output,
   repeat graph entries (to be sorted later) */
int main(int argc, char ** argv) {

	/* Validate input, but not too much. */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <chrom> <length> <glue_file> <block_file> <edge_file>\n", argv[0]);
		return -1;
	}

	string chr = argv[1];
	int    chrLen = atoi(argv[2]);


	/* Open input and output files */
	ifstream inf;
	open_file(inf, argv[3]);
	ofstream outf;
	open_file(outf, argv[4]);
	FILE * output = fopen(argv[5], "wb");
	if (output == NULL) {
		perror("opening input/output");
		return -1;
	}

	vector<Block> blocks;

	// First, convert glues to edge file for downstream graph building
	uint64_t cur_edge = 0;
	Interval intv[2];
	bool inv;
	while (inf >> intv[0].chr >> intv[0].start >> intv[0].end >> intv[1].chr >> intv[1].start >> intv[1].end >> inv) {
		Block b;
		b.intv.push_back(intv[0]);
		b.intv.push_back(intv[1]);
		b.inv.push_back(false);
		b.inv.push_back(inv);
		blocks.push_back(b);
		
		char strand;
		if (inv) strand = '-'; else strand = '+';

		assert(b.intv[0].chr == b.intv[1].chr); 

		/* Write nodes */
		pos_entry_t entry = { 0 };

		/* First node */
		strcpy(entry.pos.contigname, b.intv[0].chr.c_str());
		entry.pos.contigstart = b.intv[0].start;
		entry.pos.contigend = b.intv[0].end;
		entry.pos.strand = '+';
		entry.edgeindex = cur_edge;
		fwrite(&(entry), sizeof(pos_entry_t), 1, output);

		/* Second node */
		entry.pos.contigstart = b.intv[1].start;
		entry.pos.contigend = b.intv[1].end;
		entry.pos.strand = strand;
		entry.edgeindex = cur_edge;
		fwrite(&(entry), sizeof(pos_entry_t), 1, output);
		cur_edge++;
	}

	/* Close files */
	inf.close();
	fclose(output);

	// Second, split glues into minimal chunks.
	// Assume glues are given in terms of closed intervals.
	// Assume that there are no gaps (this case would be difficult to handle properly and would require more information in the alignment)
	vector<int> bps;
	getBlockBps(blocks, bps);
	split_blocks(blocks, bps);

	//perform transitive closure of blocks 
	//at this point, any two intervals are either disjoint or identical
	UnionFindClass uf(2 * blocks.size());
	vector<pair<Interval, int> > intervals;
	int curidx = 0;
	for (int i = 0; i < blocks.size(); i++) {
		intervals.push_back(make_pair(blocks[i].intv[0], curidx++));
		intervals.push_back(make_pair(blocks[i].intv[1], curidx++));
		uf.unionn(curidx - 1, curidx - 2, blocks[i].inv[1]); //link glued intervals
	}

	//sort intervals and then join identical intervals
	sort(intervals.begin(), intervals.end());
	for (int i = 1; i < intervals.size(); i++) {
		if (intervals[i].first == intervals[i-1].first) {
			uf.unionn(intervals[i].second, intervals[i-1].second, false); //link identical intervals
		}
	}

	//create singleton blocks to fill the gaps along the chromosome where no glues present
	vector<Interval> gaps;
	int lastEnd = 0;
	Interval gap;
	gap.chr = chr;
	for (int i = 0; i < intervals.size(); i++) {
		Interval intv = intervals[i].first;
		gap.start = lastEnd + 1;
		gap.end   = intv.start - 1;
		if (gap.end >= gap.start) {
			gaps.push_back(gap);
		}
		lastEnd = intv.end;
	}
	if (lastEnd < chrLen) { //add last gap
		gap.start = lastEnd + 1;
		gap.end   = chrLen;
		gaps.push_back(gap);
	}

	//output
	// the format is "interval invert? block_size block_index "
	vector<Interval> idx2interval(intervals.size()); 
	for (int i = 0; i < intervals.size(); i++) idx2interval.at(intervals[i].second) = intervals[i].first;

	vector<vector<pair<int, bool> > > classes;
	uf.get_classes(classes);
	int curBlock = 0;
	for (int i = 0; i < classes.size(); i++) {
		vector<Interval> outlines;
		for (int j = 0; j < classes[i].size(); j++) {
			int intv_idx  = classes[i][j].first;
			Interval intv = idx2interval.at(intv_idx);
			bool     inv  = classes[i][j].second;
			intv.label = (inv ? "1" : "0");
			outlines.push_back(intv);
		}
		sort(outlines.begin(), outlines.end());
		outlines.erase(unique(outlines.begin(), outlines.end()), outlines.end());
		for (int i = 0; i < outlines.size(); i++) outf << outlines[i] << "\t" << outlines.size() << "\t" << curBlock << endl;
		curBlock++;
	}

	for (int i = 0; i < gaps.size(); i++) {
		outf << gaps[i] << "\t0\t1\t" << curBlock++ << "\t" << endl;
	}

	outf.close();

	return 0;
}
