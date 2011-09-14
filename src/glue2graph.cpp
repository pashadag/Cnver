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

//range is a closed interval
//points should be sorted
void searchRange(vector<int> & points, Interval range, int & start, int & end) {
	vector<int>::iterator it = lower_bound(points.begin(), points.end(), range.start);
	if (it == points.end() || *it > range.end) {
		start = 0;
		end = 0;
		return;
	}
	start = it - points.begin();
	end = start + 1;
	it++;
	while (it != points.end() && *it <= range.end) {
		it++;
		end++;
	}
}

void check_glue(Glue g) {
	for (int i = 0; i < 2; i++) {
		if (g.intv[i].start > g.intv[i].end) {
			cerr << "Invalid glue: " << g << endl;
		} 
	}
	return;
}

void split_glue(Glue g, vector<int> & bps, vector<Glue> & result, int main_idx = 0) {
	Glue origGlue = g;
	check_glue(g);
	int sec_idx = abs(main_idx -1 );
	
	int bp_start, bp_end;
	Interval range = g.intv[main_idx];
	searchRange(bps, range, bp_start, bp_end);
	for (int i = bp_start; i < bp_end; i++) {
		int bp = bps[i];
		//assume alignments have no gaps
		//break off the beginning
		if (bp == g.intv[main_idx].start) continue;
		Glue newglue = g;
		g.intv[sec_idx].start += bp - g.intv[main_idx].start;
		g.intv[main_idx].start = bp;
		newglue.intv[main_idx].end = g.intv[main_idx].start - 1;
		newglue.intv[sec_idx].end = g.intv[sec_idx].start - 1;
		check_glue(newglue);
		result.push_back(newglue);
	}
	check_glue(g);
	result.push_back(g);
}



/* Takes, as input, a glue file. Produces, as output,
   repeat graph entries (to be sorted later) */
int main(int argc, char ** argv) {

	/* Validate input, but not too much. */
	if (argc != 4) {
		fprintf(stderr, "usage: %s <glue_file> <block_file> <edge_file>\n", argv[0]);
		return -1;
	}

	/* Open input and output files */
	ifstream inf;
	open_file(inf, argv[1]);
	ofstream outf;
	open_file(outf, argv[2]);
	FILE * output = fopen(argv[3], "wb");
	if (output == NULL) {
		perror("opening input/output");
		return -1;
	}

	vector<Glue> glues;
	vector<int> bps;

	// First, convert glues to edge file for downstream graph building
	vector<string> row;
	uint64_t cur_edge = 0;
	Glue g;
	while(inf >> g.intv[0].chr >> g.intv[0].start >> g.intv[0].end >> g.intv[1].chr >> g.intv[1].start >> g.intv[1].end >> g.inv) {
		glues.push_back(g);
		bps.push_back(g.intv[0].start); bps.push_back(g.intv[1].start);
		bps.push_back(g.intv[0].end + 1); bps.push_back(g.intv[1].end + 1);

		char strand;
		if (g.inv) strand = '-'; else strand = '+';

		assert(g.intv[0].chr == g.intv[1].chr); 

		/* Write nodes */
		pos_entry_t entry = { 0 };

		/* First node */
		strcpy(entry.pos.contigname, g.intv[0].chr.c_str());
		entry.pos.contigstart = g.intv[0].start;
		entry.pos.contigend = g.intv[0].end;
		entry.pos.strand = '+';
		entry.edgeindex = cur_edge;
		fwrite(&(entry), sizeof(pos_entry_t), 1, output);

		/* Second node */
		entry.pos.contigstart = g.intv[1].start;
		entry.pos.contigend = g.intv[1].end;
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
	vector<Glue> glues2;
	sort(bps.begin(), bps.end());
	for (int i = 0; i < glues.size(); i++) split_glue(glues[i], bps, glues2, 0);
	glues.clear();
	for (int i = 0; i < glues2.size(); i++) {
		check_glue(glues2[i]);
	}
	for (int i = 0; i < glues2.size(); i++) split_glue(glues2[i], bps, glues, 1);
	//for (int i = 0; i < glues.size(); i++) cout << glues[i] << endl;

	//perform transitive closure of glues
	//at this point, any two intervals are either disjoint or identical
	UnionFindClass uf(2 * glues.size());
	vector<pair<Interval, int> > intervals;
	int curidx = 0;
	for (int i = 0; i < glues.size(); i++) {
		intervals.push_back(make_pair(glues[i].intv[0], curidx++));
		intervals.push_back(make_pair(glues[i].intv[1], curidx++));
		//link glued intervals
		uf.unionn(curidx - 1, curidx - 2, glues[i].inv);
	}

	//sort intervals so that we can join identical intervals
	sort(intervals.begin(), intervals.end());
	for (int i = 1; i < intervals.size(); i++) {
		if (intervals[i].first == intervals[i-1].first) {
			//link identical intervals
			uf.unionn(intervals[i].second, intervals[i-1].second, false);
		}
	}

	//output
	vector<Interval> idx2interval(intervals.size());
	for (int i = 0; i < intervals.size(); i++) idx2interval.at(intervals[i].second) = intervals[i].first;

	vector<vector<pair<int, bool> > > classes;
	uf.get_classes(classes);
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
		for (int i = 0; i < outlines.size(); i++) outf << outlines[i] << endl;
		outf << endl;
	}
	outf.close();

	return 0;
}
