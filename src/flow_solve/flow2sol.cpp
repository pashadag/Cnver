#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "../include/defs.h"
#include "../include/interval.h"
#include "../include/block.h"

string make_key(string x, string y) {
	ostringstream o;
	o << x << "," << y;
	return o.str();
}

void updateMap(map<string, double> & store, string key, double val) {
	map<string, double>::iterator it = store.lower_bound(key);
	if (it != store.end() && it->first == key) {
		it->second += val;
	} else {
		store.insert(make_pair(key, val));
	}
}

string sbuf;
vector<string> row;
int main(int argc, char ** argv) {

	if (argc != 4) {
		fprintf(stderr, "usage: %s <block_file> <placed_link_file> <solved_graph_file>\n", argv[0]);
		return -1;
	}

	//read in blocks
	vector<Block> blocks;
	vector<IntervalFromBlock> intervals;
	ifstream inf;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	getIntervalsFromBlocks(blocks, intervals);
	inf.close();

	//read in links
	vector<Link> links;
	open_file(inf, argv[2]);
	read_links(inf, links);
	inf.close();
	map<string,int>  linkPlacement;
	for (int i = 0; i < links.size(); i++) {
		tokenize(links[i].label, row);
		string key = row[6] + "\t" + row[7] + "\t" + row[8] + "\t" + row[9];
		linkPlacement.insert(make_pair(key, i));
	}

	//find 2-cycles involving a sequence edge.  Because of the way the flow is reported,
	//this needs to be handled specially.
	vector<bool> refCycle(blocks.size(), false);
	for (int i = 1; i < intervals.size(); i++) {
		if (intervals[i-1].blockIdx == intervals[i].blockIdx && 
				blocks[intervals[i-1].blockIdx].inv[intervals[i-1].intIdx] ==
				blocks[intervals[i  ].blockIdx].inv[intervals[i  ].intIdx]) {
			refCycle[intervals[i].blockIdx] = true;
		}
	}

	//cout << "refcycles\n"; for (int i = 0; i < refCycle.size(); i++) if (refCycle[i]) cout << i << endl;

	//process flow results
	open_file(inf, argv[3]);
	int numBlocks = 0;
	map<pair<string,string>,string> arcFlow;
	while (get_row_whitespace(inf, row, sbuf)) {
		if (row[0] == "n") {
			numBlocks = max (numBlocks, atoi(row[1].c_str()) / 2);
		} else if (row[0] == "a") {
			//this does nothing if entry older exists, and that's the point
			//because the flow along parallel edges is reported as the total sum 
			arcFlow.insert(make_pair(make_pair(row[1], row[2]), row[6])); 
		}
	}

	vector<double> blockFlow(numBlocks, 0);
	map<string,double> refFlow;

	for (map<pair<string,string>,string>::iterator arcit = arcFlow.begin(); arcit != arcFlow.end(); arcit++) {
		int from    = atoi(arcit->first.first.c_str());
		int to      = atoi(arcit->first.second.c_str());
		double flow = atof(arcit->second.c_str());
		int fromBlock = (abs(from) - 1) / 2;
		int toBlock   = (abs(to)   - 1) / 2;
		int fromEnd   = (abs(from) + 1) % 2 ;
		int toEnd     = (abs(to)   + 1) % 2;

		//check if this is a link
		string key = make_string(fromBlock) + "\t" + make_string(toBlock) + "\t";
		key += make_string(fromEnd) + "\t" + make_string(toEnd);
		map<string,int>::iterator it = linkPlacement.lower_bound(key);

		if (fromBlock == toBlock && fromEnd != toEnd) { //sequence edge, something like 17-18, or a loop
			if (refCycle[fromBlock] && fromEnd == 1) {
				updateMap(refFlow, key, flow);
			} else if (it != linkPlacement.end() && it->first == key) { //this is a link 
				if (flow != 0) cout << "link\t" << flow << "\t" << links[it->second] << endl;
			} else {
				blockFlow.at(fromBlock) += (from < to) ? flow : (-1 * flow);
			}
		} else {
			if (it != linkPlacement.end() && it->first == key) { //this is a link 
				if (flow != 0) cout << "link\t" << flow << "\t" << links[it->second] << endl;
			} else { //this is a ref adj
				updateMap(refFlow, key, flow);
			}
		}
	}
	inf.close();

	for (map<string,double>::iterator it = refFlow.begin(); it != refFlow.end(); it++) {
		cout << "ref_adj\t" << it->second << "\t" << it->first << endl;
	}

	for (int i = 0; i < numBlocks; i++) {
		cout << "block\t" << blockFlow.at(i) << "\t" << i << endl;
	}
	return 0;
}
