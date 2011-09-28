#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "../include/defs.h"
#include "../include/interval.h"

string make_key(string x, string y) {
	ostringstream o;
	o << x << "," << y;
	return o.str();
}

string sbuf;
vector<string> row;
int main(int argc, char ** argv) {

	if (argc != 3) {
		fprintf(stderr, "usage: %s <placed_link_file> <solved_graph_file>\n", argv[0]);
		return -1;
	}

	//read in links
	vector<Link> links;
	ifstream inf;
	open_file(inf, argv[1]);
	read_links(inf, links);
	inf.close();
	map<string,int>  linkPlacement;
	for (int i = 0; i < links.size(); i++) {
		tokenize(links[i].label, row);
		string key = row[6] + "\t" + row[7] + "\t" + row[8] + "\t" + row[9];
		linkPlacement.insert(make_pair(key, i));
	}

	//first pass through solved_graph_file just to find the number of blocks


	//process flow results
	open_file(inf, argv[2]);
	int numBlocks = 0;
	map<pair<string,string>,string> arcFlow;
	while (get_row_whitespace(inf, row, sbuf)) {
		if (row[0] == "n") {
			numBlocks = max (numBlocks, atoi(row[1].c_str()) / 2);
		} else if (row[0] == "a") {
			//this does nothing if entry older exists, and that's the point
			arcFlow.insert(make_pair(make_pair(row[1], row[2]), row[6])); 
		}
	}

	vector<double> blockFlow(numBlocks, 0);
	for (map<pair<string,string>,string>::iterator arcit = arcFlow.begin(); arcit != arcFlow.end(); arcit++) {
		int from    = atoi(arcit->first.first.c_str());
		int to      = atoi(arcit->first.second.c_str());
		double flow = atof(arcit->second.c_str());
		int fromBlock = (abs(from) - 1) / 2;
		int toBlock   = (abs(to)   - 1) / 2;
		int fromEnd   = (abs(from) + 1) % 2 ;
		int toEnd     = (abs(to)   + 1) % 2;
		if (fromBlock == toBlock && fromEnd != toEnd) { //sequence edge, something like 17-18
			blockFlow.at(fromBlock) += (from < to) ? flow : (-1 * flow);
		} else {
			string key = make_string(fromBlock) + "\t" + make_string(toBlock) + "\t";
			key += make_string(fromEnd) + "\t" + make_string(toEnd);
			map<string,int>::iterator it = linkPlacement.lower_bound(key);
			if (it != linkPlacement.end() && it->first == key) { //this is a link 
				if (flow != 0) cout << "don_adj\t" << flow << "\t" << links[it->second] << endl;
			} else { //this is a ref adj
				cout << "ref_adj\t" << flow << "\t" << key << endl;
			}
		}
	}
	inf.close();

	for (int i = 0; i < numBlocks; i++) {
		cout << "block\t" << blockFlow.at(i) << "\t" << i << endl;
	}
	return 0;
}
