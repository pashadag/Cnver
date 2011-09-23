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

string sbuf;

vector<Link> links;

bool comp_lt(const Link &link1, const Link &link2) {
	if (link1.chr == link2.chr) {
		if (link1.from == link2.from) {
			return link1.to < link2.to;
		} 
		return link1.from < link2.from;
	} 
	return link1.chr < link2.chr;
}

int main(int argc, char ** argv) {

	/* Validate input, but not too much. */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <link_file> <max_dist>\n", argv[0]);
		return -1;
	}

	int max_dist = atoi(argv[2]);

	ifstream inf;
	open_file(inf, argv[1]);
	read_links(inf, links);
	inf.close();

	sort(links.begin(), links.end(), comp_lt);
	vector<bool> removed(links.size(), false);

	//merge links
	for (int i = 0; i < links.size(); i++) {
		if (removed[i]) continue;
		bool clusterSize = 1;
		int j = i+1;
		vector<Link> cluster;
		cluster.push_back(links[i]); //for debug
		while (j < links.size() && links[i].chr == links[j].chr && links[j].from < links[i].from + max_dist) {
			if (!removed[j] && links[i].type == links[j].type && abs(links[j].to - links[i].to) < max_dist) { //merge links
				removed[j] = true;
				links[i].from = (links[i].from * clusterSize + links[j].from) / (clusterSize + 1);
				links[i].to   = (links[i].to   * clusterSize + links[j].to  ) / (clusterSize + 1);
				clusterSize++;
				cluster.push_back(links[j]);//for debug
			}
			j++;
		}
		/*if (cluster.size() > 1) {
			for (int j = 0; j < cluster.size(); j++) {
				cout << i << "\t" << cluster.size() << "\t" << cluster[j] << endl;
			}
			cout << endl;
		}
		*/
	}

	for (int i = 0; i < links.size(); i++) {
		if (!removed[i]) {
			cout << links[i] << endl;
		}
	}

	return 0;
}
