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

map<Interval, Interval> blockMap;
vector<Link> links;

uint64_t map_pos(string chr, uint64_t pos) {
	map<Interval,Interval>::iterator it;
	for (it = blockMap.begin(); it != blockMap.end(); it++) {
		Interval intv   = it->first;
		Interval master = it->second;
		if (intv.contains(chr, pos)) {
			return pos - intv.start + master.start;
		} 
	}
	return pos;
}


int main(int argc, char ** argv) {

	/* Validate input, but not too much. */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <block_file> <link_file>\n", argv[0]);
		return -1;
	}


	// Read blocks into blockMap
	ifstream inf;
	open_file(inf, argv[1]);
	bool newblock = true;
	Interval master;
	while (getline(inf, sbuf)) {
		if (sbuf == "") {
			newblock = true;
			continue;
		} 
		Interval intv = read_interval(sbuf);
		if (newblock) {
			master = intv;
			newblock = false;
		} else {
			blockMap.insert(make_pair(intv, master));
		}
	}
	inf.close();

	//Read links
	open_file(inf, argv[2]);
	read_links(inf, links);

	//transform and output links
	for (int i = 0; i < links.size(); i++) {
		Link link = links[i];
		link.from = map_pos(link.chr, link.from);
		link.to   = map_pos(link.chr, link.to);
		cout << link << endl;
	}

	return 0;
}
