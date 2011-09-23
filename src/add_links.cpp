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

void get_bps(Link & link, int & bpfrom, int & bpto) {
	if (link.type == 0 || link.type == 2) { //+- or ++
		bpfrom = link.from + 1;
	} else {
		bpfrom = link.from;
	}
	if (link.type == 1 || link.type == 2) { //-+ or ++
		bpto = link.to + 1;
	} else {
		bpto = link.to;
	}
}


int main(int argc, char ** argv) {

	if (argc != 3) {
		fprintf(stderr, "usage: %s <block_file> <link_file>\n", argv[0]);
		exit(1);
	}

	//read in blocks
	ifstream inf;
	vector<Block> blocks;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	inf.close();

	//read in links
	vector<Link> links;
	open_file(inf, argv[2]);
	read_links(inf, links);
	inf.close();

	//create bps from links
	vector<int> bps;
	for (int i = 0; i < links.size(); i++) {
		int bpfrom, bpto;
		get_bps(links[i], bpfrom, bpto);
		bps.push_back(bpfrom);
		bps.push_back(bpto);
	}

	//split blocks up at the bps
	split_blocks(blocks, bps);

	//output blocks
	for (int i = 0; i < blocks.size(); i++) print_block(cout, blocks[i], i);

	return 0;

}
