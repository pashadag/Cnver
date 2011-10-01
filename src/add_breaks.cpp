#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/defs.h"
#include "include/interval.h"
#include "include/block.h"

int main(int argc, char ** argv) {

	if (argc != 3) {
		fprintf(stderr, "usage: %s <block_file> <bps_file>\n", argv[0]);
		return -1;
	}

	//read in blocks
	ifstream inf;
	vector<Block> blocks;
	open_file(inf, argv[1]);
	read_blocks(inf, blocks);
	inf.close();

	//read in bps 
	set<int> bps;
	int val;
	if (strcmp(argv[2],"-") == 0) {
		while (cin >> val) bps.insert(val);
	} else {
		open_file(inf, argv[2]);
		while (inf >> val) bps.insert(val);
		inf.close();
	}

	//break apart the blocks at the bps
	split_blocks(blocks, bps);

	//output
	for (int i = 0; i < blocks.size(); i++) {
		print_block(cout, blocks[i], i);
	}


	return 0;
}
