#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/graphtypes.h"
#include "include/defs.h"

/* Takes, as input, a glue file. Produces, as output,
    repeat graph entries (to be sorted later) */
int main(int argc, char ** argv) {

    /* Validate input, but not too much. */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <glue_file> <edge_file>\n", argv[0]);
		return -1;
	}

	/* Open input and output files */
	ifstream inf;
	open_file(inf, argv[1]);
	FILE * output = fopen(argv[2], "wb");
	if (output == NULL) {
        perror("opening input/output");
        return -1;
    }

	vector<string> row;
	uint64_t cur_edge = 0;
	while(get_row(inf, row)) {
		string contigname = row[0];
		uint64_t start1   = atol(row[1].c_str());
		uint64_t end1     = atol(row[2].c_str());
		string contigname2= row[3];
		uint64_t start2   = atol(row[4].c_str());
		uint64_t end2     = atol(row[5].c_str());
		string   strand   = row[6];

		assert (strand == "+" || strand == "-");
		assert(contigname == contigname2);

		/* Write nodes */
		pos_entry_t entry = { 0 };

		/* First node */
		strcpy(entry.pos.contigname, contigname.c_str());
		entry.pos.contigstart = start1;
		entry.pos.contigend = end1;
		entry.pos.strand = '+';
		entry.edgeindex = cur_edge;
		fwrite(&(entry), sizeof(pos_entry_t), 1, output);

		/* Second node */
		entry.pos.contigstart = start2;
		entry.pos.contigend = end2;
		entry.pos.strand = strand[0];
		entry.edgeindex = cur_edge;
		fwrite(&(entry), sizeof(pos_entry_t), 1, output);
		cur_edge++;
	}

    /* Close files */
    inf.close();
    fclose(output);
    
    return 0;
}
