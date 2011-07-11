#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <assert.h>

#include "include/graphtypes.h"
#include "include/dbtypes.h"


#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#define POSSTRAND(p) ((p)->strand == '+')

#define MAX_LINE_LEN        1024

using namespace std;

#define INVALID_ORIENTATION     -1
#define DEFAULT_ORIENTATION     0

#define UNKNOWN_TYPE     -1
#define FILLGRAPH_TYPE   0
#define DONOR_TYPE   1
#define FAKE_TYPE   2
#define NUM_BINS 25

class pos {
  public:
    pos(string contig, uint64_t contigpos, uint8_t strand = '+',
	bool printpos = true)
  :	contig(contig), contigpos(contigpos), strand(strand),
	printpos(printpos) {
    } bool operator<(const pos &) const;
    bool operator==(const pos &) const;
    string contig;
    uint64_t contigpos;
    uint8_t strand;
    bool printpos;
};

bool pos::operator<(const pos & other) const {
    //if (contig!=other.contig) {
	//printf("%s vs %s\n",contig.c_str(),other.contig.c_str());
    //}
    /* If contig is the same, */
    if (contig == other.contig) {
	/* Compare positions */
	return (contigpos < other.contigpos);
    } else {
	printf("ERROR: %s not %s\n",contig.c_str(),other.contig.c_str());
	exit(1);
    }

    assert(contig == other.contig);
    /* Else return contig difference */
    return (contig < other.contig);
}

bool pos::operator==(const pos & other) const {
    return !((*this < other) || (other < *this));
}

class edge_marker {
  public:
    edge_marker(uint64_t start, uint64_t end):start(start), end(end) {
    } uint64_t start;
    uint64_t end;
    bool operator<(const edge_marker & other) const {
	if (start < other.start) {
	    return true;
	}
	return end < other.end;
    }
    bool operator==(const edge_marker & other) const {
	return !((*this < other) || (other < *this));
}};

class node {
  public:
    node() {
    } set < pos > positions;
};

class edge {
  public:
    edge(int8_t orientation = DEFAULT_ORIENTATION)
  :	length(0), unmasked_length(0), num_paths(0),
	hits(0.0), is_free(false), is_cb(false) {
	switch (orientation) {
	case INVALID_ORIENTATION:
	    from_bit = to_bit = INVALID_ORIENTATION;
	    break;
	    case 0:from_bit = 0;
	    to_bit = 1;
	    break;
	    case 1:from_bit = 1;
	    to_bit = 0;
	    break;
	    case 2:from_bit = 0;
	    to_bit = 0;
	    break;
	    case 3:from_bit = 1;
	    to_bit = 1;
	    break;
	}
	expected= 0.0;
    };
    uint64_t start;
    uint64_t end;
    uint64_t length;
    uint64_t unmasked_length;
    uint64_t num_paths;
    double hits;
    bool is_free;
    bool is_cb;
    int8_t from_bit;
    int8_t to_bit;
    int8_t type;
    double expected;
};

class dg_info {
  public:
    dg_info() {
    };
    dg_info(string contig, uint64_t start, uint64_t end,
	    int8_t orientation)
  :	contig(contig), start(start), end(end),
	orientation(orientation) {
    };
    string contig;
    uint64_t start;
    uint64_t end;
    int8_t orientation;
};

typedef enum {
    SRC_NODE,
    DEST_NODE
} node_type_t;

/* Entry index comparison function */
int edgecmp(const void *v1, const void *v2)
{
    pos_entry_t *p1 = *((pos_entry_t **) v1);
    pos_entry_t *p2 = *((pos_entry_t **) v2);

    /* If the edges are not the same, */
    if (p1->edgeindex < p2->edgeindex)
	return -1;
    else if (p1->edgeindex > p2->edgeindex)
	return 1;

    /* First compare contig name, then contig position */
    int chrresult = strncmp(p1->pos.contigname,
			    p2->pos.contigname,
			    CONTIGNAME_LEN);
    if (chrresult == 0) {
	int64_t posdiff = p1->pos.contigstart - p2->pos.contigstart;
	if (posdiff > 0)
	    return 1;
	if (posdiff < 0)
	    return -1;
	return 0;
    }

    return chrresult;
}

uint64_t entry_bound(pos_entry_t * entry, node_type_t type)
{
    if (type == SRC_NODE) {
	if (entry->pos.strand == '+') {
	    return entry->pos.contigstart;
	} else {
	    return entry->pos.contigend;
	}
    } else {
	if (entry->pos.strand == '+') {
	    return entry->pos.contigend;
	} else {
	    return entry->pos.contigstart;
	}
    }
}

/*
// Hit index comparison function 
int indexcmp(const void *p1, const void *p2)
{
    // First compare contig name, then contig position 
    int chrresult = strncmp(((mapping_t_small *) p1)->contigname,
			    ((mapping_t_small *) p2)->contigname,
			    CONTIGNAME_LEN);
    if (chrresult == 0) {
	int64_t posdiff = ((mapping_t_small *) p1)->contigstart
	    - ((mapping_t_small *) p2)->contigstart;
	if (posdiff > 0)
	    return 1;
	if (posdiff < 0)
	    return -1;
	return 0;
    }

    return chrresult;
}

// Find first hit in range 
//
mapping_t_small *range_search(mapping_t_small * index, uint64_t indexsize,
			const char *contigname, uint64_t start,
			uint64_t end)
{
    mapping_t_small startkey = { 0 };
    mapping_t_small endkey = { 0 };

    // Initialize keys 
    strcpy(startkey.contigname, contigname);
    endkey = startkey;
    startkey.contigstart = start;
    endkey.contigstart = end;

    // Start binary search
    uint64_t low = 0, high = indexsize - 1, mid = high / 2;
    int res;
    while ((low <= high) && (high > 0) && mid > 0) {
	int res = indexcmp(&startkey, index + mid);
	if (res < 0) {
	    high = mid - 1;
	    mid = (low + high) / 2;
	} else if (res > 0) {
	    low = mid + 1;
	    mid = (low + high) / 2;
	} else
	    break;
    }

    // See where we're at 
    res = indexcmp(index + mid, &startkey);

    // If there's a hit for the start of the range, 
    if (res == 0) {
	// Back up to the beginning of the range 
	while (indexcmp(&startkey, index + mid) == 0) {
	    if (mid == 0)
		return index + mid;
	    mid--;
	}

	return index + mid + 1;
    }
    // Else if we're below the range start, 
    else if (res < 0) {
	// Advance until we're past the start of the range 
	while (indexcmp(index + mid, &startkey) < 0) {
	    if (++mid == indexsize)
		return NULL;
	}
    }
    // Else we're above the range start 

    // Make sure we're not past the end of the range 
    if (indexcmp(index + mid, &endkey) <= 0)
	return index + mid;
    else
	return NULL;
}
*/

double get_scov_arrivals(double * scov, uint64_t contig_length, uint64_t start, uint64_t end,double * lambdas) {
	double sum=0.0;
	//assert((startpos.contigpos + this_edgelen)<contig_length);
	if (end<=contig_length) {
		for (unsigned int j=start; j<=end; j++) {
			if (lambdas[j-1]>0.0) {
				//noting
			} else {
				printf("%f %u\n",lambdas[j-1],j-1);
			}
			assert(lambdas[j-1]>0.0);
			//only count non-repeat regions
			sum += scov[j-1];
		}      
	} else {
	        printf("##ERROR reference position %d, but contig_len is %d\n",end,contig_length);
	}
	return sum;
}


//MISKO START
void write_out(map < uint64_t, edge > &edges, map < uint64_t,
	       node > &nodes, char *edge_filename)
{
    FILE *f_edges = fopen(edge_filename, "w");
    if (f_edges == NULL) {
	printf("Opening of file, %s, failed.\n", edge_filename);
	exit(1);
    }

    int l = strlen(edge_filename);
    char edge_atrib_filename[(l + 1) + 3];
    strcpy(edge_atrib_filename, edge_filename);
    strcat(edge_atrib_filename, ".ea");
    FILE *f_ea = fopen(edge_atrib_filename, "w");
    if (f_edges == NULL) {
	printf("Opening of file, %s, failed.\n", edge_atrib_filename);
	exit(1);
    }


    char fa_filename[(l + 1) + 4];
    strcpy(fa_filename, edge_filename);
    strcat(fa_filename, ".a");
    FILE *f_fa = fopen(fa_filename, "w");
    if (f_edges == NULL) {
	printf("Opening of file, %s, failed.\n", fa_filename);
	exit(1);
    }
    char sif_filename[(l + 1) + 4];
    strcpy(sif_filename, edge_filename);
    strcat(sif_filename, ".sif");
    FILE *f_sif = fopen(sif_filename, "w");
    if (f_edges == NULL) {
	printf("Opening of file, %s, failed.\n", sif_filename);
	exit(1);
    }
    //PRINT IN GML FORMAT
    /*fprintf(f_edges,"graph [\n");

       //map<uint64_t, node>::iterator nit = nodes.begin();
       for (int i=0; i<nodes.size(); i++) {
       fprintf(f_edges,"\tnode [\n");
       fprintf(f_edges,"\t\tid %d\n",i);
       fprintf(f_edges,"\t\tlabel \"");
       set<pos>::iterator sit = nodes[i].positions.begin();
       if (sit!=nodes[i].positions.end()) {
       fprintf(f_edges,"%llu%c",sit->contigpos,sit->strand);
       }
       sit++;
       for ( ; sit!=nodes[i].positions.end(); sit++) {
       fprintf(f_edges,",%llu%c",sit->contigpos,sit->strand);
       }
       fprintf(f_edges,"\"\n");
       fprintf(f_edges,"\t]\n");
       }


       for (int i=0; i<edges.size(); i++) {
       char oreintation[3]="oo";
       if (edges[i].from_bit) {
       oreintation[0]='i';
       }
       if (edges[i].to_bit) {
       oreintation[1]='i';
       }
       oreintation[2]='\0';
       fprintf(f_edges,"\tedge [\n");
       fprintf(f_edges,"\t\tsource %llu\n",edges[i].start);
       fprintf(f_edges,"\t\ttarget %llu\n",edges[i].end);
       fprintf(f_edges,"\t\tTYPE %s\n",oreintation);
       //fprintf(f_edges,"%llu %s %llu\n",edges[i].start,orientation,edges[i].end);
       fprintf(f_edges,"\t]\n");
       }

       fprintf(f_edges,"]\n"); */


    //PRINT IN TABLE FORMAT
    fprintf(f_edges, "SOURCE\tinteraction\tDESTINATION\tTYPE\n");
    map < uint64_t, edge >::iterator it = edges.begin();
    for (; it != edges.end(); it++) {
	edge e = it->second;
	//Print first node
	fprintf(f_edges, "%d[", e.start);
	set < pos >::iterator sit = nodes[e.start].positions.begin();
	if (sit != nodes[e.start].positions.end()) {
	    fprintf(f_edges, "%llu%c", sit->contigpos, sit->strand);
	}
	sit++;
	for (; sit != nodes[e.start].positions.end(); sit++) {
	    fprintf(f_edges, ",%llu%c", sit->contigpos, sit->strand);
	}
	fprintf(f_edges, "]");
	//print orientation
	char orientation[3] = "oo";
	if (e.from_bit) {
	    orientation[0] = 'i';
	}
	if (e.to_bit) {
	    orientation[1] = 'i';
	}
	orientation[2] = '\0';
	fprintf(f_edges, "\t%s\t", orientation);
	//Print second node
	fprintf(f_edges, "%d[", e.end);
	sit = nodes[e.end].positions.begin();
	if (sit != nodes[e.end].positions.end()) {
	    fprintf(f_edges, "%llu%c", sit->contigpos, sit->strand);
	}
	sit++;
	for (; sit != nodes[e.end].positions.end(); sit++) {
	    fprintf(f_edges, ",%llu%c", sit->contigpos, sit->strand);
	}
	fprintf(f_edges, "]");
	//Print type
	fprintf(f_edges, "\t%d\t%u\t%u\t%llu\n", e.type, e.from_bit,
		e.to_bit, e.length);
    }

    fclose(f_edges);


    //PRINT IN SIF and EDGE ATTRIB FORMAT 

    fprintf(f_ea, "TYPE (class=java.lang.Integer)\n");
    fprintf(f_fa, "TYPE (class=java.lang.Integer) NODES=%llu\n",
	    nodes.size());
    it = edges.begin();
    for (; it != edges.end(); it++) {
	edge e = it->second;
	uint64_t edge_index = it->first;
	//Print first node
	fprintf(f_sif, "%d[", e.start);
	fprintf(f_ea, "%d[", e.start);
	fprintf(f_fa, "%d[", e.start);
	set < pos >::iterator sit = nodes[e.start].positions.begin();
	if (sit != nodes[e.start].positions.end()) {
	    fprintf(f_sif, "%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_ea, "%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_fa, "%llu%c", sit->contigpos, sit->strand);
	}
	sit++;
	for (; sit != nodes[e.start].positions.end(); sit++) {
	    fprintf(f_sif, ",%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_ea, ",%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_fa, ",%llu%c", sit->contigpos, sit->strand);
	}
	fprintf(f_sif, "]");
	fprintf(f_ea, "]");
	fprintf(f_fa, "]");
	//print orientation
	char orientation[3] = "oo";
	if (e.from_bit) {
	    orientation[0] = 'i';
	}
	if (e.to_bit) {
	    orientation[1] = 'i';
	}
	orientation[2] = '\0';
	fprintf(f_sif, " %s ", orientation);
	fprintf(f_ea, " (%s) ", orientation);
	fprintf(f_fa, " (%s) ", orientation);
	//Print second node
	fprintf(f_sif, "%d[", e.end);
	fprintf(f_ea, "%d[", e.end);
	fprintf(f_fa, "%d[", e.end);
	sit = nodes[e.end].positions.begin();
	if (sit != nodes[e.end].positions.end()) {
	    fprintf(f_sif, "%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_ea, "%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_fa, "%llu%c", sit->contigpos, sit->strand);
	}
	sit++;
	for (; sit != nodes[e.end].positions.end(); sit++) {
	    fprintf(f_sif, ",%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_ea, ",%llu%c", sit->contigpos, sit->strand);
	    fprintf(f_fa, ",%llu%c", sit->contigpos, sit->strand);
	}
	fprintf(f_sif, "]\n");
	fprintf(f_ea, "]");
	fprintf(f_fa, "]");
	//Print type
	fprintf(f_ea, " = %d\n", e.type);
	fprintf(f_fa, " = E%llu\n", edge_index);
    }

    fclose(f_sif);
    fclose(f_ea);
    fclose(f_fa);


}

//MISKO END



uint64_t split_edge(uint64_t edgeindex, uint64_t offset,
		    uint64_t * num_edges, map < uint64_t, edge > &edges,
		    uint64_t * num_nodes, map < uint64_t, node > &nodes,
		    map < pos, uint64_t > &pos2edge)
{
    assert(offset > 0);
    /* Create new edge */
    uint64_t newedgeindex = (*num_edges)++;
    edges[newedgeindex] = edge();
    edges[newedgeindex].type = edges[edgeindex].type;
    edges[newedgeindex].start = (*num_nodes)++;
    edges[newedgeindex].end = edges[edgeindex].end;
    edges[edgeindex].end = (*num_nodes)++;

    //printf("Created new nodes %llu and %llu\n",edges[edgeindex].end,edges[newedgeindex].start);
    /* set the new lengths of edges */
    edges[newedgeindex].length = edges[edgeindex].length - offset;
    edges[edgeindex].length = offset;


    /* Split positive and negative strand nodes  */
    set < pos > positions = nodes[edges[edgeindex].start].positions;
    set < pos >::iterator iter = positions.begin();
    while (iter != positions.end()) {
	pos oldpos = *iter;
	if (POSSTRAND(&oldpos)) {
	    pos newpos_oldedge(oldpos.contig,
			       oldpos.contigpos + edges[edgeindex].length -
			       1, oldpos.strand, false);
	    pos newpos_newedge(oldpos.contig,
			       oldpos.contigpos + edges[edgeindex].length,
			       oldpos.strand, true);
	    nodes[edges[edgeindex].end].positions.insert(newpos_oldedge);
	    nodes[edges[newedgeindex].start].positions.
		insert(newpos_newedge);
	    pos2edge[newpos_newedge] = newedgeindex;
	}
	++iter;
    }

    positions = nodes[edges[newedgeindex].end].positions;
    iter = positions.begin();
    while (iter != positions.end()) {
	pos oldpos = *iter;
	if (!POSSTRAND(&oldpos)) {
	    pos newpos_newedge(oldpos.contig,
			       oldpos.contigpos +
			       (edges[newedgeindex].length) - 1,
			       oldpos.strand, false);
	    pos newpos_oldedge(oldpos.contig,
			       oldpos.contigpos +
			       (edges[newedgeindex].length), oldpos.strand,
			       true);
	    nodes[edges[edgeindex].end].positions.insert(newpos_oldedge);
	    nodes[edges[newedgeindex].start].positions.
		insert(newpos_newedge);
	    pos2edge[newpos_oldedge] = edgeindex;
	    pos2edge[oldpos] = newedgeindex;
	}

	++iter;
    }

#ifdef DUMP
    //write_out(edges,nodes,"./splitter.edges");
#endif
    return edges[newedgeindex].start;

}

map < pos, uint64_t >::iterator make_pos(pos & check_pos,
					 uint64_t * num_edges,
					 map < uint64_t, edge > &edges,
					 uint64_t * num_nodes,
					 map < uint64_t, node > &nodes,
					 map < pos, uint64_t > &pos2edge)
{
    map < pos, uint64_t >::iterator edge_it =
	pos2edge.lower_bound(check_pos);
    if (edge_it == pos2edge.end() || !(check_pos == edge_it->first)) {
	--edge_it;
	pos p = edge_it->first;
	uint64_t edgeindex = edge_it->second;
	uint64_t offset = check_pos.contigpos - p.contigpos;
	if (!POSSTRAND(&p)) {
	    offset = edges[edgeindex].length - offset;
	}
	//printf("Splitting %llu %llu\n",edgeindex,offset);
	uint64_t new_node =
	    split_edge(edgeindex, offset, num_edges, edges, num_nodes,
		       nodes, pos2edge);
	//printf("NEW NODE %llu\n",new_node);
	edge_it = pos2edge.lower_bound(check_pos);
    }
    assert(edge_it->first == check_pos);
    return edge_it;
}


double get_expected(double* lambdas,uint64_t start,uint64_t end) {
	double result=0.0;
	while(start<=end) {
		assert(lambdas[start-1]>0.0);
		result+=lambdas[start-1];
		start++;
	}
	return result;
	
}

/* Fills in a sorted repeat graph */
int main(int argc, char **argv)
{

	/* Validate input, but not too much. */
	if (argc != 6) {
		fprintf(stderr,
				"usage: %s <graph_file> <donor_edges> <scov file> <masks> <gc file>\n",
				argv[0]);
		return -1;
	}


	/* Open input */
	int input = open(argv[1], O_RDONLY);
	if (input < 0) {
		perror("opening input");
		return -1;
	}

	/* Get size of input file */
	struct stat inputstats;


	//Set up the lambdas
	char* lambdas_filename=argv[5];
	int res=stat(lambdas_filename,&inputstats);
	assert(res==0);
	off_t filesize=inputstats.st_size;
	unsigned int chr_length=filesize/sizeof(double);
	double* lambdas=(double*)malloc(sizeof(double)*chr_length);
	if (lambdas==NULL) {
		fprintf(stderr,"failed to allocate memory for lambdas\n");
		exit(1);
	}
	FILE *f = fopen(lambdas_filename,"rb");
	res=fread(lambdas,sizeof(double),chr_length,f);
	assert(res==chr_length);
	fclose(f);


	fstat(input, &(inputstats));
	filesize = inputstats.st_size;

	/* Read input */
	uint8_t *filedata = (uint8_t *) malloc(filesize);
	if (filedata == NULL) {
		fprintf(stderr, "Out of memory\n");
		return -1;
	}
	off_t total = 0;
	while (total < filesize) {
		ssize_t count = read(input, filedata + total,
				filesize - total);
		if (count < 0) {
			perror("Reading input\n");
			return -1;
		}
		total += count;
	}

	/* Close input file */
	close(input);

	/* Calculate number of entries */
	uint64_t num_entries = filesize / sizeof(pos_entry_t);
	pos_entry_t *entries = (pos_entry_t *) filedata;

	/* Calculate number of edges */
	uint64_t num_edges = 0;
	for (uint64_t i = 0; i < num_entries; i++) {
		num_edges = MAX(entries[i].edgeindex + 1, num_edges);
	}

	/* Sort entries by edges, then by position */
	pos_entry_t **edge_index =
		(pos_entry_t **) malloc(sizeof(pos_entry_t *) * num_entries);
	for (uint64_t i = 0; i < num_entries; i++) {
		edge_index[i] = &(entries[i]);
	}
	qsort(edge_index, num_entries, sizeof(pos_entry_t *), edgecmp);

	/* Create edges */
	map < uint64_t, edge > edges;
	map < uint64_t, node > nodes;
	uint64_t num_nodes = 0;
	map < pos, uint64_t > pos2edge;
	for (uint64_t i = 0; i < num_entries; i++) {
		pos_entry_t *pe = edge_index[i];

		/* Create positions */
		pos startpos(pe->pos.contigname,
				entry_bound(pe, SRC_NODE),
				pe->pos.strand, POSSTRAND(&(pe->pos)));
		pos endpos(pe->pos.contigname,
				entry_bound(pe, DEST_NODE),
				pe->pos.strand, !POSSTRAND(&(pe->pos)));

		/* Add to edge nodes */
		if (edges.count(pe->edgeindex) == 0) {
			edges[pe->edgeindex] = edge();
			edges[pe->edgeindex].type = FILLGRAPH_TYPE;
			edges[pe->edgeindex].start = num_nodes++;
			edges[pe->edgeindex].end = num_nodes++;
			edges[pe->edgeindex].length =
				pe->pos.contigend - pe->pos.contigstart + 1;
			assert(edges[pe->edgeindex].length != 0);
		}
		nodes[edges[pe->edgeindex].start].positions.insert(startpos);
		nodes[edges[pe->edgeindex].end].positions.insert(endpos);

		/* Hook up overlap->edge mapping */
		if (POSSTRAND(&(pe->pos)))
			pos2edge[startpos] = pe->edgeindex;
		else
			pos2edge[endpos] = pe->edgeindex;
	}
#ifdef DUMP
	write_out(edges, nodes, "./fillgraph.edges");
#endif
	/* Add donor edges */
	map < uint64_t, dg_info > donor_edges;
	map < edge_marker, bool > dg_used;
	char line[MAX_LINE_LEN];
	FILE *dinput = fopen(argv[2], "r");
	while (fgets(line, MAX_LINE_LEN, dinput) != NULL) {
		char contig[CONTIGNAME_LEN];
		uint64_t start, end;
		int8_t orientation;
		if (sscanf(line, "%s %llu %llu %hhd\n",
					contig, &start, &end, &orientation) < 4)
			continue;
		if (start<=0) {
			printf("Offending edge ignored! starts at %llu\n%s\n",start,line);
			//exit(1);
			continue;
		} else if (end>chr_length) {
			printf("Offending edge ignored! ends at %llu\n%s\n",end,line);
			//exit(1);
			continue;
		}

		/* Create positions */
		pos donor_start(contig, start);
		pos donor_end(contig, end);
		if (start == end) {
			assert(donor_start == donor_end);
		}

		/* Create edge */
		edge donor_edge = edge(orientation);
		donor_edge.type = DONOR_TYPE;
		donor_edge.is_free = true;



		//Find a start edge with a start position that is
		//greater then or equal to desired 
		map < pos, uint64_t >::iterator start_edge =
			make_pos(donor_start, &num_edges, edges, &num_nodes, nodes, pos2edge);
		pos p = start_edge->first;
		assert(pos2edge.lower_bound(donor_start)->first == donor_start);
		assert(p == donor_start);
		if (!POSSTRAND(&p)) {
			donor_edge.from_bit ^= 1;
			donor_edge.start = edges[start_edge->second].end;
		} else {
			donor_edge.start = edges[start_edge->second].start;
		}


		/* Hook up donor end */
		map < pos, uint64_t >::iterator end_edge =
			make_pos(donor_end, &num_edges, edges, &num_nodes, nodes,
					pos2edge);
		p = end_edge->first;
		assert(pos2edge.lower_bound(donor_end)->first == donor_end);
		assert(p == donor_end);
		if (!POSSTRAND(&p)) {
			donor_edge.to_bit ^= 1;
			donor_edge.end = edges[end_edge->second].end;
		} else {
			donor_edge.end = edges[end_edge->second].start;
		}


		/* Add donor edge to collections */
		edge_marker current_edge_marker(donor_edge.start, donor_edge.end);
		bool already_have = (dg_used.count(current_edge_marker) != 0);
		if (!already_have) {
			donor_edges[num_edges] =
				dg_info(contig, start, end, orientation);
			edges[num_edges++] = donor_edge;
			dg_used[current_edge_marker] = true;
		}
	}
#ifdef DUMP
	write_out(edges, nodes, "./withDonor.edges");
#endif
#ifndef DUMP

	/* Read in repeat ranges */
	map < pos, uint64_t > repeats;
	FILE *rinput = fopen(argv[4], "r");
	if (rinput == NULL) {
		perror("opening ranges");
		return -1;
	}

	while (fgets(line, MAX_LINE_LEN, rinput) != NULL) {
		char contigname[CONTIGNAME_LEN];
		uint64_t contigstart;
		uint64_t contigend;

		if (sscanf(line, "%s %llu %llu\n",
					contigname, &contigstart, &contigend) == 3) {
			pos p(contigname, contigstart);
			repeats[p] = (contigend - contigstart + 1);
		}
	}
	fclose(rinput);
#endif

#ifndef DUMP
	/* Calculate hits and masked length */
	input = open(argv[3], O_RDONLY);
	if (input < 0) {
		perror("opening input");
		return -1;
	}


	/* Get size of input file */
	fstat(input, &(inputstats));
	filesize = inputstats.st_size;

	/* Memory map input */
	filedata = (uint8_t *) mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, input, 0);
	if (filedata == MAP_FAILED) {
		perror("mmap on coverage");
		return -1;
	}
	madvise(filedata, filesize, MADV_SEQUENTIAL);

	double *sum_normodds = (double*)filedata;
	uint64_t contig_length = 1+(filesize / sizeof(double));

	map < pos, uint64_t >::iterator edgeiter = pos2edge.begin();
	for (edgeiter = pos2edge.begin();
			edgeiter != pos2edge.end(); edgeiter++) {
		pos startpos = edgeiter->first;
		uint64_t edgeindex = edgeiter->second;

		uint64_t masked_edgelen = 0;
		uint64_t edgelen_remaining = edges[edgeindex].length;

		pos endpos(startpos.contig,
				startpos.contigpos + edgelen_remaining - 1);

		/* Mask out repeat regions */
		map < pos, uint64_t >::iterator repiter
			= repeats.lower_bound(startpos);
		map < pos, uint64_t >::iterator repiter_end
			= repeats.upper_bound(endpos);
		/* Check for overlap with previous repeat region */
		if (repiter != repeats.begin()) {
			repiter--;
			uint64_t prevend =
				(repiter->first).contigpos + repiter->second - 1;
			if (prevend >= startpos.contigpos) {
				uint64_t masked_len = prevend - startpos.contigpos + 1;
				if (masked_len >= edgelen_remaining) {
					edgelen_remaining = 0;
				} else {
					edgelen_remaining -= masked_len;
					startpos.contigpos += masked_len;
				}
			}
		}

		repiter = repeats.lower_bound(startpos);
		while ((edgelen_remaining > 0) && (repiter != repiter_end)) {
			/* Get coverage prior to repeat region */
			if ((repiter->first).contigpos > startpos.contigpos) {
				uint64_t this_edgelen = (repiter->first).contigpos
					- startpos.contigpos;
				if (this_edgelen > edgelen_remaining)
					this_edgelen = edgelen_remaining;
				edges[edgeindex].hits += get_scov_arrivals(sum_normodds, contig_length, startpos.contigpos, startpos.contigpos + this_edgelen - 1,lambdas);
				edges[edgeindex].expected+=get_expected(lambdas,startpos.contigpos,startpos.contigpos+this_edgelen-1);
				masked_edgelen += this_edgelen;
				edgelen_remaining -= this_edgelen;
			}

			/* Advance startpos */
			startpos.contigpos =
				((repiter->first).contigpos + repiter->second);
			if (startpos.contigpos > endpos.contigpos) {
				edgelen_remaining = 0;
				break;
			}

			if (repiter->second >= edgelen_remaining) {
				edgelen_remaining = 0;
				break;
			}
			edgelen_remaining -= repiter->second;
			repiter++;
		}

		/* Add coverage after last repeat region */
		if (edgelen_remaining > 0) {
			edges[edgeindex].hits += get_scov_arrivals(sum_normodds, contig_length, startpos.contigpos, startpos.contigpos + edgelen_remaining - 1,lambdas);
			edges[edgeindex].expected+=get_expected(lambdas,startpos.contigpos,startpos.contigpos+edgelen_remaining-1);
			masked_edgelen += edgelen_remaining;
		}


		assert(isfinite(edges[edgeindex].expected));
		edges[edgeindex].unmasked_length += masked_edgelen;
		edges[edgeindex].num_paths++;
	}
#endif

	/* Hook up free edges */
	for (edgeiter = pos2edge.begin();
			edgeiter != pos2edge.end(); ++edgeiter) {
		pos p = edgeiter->first;
		uint64_t cur_edge = edgeiter->second;

		/* Get next position/edge */
		map < pos, uint64_t >::iterator next = edgeiter;
		next++;
		if (next != pos2edge.end()) {
			pos nextpos = next->first;
			uint64_t next_edge = next->second;

			uint64_t fakeedge = num_edges++;
			edges[fakeedge] = edge();
			edges[fakeedge].type = FAKE_TYPE;
			edges[fakeedge].is_free = true;
			if (POSSTRAND(&p)) {
				edges[fakeedge].start = edges[cur_edge].end;
				edges[fakeedge].from_bit = 0;
			} else {
				edges[fakeedge].start = edges[cur_edge].start;
				edges[fakeedge].from_bit = 1;
			}

			if (POSSTRAND(&nextpos)) {
				edges[fakeedge].end = edges[next_edge].start;
				edges[fakeedge].to_bit = 1;
			} else {
				edges[fakeedge].end = edges[next_edge].end;
				edges[fakeedge].to_bit = 0;
			}
		}
	}
#ifdef DUMP
	write_out(edges, nodes, "./final.edges");
	return 1;
#endif
	/* Print out nodes */
	printf("%llu\n", num_nodes + num_edges);
	for (uint64_t i = 0; i < num_nodes; i++) {
		printf("NODE %llu 0\n", i + 1);
		set < pos >::iterator nodeiter = nodes[i].positions.begin();
		while (nodeiter != nodes[i].positions.end()) {
			pos p = *nodeiter;
			if (p.printpos) {
				fprintf(stderr, "NODE %llu %s %llu %c\n",
						i + 1, p.contig.c_str(), p.contigpos, p.strand);
			}

			++nodeiter;
		}
	}

	/* Print fake nodes */
	for (uint64_t i = 0; i < num_edges; i++) {
		printf("NODE %llu 0\n", num_nodes + i + 1);
	}

	/* Print out edges */
	uint64_t called_length = 0;
	double called_normhits = 0.0;
	uint64_t funny_edges = 0;
	uint64_t cb_edges = 0;



	for (uint64_t i = 0; i < num_edges; i++) {
		int64_t from = edges[i].start + 1;
		if (edges[i].from_bit != 0)
			from = -from;
		int64_t to = edges[i].end + 1;
		if (edges[i].to_bit == 0)
			to = -to;
		int64_t fakenode = num_nodes + i + 1;
		double edge_lambda=-1.0;
		if (edges[i].is_cb) {
			/*if(edges[i].num_paths > 1)
			  {
			  fprintf(stderr, "Contig break with more than two paths!\n");
			  fprintf(stderr, "ARC %llu from %ld to %ld!\n",
			  i, from, to);
			  return -1;
			  } */

			printf("ARC %ld %ld 0.0 0 0.0\n", from, fakenode);
			cb_edges++;
		} else {
			edge_lambda = (edges[i].is_free || edges[i].expected==0 || edges[i].unmasked_length==0 ?
					0.0 : (edges[i].expected /
						edges[i].unmasked_length));
			assert(isfinite(edge_lambda)); 
			//if (edges[i].is_free!=true && edges[i].expected==0 && edges[i].unmasked_length>0) {
			//		fprintf(stderr,"#EDGE %llu , %f , %llu\n",i,edges[i].expected,edges[i].unmasked_length);
			//  }
			uint64_t edge_called_length = (edges[i].is_free || edges[i].expected==0 ?
					0 : (edges[i].unmasked_length /
						edges[i].num_paths));
			printf("ARC %ld %ld %f %lu %f\n",
					from,
					fakenode,
					edges[i].hits, edge_called_length, edge_lambda);
			if ((edge_called_length == 0)
					|| (edges[i].unmasked_length == 0))
				funny_edges++;

			//called_length += edge_called_length;
			called_length +=
				(edges[i].is_free ? 0 : edges[i].unmasked_length);
			called_normhits += edges[i].hits;
		}
		printf("ARC %ld %ld %f %lu\n", fakenode, to, 0.0, 0, 0.0);

		fprintf(stderr, "ARC %lu %ld %ld %ld %lu %f\n",
				i,
				labs(from),
				fakenode,
				labs(to), edges[i].length, edge_lambda);
	}

	/* Add last (fake) edge */
	printf("ARC %lu %lu 0.0 -1 0\n",
			edges[(pos2edge.rbegin()->second)].end + 1,
			edges[(pos2edge.begin()->second)].start + 1);

	/* Print outgoing edge positions */
	edgeiter = pos2edge.begin();
	while (edgeiter != pos2edge.end()) {
		pos from_pos = edgeiter->first;
		uint64_t edgeindex = edgeiter->second;
		uint64_t from_node;
		if (POSSTRAND(&from_pos))
			from_node = edges[edgeindex].start;
		else
			from_node = edges[edgeindex].end;
		fprintf(stderr, "PATH %s %lu %c %lu %lu\n",
				from_pos.contig.c_str(),
				from_pos.contigpos,
				from_pos.strand, from_node + 1, edgeindex);
		edgeiter++;
	}

	/* Print donor edge info */
	for (map < uint64_t, dg_info >::iterator dgiter = donor_edges.begin();
			dgiter != donor_edges.end(); dgiter++) {
		uint64_t edgeindex = dgiter->first;
		dg_info info = dgiter->second;
		fprintf(stderr, "DG: %lu %s %lu %lu %d\n",
				edgeindex,
				info.contig.c_str(),
				info.start, info.end, info.orientation);
	}

	fprintf(stderr, "GRAPHINFO: %llu %llu %llu %llu %llu\n", called_length,
			(uint64_t) called_normhits, funny_edges, cb_edges, 0);



	return 0;
}
