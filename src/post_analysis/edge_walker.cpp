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
#include <stack>
#include <assert.h>

#include "../include/graphtypes.h"
#include "../include/dbtypes.h"

#ifdef SINGLE
#define CUTOFF      150
#else
#define CUTOFF 1000
#endif

#define GAIN 0.5
#define LOSS -0.5
#define SAME 0.0

#define MAX(x, y)   ((x) > (y) ? (x) : (y))
#define POSSTRAND(p) (p->strand == '+')

#define MAX_LINE_LEN        1024

using namespace std;

const char * STR_NODE   = "NODE";
const char * STR_ARC    = "ARC";
const char * STR_PATH   = "PATH";

/* Position */
typedef struct {
    char            contigname[CONTIGNAME_LEN];
    uint64_t        contigpos;
    char            strand;
} pos_t;

struct ltpos
{
    bool operator()(const pos_t * p1, const pos_t * p2) const
    {
        int res = strcmp(p1->contigname, p2->contigname);
        if (res == 0)
            return p1->contigpos < p2->contigpos;
        else
            return res < 0;
    }
};


class node {
public:
    node() {};
    void add(pos_t * pos, bool override = false);
    void merge(node * other);
    bool is_pos_positive(pos_t * pos);
    set<pos_t *, ltpos> positions();
    set<pos_t *, ltpos> immutable_positions;
    set<pos_t *, ltpos> mutable_positions;
};

set<pos_t *, ltpos> node::positions()
{
    set<pos_t *, ltpos> result;
    set_union(immutable_positions.begin(),
                immutable_positions.end(),
                mutable_positions.begin(),
                mutable_positions.end(),
                insert_iterator<set<pos_t *, ltpos> >(result, result.begin()));
    return result;
}

void node::add(pos_t * pos, bool override)
{
    /* Erase from mutable positions */
    mutable_positions.erase(pos);
    
    /* Add */
    if (override)
        immutable_positions.insert(pos);
    else if (immutable_positions.find(pos) == immutable_positions.end())
        mutable_positions.insert(pos);
}

void node::merge(node * other)
{
    /* Merge mutable positions */
    set<pos_t *, ltpos>::iterator iter;
    iter = mutable_positions.begin();
    while (iter != mutable_positions.end())
    {
        add(*iter, false);
        iter++;
    }
    iter = immutable_positions.begin();
    while (iter != immutable_positions.end())
    {
        add(*iter, true);
        iter++;
    }
}

bool node::is_pos_positive(pos_t * pos)
{
    /* Find node */
    set<pos_t *, ltpos>::iterator iter;
    
    /* Check immutable positions */
    iter = immutable_positions.find(pos);
    if (iter != immutable_positions.end())
    {
        return ((*iter)->strand == '+');
    }
    
    /* Check mutable positions */
    iter = mutable_positions.find(pos);
    return ((*iter)->strand == '+');
}

class edge {
public:
    edge(uint64_t from, uint64_t to, uint64_t toalias, uint64_t edgelen,
        bool is_donor = false, int8_t orientation = -1) 
        : from(from), to(to), toalias(toalias), edgelen(edgelen),
            orientation(orientation), is_donor(is_donor), actual_flow(0.0),
            solved_flow(0.0) {
	saved_solved_flow=solved_flow;
	hits=0.0;
        expected=0.0;
	unmasked_length=edgelen;
	dirty=false;
};
    bool	dirty;
    uint64_t    from;
    uint64_t    to;
    uint64_t    toalias;
    uint64_t    edgelen;
    int8_t      orientation;
    bool        is_donor;
    double      actual_flow;
    double      solved_flow;
    double	hits;
    double	saved_solved_flow;
    double expected;
    uint64_t	unmasked_length;
};

typedef enum {
    SRC_NODE,
    DEST_NODE
} node_type_t;

struct ltnode
{
    /* Just do a pointer comparison... */
    bool operator()(const node * n1, const node * n2) const
    {
        return ((uint64_t)n1) < ((uint64_t)n2);
    }
};

struct ltedge
{
    /* Just do a pointer comparison... */
    bool operator()(const edge * n1, const edge * n2) const
    {
        return ((uint64_t)n1) < ((uint64_t)n2);
    }
};

class dg_info {
public:
    dg_info() {};
    dg_info(string contig, uint64_t start, uint64_t end, int8_t orientation)
        : contig(contig), start(start), end(end), orientation(orientation) {};
    string contig;
    uint64_t start;
    uint64_t end;
    int8_t orientation;
};


#ifdef DUMP
int has_opened=0;
void write_out(map<pos_t *, uint64_t, ltpos>::iterator walker,map<uint64_t, edge *> &edges, double max_path_diff, int max_path_len, uint64_t max_edges_used) {
    /* Write out graph */
    stack<uint64_t> path_taken;
    char filename[200];
    sprintf(filename,"dumpcnvs.paths");
    FILE* f_output=NULL;
    if (has_opened==0) {
	f_output=fopen(filename,"w");
	has_opened=1; 
   } else {	
	f_output = fopen(filename, "a");
    }
    if (f_output==NULL) {
	printf("ERROR OPENING LOG FILE\n");
	exit(1);
    }
    fprintf(f_output, "MAX_PATH_DIFF %f , MAX_PATH_LEN %d\n",max_path_diff,max_path_len); 
    while (max_edges_used > 0 )
    {
	max_edges_used--;
	pos_t * p = walker->first;
	edge * e = edges[walker->second];
	path_taken.push(p->contigpos);
	//e->solved_flow -= max_path_diff;
	double ratio=e->hits/(e->unmasked_length*e->actual_flow);
	max_path_len -= e->unmasked_length;
	fprintf(f_output,"E:%10d %12llu %12llu %c  SS:%2.2f S:%2.2f ACT:%2.2f ARR:%12.4f  LEN:%12llu UN-LEN:%12llu  ARR/(UN-LEN*ACT):%2.2f\n",walker->second,p->contigpos,p->contigpos+e->edgelen-1,p->strand,e->saved_solved_flow,e->solved_flow,e->actual_flow,e->hits,e->edgelen,e->unmasked_length,ratio);

	walker++;
    }
    fprintf(f_output,"FULL PATH: ");
    while (!path_taken.empty()) {
	fprintf(f_output, "%llu ", path_taken.top());
	path_taken.pop();
    }
    fprintf(f_output,"\n");
    fclose(f_output);
}   
#endif

void update_path(stack<uint64_t> &used_edges, map<uint64_t, edge *> &edges,double flow) {
        stack<uint64_t> temp;
        while(!used_edges.empty()) {
                uint64_t edge_num=used_edges.top();
                edge * e = edges[edge_num];
                e->solved_flow += flow;
                temp.push(edge_num);
                used_edges.pop();
        }
        while(!temp.empty()) {
                used_edges.push(temp.top());
                temp.pop();
        }

}

int avg=0;

/* Calls CNVs in a solved flow graph */
int main(int argc, char ** argv)
{
    /* Validate input, but not too much. */
    if (argc < 5)
    {
		fprintf(stderr, "usage: %s <graph_info> <solved_flow> <start_pos> <end_pos> <AVG_CALL(Y?)>\n", argv[0]);
		return -1;
	} else if (argc==6) {
		avg=1;
	}

	uint64_t start_pos=atoll(argv[3]);
	uint64_t end_pos=atoll(argv[4]);    

	/* Open graph info */
	FILE * input = fopen(argv[1], "r");
	if (input == NULL)
	{
		perror("opening input");
		return -1;
	}

	/* Read in graph info */
	map<uint64_t, node *> nodes;
	map<pos_t *, uint64_t, ltpos> pos2node;
	map<uint64_t, edge *> edges;
	map<uint64_t, uint64_t> edgealias;
	map<pos_t *, uint64_t, ltpos> pos2edge;
	map<uint64_t, dg_info> donor_edges;

	char line[MAX_LINE_LEN];
	while (fgets(line, MAX_LINE_LEN, input) != NULL)
	{
		/* Cheap */
		switch (line[0])
		{
			case 'N':
				{
					/* Node/pos entry */
					char linetype[MAX_LINE_LEN];
					uint64_t nodenum;
					char contigname[CONTIGNAME_LEN];
					uint64_t contigpos;
					char strand;

					sscanf(line, "%s %llu %s %llu %c\n",
							linetype,
							&nodenum,
							contigname,
							&contigpos,
							&strand);

					/* Create new node if it doesn't exist */
					node * n;
					if (nodes.find(nodenum) == nodes.end())
					{
						n = new node();
						nodes[nodenum] = n;
					}
					else
					{
						n = nodes[nodenum];
					}

					/* Add position to node */
					pos_t * pos = new pos_t();
					strcpy(pos->contigname, contigname);
					pos->contigpos = contigpos;
					pos->strand = strand;
					n->add(pos);
					pos2node[pos] = nodenum;
				}
				break;
			case 'A':
				{
					/* Arc entry */
					char linetype[MAX_LINE_LEN];
					uint64_t edgenum;
					uint64_t from;
					uint64_t toalias;
					uint64_t to;
					uint64_t len;

					sscanf(line, "%s %llu %llu %llu %llu %llu\n",
							linetype,
							&edgenum,
							&from,
							&toalias,
							&to,
							&len);

					/* Create edge */
					edge * e = new edge(from, to, toalias, len);
					edges[edgenum] = e;
					edgealias[toalias] = to;
				}
				break;
			case 'P':
				{
					/* Path entry */
					char linetype[MAX_LINE_LEN];
					char contigname[CONTIGNAME_LEN];
					uint64_t contigpos;
					char strand;
					uint64_t nodenum;
					uint64_t edgenum;

					sscanf(line, "%s %s %llu %c %llu %llu\n",
							linetype,
							contigname,
							&contigpos,
							&strand,
							&nodenum,
							&edgenum);

					/* Add pos2edge entry */
					pos_t * pos = new pos_t();
					strcpy(pos->contigname, contigname);
					pos->contigpos = contigpos;
					pos->strand = strand;
					pos2edge[pos] = edgenum;
				}
				break;
			case 'D':
				{
					/* Donor edge */
					char linetype[MAX_LINE_LEN];
					uint64_t edgenum;
					char contigname[CONTIGNAME_LEN];
					uint64_t start;
					uint64_t end;
					int8_t orientation;

					sscanf(line, "%s %llu %s %llu %llu %d\n",
							linetype,
							&edgenum,
							contigname,
							&start,
							&end,
							&orientation);

					/* Add dg_info entry */
					donor_edges[edgenum] = 
						dg_info(contigname, start, end, orientation);
				}
				break;
		}
	}

	fclose(input);

	/* Read in flow */
	input = fopen(argv[2], "r");
	if (input == NULL)
	{
		perror("opening input");
		return -1;
	}

	uint64_t cur_edge = 0;
	while (fgets(line, MAX_LINE_LEN, input) != NULL)
	{
		/* Cheap again */
		if (line[0] == 'A')
		{
			/* Arc entry */
			char linetype[MAX_LINE_LEN];
			int64_t from;
			int64_t to;
			double hits;
			int64_t len;
			double flow_count = 0.0;
			double edge_lambda=0.0;
			int scanned=sscanf(line, "%s %lld %lld %lf %lld %lf %lf\n",
					linetype,
					&from,
					&to,
					&hits,
					&len,
					&edge_lambda,
					&flow_count);

			if (scanned==6) {
				flow_count=edge_lambda;
				edge_lambda=-1.0;
			}
			//printf("%s,%lld,%lld,%lf,%lld,%lfhits\n%s",linetype,from,to,hits,len,flow_count,line);
			/* Skip dummy edges */
			if ((edgealias.find(from) != edgealias.end()))
				continue;



			/* Assumes edges are in the same order as we gave them */
			edge * e = edges[cur_edge];
			if (e == NULL)
			{
				cur_edge++;
				continue;
			}

			/* Sanity check */
			if ((e->from != labs(from)) || (e->to != edgealias[labs(to)]))
			{
				fprintf(stderr, "Invalid edge line: \"%s\"\n", line);
				fprintf(stderr, "Was expecting %llu->%llu\n", e->from, e->to);
				return -1;
			}


			/* Set flow count for edge */
			//printf("READ IN %f\n%s\n\n",edge_lambda,line);
			e->expected=edge_lambda*len/2;
			e->solved_flow = flow_count;
			e->unmasked_length=len;
			e->saved_solved_flow = flow_count;
			if (e->hits<hits) {
				//printf("Setting hits %f\n",hits);
				e->hits=hits;
			}
			/* If this is a donor edge, and there's flow on it, */

			cur_edge++;
		}
	}
	fclose(input);


	/* Walk graph, increasing actual_flow */
	map<pos_t *, uint64_t, ltpos>::iterator walker = pos2edge.begin();
	while (walker != pos2edge.end())
	{
		/* Increment edge */
#ifndef SINGLE
		edges[walker->second]->actual_flow += 2.0;
#else
		edges[walker->second]->actual_flow += 1.0;
#endif
		walker++;
	}

	walker=pos2edge.begin();
	/* Walk graph, increasing actual_flow */
	while (walker != pos2edge.end())
	{
		uint64_t edgenum=walker->second;
		double solved_flow=edges[edgenum]->solved_flow;
		double actual_flow=edges[edgenum]->actual_flow;
		if ((solved_flow-int(solved_flow))>0.0) {
			double new_solved_flow;
			if (solved_flow>actual_flow) {
				new_solved_flow=solved_flow-0.5;
			} else {
				new_solved_flow=solved_flow+0.5;
			}
			edges[edgenum]->solved_flow=new_solved_flow;
		}
		walker++;
	}

	walker=pos2edge.begin();
	pos_t tmp;
	tmp = *(walker->first);
	tmp.contigpos=start_pos;

	walker=pos2edge.lower_bound(&tmp);
	if (walker!=pos2edge.begin()) {
		walker--;
	}
	tmp.contigpos=end_pos;
	map<pos_t *, uint64_t, ltpos>::iterator walker_end = pos2edge.lower_bound(&tmp);
	if (walker_end!=pos2edge.end()) {
		if ( walker_end->first->contigpos==end_pos) {
			walker_end++;
		}
	}

	double wavg_length=0;
	double wavg_weight=0;

	double wavg_weight_act=0;


#define BUCKETS 2000   
	double counter[BUCKETS*2];
	for (int i=0; i<BUCKETS*2; i++) {
		counter[i]=0;
	}

	//get the count of each edge
	map<uint64_t,uint64_t> edge_usage; 
	while (walker!=walker_end) {
		uint64_t edgenum=walker->second;
		edge * e=edges[edgenum];
		pos_t * p=walker->first;
		//printf("Walking %llu\n",edgenum);
		if (edge_usage.count(edgenum)>0) {
			edge_usage[edgenum]++;
		} else {
			edge_usage[edgenum]=1;
		}
		walker++;
	} 

	map<uint64_t,uint64_t> edge_remainder; 
	map<uint64_t,uint64_t> edge_flow; 
	map<uint64_t,uint64_t>::iterator it;
	for ( it=edge_usage.begin() ; it != edge_usage.end(); it++ ) {
		uint64_t edgenum = (*it).first;
		edge_remainder[edgenum]=int(edges[edgenum]->saved_solved_flow)%edge_usage[edgenum];
		edge_flow[edgenum]=int(edges[edgenum]->saved_solved_flow)/edge_usage[edgenum];
	}

	//go to the begining again 
	walker=pos2edge.begin();
	tmp = *(walker->first);
	tmp.contigpos=start_pos;

	walker=pos2edge.lower_bound(&tmp);
	if (walker!=pos2edge.begin()) {
		walker--;
	}
	while (walker!=walker_end) {
		uint64_t edgenum=walker->second;
		edge * e=edges[edgenum];
		pos_t * p=walker->first;
		double ratio=e->hits/(e->unmasked_length*e->actual_flow);
		int64_t edge_start=(int64_t)p->contigpos;
		int64_t edge_end=(int64_t)p->contigpos+(int64_t)e->edgelen-1;
		int64_t relative_edge_start=edge_start-(int64_t)start_pos;
		int64_t relative_edge_end=edge_end-(int64_t)start_pos;

		//weighted average
		int64_t wavg_start=edge_start;
		if ((int64_t)start_pos>wavg_start) {
			wavg_start=(int64_t)start_pos;
		}

		int64_t wavg_end=edge_end;
		if (end_pos<wavg_end) {
			wavg_end=end_pos;
		}

		//if (!(e->dirty)) {
		wavg_length+=wavg_end-wavg_start+1;
		wavg_weight+=(wavg_end-wavg_start+1)*(e->saved_solved_flow-e->actual_flow);
		wavg_weight_act+=(wavg_end-wavg_start+1)*(e->actual_flow);
		assert(abs(int(e->saved_solved_flow))<(BUCKETS-1));
		//at BUCKETS = 0, less is negative, more is positive
		uint64_t flow = edge_flow[edgenum];
		if (edge_remainder[edgenum]>0) {
			flow++;
			edge_remainder[edgenum]--;
		}	
		counter[flow+BUCKETS]+=wavg_end-wavg_start+1;
		e->dirty=true;
		//printf("U: %llu, %llu\n",edge_usage[edgenum],flow);
		//printf("E:%10d %12lld %12lld, %8lld %8lld, %c  SS:%2.2f S:%2.2f ACT:%2.2f ARR:%12.4f  LEN:%12llu UN-LEN:%12llu EX:%f EXBP:%f,\t%d\n",walker->second,edge_start,edge_end,relative_edge_start,relative_edge_end,p->strand,e->saved_solved_flow,e->solved_flow,e->actual_flow,e->hits,e->edgelen,e->unmasked_length,e->expected*e->actual_flow,e->expected,edge_usage[edgenum]);
		//}	
		if (avg==0) {
			printf("E:%10d %12lld %12lld, %8lld %8lld, %c  SS:%2.2f S:%2.2f ACT:%2.2f ARR:%12.4f  LEN:%12llu UN-LEN:%12llu EX:%f EXBP:%f\n",walker->second,edge_start,edge_end,relative_edge_start,relative_edge_end,p->strand,e->saved_solved_flow,e->solved_flow,e->actual_flow,e->hits,e->edgelen,e->unmasked_length,e->expected*e->actual_flow,e->expected);
		}
		walker++;

	}

	//Compute statistics
	double mean=0;
	double mode=0;
	double median=0;
	double total=0;
	for (int i=0; i<BUCKETS*2; i++) {
		total+=counter[i];
	}
	median=total/2;
	for (int i=0; i<BUCKETS*2; i++) {
		if (counter[i]!=0) {
			printf("BUCKET: %d\t%f\n",i-BUCKETS,counter[i]);
		}
	}
	double mode_freq=0;
	bool found_median=false;
	for (int i=0; i<BUCKETS*2; i++) {
		if (!found_median) {
			median-=counter[i];
			if (median<=0 ) {
				median=i-BUCKETS;
				//printf("Median is %f\n",median);
				found_median=true;
			}
		}
		if (counter[i]>mode_freq) {
			mode=i-BUCKETS;
			mode_freq=counter[i];
		}
		mean+=counter[i]*(i-BUCKETS);
	}
	mean=mean/total;

	printf("START\tEND\tLENGTH\tWAVG\tWAVG_ACTUAL\tWEIGHTED_MEAN\tWEIGHTED_MODE\tWEIGHTED_MEDIAN\n");
	if (avg==1 && wavg_length!=0) {
		printf("%llu\t%llu\t%.0f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",start_pos,end_pos,wavg_length,wavg_weight/wavg_length,wavg_weight_act/wavg_length,mean,mode,median);
	} else if (avg==1) {
		printf("%llu\t%llu\t%.0f\t%.4f\t%.4f\t%.4f\t%.4f\t%.f\n",start_pos,end_pos,wavg_length,0,wavg_weight_act/wavg_length,0,0,0);
	}
	return 0;
}
