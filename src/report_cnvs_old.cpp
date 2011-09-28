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

#include "include/graphtypes.h"
#include "include/dbtypes.h"

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
				dirty=false;
				expected=0.0;
				unmasked_length=edgelen;
				decide_length=edgelen;
			};
		bool dirty;
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
		uint64_t	unmasked_length;
		double expected;
		uint64_t    decide_length;
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
void write_out(map<pos_t *, uint64_t, ltpos>::iterator walker,map<uint64_t, edge *> &edges, double max_path_diff, int max_path_len, uint64_t max_edges_used, uint64_t cut_off, uint8_t length_type) {
	/* Write out graph */
	stack<uint64_t> path_taken;
	char filename[200];
	sprintf(filename,"dump.c%llu.t%u.paths",cut_off,length_type);
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
		fprintf(f_output,"E:%10d %12llu %12llu %c  SS:%2.2f S:%2.2f ACT:%2.2f ARR:%12.4f  LEN:%12llu UN-LEN:%12llu DLEN:%12llu ARR/(ULEN*ACT):%2.2f\n",walker->second,p->contigpos,p->contigpos+e->edgelen-1,p->strand,e->saved_solved_flow,e->solved_flow,e->actual_flow,e->hits,e->edgelen,e->unmasked_length,e->decide_length,ratio);

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


/* Calls CNVs in a solved flow graph */
int main(int argc, char ** argv)
{
	/* Validate enput, but not too much. */
	if (argc != 5)
	{
		fprintf(stderr, "usage: %s <graph_info> <solved_flow> <cutoff> <[0-unmaksed_length|1-total_length]>\n", argv[0]);
		return -1;
	}

	char* graph_info_filename=argv[1];
	char* solved_flow_filename=argv[2];
	uint64_t cut_off=atol(argv[3]);
	uint8_t length_type=atoi(argv[4]);

	/* Open graph info */
	FILE * input = fopen(graph_info_filename, "r");
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
	input = fopen(solved_flow_filename, "r");
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
			e->solved_flow = flow_count;
			if (length_type==0) {
				e->decide_length=len;
			}
			assert(e->expected==0.0);
			e->expected=len*edge_lambda/2;
			e->unmasked_length=len;
			e->saved_solved_flow = flow_count;
			if (e->hits<hits) {
				//printf("Setting hits %f\n",hits);
				e->hits=hits;
			}
			/* If this is a donor edge, and there's flow on it, */
			if (donor_edges.count(cur_edge) != 0)
			{
				if (flow_count != 0.0)
				{
					/* Print out the edge info */
					printf("%s %llu %llu %d %f\n",
							donor_edges[cur_edge].contig.c_str(),
							donor_edges[cur_edge].start,
							donor_edges[cur_edge].end,
							donor_edges[cur_edge].orientation,
							flow_count);
				}
			}

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

	fprintf(stderr,"#contig_name start_position end_position length flow_diff non_dirty_unmasked non_dirty_total\n");
	/* Walk again, looking for differences */
	bool no_diffs;
	do {
		stack<uint64_t> used_edges;
		no_diffs = true;
		double last_diff = 0.0;
		double last_sign = 0.0;
		char last_contigname[CONTIGNAME_LEN];
		uint64_t last_start = 0;
		uint64_t last_end = 0;
		uint64_t path_len = 0;
		pos_t max_path_start = { 0 };
		uint64_t max_path_len = 0;
		double max_path_diff = 0.0;
		uint64_t this_edges_used=0;
		uint64_t max_edges_used=0;
		pos_t this_path_start = { 0 };
		uint64_t this_path_len = 0;
		double this_path_mindiff = 0.0;
		double max_path_flow_sign=0.0;
		pos_t *  look_back_pos=NULL;
		edge * look_back_edge=NULL; 
		walker = pos2edge.begin();
		while (walker != pos2edge.end()) {
			uint64_t edgenum = walker->second;
			edge * e = edges[edgenum];
			pos_t * p = walker->first;

			double flow_diff = (e->solved_flow - e->actual_flow);
			assert(flow_diff!=0.5 || flow_diff==-0.5);
			double flow_sign = fmin(fmax(flow_diff, -0.5), 0.5);
			assert(flow_sign == GAIN || flow_sign == LOSS || flow_sign == SAME );

			/*if (edgenum==12994 || edgenum==12995) {
			  printf("X: %f %f\n",last_sign,flow_sign);
			  }
			  printf("Used edges in Stack %d, used edges %d\n",used_edges.size(),this_edges_used);*/
			assert(used_edges.size()==this_edges_used);

			if (flow_sign != last_sign) {
				/* If we're finishing a path */
				if (last_sign == GAIN || last_sign == LOSS) {
					/* Check path len against max */
					if (this_path_len > max_path_len) {
						strcpy(max_path_start.contigname, 
								this_path_start.contigname);
						max_path_start.contigpos = this_path_start.contigpos;
						max_path_len = this_path_len;
						max_path_flow_sign = last_sign;
						max_path_diff = this_path_mindiff;
						max_edges_used=this_edges_used;
					}
					update_path(used_edges,edges,this_path_mindiff);
					while(!used_edges.empty()) {
						used_edges.pop();
					}
					this_edges_used=0;       
				}
				//printf("Used edges size is %d\n",used_edges.size());
				assert(used_edges.size()==0);
				if (look_back_pos!=NULL) {
					//printf("Reseting to lookback\n");
					walker=pos2edge.find(look_back_pos);
					look_back_pos=NULL;
					look_back_edge=NULL;
					last_sign=SAME;
					assert(walker!=pos2edge.end());
					continue;
					/*assert(walker!=pos2edge.end());
					  p=look_back_pos;
					  assert(walker->first==p);
					  e=look_back_edge;
					  look_back_pos=NULL;
					  look_back_edge=NULL;
					  flow_diff = (e->solved_flow - e->actual_flow);
					  assert(flow_diff!=0.5 || flow_diff==-0.5);
					  flow_sign = fmin(fmax(flow_diff, -0.5), 0.5);
					//printf("%f\n",flow_sign);
					assert(flow_sign == GAIN || flow_sign == LOSS || flow_sign == SAME );*/
				}

				if (flow_diff != 0.0) {
					/* Start a new path */
					strcpy(this_path_start.contigname, p->contigname);
					this_path_start.contigpos = p->contigpos;
					this_path_len = e->decide_length;
					this_path_mindiff = flow_diff;
					this_edges_used=1;
					e->solved_flow-=this_path_mindiff;
					used_edges.push(edgenum);
					no_diffs = false;
				}

				last_sign = flow_sign;
			} else {
				/* If we're continuing a path, */
				if (flow_sign == GAIN || flow_sign == LOSS) {
					/* Update end of path */
					//update_path(used_edges,edges,this_path_mindiff);
					double old_flow=this_path_mindiff,new_flow,change_flow;;
					this_path_len += e->decide_length;
					assert(last_sign==flow_sign);
					this_edges_used += 1;
					/* Update minimum difference */
					if (flow_sign == GAIN ) {
						new_flow = fmin(this_path_mindiff, flow_diff);
					} else {
						assert(flow_sign==LOSS);
						new_flow = fmax(this_path_mindiff, flow_diff);
					}
					used_edges.push(edgenum);
					e->solved_flow-=old_flow;
					change_flow=old_flow-new_flow;
					if (change_flow!=0.0) {
						update_path(used_edges,edges,change_flow);
						this_path_mindiff=new_flow;
					}
				}
			}
			walker++;
		}

		/* Check last path if we're on one... */
		if (last_sign == GAIN || last_sign == LOSS)
		{
			//if (this_unmasked_path_len > max_unmasked_path_len)
			if (this_path_len > max_path_len)
			{
				strcpy(max_path_start.contigname, 
						this_path_start.contigname);
				max_path_start.contigpos = this_path_start.contigpos;
				max_path_len = this_path_len;
				max_path_flow_sign = last_sign;
				max_path_diff = this_path_mindiff;
				max_edges_used=this_edges_used;
			}
			update_path(used_edges,edges,this_path_mindiff);
			while(!used_edges.empty()) {
				used_edges.pop();
			}
		}

		assert(used_edges.empty());

		/* Check if we've reached our cutoff */
		if (max_path_len < cut_off) {
			break;
		}

		/* Print out maximal path and correct flow */
		if (max_edges_used != 0)
		{

#ifdef DUMP		    
			walker = pos2edge.find(&max_path_start);
			write_out(walker,edges,max_path_diff,max_path_len,max_edges_used,cut_off,length_type);
#endif
			uint64_t non_dirty_unmasked_length=0;
			uint64_t non_dirty_total_length=0;
			uint64_t final_total_length=0;
			walker = pos2edge.find(&max_path_start);
			while (max_edges_used > 0)
			{
				assert(walker!=pos2edge.end());
				edge * e = edges[walker->second];
				double flow_diff = (e->solved_flow - e->actual_flow);
				assert(flow_diff!=0.5 || flow_diff==-0.5);
				double flow_sign = fmin(fmax(flow_diff, -0.5), 0.5);
				assert(flow_sign == GAIN || flow_sign == LOSS || flow_sign == SAME );
				//printf("%llu: Max flow sign %f, current %f\n",walker->second,max_path_flow_sign,flow_sign);
				if (!e->dirty) {
					non_dirty_unmasked_length+=e->unmasked_length;
					non_dirty_total_length+=e->edgelen;
					e->dirty=true;
				}
				final_total_length+=e->edgelen;



				if (flow_sign==max_path_flow_sign) {
					e->solved_flow -= max_path_diff;
					max_edges_used--;
				}
				walker++;
			}
			fprintf(stderr, "%s %llu %llu %llu %f %llu %llu\n",
					max_path_start.contigname,
					max_path_start.contigpos,
					max_path_start.contigpos + final_total_length-1,
					max_path_len,
					max_path_diff,non_dirty_unmasked_length,non_dirty_total_length);
		}        
	} while (!no_diffs);

	return 0;
}
