#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<cstdio>
#include<algorithm>

#include "../include/defs.h"

using namespace std;
#define SCALE_FACTOR 10000
#define NUM_EXTRA_SEGMENTS 0 
#define MAX_SEQ_EDGES 1189999
#define ALLOW_NONE 2

//scales the value of a slope so that it is an integer.
//most of the function is doing error processing
int make_int(double cost) {
	int val;
	val = int (SCALE_FACTOR * cost);
	return val;
	//if (val != 0 || cost == 0) return val;
	//cerr << "Scaled integer edge cost is 0 -- the float value is " << cost << endl;
	//cerr << "THE SCALE_FACTOR is probably too small\n";
	//exit(1);
	return -1;
}

//l is the average length of the edge in the REF
//k is the total number of normalized reads mapping to the region (DOC)
//lambda is the expected number of reads mapping to a single position
//f is the number of edges in the donor (the flow)
inline double eval_poisson (double l, double k, double lambda, double f) {
	//log in c++ is the natural log (ln)
	double x = lambda * f * l;
	return x - k * log(x);
}

void define_arc(int from, int to, int low, int cap, double cost){
	cout << "a " << from << " " << to << " " << low << " " << cap << " " << make_int(cost) << endl;
}

void process_convex_arc(vector<string> & row) {
	int from = atoi(row[1].c_str());
	int to   = atoi(row[2].c_str());

	//k, l, and lambda are defined as in the GR paper:
	//k is the total number of normalized reads mapping to the region (DOC)
	//l is the average length of the edge in the REF
	//lambda is the expected number of reads mapping to a single position
	double k = atof(row[3].c_str());
	double l = atof(row[4].c_str());
	double lambda = atof(row[5].c_str());

	if (from == to) {
		cerr << "Found loop at " << from << endl;
		return;
		//TODO should this be a special case?
	}

	if (l == 0) { //fully masked region
		define_arc(from, to, 0, MAX_SEQ_EDGES, 1);
		return;
	}


	double center = k / (lambda * l); 
	double temp;
	double slope, lastSlope;
	int low  = max (1, (int) floor(center) - 1);			// we need at least two/three segments.
	int high = (int) ceil(center) + 1; 
	low -= NUM_EXTRA_SEGMENTS;
	high += NUM_EXTRA_SEGMENTS;
	low = max(low, 1);
	high = max(high, 2); 
	if (high > MAX_SEQ_EDGES - 1000) {
		cerr << "warning: very high copy count segment(center is " << center << "), capacity maxed out  (edge " << from << "->" << to <<  ")";
		high = MAX_SEQ_EDGES - 1000;
	}

	//add first edge, and an edge to allow zero flow
	double costofzero = eval_poisson(l, k, lambda, 1) - eval_poisson(l, k, lambda, 1.0 / (1 + ALLOW_NONE)); 
	slope = eval_poisson(l, k, lambda, low + 1) - eval_poisson(l, k, lambda, low);
	if (high - low == 1) {  //there is only one segment	
		assert(center <= 1);
		if (costofzero > 0) {
			define_arc(from, to, 0, 1, costofzero);
			define_arc(from, to, 0, MAX_SEQ_EDGES - high, slope);
		} else if (costofzero < 0) {
			define_arc(to, from, 0, 1, -1 * costofzero);
			define_arc(from, to, 1, MAX_SEQ_EDGES - high, slope);
		}
	} else if (slope >= 0) {
		assert (low == 1); //why? ah ok, think about the possible cases and you'll get it
		assert (center >= 1);
		assert (costofzero <= 0);
		define_arc(from, to, low, low + 1, slope);
		define_arc(to, from, 0, 1, -1 * costofzero);
	} else {
		assert (costofzero <= 0);
		define_arc(to, from, 0, 1, -1 * costofzero);
		define_arc(to, from, 0, low, -1 * slope);
	}

	//middle edges
	lastSlope = slope;
	for (int x = low + 1; x < high - 1; x++) {
		slope = eval_poisson(l, k, lambda, x + 1) - eval_poisson(l, k, lambda, x);
		if (slope >= 0 && lastSlope < 0) {
			define_arc(from, to, x, x + 1, slope);
		} else if (slope >= 0) {
			define_arc(from, to, 0, 1, slope);
		} else {	
			define_arc(to, from, 0, 1, -1 * slope);
		}
		lastSlope = slope;
	}

	//last edge
	if (high - low > 1) {
		slope = eval_poisson(l, k, lambda, high) - eval_poisson(l, k, lambda, high - 1); 
		assert (slope > 0);  //if crash here try increasing scalefactor
		if (lastSlope < 0) {
			define_arc(from, to, high - 1, MAX_SEQ_EDGES - high, slope);
		} else {
			define_arc(from, to, 0, MAX_SEQ_EDGES - high, slope);
		}
	}

	//TODO read the solution back?
	//TODO are loops going to work?
}

int main(int argc, char ** argv) {
	if (argc != 1) {
		fprintf(stderr, "usage: %s \n", argv[0]);
		return -1;
	}

	vector<string> row;
	string sbuf;
	while(get_row_whitespace(cin, row, sbuf)){
		if (row[0] != "convex_arc") {
			cout << sbuf << endl;
		} else {
			process_convex_arc(row);
		} 
	}
}
