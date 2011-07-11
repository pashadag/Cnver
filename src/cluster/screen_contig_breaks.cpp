#include<cassert>
#include<cmath>
#include<cstring>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
#include "getopt.h"

using namespace std;
//#include "defs.h"

class interval {
public:
	string chr;
	int start;
	int end;
	bool operator< (const interval & i) const {
		if (chr == i.chr) return end < i.end;
		return chr < i.chr;
	}
	

	bool contains(const interval & c1) {
		if (c1.chr != chr) return false;
		if ((start <= c1.start) && (c1.end <= end)) return true;
		return false;
	}
};

bool mycomp (const interval & i1, const interval & i2) {
	if (i1.chr != i2.chr) return (i1.chr < i2.chr);
	return (i1.end <= i2.start);
	//assuming intervals are non-overlapping
}



int main(int argc, char * argv[]) {

	int tolerance, mean, sd;
	int option_index = 0; 
	char ch;
	vector<interval> contigBreaks;
	string sbuf, wholeLine;
	string filename;

	static struct option long_options[] =
	{
		{"tolerance", required_argument, 0, 'a'},
		{"mean", required_argument, 0, 'a'},
		{"stdev", required_argument, 0, 'a'},
		{"breaksFile", required_argument, 0, 'a'},
		{0, 0, 0, 0}
	};

	while ((ch = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
		switch (ch) {
			case 0:
				//printf ("option %s", long_options[option_index].name);
				//if (optarg) printf (" with arg %s", optarg);
				if (strcmp(long_options[option_index].name, "tolerance")==0) tolerance = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "breaksFile")==0) filename = optarg;
				else if (strcmp(long_options[option_index].name, "mean")==0) mean = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "stdev")==0) sd = atoi(optarg);
				break;
		}
		option_index=0;
	}

	ifstream breaks_file(filename.c_str());

	while (getline(breaks_file, sbuf)) {
		istringstream line(sbuf);
		interval i;
		line >> i.chr >> i.start >> i.end;
		contigBreaks.push_back(i);
	}

	//sort contig breaks
	sort(contigBreaks.begin(), contigBreaks.end());


	while (getline(cin, wholeLine)) {
		istringstream line(wholeLine);
		interval origKey, startKey, endKey;
		int id, clusType ;
		line >> origKey.chr >> origKey.start >> origKey.end >> clusType >> id;


		//searching for left point
		startKey = origKey;
		if (clusType != 3) {
			startKey.start += tolerance;
		} else {
			startKey.start -= tolerance;
		}
		startKey.end = startKey.start;
		vector<interval>::iterator hit = lower_bound(contigBreaks.begin(), contigBreaks.end(), startKey);
		if (hit != contigBreaks.end()) {
			if (hit->contains(startKey)) {
				//cout << "DISQUALIFY LEFT \t";
				//cout << "KEY " << origKey.chr << " " << add_commas(origKey.start) <<  " " << add_commas(origKey.end) << " " << clusType << " " << id << "\t\t\t";
				//cout << " HIT " << hit->chr << " " << add_commas(hit->start) <<  " " << add_commas(hit->end) << endl;
				continue;
			}
		} 

		//searching for right point
		endKey = origKey;
		if (clusType != 2) {
			endKey.end -= tolerance;
		} else {
			endKey.end += tolerance;
		}
		endKey.start = endKey.end;
		hit = lower_bound(contigBreaks.begin(), contigBreaks.end(), endKey);
		if (hit != contigBreaks.end()) {
			if (hit->contains(endKey)) {
				//cout << "DISQUALIFY RIGHT\t";
				//cout << "KEY " << origKey.chr << " " << add_commas(origKey.start) <<  " " << add_commas(origKey.end) << " " << clusType << " " << id << "\t\t\t";
				//cout << " HIT " << hit->chr << " " << add_commas(hit->start) <<  " " << add_commas(hit->end) << endl;
				continue;
			}
		} 

		cout << wholeLine << endl;

	}
}

