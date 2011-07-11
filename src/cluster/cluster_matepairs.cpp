#include<cassert>
#include<cmath>
#include<cstring>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<set>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
#include<algorithm>
#include "getopt.h"
#include "limits.h"

using namespace std;
ostringstream oMsg;
string sbuf;

int colDist, colChr, colLeft, colRight, colID, colTemplate, baseLenFactor, mean, sd, mdJoinTolerance, avg_md_cuttoff;
int mdSignifTol;
int CONCISE;

template<class T>
class winScan {  //read in input file using a sliding window so as not to keep everything in memory.
public:
	
	deque<T> buf;	
	deque<bool> visited;

	winScan () {}
	winScan (istream * _in) : in(_in) {
		load_line();
	}
	void init (istream * _in) {
		in = _in;
		load_line();
	}
	void slideWin () {
		if (!buf.empty()){
			buf.pop_front();
			visited.pop_front();
		}
		if (buf.empty()) {
			load_line();
		}
	}	
	void incr(int & index) {
		assert (index >= 0);
		assert (!buf.empty());
		if (index == (buf.size() - 1)) {
			load_line();
		}
		index++;
	}

	istream * in;

private:
	bool load_line () {
		string s;
		if (getline(*in, s)) {
			//cout << "read line " << s << endl << flush;
			buf.push_back(make_T(s, dummyT));
			visited.push_back(false);
			return true;
		} else {
			return false;
		}
	}

	T dummyT;
	
};


string & make_T(string & s, string & dummy) {
	return s;
}

vector<string> make_T (string & s, vector<string> & dummy) {
	//tokenize the string
	vector<string> row;
	istringstream istrm(s);
	string val;
	while (istrm >> val) {
	       row.push_back(val);
	}
	return row;
}

winScan<vector<string> > in;

int make_int(const string & s) {
	return atoi(s.c_str());
}


void dump_line (vector<string> * line) {

	if (colID != -1) cout << line->at(colID);  
	cout << "\tchr = " << line->at(colChr) << "\tleft = " << line->at(colLeft) << "\tright = " << line->at(colRight);
	cout << "\tmd = " << line->at(colDist) << "\ttemplate_id = " << line->at(colTemplate);
}

void dump_cluster(vector<vector< string> > curCluster, int & clusterIndex, int clusterType) {
	if (curCluster.size() == 1) return; //single hit is not a cluster
	int i;

	//get the ranges of the matepair, as well as avg_md and sd

	int valLeft, valRight;
	double avg_md = 0;
	int x1, x2, x3, x4;
	//x1 is the leftmost start point of the left side of the generalized matepair
	//x2 is the rightmost start point of the left side of the generalized matepair
	//x3 is the leftmost start point of the right side of the generalized matepair
	//x4 is the rightmost start point of the right side of the generalized matepair

	//starting values	
	x1 = make_int(curCluster[0].at(colLeft)); 
	x2 = 0;
	x3 = make_int(curCluster[0].at(colRight));
	x4 = 0;
	for (i = 0; i < curCluster.size(); i++) {
		valRight = make_int(curCluster[i].at(colRight));
		valLeft  = make_int(curCluster[i].at(colLeft));
		if (valRight < x3) x3 = valRight;
		if (valRight > x4) x4 = valRight;
		if (valLeft  > x2) x2 = valLeft;
		if (valLeft  < x1) x1 = valLeft;

		avg_md += make_int(curCluster[i].at(colDist)) / double(curCluster.size());
	}

	//throw out if avg_md is too small
	if (avg_md < (avg_md_cuttoff)) return;

	//throw away identical mappings, even if its different templates.  this is probably an artifact.
	vector<bool> skip(curCluster.size(), false);
	string lastMd = "0";
	string lastLeft = "0";
	for (i = 0; i < curCluster.size(); i++) {
		if ((curCluster[i].at(colLeft) != lastLeft)) {
			lastMd = "0";
		}
		if (curCluster[i].at(colDist) == lastMd) skip[i] = true; 
		lastMd = curCluster[i].at(colDist);
		lastLeft = curCluster[i].at(colLeft);
	}

	//throw away all but the leftmost mapping of any template.
	set<string> templates;
	for (i = 0; i < curCluster.size(); i++) {
		if (!skip[i]) {
			pair<set<string>::iterator,bool> ret = templates.insert(curCluster[i].at(colTemplate)); //returns false if such an element already exists
			if (ret.second == false) skip[i] = true;
		}
	}

	//find standard deviation
	int total = 0;
	double sd = 0;
	for (i = 0; i < curCluster.size(); i++) {
		if (!skip[i]) {
			sd += pow((atoi(curCluster[i].at(colDist).c_str()) - avg_md),  2.0) / double(curCluster.size());
			total++;
		}
	}
	if (total == 1) return;
	sd = sqrt(sd);
	
	if (!CONCISE) {
		cout << endl;
		cout << "HEADER1 " << clusterIndex << "\t" << curCluster.at(0).at(colChr) << "\t leftBase = |(" << x1 << ", " << x2 << ")| = " << x2 - x1;
		cout<< "\t rightBase = |(" << x3 << ", " << x4 << ")| = " << x4 - x3 << "\tmd = " << avg_md << " +- " << sd << endl;
		cout << "HEADER2 "; 
		for (int i = 0; i < curCluster.size(); i++) {
			if (skip[i]) continue;
			//cout << "CLUSTER" << clusterType << " " << clusterIndex<< "\t";
			//dump_line(&curCluster[i]);
			//cout << endl << flush;
			if (colID != -1) {
				cout << curCluster[i].at(colID) << " " ;
			} else {
				cout << curCluster[i].at(colTemplate) << " ";
			}
		}
		cout << endl;
	}

	int from = -1;
	int to = 1;
	int l_from = x2 - x1; 
	int l_to = x4 - x3;

	if (clusterType == 0) {
		from = x2; to = x3; 
	} else if (clusterType == 1) {
		from = x1; to = x4; 
	} else if (clusterType == 2) {
		from = x2; to = x4;
	} else if (clusterType == 3) {
		from = x1; to = x3; 
	}

	cout << curCluster.at(0).at(colChr) << "\t";
	cout << from << "\t" << to << "\t" << clusterType << "\t" << total << "\t" << l_from << "\t" << l_to << "\t" << avg_md << "\t"  << clusterIndex << "\tEDGE\n";

	clusterIndex++;
	return;
}

void cluster_indels23(int clusterType) {

	int curIndex;
	vector<string> * outerLine;
	vector<string> * innerLine;
	vector<vector<string> >  curCluster;  
	int outerMd, innerMd;
	int curClusterIndex = 1;
	assert(!in.buf.empty());

	int colMain, colSec;
	if (clusterType == 2) {
		colMain = colLeft;
		colSec  = colRight;
	} else {
		colMain = colRight;
		colSec  = colLeft;
	}	

	//Assume input is sored by colMain
	//For every unclustered line
	//	start a potential cluster with this line, topLine
	//	nextline = curline++
	//	while nextline.colMain is within mean + baseLenFactor*sd from topline.colMain
	//		goto nextline if this line has already been clustered
	//		if the difference between the main coordinates is within mdJoinToler*sd from the dif in sec coord
	//		then add this line to cluster

	do {
		if (!in.visited.front()) { //visited really means clustered
			outerLine = &in.buf.at(0);
			curCluster.push_back(*outerLine);
			//cout << "outer " << flush << outerLine->at(0) << " " << outerLine->at(1) << " " << outerLine->at(2) << " " << outerLine->at(3)  << endl << flush;
			outerMd = make_int(outerLine->at(colDist));
			curIndex = 0;
			in.incr(curIndex);
			while (curIndex < in.buf.size()) {   
				innerLine = &in.buf.at(curIndex);
				if ((innerLine->at(colChr) != outerLine->at(colChr)) || ((make_int(innerLine->at(colMain)) - make_int(outerLine->at(colMain))) > mean + baseLenFactor * sd)) 
					break;
				if (!in.visited.at(curIndex)) {
					innerMd = make_int(innerLine->at(colDist));
					int mainDiff = make_int(innerLine->at(colMain)) - make_int(outerLine->at(colMain));
					int secDiff = make_int(outerLine->at(colSec)) - make_int(innerLine->at(colSec));

					if (abs(mainDiff - secDiff) < mdJoinTolerance*sd) {
						in.visited.at(curIndex)  = true;
						curCluster.push_back(*innerLine);
					}
				}
				in.incr(curIndex); 
			}
			dump_cluster(curCluster, curClusterIndex, clusterType);
			curCluster.clear();
		}
		in.slideWin();

	} while (!in.buf.empty()); 



}

void cluster_indels01(int clusterType) {

	int curIndex;
	vector<string> * outerLine;
	vector<string> * innerLine;
	vector<vector<string> >  curCluster;  
	int outerMd, innerMd;
	int curClusterIndex = 1;
	assert(!in.buf.empty());

	int colMain;
	if (clusterType == 0) {
		colMain = colLeft;
	} else {
		colMain = colRight;
	}	

	do {
		if (!in.visited.front()) { //visited really means clustered
			outerLine = &in.buf.at(0);
			curCluster.push_back(*outerLine);
			//cout << "outer " << flush << outerLine->at(0) << " " << outerLine->at(1) << " " << outerLine->at(2) << " " << outerLine->at(3)  << endl << flush;
			outerMd = make_int(outerLine->at(colDist));
			curIndex = 0;
			in.incr(curIndex);
			while (curIndex < in.buf.size()) {   
				innerLine = &in.buf.at(curIndex);
				if ((innerLine->at(colChr) != outerLine->at(colChr)) || ((make_int(innerLine->at(colMain)) - make_int(outerLine->at(colMain))) > mean + baseLenFactor * sd)) 
					break;
				if (!in.visited.at(curIndex)) {
					innerMd = make_int(innerLine->at(colDist));
					if (abs(innerMd - outerMd) < mdJoinTolerance*sd) {
						in.visited.at(curIndex)  = true;
						curCluster.push_back(*innerLine);
					}
				}
				in.incr(curIndex); 
			}
			dump_cluster(curCluster, curClusterIndex, clusterType);
			curCluster.clear();
		}
		in.slideWin();

	} while (!in.buf.empty()); 

}


int main(int argc, char * argv[]) {

	in.init(&cin);
	int type;
	char ch;
	colID = -1;
	ifstream ifs;
	int option_index = 0; // getopt_long stores the option index here.

	static struct option long_options[] =
	{
		{"colDist", required_argument, 0, 0}, 
		{"colChr", required_argument, 0, 0},
		{"colLeft", required_argument, 0, 0},
		{"colRight", required_argument, 0, 0},
		{"colTemplate", required_argument, 0, 0}, 
		{"colID", required_argument, 0, 0},
		{"mean", required_argument, 0, 0},
		{"stdev", required_argument, 0, 0},
		{"type", required_argument, 0, 0},
		{"baseLenFactor", required_argument, 0, 0},
		{"mdJoinTolerance", required_argument, 0, 0},
		{"concise", required_argument, 0, 0},
		{0, 0, 0, 0}
	};


	while ((ch = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
		switch (ch) {
			case 0:
				//printf ("option %s", long_options[option_index].name);
				//if (optarg) printf (" with arg %s", optarg);
				if (strcmp(long_options[option_index].name, "colDist")==0) colDist = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "colChr")==0) colChr = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "colLeft")==0) colLeft = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "colRight")==0) colRight = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "colTemplate")==0) colTemplate = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "colID")==0) colID = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "mean")==0) mean = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "stdev")==0) sd = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "baseLenFactor")==0) baseLenFactor = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "mdJoinTolerance")==0) mdJoinTolerance = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "type")==0) type = atoi(optarg);
				else if (strcmp(long_options[option_index].name, "concise")==0) CONCISE = atoi(optarg);
				break;
		}
		option_index=0;
	}
	avg_md_cuttoff = 2 * mean;

	/*if (type < 2) 
		cerr << "PARAMS Things will be clustered with a leader if  their mds are within " << mdJoinTolerance << " * sd = " << mdJoinTolerance *sd << endl; 
	else 
		cerr << "PARAMS Things will be clustered with a leader if  the differences between their endpoints are within " << mdJoinTolerance << " * sd = " << mdJoinTolerance *sd << endl; 

	cout <<  "PARAMS The cluster will be extended at most mean + " << baseLenFactor << " * sd = " << mean + baseLenFactor * sd << " from the leader.\n";
	cout <<  "PARAMS All clusters with avg_md < " << avg_md_cuttoff << " are dropped.\n";
	cout << "CHR\tFROM\tTO\tTYPE\t#MP\tL_FROM\tL_TO\tDIST\tID\tEDGE\n";
*/ 

	if (in.buf.empty()) {
		cerr << "No matepairs to cluster.\n";
		return 0;
	}

	assert ((type >= 0) && (type <= 3));
	if (type == 0 || type == 1) {
		cluster_indels01(type);
	} else if (type == 2 || type == 3) {
		cluster_indels23(type);
	}
}

