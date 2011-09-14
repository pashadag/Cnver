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
#include<algorithm>

using namespace std;
ostringstream oMsg;
string sbuf;
#include "../include/defs.h"
#include "../include/interval.h"



string filename;
vector<Link> links; 
vector<Interval> calls;
int looseness = 100;
int maxdelsize = 10000;
double minDelPercentage = 0.2;
typedef pair<int, Link> Hit;
FILE *scov_file;
FILE *gc_file;

int getChange(const Interval & i) { //gets the +1 or -2 from the label of the interval, assuming its the first thing
	istringstream line(i.label);
	int change = 0;
	line >> change;
	return change;
}


double percDeletion(Interval query, Interval & biggestHit) {
	if (query.end - query.start <= 0) return false;

	vector<Interval>::iterator call = lower_bound(calls.begin(), calls.end(), query);

	
	int totalOverlap = 0;
	int biggestOverlap = 0;
	if (call != calls.begin() && (call - 1)->overlaps(query)) call--;
	while (call != calls.end() && !call->fullyRightOf(query)) {
		if (getChange(*call) < 0) {
			if (call->amountThatOverlaps(query) > biggestOverlap) biggestHit = *call;
			totalOverlap += call->amountThatOverlaps(query);
		}
		call++;
	}

	return (totalOverlap / (query.end - query.start));

}

double get_doc_ratio(Interval i) {
	long range = i.end - i.start + 1;
	if (range <= 0) return false;

	double * fs= (double*)malloc(sizeof(double) * range);
	if (fs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}
	double * gs= (double*)malloc(sizeof(double) * range);
	if (gs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}


	long go_to = (i.start - 1) * sizeof(double);
	int res = fseek(scov_file, go_to, SEEK_SET);
	if (res != 0) {
		cerr << "error at interval " << i << endl;
		exit(1);
	}
	res = fseek(gc_file, go_to, SEEK_SET);
	if (res != 0) {
		cerr << "error at interval " << i << endl;
		exit(1);
	}
	res = fread(fs, sizeof(double) * range, 1, scov_file);
	if (res != 1) {
		cerr << "error at interval " << i << endl;
		exit(1);
	}
	res = fread(gs, sizeof(double) * range, 1, gc_file);
	if (res != 1) {
		cerr << "error at interval " << i << endl;
		exit(1);
	}
	double sum_masked = 0.0;
	double sum_expected = 0.0;
	int index = 0;
	while (index < range) {
		//printf("%d\t%.4f\t%.4f\n",index+start,fs[index],gs[index]);
		if (gs[index] > 0.0) {
			sum_masked += fs[index];
			sum_expected += gs[index];
		} 
		index++;
	}
	double doc_ratio = (double) sum_masked / (double) sum_expected;
	delete fs;
	delete gs;
	return doc_ratio;


}
vector<Hit> rangeSearch(Interval range) {
	vector<Hit> hits;
	for (int i = 0; i < links.size(); i++) {
		if (links[i].from <= range.end && links[i].from >= range.start) {
			hits.push_back(make_pair(links[i].to, links[i]));
		}
		//note that we might add two hits per single link, but that's good
		if (links[i].to <= range.end && links[i].to >= range.start) {
			hits.push_back(make_pair(links[i].from, links[i]));
		}
	}
	return hits;
}


void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " work_dir ref_name [-u used_link_file -c calls_file -l looseness -d maxdelsize -p minDelPercentage]" << endl;
	cerr << "\t-l: when searching for a link at an endpoint of a gain, allow this much +/- bp on either side.\n";
	cerr << "\t-d: only consider pairs of links that lead to a range of this much.\n";
	cerr << "\t-p: only report pairs of links where there is at least this much percentage of the range being in a deletion.\n";

	exit(1);
}


int main(int argc, char * argv[]) {


	if (argc == 1) usage(argc, argv);
	char ch;
	string  links_filename, calls_filename;
	while ((ch = getopt(argc, argv, "u:c:l:d:p:")) != -1) {
		switch (ch) {
			case 'u':
				links_filename = optarg;
				break;
			case 'c':
				calls_filename = optarg;
				break;
			case 'l':
				looseness = atoi(optarg);
				cerr << argv[0] << ": +/- at gain endpoints is " << looseness << ".\n";
				break;
			case 'd':
				maxdelsize = atoi(optarg);
				cerr << argv[0] << ": max range for link destinations is " << maxdelsize << ".\n";
				break;
			case 'p':
				minDelPercentage = atof(optarg);
				cerr << argv[0] << ": min percentage of range in deletion is " << minDelPercentage << ".\n";
				break;

		}
	}

	if ( optind != argc - 2 ) { //should be two other parameter
		cerr << "Invalid parameters.\n";
		usage(argc, argv);
	}

	string base_dir = argv[optind];
	string ref_name = argv[optind+1];
	// cerr << "Using base_dir " << base_dir << endl; cerr << "Using ref_name " << ref_name << endl;
	if (links_filename == "") links_filename = base_dir + "/" + ref_name + "/" + ref_name + ".used_dgs";
	if (calls_filename == "") calls_filename = base_dir + "/" + ref_name + "/" + ref_name + ".cnvs.raw.annot";

	ifstream ifile;
	int num_loaded;
	open_file(ifile, links_filename);
	num_loaded = read_links(ifile, links);
	cerr << argv[0] << ": loaded " << num_loaded << " used links.\n";
	ifile.close();
	open_file(ifile, calls_filename);
	num_loaded = read_intervals(ifile, calls);
	sort(calls.begin(), calls.end());
	cerr << argv[0] << ": loaded " << num_loaded << " cnv calls.\n";
	ifile.close();


	filename = base_dir + "/" + ref_name + "/" + ref_name + ".scov";
	scov_file = fopen(filename.c_str(),"rb");
	if (scov_file == NULL){
		printf("Error opening scov file\n");
		exit(1);
	}
	filename = base_dir + "/" + ref_name + "/" + ref_name + ".gc";
	gc_file = fopen(filename.c_str(),"rb");
	if (gc_file == NULL){
		printf("Error opening gc file\n");
		exit(1);
	}




	for (int i = 0; i < calls.size(); i++) {
		int change = getChange(calls[i]);
		if (change <= 0) continue;
		vector<Hit> lhits, rhits;
		lhits = rangeSearch(Interval(calls[i].chr, calls[i].start - looseness, calls[i].start + looseness));
		rhits = rangeSearch(Interval(calls[i].chr, calls[i].end   - looseness, calls[i].end   + looseness));
		for (int j = 0; j < lhits.size(); j++) {
			int pos1  = lhits[j].first;
			Link link1 = lhits[j].second;
			for (int k = 0; k < rhits.size(); k++) {
				int pos2  = rhits[k].first;
				Link link2 = rhits[k].second;
				if (pos1 > pos2) swap(pos1, pos2);
				Interval query(link1.chr, pos1, pos2);
				Interval bestDel;
				double perc = percDeletion(query, bestDel);
				if (query.end - query.start < maxdelsize && perc >= minDelPercentage) {
					/*cout << calls[i] << "\tInsertion" << endl;
					  cout << query << "\t" << "(size " << query.end - query.start << ")" << "\tLinked_range" << endl;
					  cout << link1 << "\tLeft_link" << endl;
					  cout << link2 << "\tRight_link" << endl;
					  if (perc > 0) {
					  cout << bestDel << "\tBest_deletion" << endl;
					  } else {
					  cout << "none\tnone\tnone\tBest_deletion" << endl;
					  }
					  cout << perc << "\tPercent_del" << endl;
					  cout << endl;
					 */

					cout << "Insertion\t" << calls[i] << "\tSize = " << calls[i].end - calls[i].start << endl;
					cout << "Linked_range\t" << query << "\tSize = " << query.end - query.start << "\tPerc = " << perc << "\tDoc = " << get_doc_ratio(query)  << endl;
					cout << "Left_link\t" << link1 << endl;
					cout << "Right_link\t" << link2 << endl;
					if (perc > 0) {
						cout << "Best_deletion\t" << bestDel << endl;
					} else {
						cout << "Best_deletion" << endl;
					}

					cout << endl;

				}
			}
		}
	}

	return 0;
}


