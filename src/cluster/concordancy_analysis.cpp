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

using namespace std;
ostringstream oMsg;
string sbuf;
string filename, baseFilename;
#include "../include/defs.h"
#include "../include/dbtypes.h"
#include "../include/InputReader.h"

typedef	enum dirEnum { FOR = 0, BACK = 1};
typedef list<Map_t> mapsType;
bool CONCISE_OUTPUT = false;
bool PREPEND_OUTPUT = false;
int upperd = -2;
int lowerd = -1;
int matepair_mappings_counter = 0;
int maxDist = 1000000000; //1G
int readLen= -1;
vector< ofstream * > outf; 
vector< ostringstream * > tmpStm;
typedef map<string,int> ref_map_t;
ref_map_t ref_names_map;


	/*typedef struct {
		uint64_t    read_id;
		char        orientation;
		char        ref_name[CONTIGNAME_LEN_SMALL];
		uint64_t    ref_pos;
		double      normodds;
	} Map_t;
	*/

void parse_map(Map_t & mapt, uint64_t & matepair, dirEnum &type) {
	int parity = mapt.read_id % 2;
	matepair = mapt.read_id - parity;
	if (parity== 0) 
		type = FOR;
	else 
		type = BACK;
	return;
}

int get_ref_index(string ref_name) {
	ref_map_t::iterator it = ref_names_map.find(ref_name); //Identify reference index
	if (it != ref_names_map.end()) {
		return it->second; 
	} 
	return -1;
}

int mapped_dist(const Map_t & fmap, const Map_t & bmap) {
	const Map_t * lmap;
	const Map_t * rmap;
	if (fmap.ref_pos < bmap.ref_pos) {
		lmap = &fmap;
		rmap = &bmap;
	} else {
		lmap = &bmap;
		rmap = &fmap;
	}

	if ((lmap->orientation == '+') && (rmap->orientation == '-')) {
		return rmap->ref_pos - lmap->ref_pos + readLen;
	} else if ((lmap->orientation == '-') && (rmap->orientation == '+')) {
		return rmap->ref_pos - lmap->ref_pos +  3*readLen;
	} else if ((lmap->orientation == '+') && (rmap->orientation == '+')) {
		return rmap->ref_pos - lmap->ref_pos +  2*readLen;
	} else if ((lmap->orientation == '-') && (rmap->orientation == '-')) {
		return rmap->ref_pos - lmap->ref_pos +  2*readLen;
	}

	return 0; //should never get here
}


bool is_concordant(const Map_t & fmap, const Map_t & bmap, int & dist) {
	dist =  mapped_dist(fmap, bmap);
	bool correct_orientation = false;
	correct_orientation = 
		(fmap.orientation  == '+') && (bmap.orientation == '-') && (fmap.ref_pos < bmap.ref_pos) || 
		(fmap.orientation  == '-') && (bmap.orientation == '+') && (fmap.ref_pos > bmap.ref_pos);
	if (correct_orientation && (dist >= lowerd) && (dist <= upperd)) return true;
	return false;
}

void process_matepair_disc(mapsType & fmaps, mapsType & bmaps, const int & totalMaps) {
	list<Map_t>::const_iterator itf, itb;
	int dist, chr;
	int counterRecall = matepair_mappings_counter;
	bool foundDiscordantMapping = false;
	uint64_t read_ids[2];

	for (int i = 0; i < tmpStm.size(); i++) tmpStm[i]->str("");
	for (itf = fmaps.begin(); itf != fmaps.end(); itf++) {
		for (itb = bmaps.begin(); itb != bmaps.end(); itb++) {
			if (is_concordant(*itf, *itb, dist)) {
				matepair_mappings_counter = counterRecall;
				return;
			} 
			if ((0 == strcmp(itf->ref_name,itb->ref_name)) && (dist <= maxDist)) {
				if (!CONCISE_OUTPUT) { 
					const Map_t * lmap;
					const Map_t * rmap;
					int type = 0; // 0 not necessary, its just to get rid of warning
					if (itf->ref_pos < itb->ref_pos) {
						lmap = &(*itf);
						rmap = &(*itb);
					} else {
						lmap = &(*itb);
						rmap = &(*itf);
					}
					if ((lmap->orientation == '+') && (rmap->orientation == '-')) {
						type = 0;
					} else if ((lmap->orientation == '-') && (rmap->orientation == '+')) {
						type = 1;
					} else if ((lmap->orientation == '+') && (rmap->orientation == '+')) {
						type = 2;
					} else if ((lmap->orientation == '-') && (rmap->orientation == '-')) {
						type = 3;
					}

					int file_index = get_ref_index(itf->ref_name);
					if (PREPEND_OUTPUT) *(tmpStm[file_index]) << ++matepair_mappings_counter << "\t";
					Mmap_t mmap;
					mmap.dist = dist;
					strcpy( mmap.ref_name, itf->ref_name);
					mmap.left_pos = lmap->ref_pos;
					mmap.right_pos = rmap->ref_pos;
					mmap.read_id = itf->read_id;
					mmap.type = type;
					mmap.l_orientation = lmap->orientation;
					mmap.r_orientation = rmap->orientation;

					*(tmpStm[file_index]) << mmap << endl;

					/*
					*(tmpStm[file_index]) << dist << "\t" << itf->ref_name << "\t";
					*(tmpStm[file_index]) << lmap->ref_pos << "\t" << rmap->ref_pos << "\t";
					*(tmpStm[file_index]) << itf->read_id << "\t" << type << "\t" << lmap->orientation << "\t" << rmap->orientation;
					*(tmpStm[file_index]) << endl;
					*/
					//[id] dist, chr, left, right, template, type, left-orientation, right-orientation
				} else {
					foundDiscordantMapping = true; //this is for the concise output
					read_ids[0] = itf->read_id;
					read_ids[1] = itb->read_id;
				}
			} 
		}
	}

	//Write output
	if (!CONCISE_OUTPUT) {
		//cout << "chr = " << chr << ", and length tmpStm = " << tmpStm.str().size() << endl;
		for (int i = 0; i < tmpStm.size(); i++) {
			*(outf[i]) << tmpStm[i]->str();
			//outf[i]->flush();
		}
	} else {
		//just write to cout the template id
		if (foundDiscordantMapping) {
			cout << read_ids[0] << endl << read_ids[1] << endl;
		}
	}

	return;
}


void usage(int argc, char * argv[]) {
	cout << "Usage: " << argv[0] << " <map_file> <> {<optparam3>} " << endl;
	exit(1);
	/*cout << "\tmap_list     : file with names of mapping files.\n";
	  cout << "\tref_names    : file with names of reference sequences.\n";
	  cout << "\tinput_format : format of mapping files, either \"bam\" or \"text\".\n";
	 */
}

int main(int argc, char * argv[]) {
	//if (argc != 1+0) usage(argc, argv);

	vector<string> map_list;
	vector<string> ref_names;
	string input_format;
	ifstream inf;

	char ch;

	while ((ch = getopt(argc, argv, "pcf:b:u:l:t:d:r:m:n:i:")) != -1) {
		switch (ch) {
			case 'p':
				PREPEND_OUTPUT = true;
				break;
			case 'c':
				CONCISE_OUTPUT = true;
				break;
			case 'u':
				upperd = atoi(optarg);
				break;
			case 'l':
				lowerd = atoi(optarg);
				break;
			case 'd':
				maxDist = atoi(optarg);
				break;
			case 'r':
				readLen = atoi(optarg);
				break;
			case 'm':
				filename = optarg;
				read_file_into_vector(filename, map_list);
				break;
			case 'n':
				filename = optarg;
				read_file_into_vector(filename, ref_names);
				break;
			case 'i':
				input_format = optarg;
				break;
		}
	}

	//sanity checks
	if (!((upperd >= lowerd) && (lowerd >= 0) && (readLen > 0) && (optind == argc ) )) {//should be no other parameter
		cerr << "Invalid options.\n";
		usage(argc, argv);
	}


	//open outputfiles and initialize ref_names_map
	outf.resize(ref_names.size());
	tmpStm.resize(ref_names.size());
	for (int i = 0; i < ref_names.size(); i++)  {
		ref_names_map.insert(make_pair(ref_names[i], i));
		outf[i] = new ofstream;
		open_file(*outf[i], ref_names[i]  + ".mmap");
		tmpStm[i] = new ostringstream;
	}

	//main loop
	mapsType fmaps;
	mapsType bmaps;
	uint64_t lastReadId = -1;
	for (int map_file_index = 0; map_file_index < map_list.size(); map_file_index++) {
		InputReader reader(input_format, map_list[map_file_index], lastReadId + 1);
		Map_t curMap;
		uint64_t lastmate = -1;
		fmaps.clear();
		bmaps.clear();
		int totalMaps = 0;
		while (reader.getNext(curMap)) {
			dirEnum dir;
			uint64_t mate;
			parse_map(curMap, mate, dir);
			if (mate != lastmate) { //finished reading in all mappings of current matepair, so process them.
				//cout << "Processing fmaps: " << endl; for (mapsType::iterator itf = fmaps.begin(); itf != fmaps.end(); itf++) cout << *itf << endl; cout << "Processing bmaps: " << endl; for (mapsType::iterator itf = bmaps.begin(); itf != bmaps.end(); itf++) cout << *itf << endl; 
				process_matepair_disc(fmaps, bmaps, totalMaps);
				fmaps.clear();
				bmaps.clear();
				totalMaps = 0;
			}
			if (get_ref_index(curMap.ref_name) != -1) {
				totalMaps++;
				if (dir == FOR) fmaps.push_back(curMap);
				else if (dir == BACK) bmaps.push_back(curMap);
			}
			lastmate = mate;
			lastReadId = curMap.read_id;
		}
	}




	//close files
	for (int i = 0; i < outf.size(); i++) {
		outf[i]->close();
		delete outf[i];
		delete tmpStm[i];
	}
	return 0;
}

