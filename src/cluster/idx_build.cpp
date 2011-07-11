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

#include "../include/dbtypes.h"
#include "../include/defs.h"

typedef uint64_t pType;

pType get_bowtie_rid(const string & buf) {
	// example: 3+:<24,11910070,2>
	istringstream line(buf);
	int i;
	char c;
	for (i = 0; i < buf.length(); i++) {
		line.get(c);
		if ((c == '+')  || (c == '-')) {
			return atol(buf.substr(0, i).c_str());
		}
	}
	cerr << "Encountered bad readname: " << buf << endl;
	return 0;

}


int main(int argc, char * argv[]) {
	if (argc < 3) {
		cerr << "Usage: idx_build <index_filename> <data_filename> {<mode> {<delimeter>}}.\n";
		exit(1);
	}

	string mode = "basic";
	string delim = "";

	if (argc > 3) {
		mode = argv[3];
		if (mode == "delim") delim = argv[4];
	}

	//open data file
	ofstream indexFile(argv[1], ios::out | ios::binary);

	string sbuf;
	if (mode == "basic") {
		ifstream dataFile(argv[2]);
		pType pos = dataFile.tellg(); 
		while (getline(dataFile, sbuf)) {
			indexFile.write((char*) &pos, sizeof(pType));
			pos = dataFile.tellg(); 
		}
	} else if (mode == "delim") {
		ifstream dataFile(argv[2]);
		pType pos = dataFile.tellg(); 
		while (dataFile >> sbuf) {
			if (sbuf == delim) {
				indexFile.write((char*) &pos, sizeof(pType));
			}
			pos = dataFile.tellg();
		}
	} else if (mode == "bowtie_concise") {
		ifstream dataFile(argv[2]);
		pType pos = dataFile.tellg(); 
		pType last = 948238;
		while (getline(dataFile,sbuf)) {
			pType rid = get_bowtie_rid(sbuf);
			if (rid != last) {
				indexFile.write((char*) &rid, sizeof(pType));
				indexFile.write((char*) &pos, sizeof(pType));
				last = rid;
			}
			pos = dataFile.tellg();
		}
	} else if (mode == "mmap") {
		ifstream dataFile(argv[2]);
		pType pos = dataFile.tellg(); 
		pType last_id = 948238;
		Mmap_t cur;
		while (dataFile >> cur) {
			if (cur.read_id != last_id) {
				indexFile.write((char*) &cur.read_id, sizeof(pType));
				indexFile.write((char*) &pos, sizeof(pType));
				last_id = cur.read_id;
			}
			pos = dataFile.tellg();
		}
	} else if (mode == "rmap") {
		ifstream dataFile;
		open_file_binary(dataFile, argv[2]);
		pType pos = dataFile.tellg(); 
		pType last_id = 948238;
		Map_t cur;
		while (dataFile.read((char *) &cur, sizeof(Map_t))) { 
			//cout << "Read in " << cur << endl;
			if (cur.read_id != last_id) {
				//cout << "Processing...\n";
				indexFile.write((char*) &cur.read_id, sizeof(pType));
				indexFile.write((char*) &pos, sizeof(pType));
				last_id = cur.read_id;
			}
			pos = dataFile.tellg();
		}
	} else {
		cerr << "Unknown mode: " << mode << endl;
	}


}


