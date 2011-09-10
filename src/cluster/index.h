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

#include "../include/defs.h"

typedef unsigned long int pType;

class indexWriterClass {  
	public:
		void init (string indexFilename, ofstream * _dataFile) {
			indexFile.open(indexFilename.c_str(), ios::out | ios::binary);
			dataFile = dataFile;
		}

		void init(string indexFilename ) {
			indexFile.open(indexFilename.c_str(), ios::out | ios::binary);
		}

		void mark () {
			pType pos = dataFile->tellp(); 
			indexFile.write((char*) &pos, sizeof(pType));
		}

		void mark (pType pos) { 
			indexFile.write((char*) &pos, sizeof(pType));
		}

	private:

		ofstream indexFile;
		ofstream * dataFile;
};

long binary_search_via_index(ifstream &indexFile, const pType & val, pType first, pType last, pType & nextOne) {
	//cerr << "calling binary search with " << first << " " << last << endl;
	pType rid, data;
	if (last < first) {
		return -1;
	} else if (first == last) {
		indexFile.seekg(2 * sizeof(pType) * first);
		indexFile.read((char*) &rid, sizeof(pType));
		if (rid == val) {
			indexFile.read((char*) &data, sizeof(pType)); 
			indexFile.read((char*) &nextOne, sizeof(pType));  //dummy
			indexFile.read((char*) &nextOne, sizeof(pType)); 
			return data;
		} else {
			return -1;
		}
	} else {
		int mid = int ((last + first) / 2);
		//read the mid val;
		indexFile.seekg(2 * sizeof(pType) * mid);
		indexFile.read((char*) &rid, sizeof(pType));
		//cerr << "\tread mid = " << mid << " and got read id " << rid << endl;
		if (rid == val) {
			indexFile.read((char*) &data, sizeof(pType)); 
			indexFile.read((char*) &nextOne, sizeof(pType));  //dummy
			indexFile.read((char*) &nextOne, sizeof(pType)); 
			return data;
		} else if (rid > val) {
			return binary_search_via_index(indexFile, val, first, mid - 1, nextOne);
		} else {
			return binary_search_via_index(indexFile, val, mid + 1, last, nextOne);
		}
	}
}

string binary_search_lookup(char * indexFilename, char * dataFilename, pType val) {
	ifstream indexFile(indexFilename, ios::in | ios::binary);
	ifstream dataFile(dataFilename);


	//get file length 
	indexFile.seekg (0, ios::end);
	pType numberOfRecords = ((long) indexFile.tellg() + 1) / (2*sizeof(pType));
	pType nextOne;
	long loc = binary_search_via_index(indexFile, val, 0, numberOfRecords - 1, nextOne);
	if (loc == -1) {
		cerr << "Cannot find val=" << val << " in index file.\n";
		exit(1);
	} else {
		dataFile.seekg(loc);
		string sbuf;
		ostringstream o; 
		while (dataFile.tellg() < nextOne) {
			getline(dataFile, sbuf);
			o << sbuf << endl;
		}
		return o.str();
	}
}

template<class T>
pType rewind(ifstream &dataFile, pType index, const T & val) {
	T cur;
	dataFile.seekg(sizeof(val) * index);
	dataFile.read((char*) &cur, sizeof(cur));
	assert(cur == val);
	while ( index > 0 ) {
		dataFile.seekg(sizeof(val) * (index - 1));
		dataFile.read((char*) &cur, sizeof(cur));
		if (!(cur == val)) return index;
		index --;
	}
	return index;
}




template<class T>
long binary_search_without_index(ifstream &dataFile, const T & val, int64_t first, int64_t last) {
	//cout << "calling binary search with " << first << " " << last << endl;
	T cur;
	if (last < first) {
		return -1;
	} else if (first == last) {
		dataFile.seekg(sizeof(val) * first);
		dataFile.read((char*) &cur, sizeof(cur));
		if (cur == val) {
			return sizeof(val) * rewind(dataFile, first, val);
		} else {
			return -1;
		}
	} else {
		int64_t mid = int ((last + first) / 2);
		//read the mid val;
		dataFile.seekg(sizeof(val) * mid);
		dataFile.read((char*) &cur, sizeof(cur));
		//cout << "\tread mid = " << mid << " and got read " << cur << endl;
		if (cur == val) {
			//goto rewind
			return sizeof(val) * rewind(dataFile, mid, val);
		} else if (cur < val) {
			return binary_search_without_index(dataFile, val, mid + 1, last);
		} else {
			return binary_search_without_index(dataFile, val, first, mid - 1);
		}
	}
}


template<class T>
string binary_search_lookup_in_binary_file(char * dataFilename, T val) { 
	ifstream dataFile;
	open_file_binary(dataFile, dataFilename);

	//get dataFile length 
	dataFile.seekg (0, ios::end);
	pType filesize =  dataFile.tellg();
	pType numberOfRecords = ((long) dataFile.tellg() + 1) / sizeof(T);
	pType nextOne;
	int64_t loc = binary_search_without_index(dataFile, val, 0, numberOfRecords - 1);
	if (loc == -1) {
		cerr << "Cannot find val=" << val << " in index file.\n";
		exit(1);
	} else {
		ostringstream o; 
		while (loc < filesize) {
			T cur;
			dataFile.seekg(loc);
			dataFile.read((char *) &cur, sizeof(T));
			//cout << "Processing " << cur << endl;
			if (!(cur == val)) break;
			o << cur << endl;
			loc += sizeof(T);
		}
		return o.str();
	}
}


string lookup(char * indexFilename, char * dataFilename, int val) {
	//cerr << "In lookup.\n";
	ifstream indexFile(indexFilename, ios::in | ios::binary);
	ifstream dataFile(dataFilename);
	pType pVal1, pVal2;

	//get file length for error checking	
	indexFile.seekg (0, ios::end);
	int indexLen = indexFile.tellg();
	if (val >= indexLen) {
		cerr << "Out of bounds error in lookup(): seeking " << val << ", file size is " << indexLen << endl;
		exit(1);
	}

	indexFile.seekg(sizeof(pType) * val);
	indexFile.read((char*) &pVal1, sizeof(pType));
	indexFile.read((char*) &pVal2, sizeof(pType)); // this will crash if its the last value, but its a degenerate case and is ignored for now.

	dataFile.seekg(pVal1);
	//cout << "pVal1\t" << pVal1 <<  "pVal2\t" << pVal2 << endl;

	string sbuf;
	ostringstream o; 

	while (dataFile.tellg() < pVal2) {
		getline(dataFile, sbuf);
		o << sbuf << endl;
	}
	return o.str();
}


