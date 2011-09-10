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
#include<cstdarg>
#include<algorithm>
#include<sys/stat.h>
#include<sys/types.h>
#include<fcntl.h>

using namespace std;
#include "include/defs.h"
#include "include/dbtypes.h"


long offset = 0;
string in_filename;

void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " binfile type(rmap/int/double/bamnames/char) offset (default is 0, can be negative)\n";
	exit(1);
}


template<class T>
void process_file(T data) {
	//struct stat statinfo; stat(in_filename.c_str(), &statinfo); long filesize = statinfo.st_size; assert (filesize % sizeof(data) == 0); cerr << "Filesize is " << filesize << endl; long seekloc;
	ifstream inf;
	open_file_binary(inf, in_filename);
	if (offset >= 0) {
		inf.seekg(offset * sizeof(data), ios_base::beg);
	} else {
		inf.seekg(offset * sizeof(data), ios_base::end);
	}

	while (inf.read((char*) & data, sizeof(data))) {
		cout << data << endl;
	}
}

struct BamReadName {
	char name[MAX_BAM_READNAME];
};

ostream & operator<< (ostream & out, const BamReadName & bamreadname) {
    out << bamreadname.name;
    return out;
}


int main(int argc, char * argv[]) {

	if (argc < 3) usage(argc, argv);
	in_filename = argv[1];
	string datatype = argv[2];
	if (argc > 3) {
		offset = atol(argv[3]);
	}
	if (datatype == "rmap") {
		Map_t data;
		process_file(data);
	} else if (datatype == "int") {
		uint64_t data;
		process_file(data);
	} else if (datatype == "double") {
		double data;
		process_file(data);
	} else if (datatype == "bamnames") {
		BamReadName data;
		process_file(data);
	} else if (datatype == "char") {
		char data;
		process_file(data);
	}
	return 0 ;
}

