#include<cassert>
#include<cmath>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<cstdarg>
#include<cstdio>

#include "../include/defs.h"

using namespace std;

//takes a cs2 definition file and adds a header "p" line with the number of nodes and vertices
int main(int argc, char ** argv) {

	if (argc != 1) {
		fprintf(stderr, "usage: %s \n", argv[0]);
		return -1;
	}

	ostringstream lines;
	vector<string> row;
	string sbuf;
	int totArcs  = 0;
	int totNodes = 0;
	while(get_row_whitespace(cin, row, sbuf)){
		if (row[0] == "n") {
			totNodes++;
			lines << sbuf << endl;
		} else if (row[0] == "a") {
			totArcs++;
			lines << sbuf << endl;
		} else if (row[0] != "p") {
			lines << sbuf << endl;
		} 
	}

	cout << "p min " << totNodes << " " << totArcs << endl;
	cout << lines.str();
}
