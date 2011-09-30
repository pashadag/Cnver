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



int main(int argc, char * argv[]) {

	string labels = "";
	if (argc == 2) labels = argv[1];
	cout << "digraph cs2_Graph {" << endl;
	vector<string> row;
	while (get_row_whitespace(cin, row)) {
		if (row[0] == "a" || row[0] == "convex_arc") {
			int from = atoi(row[1].c_str());
			int to   = atoi(row[2].c_str());
			if (from < 0 && to < 0) {
				swap(from, to);
			}
			cout << abs(from) << " -> " << abs(to);
			if (from < 0 && to > 0) {
				cout << " [dir=none] ";
			} else if (from > 0 && to < 0) {
				cout << " [dir=both] ";
			} else if (row[0] == "convex_arc") {
				cout << " [color=green] ";
			}
			if (labels == "label") {
				vector<string> rest(row.begin() + 3, row.end());
				cout << " [label=\"" << rest << "\"] ";
			}
			cout << ";" << endl;
		}
	}
	cout << "}" << endl;
	return 0;
}


