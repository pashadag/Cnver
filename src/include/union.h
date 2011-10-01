#ifndef CPCOUNT_UNION_H
#define CPCOUNT_UNION_H

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


bool pred_empty_inv(const vector< pair< int, bool> > & v) {return v.size() == 0;}
bool pred_empty(const vector< int>  & v) {return v.size() == 0;}

//Straight out of Cormen
class UnionFindClass {
	public:
		UnionFindClass(int size) : data(size, -1), rank(size, 0) {}

		void unionn( const int & x, const int & y) {
			int fsx = find_set(x);
			int fsy = find_set(y);
			link(fsx,fsy);
		}

		int find_set(const int & x) {
			if (data.at(x) == -1) {
				data.at(x) = x;
			} else if (data.at(x) != x) {
				data[x] = find_set(data[x]);
			}
			return data.at(x);
		}

	private:
		void link (const int & x, const int & y) {
			if (rank.at(x) > rank.at(y)) {
				data.at(y) = x;
			} else {
				data.at(x) = y;
				if (rank.at(x) == rank.at(y)) rank.at(y)++;
			}
		}

		vector<int> data;
		vector<int> rank;


	public:

		int num_classes() {
			int count = 0;
			for (int i = 0; i < data.size(); i++) {
				if (data[i] == i) count++;
			}
			return count; //return *max_element(data.begin(), data.end());
		}

		int num_elements() {
			return count_if(data.begin(), data.end(), bind2nd(not_equal_to<int>(), -1));
		}

		int size() {
			return data.size();
		}

		int operator[](int i) {
			return data.at(i);
		}

		string dump() {
			ostringstream o;
			int i;
			for (i = 0; i < data.size(); i++) {
				o << i << "\t" << data[i] << "\t" << rank[i] <<  endl;
			}
			return o.str();
		}

		void get_classes( vector<vector<int> > & otherWay) {
			otherWay.resize(data.size());
			for (int i = 0; i < data.size(); i++) {
				if (data[i] != -1) {
					find_set(i);
					otherWay.at(data.at(i)).push_back(i);
				}
			}
			otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), pred_empty), otherWay.end());
		}


		void get_classes (vector<int> & consolidatedData, vector<vector< int > > & otherWay) {
			get_classes(otherWay);
			consolidatedData.resize(data.size(), -1);
			for (int i = 0; i < otherWay.size(); i++) 
				for (int j = 0; j < otherWay[i].size(); j++) 
					consolidatedData.at(otherWay[i][j]) = i;		
		}


};


//Straight out of Cormen
class UnionFindInvClass {
public:
	UnionFindInvClass(int size) : data(size, -1), dataRev(size,false), rank(size, 0) {}

	void unionn( const int & x, const int & y, bool rev) {
		int fsx = find_set(x);
		int fsy = find_set(y);
		rev = rev ^ dataRev.at(x) ^ dataRev.at(y);
		link(fsx,fsy, rev);
	}

	int find_set(const int & x) {
		if (data.at(x) == -1) {
			data.at(x) = x;
			dataRev.at(x) = false;
		} else if (data.at(x) != x) {
			dataRev[x] = dataRev[x] ^ dataRev[data[x]];
			data[x] = find_set(data[x]);
		}
		return data.at(x);
	}

private:
	void link (const int & x, const int & y, bool rev) {
		if (rank.at(x) > rank.at(y)) {
			data.at(y) = x;
			dataRev.at(y) = rev;
		} else {
			data.at(x) = y;
			dataRev.at(x) = rev;
			if (rank.at(x) == rank.at(y)) rank.at(y)++;
		}
	}

	vector<int> data;
	vector<bool> dataRev;
	vector<int> rank;


public:

	int num_classes() {
		int count = 0;
		for (int i = 0; i < data.size(); i++) {
			if (data[i] == i) count++;
		}
		return count; //return *max_element(data.begin(), data.end());
	}

	int num_elements() {
		return count_if(data.begin(), data.end(), bind2nd(not_equal_to<int>(), -1));
	}

	int size() {
		return data.size();
	}

	int operator[](int i) {
		return data.at(i);
	}

	string dump() {
		ostringstream o;
		int i;
		for (i = 0; i < data.size(); i++) {
			o << i << "\t" << data[i] << "\t" << rank[i] <<  endl;
		}
		return o.str();
	}

	void get_classes( vector<vector< pair<int, bool> > > & otherWay) {
		otherWay.resize(data.size());
		for (int i = 0; i < data.size(); i++) {
			if (data[i] != -1) {
				find_set(i);
				otherWay.at(data.at(i)).push_back(make_pair(i,dataRev.at(i)));
			}
		}
		otherWay.erase(remove_if(otherWay.begin(), otherWay.end(), pred_empty_inv), otherWay.end());
	}


	void get_classes (vector<int> & consolidatedData, vector<vector< pair<int, bool> > > & otherWay) {
		get_classes(otherWay);
		consolidatedData.resize(data.size(), -1);
		for (int i = 0; i < otherWay.size(); i++) 
			for (int j = 0; j < otherWay[i].size(); j++) 
				consolidatedData.at(otherWay[i][j].first) = i;		
	}


};


#endif 
