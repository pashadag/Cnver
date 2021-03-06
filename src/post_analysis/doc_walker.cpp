#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include "../include/interval.h"

string sbuf;
FILE *f;
FILE *g;

void process_interval(Interval i) {
	double * fs= (double*)malloc(sizeof(double) * (i.end - i.start + 1));
	if (fs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}
	double * gs= (double*)malloc(sizeof(double) * (i.end - i.start + 1));
	if (gs==NULL) {
		printf("Cannot malloc mem\n");
		exit(1);
	}


	long go_to = (i.start - 1) * sizeof(double);
	int res = fseek(f, go_to, SEEK_SET);
	assert(res == 0);
	res = fseek(g, go_to, SEEK_SET);
	assert(res == 0);
	res = fread(fs, sizeof(double) * (i.end - i.start + 1), 1, f);
	assert(res == 1);
	res = fread(gs, sizeof(double) * (i.end - i.start + 1), 1, g);
	assert(res == 1);
	double sum_unmasked = 0.0;
	double sum_masked = 0.0;
	double sum_expected = 0.0;
	int masked = 0;
	int index = 0;
	while (index < i.end - i.start + 1) {
		//printf("%d\t%.4f\t%.4f\n",index+start,fs[index],gs[index]);
		sum_unmasked += fs[index];
		if (gs[index] > 0.0) {
			sum_masked += fs[index];
			sum_expected += gs[index];
		} else {
			masked++;
		}
		index++;
	}
	double doc_ratio = (double) sum_masked / (double) sum_expected;
	cout << i << "\t";
	printf ("DOC_ratio:\t%.4f\tObserved_coverage:\t%.4f\tExpected_coverage:\t%.4f\tNum_masked_positions:\t%d\tObserved_coverage_including_masked_regions:\t%.4f\n", doc_ratio, sum_masked, sum_expected, masked, sum_unmasked);
	delete fs;
	delete gs;

}

int main(int argc, char** argv) {
	unsigned int start = -1;
	unsigned int end = -1;
	if (argc == 5) {
		start = atoi(argv[3]);
		end   = atoi(argv[4]);
	} else if (argc != 3) {
		printf("%s scov_file gc_file start end\n", argv[0]);
		printf("The .scov and .gc files can be found in the work_dir.\n");
		printf("If start/end is not specified, an interval file is taken as an input.\n");
		exit(1);
	}

	f = fopen(argv[1],"rb");
	if (f==NULL){
		printf("Error opening file\n");
		exit(1);
	}
	g = fopen(argv[2],"rb");
	if (g==NULL){
		printf("Error opening file\n");
		exit(1);
	}

	if (start == -1) {
		while (getline(cin, sbuf)) {
			Interval i = read_interval(sbuf);
			process_interval(i);
		}
	} else {
		process_interval(Interval("-1", start, end));
	}
	fclose(f);
	fclose(g);
	return 0;
}
