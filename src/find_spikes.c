#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include<sys/stat.h>
#include <limits.h>

struct coverage_s {
	double normodds;
};

typedef struct coverage_s coverage_t;


#define LOW 0
#define HIGH 1
#define NORMAL 2

void usage() {
	fprintf(stderr,"usage: spikes <covmap> <lower> <upper>\n");
	exit(1);
}


void extractChr(char * tochr, char *fromstr) {
	char * t;
	char *pch=fromstr;
	while ((pch = strstr(pch+1,"chr"))!=NULL) {
		t=pch;
	}
	pch=t;

	int pos = 0;
	while (pch[pos] != '.') {
		tochr[pos] = pch[pos];
		pos++;
	}

	tochr[pos] = '\0';

}

//1-based, closed interval reporting
int main (int argc, char *argv[]) {

	if (argc != 1+3) { usage(); }

	int LOWERBOUND = atoi(argv[2]);
	int UPPERBOUND = atoi(argv[3]);

	char chr[6];
	extractChr(chr,argv[1]);

	//fprintf(stderr,"chromosome: %s\n",chr);

	// pointer to covmap file
	FILE *cmfp;

	// get length of contig
	struct stat st;
	stat(argv[1], &st);
	int num_recs = st.st_size / sizeof(coverage_t);

	// open the covmap
	if ( (cmfp = fopen(argv[1],"r")) == NULL ) {
		perror("fopen");
		exit(1);
	}

	// pointer to a struct (read line by line from coverage map)
	coverage_t * cov = (coverage_t*) malloc (sizeof(coverage_t));

	int last_state = NORMAL;
	int current_state;

	unsigned int pos = 0;

	//fprintf(stderr, "num_recs = %d\n", num_recs);
	while (!feof(cmfp)) {
		// read an element from the coverage map into cov
		fread (cov,sizeof(coverage_t),1,cmfp);
		pos++;
		current_state = NORMAL;

		if (cov->normodds <= LOWERBOUND) { current_state = LOW; }
		else if (cov->normodds >= UPPERBOUND) { current_state = HIGH; }
		//if (cov->normodds > 1000000) fprintf(stderr, "high normodds: %f\n", cov->normodds);
		//fprintf(stderr, "reg normodds: %f\n", cov->normodds);
		if (last_state != current_state) {
			if (last_state != NORMAL) {
				char *state;
				if ( last_state == LOW ) {
					state = "low";
				} else if (last_state == HIGH) {
					state = "high";
				}
				printf("%u \t %s\n",pos - 1 ,state);
			}
			if (current_state != NORMAL)  printf("%s \t %u \t ",chr, pos);
		}
		last_state = current_state;

	}

	if (last_state == LOW || last_state == HIGH) {
		char *state;
		if ( last_state == LOW ) {
			state = "low";
		} else if (last_state == HIGH) {
			state = "high";
		}
		printf("%u \t %s\n",pos ,state);
	}

	fclose(cmfp);
	//fprintf(stderr,"\nRead in a total of %u positions\n",pos);
	return 0;
}
