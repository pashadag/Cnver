#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "include/dbtypes.h"
#include "include/defs.h"




int window_size=35;
int num_bins=35;


int whats_my_bin(int gc_count,int at_count) {
	assert(window_size%num_bins==0);
	if (gc_count==0 && at_count==0) {
		return -1;
	}
	int bin = ((gc_count*100)/(at_count+gc_count))/2;
	if (bin==num_bins) {
		bin--;
	}
	if (bin>=num_bins || bin<0) {
		fprintf(stderr,"Window error!");
	}
	assert(bin>=0 && bin<num_bins);
	return bin;
}


void read_in_masks(string repeat_filename, short * repeats_mask) {
	ifstream inf;
	open_file(inf, repeat_filename);

	vector<string> line;
	while (get_row_whitespace(inf, line)) {
		int from = atoi(line[1].c_str());
		int to   = atoi(line[2].c_str());
		assert (from >= 0);
		assert (from <= to);
		while (from <= to) {
			repeats_mask[from-1] = 1; //make the repeats 0 based
			from++;
		}
	}
	inf.close();
	return;
}

void usage(int argc, char** argv) {
	cout << "Usage: " << argv[0] << " fasta_filename masks_filename scov_filename output_gc_filename contig_name window_size(default is 35) num_bins(default is 35)\n";
	exit(1);
}

int main (int argc, char** argv) {
	if (argc!=6 && argc!=7 && argc!=8) usage(argc, argv);

	if (argc>=7) {
		window_size=atoi(argv[6]);
		if (argc==8) num_bins=atoi(argv[7]);
	}
	assert(window_size%num_bins == 0);

	//fprintf(stderr,"Using window_size %d, num_bins %d\n",window_size,num_bins);
	string fasta_filename    = argv[1];
	string masks_filename    = argv[2];
	string scov_filename     = argv[3];
	string gc_filename       = argv[4];
	string contig_name       = argv[5];

	string sequence = read_genome(fasta_filename); //read in fasta file
	//printf("Read in %d chars from %s\n",sequence_length,fasta_filename);

	ifstream inf;
	open_file_binary(inf, scov_filename);
	double * scov = new double[sequence.length()];
	inf.read((char *) scov, sizeof(double) * sequence.length());
	inf.close();

	//read in the masks
	short * masks = (short*)malloc(sizeof(short)*sequence.length());
	if (masks ==NULL) {
		fprintf(stderr,"Problem with allocating memory!\n");
		exit(1);
	}
	memset(masks, 0, sizeof(short)*sequence.length());	
	read_in_masks(masks_filename, masks);

	
	int i;
	
	//The gc content percentages are binned, where each bin represents a range of gc content percentage
	//Allocate the bins
	double bins_hits[num_bins];
	double bins_length[num_bins];
	double bins_lambdas[num_bins];
	for (i=0; i<num_bins; i++) {
		bins_hits[i]=0.0;
		bins_length[i]=0.0;
		bins_lambdas[i]=0.0;
	}

	//calculate pos2bins, i.e. the gc content at every position.
	int * pos2bin=(int*)malloc(sizeof(int)*sequence.length());
	for (i=0; i<sequence.length(); i++) pos2bin[i]=-1;

	for (i = 0; i < sequence.length(); ) {
		int left = sequence.length() - i;
		if (left > window_size) left = window_size;

		int gc_content=0;
		int at_content=0;
		int repeats_inside=0;

		int j;
		for (j=0; j<left; j++) {
			char c=sequence[i+j];
			if (c=='c' || c=='C' || c=='g' || c=='G') {
				gc_content+=1;
			}
			if (c=='a' || c=='A' || c=='t' || c=='T') {
				at_content+=1;
			}
			if (masks[i+j] || c=='N' || c=='n') {
				repeats_inside++;
			}
		}	

		int bin = whats_my_bin(gc_content,at_content);
		for (j=0; j<left; j++) {
			char c=sequence[i+j];
			if (masks[i+j] || c=='N' || c=='n') {
				pos2bin[i+j] = -1;
			} else {
				pos2bin[i+j] = bin;
			}
		}

		if (repeats_inside==0) {
			for (j=0; j<left; j++) bins_hits[bin] += scov[i+j];
			bins_length[bin] += window_size;
		}
		i += window_size;
	}

	//compute bins_lambdas:  the arrival rate for each bin
	for (i=0; i<num_bins; i++) {
		if (bins_length[i] == 0) {
			bins_lambdas[i] = -2.0;
		} else {
			bins_lambdas[i] = bins_hits[i] / bins_length[i];
		}
		//printf("B:%d, %f, H: %f , L: %f\n",i,bins_lambdas[i],bins_hits[i],bins_length[i]);
	}

	//calculate pos2lambda : expected arrivals per position
	free(masks);
	delete[] scov;
	double* pos2lambda = (double*)malloc(sizeof(double)*(sequence.length()));
	if (pos2lambda == NULL) {
		fprintf(stderr,"Trouble allocating pos2lambda\n");
		exit(1);
	}
	unsigned int repeats = 0;
	unsigned int not_repeats=0;
	unsigned int shifted_positions=0;
	for(i = 0; i < sequence.length(); i++) {
		int bin = pos2bin[i];
		//cout << "Looking at position " << i << endl;
		if (bin >= 0) {
			not_repeats++;
			if (bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0) {
				shifted_positions++;
			}
			//printf("Thinking about bin %d, %f\n",bin,bins_lambdas[bin]);
			if (bin > (num_bins/2)) {
				while(bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0)  bin--;
			} else {
				while(bins_lambdas[bin] == 0.0 || bins_lambdas[bin] == -2.0)  bin++;
			}

			if ( bin >= num_bins  || bin < 0 ) { 
				pos2lambda[i] = 0.0;
			} else {
				pos2lambda[i] = bins_lambdas[bin];
				assert(pos2lambda[i] > 0.0);
			}
		} else {
			assert(bin == -1);
			repeats++;
			pos2lambda[i] = -1.0;
		}
	}

	//write out gc file
	ofstream outf;
	open_file_binary(outf, gc_filename);
	outf.write((char *) pos2lambda, sizeof(double) * sequence.length());
	outf.close();

	//printf("Repeats: %d, Un-Repeats: %d\n",repeats, not_repeats);
	//printf("Shifted positions: %d\n",shifted_positions);

	free(pos2lambda);
	return 0;

}
