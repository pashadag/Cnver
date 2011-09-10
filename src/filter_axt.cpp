#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include "include/defs.h"
#include "include/dbtypes.h"

#define MAX_LINE_LEN    1000 /* Only for header lines - not for alignments */
#define MAX_SEQ_LEN     1048576 /* Let's hope so... */

#define SIMILARITY      99      /* TODO: make tunable */
#define WINDOW          100
#define MATCH_SCORE     100
#define RMATCH_SCORE    100
#define MISMATCH_SCORE  0
#define GAP_SCORE       (WINDOW * -MATCH_SCORE) /* No Gaps */
#define SCORE_SCALING   WINDOW * MATCH_SCORE

/* Produces a score for a set of inputs */
int align_score(char * s1, char * s2, int index)
{
    if ((s1[index] == '-') || (s2[index] == '-'))
        return GAP_SCORE;
    else if (tolower(s1[index]) == tolower(s2[index])) {
        if ((islower(s1[index])) || (islower(s2[index])))
            return RMATCH_SCORE;
        else
            return MATCH_SCORE;
    } else
        return MISMATCH_SCORE;
}

/* Takes, as input, AXT-format alignment results. Produces, as output,
    repeat graph entries (to be sorted later) */
int main(int argc, char ** argv)
{
    FILE * input;

    /* Validate input, but not too much. */
	if (argc != 4) {
		fprintf(stderr,"Compiled with:\nSIMILARITY:\t%3d\nWINDOW:\t\t%3d\nMATCH_SCORE:\t%d\nRMATCH_SCORE:\t%3d\nMISMATCH_SCORE:\t%3d\n",
				SIMILARITY,WINDOW,MATCH_SCORE,RMATCH_SCORE,MISMATCH_SCORE);
		fprintf(stderr, "usage: %s <axt_file> <contig_length> <output_file>\n", argv[0]);
		return -1;
	}


	/* Open input and output files */
	input = fopen(argv[1], "r");
	if (input == NULL) {
        perror("opening input/output");
        return -1;
    }

	ofstream outf;
	open_file(outf, argv[3]);
    
    uint64_t contiglen = atol(argv[2]);

    /* Read input and write to output */    
    char line[MAX_LINE_LEN];
    char s1[MAX_SEQ_LEN];
    char s2[MAX_SEQ_LEN];
    uint64_t cur_edge = 0;

    for (;;) {
        if (fgets(line, MAX_LINE_LEN, input) == NULL)
            break;

        if (line[0] == '#')
            continue;
        

        uint32_t alignnum;
		char contigname1[CONTIGNAME_LEN_SMALL];
        uint64_t contigstart1;
        uint64_t contigend1;
		char contigname2[CONTIGNAME_LEN_SMALL];
        uint64_t contigstart2;
        uint64_t contigend2;
        char strand;
        uint32_t score;

        /* Read in AXT header line */
        sscanf(line,
                "%u %s %llu %llu %s %llu %llu %c %u",
                &(alignnum),
                contigname1,
                &(contigstart1),
                &(contigend1),
                contigname2,
                &(contigstart2),
                &(contigend2),
                &(strand),
                &(score));
        
        /* Read in sequences */
        fgets(s1, MAX_SEQ_LEN, input);
        uint32_t s1len = strlen(s1) - 1;
        fgets(s2, MAX_SEQ_LEN, input);
        uint32_t s2len = strlen(s2) - 1;
        uint32_t minlen = (s1len < s2len ? s1len : s2len);

        /* Skip short alignments and different contigs */
        if ((minlen < WINDOW) || (strcmp(contigname1, contigname2))) {
            /* Skip extra newline) */
            while (fgetc(input) != '\n')
                ;
            
            /* Go to next alignment */
            continue;
        }
        
        /* Find similar substrings */
        int32_t score_total = 0;
        int32_t window_start = -1;
        int32_t window_end = -1;
        uint32_t s1gaps = 0;
        uint32_t s2gaps = 0;
        uint32_t skip = 0;
        
        /* Prime score total */
        for (int i = 0; i < WINDOW; i++) {
            /* Shift in new scores */
            score_total += align_score(s1, s2, i);
        }

        /* Shift window along */
        int i = 0;
        do {
            /* If we're in a similar window, */
            if (skip > 0) {
                skip--;
            } else if (score_total >= ((SIMILARITY * SCORE_SCALING) / WINDOW)) {
                /* If this is the start of the window, */
                if (window_start == -1) {
                    /* Set the window start */
                    window_start = i;
                }
                
                /* Update window end */
                window_end = i + WINDOW - 1;
            } else if (window_end != -1) {
                /* Trim window */
                while (align_score(s1, s2, window_start) != MATCH_SCORE) window_start++;
                while (align_score(s1, s2, window_end) != MATCH_SCORE) window_end--;

                /* Write out axt*/
				uint64_t start1 = contigstart1 + window_start - s1gaps;
				uint64_t end1   = contigstart1 + window_end - s1gaps;
				uint64_t start2, end2;
                if (strand == '+') {
                    start2 = contigstart2 + window_start - s2gaps;
                    end2   = contigstart2 + window_end - s2gaps;
                } else {
                    start2 = contiglen - (contigstart2 + window_end - s2gaps) + 1;
                    end2   = contiglen - (contigstart2 + window_start - s2gaps) + 1;
                }
				outf << contigname1 << "\t" << start1 << "\t" << end1 << "\t" << contigname2 << "\t" << start2 << "\t" << end2 << "\t" << strand << endl;

                skip = window_end - i;
                window_start = -1;
                window_end = -1;
            }
            
            /* Count gaps */
            if (s1[i] == '-')
                s1gaps++;
            if (s2[i] == '-')
                s2gaps++;

            /* Shift out old score and shift in new score*/
			score_total += align_score(s1, s2, i + WINDOW) - align_score(s1, s2, i);
			//printf("Window start %d window end %d\n",window_start,window_start);
            
        } while (i++ < (minlen - WINDOW));
        
        /* If we've got a left over window, */
        if (window_start != -1) {
            /* Trim window */
            while (align_score(s1, s2, window_start) != MATCH_SCORE) window_start++;
            while (align_score(s1, s2, window_end) != MATCH_SCORE) window_end--;

			/* Write out axt*/
			uint64_t start1 = contigstart1 + window_start - s1gaps;
			uint64_t end1   = contigstart1 + window_end - s1gaps;
			uint64_t start2, end2;
			if (strand == '+') {
				start2 = contigstart2 + window_start - s2gaps;
				end2   = contigstart2 + window_end - s2gaps;
			} else {
				start2 = contiglen - (contigstart2 + window_end - s2gaps) + 1;
				end2   = contiglen - (contigstart2 + window_start - s2gaps) + 1;
			}
			outf << contigname1 << "\t" << start1 << "\t" << end1 << "\t" << contigname2 << "\t" << start2 << "\t" << end2 << "\t" << strand << endl;
        }    
        
        /* Skip extra newline) */
        while (fgetc(input) != '\n')
            ;
    }
    
    /* Close files */
    fclose(input);

	outf.close();
    
    return 0;
}
