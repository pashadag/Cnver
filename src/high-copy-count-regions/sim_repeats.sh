#!/bin/sh
CHR=$1

#This is the main script for generating high-copy count regions of the genome.  
#It must be run separately for each chromosome, which is passed as a parameter to the script
#Please set the variables below to the appropriate values first.

### PARAMETERS ####

READLEN=36
OUTDIR=. #this can be any directory in which the resulting annotations will be stored
BOWTIE=~/bowtie-0.10.1/bowtie #this is the bowtie command.  Needs to be version  0.10.1
INDEX=${CHR} #bowtie index for the chromosome
REFGEN=${CHR}.fa  #reference chromosome
CHROM_LENGTH=chrom_length #this file contains a line for every chromosome, with the chromosome name (e.g. chr1) and its length.
THREADS=1 #number of threads for parallel bowtie execution



### MAIN SCRIPT ####

mkdir -p $OUTDIR
CHR_LEN=`cat $CHROM_LENGTH | awk -v chr=${CHR} '{ if ($1 == chr) print $2; }'` 
START=0
END=`expr ${CHR_LEN} - 1`

#This is the main command
sequence $REFGEN $START $END $READLEN | $BOWTIE --concise -f -p $THREADS -v 2 -k 1 -m 400  ${INDEX} - | findRepeats ${CHR_LEN} | awk -v chr=${CHR} '{ print chr , $0 }' > ${OUTDIR}/${CHR}.repeats

