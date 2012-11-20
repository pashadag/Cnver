#Build companion package based in current directory
#Need to have fasta files ready already in the fasta_files_folder, and a chrom_lenghths file in the cur dir.
#Chromosome names are hardocoded into the script, they must be changed to refled what's relevant for your organism
#Most of the script doesn't execute anyting but instead prints statements you should then feed into sh, probably paralellized to speed things up
#Though each of the commands have been used, this glue script hasn't been run or debugged, so expect to make some tweaks

#build cbs files
mkdir -p contig_breaks_folder
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y PUC12; do
	echo "chr$chr 1 1" > contig_breaks_folder/chr$chr.cbs
done

#build bowtie index for each chromosome
mkdir -p bowtie_index
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y PUC12; do
	echo "bowtie-build fasta_files_folder/chr$chr.fa bowtie_index/chr$chr"
done

mkdir -p self_alignments_folder
#build self-alignment for each chromosome
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y PUC12; do
	echo "lastz fasta_files_folder/chr$chr.fa --self --nomirror --format=axt > self_alignments_folder/chr$chr.axt"
done

#build repeat regions by doing alignment
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y PUC12; do
	CHR=chr$chr

	### PARAMETERS ####
	READLEN=100
	OUTDIR=./repeat_regions_folder #this can be any directory in which the resulting annotations will be stored
	BOWTIE=~/bowtie-0.12.8/bowtie #this is the bowtie command.  Needs to be version  0.10.1
	INDEX=./bowtie_index/chr$chr #bowtie index for the chromosome
	REFGEN=./fasta_files_folder/chr$chr.fa  #reference chromosome
	CHROM_LENGTH=./chrom_lengths #this file contains a line for every chromosome, with the chromosome name (e.g. chr1) and its length.
	THREADS=10 #number of threads for parallel bowtie execution

	mkdir -p $OUTDIR
	CHR_LEN=`cat $CHROM_LENGTH | awk -v chr=${CHR} '{ if ($1 == chr) print $2; }'` 
	START=0
	END=`expr ${CHR_LEN} - 1`

	#This is the main command
	echo "~/Cnver/src/high-copy-count-regions/sequence $REFGEN $START $END $READLEN | $BOWTIE --concise -f -p $THREADS -v 2 -k 1 -m 400  ${INDEX} - | ~/Cnver/src/high-copy-count-regions/findRepeats ${CHR_LEN} | awk -v chr=${CHR} '{ print chr , \$0 }' > ${OUTDIR}/${CHR}.repeats"

done
