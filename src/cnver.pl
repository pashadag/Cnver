#!/usr/bin/perl
use Getopt::Long;
use strict;
use POSIX;
use File::Basename;
use File::Path;
use File::Copy;
use Env;

sub getTime{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return $theTime;

}

sub execCommand {
	my $retval;
	my ($command) = @_;
	print LOGFILE getTime() ."\t\t". $command . "\n";
	print "Exec at " . getTime() . " $command...\n";
	$retval  = `$command`;
	if ( $? != 0 ) {
		print "command failed: $!\n";
		exit;
	}
	chomp $retval;
	print $retval . "\n";
	return $retval;
}

sub fasta_length {
	my ($filename) = @_;
	my $lines  = `cat $filename | grep -v ">" | wc -l`;
	my $chars  = `cat $filename | grep -v ">" | wc -c`;
	my $length = $chars - $lines;
	return $length;
}

sub usage{
print ( "Usage: 
	--map_list <filename>      : A short text file containing the absolute filenames (one per line) of all the bam files you would like CNVer to use.
	--ref_folder <path>        : Location of downloaded hg18 human companion package
	--work_dir <path>          : Working directory for CNVer (will be created if it does not exist)
	--read_len <int>           : Length of reads
	--mean_insert <int>        : Mean insert size
	--stdev_insert <int>       : Standard error of insert size
	--min_mps <int>            : Minimum number of matepairs required for a cluster 

	Optional:
	--mapping_format <bam|txt> : format in which alignments are stored (default is bam)
	--mem_lim <int>            : maximum amount of memory to use during sorting stage (in gigabytes)
	--ref_names <filename>     : file with names of chromosomes of interest (default: autosomes.txt in ref_folder directory, with all the autosomes)
	--ref_single <string>      : used to limit the execution of the second stage to one reference sequence. This will override the ref_names option during the second stage.
	--mode <string>            : you can specify either stage1 or stage2 in order to only run those stages of CNVer.
	--dump                     : forces CNVer to run in dump mode, creating auxiliary files usefule for debugging but slowing down the execution.

	Note:
	1. Please give all filenames in absolute paths\n");
	exit;
}

my $mapping_format = "bam";
my $map_list;
my $ref_folder;
my $mean_insert;
my $stdev_insert;
my $min_mps;
my $command;
my $ref_name_file_of_interest = "nothing";
my $ref_single    = "all";
my $ploidy = 2; # note ploidy is not really supported.  I think report_cnvs has to change to support it.
my $temp_folder; 
my $read_len;
my $mode = "all";
my $mem_lim = "";
my $dump = 0;
my $gc_win_size = 500;
#my $gc_win_size = 200;

open ( LOGFILE, ">>log.txt") or die ("$0 : failed to open log file for output: $!\n");
my $paramLine = "";
foreach (@ARGV) { $paramLine .= "$_ " }
print LOGFILE getTime() . "\t cnver called with $paramLine\n";

GetOptions ('map_list=s' => \$map_list, 'ref_folder=s' => \$ref_folder, 'ref_names=s' => \$ref_name_file_of_interest, 'read_len=i' => \$read_len, 'mean_insert=i' => \$mean_insert, 
	'stdev_insert=i' => \$stdev_insert, 'min_mps=i' => \$min_mps, 'mapping_format=s' => \$mapping_format, 'work_dir=s' => \$temp_folder, 'mode=s' => \$mode, 
	'mem_lim=i' => \$mem_lim, 'ref_single=s' => \$ref_single, 'dump' => \$dump, 'ploidy=i' => \$ploidy ) or usage();

if (!defined($map_list) or !defined($ref_folder) or $mean_insert < 0 or $stdev_insert < 0 or $min_mps < 1 or $ploidy != 2) {
	usage();
}




my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};
defined($CNVER_FOLDER) or die ("Error: please set CNVER_FOLDER enviornment variable to the correct path");
mkpath($temp_folder);

my $contig_breaks_folder      = "$ref_folder/contig_breaks_folder";
my $repeat_regions_folder     = "$ref_folder/repeat_regions_folder";
my $self_alignments_folder    = "$ref_folder/self_alignments_folder";
my $fasta_files_folder        = "$ref_folder/fasta_files_folder";
   $ref_name_file_of_interest = "$ref_folder/autosomes.txt" if ($ref_name_file_of_interest eq "nothing");
my $ref_name_file_full        = "$ref_folder/allchr.txt";
my $mapping_files_folder      = "$temp_folder/mapping_files";
my $results_folder            = "$temp_folder/calls/";

(-d $contig_breaks_folder)       or die( "Folder $contig_breaks_folder does not exist!");
(-d $repeat_regions_folder)      or die ("Folder $repeat_regions_folder does not exist!");
(-d $self_alignments_folder)     or die("Folder $self_alignments_folder does not exist!");
(-d $fasta_files_folder)         or die("Folder $fasta_files_folder does not exist!");
(-d $temp_folder)                or die("Folder $temp_folder does not exist!");
(-e $map_list)                   or die("File $map_list does not exist");


# CNVer will not make calls shorter than this many bases (unmasked).
my $min_cnv_length = 1000;

# The minimum and maximum mapped distance which is considered concordant by the algorithm
my $min_concordant_mapped_dist = $mean_insert - 3 * $stdev_insert;
my $max_concordant_mapped_dist = $mean_insert + 3 * $stdev_insert;
$min_concordant_mapped_dist = 1 if ($min_concordant_mapped_dist < 1);

# The maximum mapped distance allowed for links.  All clusters with a bigger span then this are discarded.
my $link_len_cutoff = 10000000;

# The memory limit for sort command.
$mem_lim = "-S$mem_lim" . "g" unless ($mem_lim eq "");

# All links with at least one endpoint that falls within this distance from a contig break are discarded.
my $tolerance_around_contig_breaks = $mean_insert + 3 * $stdev_insert;

# The amount of standard deviations to allow in the width of a cluster is mean + base_len_factor.
my $base_len_factor = 1;

# The maximum difference between the main coordinates must be within MD_JOIN_TOLERANCE * STDEV from the dif in sec coord
my $md_join_tolerance = 6;



open(REF_NAMES, $ref_name_file_of_interest)            or die("File $ref_name_file_of_interest does not exist");
while (<REF_NAMES>) {
	chomp;
	(-e "$contig_breaks_folder/$_.cbs")    or die ("Missing $contig_breaks_folder/$_.cbs");
	(-e "$repeat_regions_folder/$_.rep")   or die ("Missing $repeat_regions_folder/$_.rep");
	(-e "$self_alignments_folder/$_.axt")  or die ("Missing $self_alignments_folder/$_.axt");
	(-e "$fasta_files_folder/$_.fa")       or die ("Missing $fasta_files_folder/$_.fa");
}
close(REF_NAMES);

if ($mode eq "stage1") {
	goto STAGE1;
} elsif ($mode eq "stage2") {
	goto STAGE2;
} elsif ($mode eq "all") {
	goto STAGE1;
} elsif ($mode eq "misc") {
	goto STAGE2;
} else {
	print stderr "Unknown mode: $mode.\n";
	exit;
}



STAGE1:

#Create the binary read mapping files (rmap) and the matepair mapping files (mmap)
mkpath($mapping_files_folder);
chdir $mapping_files_folder;
execCommand ("($CNVER_FOLDER/src/maps2bin $map_list $ref_name_file_full $mapping_format &);  $CNVER_FOLDER/src/cluster/concordancy_analysis  -r$read_len -d$link_len_cutoff -l$min_concordant_mapped_dist -u$max_concordant_mapped_dist -m$map_list -n$ref_name_file_full -i$mapping_format ");
chdir "..";

if ($mode eq "stage1") {
	exit;
}

STAGE2:
open(REF_NAMES, $ref_name_file_of_interest) or die("File $ref_name_file_of_interest does not exist");
while (<REF_NAMES>) {
	chomp;
	my $ref = $_;
	next if ($ref_single ne "all" and $ref_single ne $ref); 
	my $work_fol = "$temp_folder/$ref/";
	my $ref_length = fasta_length ("$fasta_files_folder/$ref.fa");
	if ($mode eq "misc") {
		goto MISC;
	}
	mkpath($work_fol);


	#make scov files
	execCommand("$CNVER_FOLDER/src/make_simple_coverage $temp_folder/mapping_files/$ref.rmap $work_fol/$ref.scov $ref $ref_length");

	#make spikes
	execCommand ("$CNVER_FOLDER/src/find_spikes $work_fol/$ref.scov  0 15 | grep -v low | awk '{ s=\$2; if (s>1) {s=s-2} else {s=0} print \$1,s,\$3+2,\$4 }' > $work_fol/$ref.spikes");

	#find N regions
	execCommand ("$CNVER_FOLDER/src/make_ns $fasta_files_folder/$ref.fa $ref > $work_fol/$ref.ns");

	#make the hard masks
	execCommand ("cat $contig_breaks_folder/$ref.cbs $repeat_regions_folder/$ref.rep $work_fol/$ref.ns | sort -k 2n | $CNVER_FOLDER/src/interval_join > $work_fol/$ref.breaks");

	#make the soft masks
	execCommand ("cat $contig_breaks_folder/$ref.cbs $repeat_regions_folder/$ref.rep $work_fol/$ref.ns $work_fol/$ref.spikes | sort -k 2n | $CNVER_FOLDER/src/interval_join > $work_fol/$ref.softmasks");

	#make GC map
	execCommand ("$CNVER_FOLDER/src/make_gc_map $fasta_files_folder/$ref.fa $work_fol/$ref.softmasks $work_fol/$ref.scov $work_fol/$ref.gc $ref $gc_win_size 50");


	#now we work separately for each type of cluster	
	my @sort_opts              = ("-k3n,3 -k4n,4", "-k4n,4 -k3n,3", "-k3n,3 -k4n,4", "-k4n,4 -k3n,3");
	# dist(1) chr(2) left(3) right(4) template(5) type(6) lstrand(7) rstrand(8)
	for (my $type = 0; $type < 4; $type++) {
		#sort disc files and cluster
		my $cluster_opts  = "--mean=$mean_insert --stdev=$stdev_insert --colID=-1 --colDist=0 --colChr=1 --colLeft=2 --colRight=3 --colTemplate=4 ";
		my $concise = 1 - $dump;
		$cluster_opts    .= "--baseLenFactor=$base_len_factor --mdJoinTolerance=$md_join_tolerance --type=$type --concise=$concise";
		execCommand ("awk '{ if (\$6 == $type) print \$0 }'  $mapping_files_folder/$ref.mmap  | sort $sort_opts[$type] -k1n,1 $mem_lim | $CNVER_FOLDER/src/cluster/cluster_matepairs $cluster_opts > $work_fol/$ref.t$type");
	}

	if ($dump) {
		for (my $type = 0; $type < 4; $type++) {
			execCommand ("$CNVER_FOLDER/src/cluster/idx_build $work_fol/$ref.t$type.idx $work_fol/$ref.t$type delim HEADER1");
			#for lookup, use ~/cnver/src/cluster/idx_lookup normal chr6.t0.idx chr6.t0 5
		}
		execCommand("$CNVER_FOLDER/src/cluster/idx_build $mapping_files_folder/$ref.mmap.idx $mapping_files_folder/$ref.mmap mmap");
		#for lookup, use ~/cnver/src/cluster/idx_lookup double chr6.mmap.idx chr6.mmap 5
	}


	#create link files 
	my $screen_opts="--mean=$mean_insert --stdev=$stdev_insert --tolerance=$tolerance_around_contig_breaks --breaksFile = $contig_breaks_folder/$ref.cbs";
	execCommand ("grep -h EDGE $work_fol/$ref.t[0123] | $CNVER_FOLDER/src/cluster/screen_contig_breaks $screen_opts | awk '{ if ((\$5 >= $min_mps) && (\$8 < $link_len_cutoff)) print \$0 }' > $work_fol/$ref.links");
	
MISC: 

	#NOTE: Clusters are now switched to one based instead of zero-based as they were before....need to correct for this later.

	#Create the graph
	execCommand ("$CNVER_FOLDER/src/filter_axt $self_alignments_folder/$ref.axt $ref_length $work_fol/$ref.glue");
	execCommand ("$CNVER_FOLDER/src/glue2block $ref $ref_length $work_fol/$ref.glue $work_fol/$ref.blocks $work_fol/$ref.edges");
	execCommand ("cat $work_fol/$ref.breaks | awk '{ print \$2; print \$3 + 1 }' | $CNVER_FOLDER/src/add_breaks $work_fol/$ref.blocks - > $work_fol/$ref.blocks.breaks");
	execCommand ("$CNVER_FOLDER/src/add_links $work_fol/$ref.blocks.breaks $work_fol/$ref.links > $work_fol/$ref.blocks.breaks.links");
	execCommand ("cat $work_fol/$ref.blocks.breaks.links | $CNVER_FOLDER/src/post_analysis/doc_walker $work_fol/$ref.scov $work_fol/$ref.gc internal  > $work_fol/$ref.blocks.breaks.links.cov");
	execCommand ("$CNVER_FOLDER/src/block2graph $work_fol/$ref.blocks.breaks.links.cov $work_fol/$ref.links $work_fol/$ref.graph $ploidy");

	#solve the flow
	execCommand ("cat $work_fol/$ref.graph | $CNVER_FOLDER/src/flow_solve/linearize > $work_fol/$ref.graph.lin");
	execCommand ("$CNVER_FOLDER/src/flow_solve/biflow_solve.pl $work_fol/$ref.graph.lin $work_fol/$ref.graph.lin.mon");
	execCommand ("$CNVER_FOLDER/src/flow_solve/flow2sol $work_fol/$ref.links.placed $work_fol/$ref.graph.lin.sol > $work_fol/$ref.sol");
	execCommand ("$CNVER_FOLDER/src/report_cnvs $work_fol/$ref.blocks.breaks.links.cov $work_fol/$ref.sol $ploidy $min_cnv_length > $work_fol/$ref.cnvs.raw");

	# remove calls that overlap contig breaks
	execCommand ("cat $work_fol/$ref.cnvs.raw | $CNVER_FOLDER/src/intRemove2  $contig_breaks_folder/$ref.cbs > $work_fol/$ref.cnvs.raw.screened");

	#smooth results
	chdir ($work_fol);
	execCommand ("$CNVER_FOLDER/src/smoother.sh 5 $min_cnv_length $ref.cnvs.raw.screened > $work_fol/$ref.cnvs.smoothed");
	chdir "..";

	next;

	#old graph code
	execCommand ("$CNVER_FOLDER/src/make_reference_graph $work_fol/$ref.edges $work_fol/$ref.graph.old");
	execCommand ("$CNVER_FOLDER/src/fill_reference_graph $work_fol/$ref.graph.old $ref_length");
	#make the donor graph
	execCommand("$CNVER_FOLDER/src/make_donor_graph $work_fol/$ref.graph.old $work_fol/$ref.links $work_fol/$ref.scov $work_fol/$ref.breaks $work_fol/$ref.gc 1> $work_fol/$ref.problem.old 2> $work_fol/$ref.graphinfo");

	my $contiglen = `tail -n 1 $work_fol/$ref.graphinfo | awk '{print \$2 }' `;
	my $numreads  = `tail -n 1 $work_fol/$ref.graphinfo | awk '{print \$3 }' `;
	chomp $contiglen;
	chomp $numreads;

	chdir ($work_fol);
	#prepare graph for cs2 solver and solve it
	execCommand ("perl $CNVER_FOLDER/src/flow_solve/make_cs2_graph.pl $contiglen $numreads 10 $ref.problem $ploidy 0 20 1 0 $CNVER_FOLDER/src/cs2-4.6/");

	#get the output cnvs
	execCommand ("$CNVER_FOLDER/src/report_cnvs_old $work_fol/$ref.graphinfo $work_fol/$ref.problem.out 100 0 2> $work_fol/$ref.cnvs.oldraw 1> $work_fol/$ref.used_dgs");




	#add doc ratios and copy cnvs to results folder
	mkpath($results_folder);
	execCommand("cat $work_fol/$ref.cnvs.smoothed | awk -v OFS=\"\\t\" '{ if (\$4 > 0) print \$1, \$2, \$3, \"gain\"; else print \$1, \$2, \$3, \"loss\" }' | $CNVER_FOLDER/src/post_analysis/doc_walker $work_fol/$ref.scov $work_fol/$ref.gc | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$4, \$6 }' > $results_folder/$ref.cnvs"); 
	execCommand("cat $work_fol/$ref.cnvs.smoothed | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$4}' | $CNVER_FOLDER/src/post_analysis/doc_walker $work_fol/$ref.scov $work_fol/$ref.gc | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$4, \$6 }' | $CNVER_FOLDER/src/post_analysis/ref_walker.pl $temp_folder | > $work_fol/$ref.cnvs.smoothed.annot"); 
	execCommand("cat $work_fol/$ref.cnvs.raw.screened | grep -v \"#\" | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$5}' | sort -k2n,2 | $CNVER_FOLDER/src/post_analysis/doc_walker $work_fol/$ref.scov $work_fol/$ref.gc | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$4, \$6 }' | $CNVER_FOLDER/src/post_analysis/ref_walker.pl $temp_folder > $work_fol/$ref.cnvs.raw.screened.annot"); 
	execCommand("cat $work_fol/$ref.cnvs.raw | grep -v \"#\" | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$5}' | sort -k2n,2 | $CNVER_FOLDER/src/post_analysis/doc_walker $work_fol/$ref.scov $work_fol/$ref.gc | awk -v OFS=\"\\t\" '{ print \$1, \$2, \$3, \$4, \$6 }' | $CNVER_FOLDER/src/post_analysis/ref_walker.pl $temp_folder > $work_fol/$ref.cnvs.raw.annot"); 
	

}

close(REF_NAMES);

