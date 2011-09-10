#!/usr/bin/perl
use Getopt::Long;
use strict;
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
	$retval  = `$command`;
	if ( $? == -1 ) {
		print "command failed: $!";
	} 
	#print $retval;
	chomp $retval;
	return $retval;
}

sub usage{
	print ( "Usage: ref_walker.pl work_dir ref_name from to.
	work_dir <path>      : Location of CNVer work_dir, as specificied for original CNVer run.\n");
	print ("ref_walker can be used to report the absolute copy counts of a given region in the reference (taken from stdin). The program outputs the sequence edges of the graph that correspond to the walk of that region. For each edge, it outputs the amount of times it appears in the reference, the amount of times it appears in the donor, its total length, the length of the unmasked part, the DOC ratio along that egde, and the right point of the edge minus the start point of the region (the EndOffset).\n\n");
	exit;
}



my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};
defined($CNVER_FOLDER) or die ("Error: please set CNVER_FOLDER enviornment variable to the correct path");


if (scalar(@ARGV) != 1) {
	usage();
}

my $work_dir           = $ARGV[0];
(-d "$work_dir")            or die ("Missing $work_dir folder.");

my $callnum = 0;
while (my $line = <STDIN>) {
	my ($ref, $from, $to) = split (/\s+/, $line);
	chomp $line;

	my $graphinfo_file     = "$work_dir/$ref/$ref.graphinfo";
	my $problemout_file    = "$work_dir/$ref/$ref.problem.out";
	my $scov_file          = "$work_dir/$ref/$ref.scov";
	my $gc_file            = "$work_dir/$ref/$ref.gc";

	my $retval;


	#print "Don\tRef\tTotLen\tUnmaskedLen\tDOC_Ratio\tEndOffset\n";

	my $intlen = $to - $from + 1;

	# $13 is the ref flow, $17 is the length, $6 is the end position of the segment relative to the query beginning
	my $count = execCommand("$CNVER_FOLDER/src/post_analysis/edge_walker $graphinfo_file $problemout_file $from $to | grep E: | tr ':,' ' ' | awk '{ if (\$6 > -1) tot = tot + \$13 * \$17 } END { print tot / $intlen; }' ");
	print "$line\t$count\n";

	#$retval = execCommand("$CNVER_FOLDER/src/post_analysis/edge_walker $graphinfo_file $problemout_file $from $to | grep E: | tr ':,' ' ' | awk '{ if (\$21 > 0) rat=int(100*\$15/\$21)/100; else rat = 0; print \$9, \$13, \$17, \$19,  rat, \$6, \$3, \$4 }' | awk '{ if (\$6 != -1) { tot = tot + \$2 * ((\$8 - \$7 + 1)/$intlen); } } END { print tot; }' ");


	$callnum++;
}



