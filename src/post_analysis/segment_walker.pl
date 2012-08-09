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
		print "command failed: $!\n";
	} 
	print $retval;
	chomp $retval;
	return $retval;
}

sub usage{
	print ( "Usage: segment_walker.pl work_dir ref_name from to.
	work_dir <path>      : Location of CNVer work_dir, as specificied for original CNVer run.
	ref_name <string>    : Name of reference sequence.
	from     <int>       :
	to       <int>       : Queried segment of the reference (optional).\n");
	print ("segment_walker can be used to report the absolute copy counts found by CNVer for a given region. The program outputs the sequence edges of the graph that correspond to the walk of that region. For each edge, it outputs the amount of times it appears in the reference, the amount of times it appears in the donor, its total length, the length of the unmasked part, the DOC ratio along that egde, and the right point of the edge minus the start point of the region (the Offset).\n\n");
	exit;
}



my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};
defined($CNVER_FOLDER) or die ("Error: please set CNVER_FOLDER enviornment variable to the correct path");


if (scalar(@ARGV) != 4 && scalar(@ARGV) != 3) {
	usage();
}

my $work_dir           = $ARGV[0];
my $ref                = $ARGV[1];
my $from               = $ARGV[2];
my $to;
my $block_file         = "$work_dir/$ref/$ref.blocks.breaks.links.cov";
my $sol_file           = "$work_dir/$ref/$ref.sol";
my $scov_file          = "$work_dir/$ref/$ref.scov";
my $gc_file            = "$work_dir/$ref/$ref.gc";
my $ploidy             = 2; # this needs to be taken out.

(-d "$work_dir")            or die ("Missing $work_dir folder.");
(-e "$scov_file")           or die ("Missing $scov_file file.");
(-e "$sol_file")           or die ("Missing $sol_file file.");
(-e "$block_file")           or die ("Missing $block_file file.");
(-e "$gc_file")             or die ("Missing $gc_file file.");

if (scalar(@ARGV) == 4) {
	$to                 = $ARGV[3];
} else {
	$to = "";
}
print "From\tTo\tDon\tRef\tTotLen\tUnmaskedLen\tDOC_Ratio\tOffset\tBlock#\n";
execCommand("$CNVER_FOLDER/src/post_analysis/ppwalker $block_file $sol_file $ploidy $ref $from $to");

