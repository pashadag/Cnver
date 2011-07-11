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
	print $command;
	$retval  = `$command`;
	if ( $? == -1 ) {
		print "command failed: $!";
	}
	print $retval;
	chomp $retval;
	return $retval;
}

sub usage{
	print ( "Usage: find_support.pl work_dir
		work_dir <path>      : Location of CNVer work_dir, as specificied for original CNVer run.
		stdin <call file>    : List of calls for which to find doc and link support\n");

	exit;
}



my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};

defined($CNVER_FOLDER) or die ("Error: please set CNVER_FOLDER enviornment variable to the correct path");


if (scalar(@ARGV) != 1) {
	usage();
}

my $work_dir           = $ARGV[0];

while (my $line = <STDIN>) {
	my ($ref, $from, $to) = split (/\s+/, $line);
	chomp $line;

	my $graphinfo_file     = "$work_dir/$ref/$ref.graphinfo";
	my $problemout_file    = "$work_dir/$ref/$ref.problem.out";
	my $scov_file          = "$work_dir/$ref/$ref.scov";
	my $gc_file            = "$work_dir/$ref/$ref.gc";
	my $link_file          = "$work_dir/$ref/$ref.links";

	my $docratio = `$CNVER_FOLDER/src/post_analysis/doc_walker $scov_file $gc_file $from $to | awk '{ print \$5 }'`;
	my $from_loose = $from - 3000;
	my $to_loose = $to + 3000;

	my $linkevid = `cat $link_file | intContain2 $ref $from_loose $to_loose | grep EDGE | tr '\\n' '\\t' `;
	chomp $docratio;
	chomp $linkevid;

	print "$line\t$docratio\t$linkevid\n";


}

