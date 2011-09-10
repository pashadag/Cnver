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
	#print "Exec $command\n";
	$retval  = `$command`;
	#print "\n";
	if ( $? == -1 ) {
		print "command \"$command\" failed: $!";
	} 
	#print $retval;
	chomp $retval;
	return $retval;
}

sub usage{
	print ( "Usage: cluster_info.pl category work_dir [{ref_name category id type} | {ref_names id}]
	category <cluster|read> : Print cluster or read info
	work_dir <path>      : Location of CNVer work_dir, as specificied for original CNVer run.
	if cluster:
		ref_name <string>    : Name of reference sequence.
		id <int>             : Index of cluster
		type <0-3>           : Type of cluster.
	if read:
		ref_names_file <string>
		id <int> \n");
	
	exit;
}



my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};
defined($CNVER_FOLDER) or die ("Error: please set CNVER_FOLDER enviornment variable to the correct path");


if (scalar(@ARGV) < 4) {
	usage();
}

my $category           = $ARGV[0];
my $work_dir           = $ARGV[1];
my $mapping_files_folder = $work_dir . "/mapping_files/";

(-d "$mapping_files_folder")            or die ("Missing $mapping_files_folder folder.");



if ($category eq "cluster") {
	my $ref                = $ARGV[2];
	my $cl_idx             = $ARGV[3] - 1;
	my $cl_type            = $ARGV[4];

	(-d "$work_dir/$ref")            or die ("Missing $work_dir/$ref folder.");

	my $retval = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup normal $work_dir/$ref/$ref.t$cl_type.idx $work_dir/$ref/$ref.t$cl_type $cl_idx");
	my $head1;
	my $head2;
	($head1, $head2) = $retval =~ /\nHEADER1(.*)\nHEADER2(.*)\n/;
	print "CLUSTER HEADER: $head1\n\n";
	print "TEMPLATE MEMBERS: $head2\n\n";

	my @list_of_idx = split / /,  $head2;

	for (my $i = 1; $i < scalar(@list_of_idx); $i++) {
		my $rid = $list_of_idx[$i];
		my $bam_name1 = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup bamnames $mapping_files_folder/readNames.idx $rid");
		my $bam_name2 = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup bamnames $mapping_files_folder/readNames.idx " . ($rid + 1));
		print "MAPPINGS\t$rid\t$bam_name1\t$bam_name2\n";
		my $retval = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup double $mapping_files_folder/$ref.mmap.idx $mapping_files_folder/$ref.mmap $rid");
		print $retval;
		print "\n\n";
	}
} elsif ($category eq "read") {
	my $ref_name_file = $ARGV[2];
	my $read_id = $ARGV[3];

	open(REF_NAMES, $ref_name_file) or die("File $ref_name_file does not exist");
	my $bam_name = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup bamnames $mapping_files_folder/readNames.idx $read_id");
	print "BAM_READNAME\t$bam_name\n";
	while (<REF_NAMES>) {
		chomp;
		my $ref = $_;
		my $rmap_file = $mapping_files_folder . "/$ref.rmap";
		(-e $rmap_file) or die ("Missing $rmap_file file.");

		my $retval = execCommand("$CNVER_FOLDER/src/cluster/idx_lookup rmap $rmap_file $read_id 2> /dev/null");
		if ($retval ne "") {
			print $retval . "\n";
		}
		
	}
	close(REF_NAMES);

}



