#!/usr/bin/perl
use strict;
sub execCommand {
	my $retval;
	my ($command) = @_;
	print "Exec $command...\n";
	$retval  = `$command`;
	if ( $? != 0 ) {
		print "command failed: $!\n";
		exit;
	}
	chomp $retval;
	print $retval . "\n";
	return $retval;
}


sub print_arc {
	my ($from, $to, $rest) = @_;
	print MONFLOW "a $from $to $rest\n";
}

if ($#ARGV + 1 != 2) {
	print STDERR "Usage: $0 <flow_def> <mon_flow_name>\n";
	exit(1);
}

my $origFlow = $ARGV[0];
my $monFlow = $ARGV[1];
my $CNVER_FOLDER=$ENV{"CNVER_FOLDER"};

open(ORIGFLOW, $origFlow)   || die("Could not open file!");
open(MONFLOW, ">$monFlow") || die("Could not open file!");
my $line;
my @row;
while ($line = <ORIGFLOW>) {
	chomp $line;
	if ($line =~ /^n\s+([0-9]+)\s+([0-9]+)/) {
		my $idx = $1;
		my $bal = $2;
		print MONFLOW "n " . 2 * $idx     . " $bal\n";
		print MONFLOW "n " . (2 * $idx - 1) . " $bal\n";
	} elsif ($line =~ /^a\s+(-*[0-9]+)\s+(-*[0-9]+)\s+(.*)/) {
		my $from = $1;
		my $to   = $2;
		my $rest = $3;
		if ($from > 0 && $to > 0) {
			print_arc(2 * $from - 1, 2 * $to - 1, $rest);
			print_arc(2 * $to, 2 * $from, $rest);
		} elsif ($from < 0 && $to < 0) {
			print_arc(-2 * $to - 1, -2 * $from -1, $rest);
			print_arc(-2 * $from, -2 * $to, $rest);
		} elsif ($from > 0 && $to < 0) {
			print_arc(2 * $from - 1, -2 * $to, $rest);
			print_arc(-2 * $to - 1, 2 * $from, $rest);
		} elsif ($from < 0 && $to > 0) {
			print_arc(-2 * $from, 2 * $to - 1, $rest);
			print_arc(2 * $to, -2 * $from - 1, $rest);
		}
	} else {
		print MONFLOW "$line\n";
	}
}

close(ORIGFLOW);
close(MONFLOW);

execCommand("cat $monFlow | $CNVER_FOLDER/src/flow_solve/add_cs2_header | $CNVER_FOLDER/src/cs2-4.6/cs2 > $monFlow.sol");

# read in flow results into hash
my %hash = ();
open(MONFLOWSOL, "$monFlow.sol") || die("Could not open file!");
while ($line = <MONFLOWSOL>) {
	if ($line =~ /f\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)/) {
		my $from = $1;
		my $to = $2;
		my $flow = $3;

		#unmonotonize node indices
		my $fromb = int(($from + 1)/ 2);
		my $tob = int(($to + 1)/ 2);
		if (($from % 2) == 0) {
			$fromb *= -1;
		}
		if (($to % 2) == 0) {
			$tob *= -1;
		}

		my $key = $fromb . "," . $tob;
		if (defined($hash{$key})) { #this is for the special case of self-loops (direction switching)
			$hash{$key} += $flow / 2;
		} else {
			$hash{$key} = $flow / 2;
		}
	}
}
close(MONFLOWSOL);

#print "Hash Contents: \n"; print "$_\t$hash{$_}\n" for (sort keys %hash); print "done\n\n";

# read in flow definition file again and append it with the hash results
open(ORIGFLOW, $origFlow)   || die("Could not open file!");
open(ORIGFLOWSOL, ">$origFlow.sol")   || die("Could not open file!");
while ($line = <ORIGFLOW>) {
	if ($line =~ /^a\s+(-*[0-9]+)\s+(-*[0-9]+)\s+(.*)/) {
		my $from = $1;
		my $to   = $2;
		my $rest = $3;
		my $flow;
		if ($from != -$to) { #if not a self-loop
			$flow = $hash{"$from,$to"} + $hash{-1*$to . "," . -1*$from};
		} else { #direction switching self-loop 
			$flow = $hash{"$from,$to"};
		}
		print ORIGFLOWSOL "a $from $to $rest $flow\n";
	} else {
		print ORIGFLOWSOL $line;
	}
}

close(ORIGFLOW);
close(ORIGFLOWSOL);


