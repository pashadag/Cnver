use FindBin qw($Bin);
$affiner_dir = "$Bin/";
if ($#ARGV+1!=10) {
	print "Usage: $0 genome_length num_reads maxlines infile baseflow freelimit allownone cs2_dir";
	exit;
}

sub unmap{
    my($number) = shift;
    $preoff  = $number - 2*$numnodes -1;
    if ($preoff > 0) {
	return int(($preoff + 1)/2); 
    }
    else {
	return int(($preoff - 1)/2); 
    }
}

sub leftmap {
    my($number) = shift;
    return 2*$number + 2*$numnodes+1;
}
sub rightmap {
    my($number) = shift;
    if ($number > 0) {
	return 2*$number + 1  + 2*$numnodes+1;
    }
    else {
	return 2*$number - 1  + 2*$numnodes+1;
    }
}

sub isnode{
    my($num1) = shift;
    my($num2) = shift;
    
}

sub round{
    my($number) = shift;
    return int($number+0.5);
}


$genome_length = $ARGV[0];
$numreads = $ARGV[1];
$maxlines = $ARGV[2];
$infile = $ARGV[3];
$baseflow = $ARGV[4];
$freelimit = $ARGV[5];
$allownone = $ARGV[6];
$cs2_dir = $ARGV[9];
$genome_length = $genome_length * $baseflow;

my %nodehash = ();
my %archash = ();

open(OUTFILE, ">outfile.txt");
open(INFILE, "<$infile");

$arccntr = 0;

while ($line = <INFILE>) {
	chomp $line;
	if ($line =~ /^\/\//) {

	}
	elsif ($line =~ /^([0-9]+)/) {
		$numnodes = $1;       
		$read_arrival = $genome_length/$numreads;
	}
	elsif ($line =~ /^NODE\s+(-*[0-9]+)\s+(-*[0-9]+)/) {
		$node = $1;
		$karrivals = $2;
		$nodehash{$node} = 0;
		$from = 2*$node + 2*$numnodes+1;
		$to = 2*$node -1 + 2*$numnodes+1;
		print OUTFILE "$from $to $karrivals 0 k\n";
		$from = -2*$node + 2*$numnodes+1;
		$to = -2*$node +1 + 2*$numnodes+1;
		print OUTFILE "$from $to $karrivals 0 k\n";
	}
	elsif ($line =~ /^ARC/) {
		@arr = split(/ /, $line);

		$arccntr = $arccntr+1;
		$frnode = $arr[1];
		$tonode = $arr[2];
		$read_cnt = $arr[3];
		$len = $arr[4];

		if (defined $arr[5]) {
			$suffix = "$arr[5]";
		}
		else {
			$suffix = "";
		} 

		if ($frnode == $tonode) {
			if ($len == 0 || $read_arrival == 0) {
				$flow = 0;
			}
			else {
				$flow = round($read_cnt/($len/$read_arrival));
			}
			if ($flow == 0) { $flow = 1; }
			$archash{"$frnode,$tonode" } = 2*$flow;
			next;
		}

		$archash{"$frnode,$tonode" } = 0;
		if (defined($archash{(-$tonode).",".(-$frnode)})) {
			if (-$tonode != $frnode) {
				print STDERR "INPUT doubledge $frnode $tonode\n";
				exit(1);
			}
		}


		if ($frnode > 0) {
			$from = 2*$frnode - 1 + 2*$numnodes+1;
		} else {
			$from = 2*$frnode + 1 + 2*$numnodes+1;
		}
		$to = 2*$tonode + 2*$numnodes+1;
		print OUTFILE  "$from $to $read_cnt $len r $suffix\n";
		$tt = $frnode;
		$frnode = -1* $tonode;
		$tonode = -1* $tt;
		if ($frnode > 0) {
			$from = 2*$frnode - 1 + 2*$numnodes+1;
		} else {
			$from = 2*$frnode + 1 + 2*$numnodes+1;
		}
		$to = 2*$tonode + 2*$numnodes+1;
		print OUTFILE "$from $to $read_cnt $len r $suffix\n";
#	print STDERR "arc\n";
	}
	else {
		print STDERR "nowhere $line\n";

	}
}
$nodecnt = $numnodes * 2 + 2*$numnodes+1;
$totarccnt = 2*$numnodes + 2* $arccntr;
close(OUTFILE);
close (INFILE);
print STDERR "0 read_arrival=$read_arrival nodecnt=$nodecnt totalarccnt=$totarccnt\n";

print "Running perl $affiner_dir/affiner.pl $baseflow $read_arrival $nodecnt $maxlines $freelimit $allownone < outfile.txt > nf\n";
system("perl $affiner_dir/affiner.pl $baseflow $read_arrival $nodecnt $maxlines $freelimit $allownone < outfile.txt > nf");
$arccnt = `wc -l nf`;
$arccnt =~ /([0-9]+)\s(.*)/;
$arccnt = $1;
print STDERR "Arc Count $arccnt\n";
system("echo \"p min $nodecnt $arccnt\" > zf" );
system("echo \"n 1 0\" >> zf" );
system("cat nf >> zf");
print "Running $cs2_dir/cs2 < zf > whole_solution\n";
system("$cs2_dir/cs2 < zf > whole_solution");
print "Done with cs2.\n";
system("egrep \"^f\" whole_solution > flow_solution");

#parse flow results
open(FLOWSOLV, "<flow_solution");

while ($line = <FLOWSOLV>) {
	chomp $line;
	if ($line =~ /f\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)/) {
		$from = $1;
		$to = $2;
		$flow = $3;
		if (unmap($from) == unmap($to)) { # node
			$base = unmap($from);
			if ($base < 0) {
				$base = -1*$base;
			}
			if (!defined($nodehash{$base})) {
				print STDERR "INVALID NODE $line $base\n";
				exit(1);
			}
			if ($from%2 ==1) {
				$nodehash{$base} = $nodehash{$base} + $flow;
			}
			else {
				$nodehash{$base} = $nodehash{$base} - $flow;
			}
		}
		else { #arc
			$sign = 1;
			if ($to % 2 == 0) {
				$tt = $from;
				$from = $to;
				$to = $tt;
				$sign = -1;
			}
			$basefr = unmap($from);
			$baseto = unmap($to);
			if (defined($archash{$basefr.",".$baseto})) {
				if (defined($archash{(-$baseto).",".(-$basefr)})) {
					if (-$basefr != $baseto) {
						print STDERR "ARC valid twice $line\n";
						exit(1);
					}
				}
				$archash{$basefr.",".$baseto} =  $archash{$basefr.",".$baseto} +$sign*$flow;


			} elsif (defined($archash{(-$baseto).",".(-$basefr)})) {
				if (defined($archash{$basefr.",".$baseto})) {
					if (-$basefr != $baseto) {
						print STDERR "ARC valid twice $line\n";
						exit(1);
					}
				}
				$archash{-$baseto.",".-$basefr} =  $archash{-$baseto.",".-$basefr} + $sign*$flow;
			}
			else {
				if ($basefr == 0 || $baseto == 0) {
					if ($flow > 0) {
						print STDERR "FLOW START/END $basefr $baseto $flow\n";
					}
				}
				else {
					print STDERR "INVALID ARC $line\n $basefr $baseto\n";
					exit(1);
				}
			}

		}
	}
	else {
		print STDERR "PARSE ERROR";
		exit 1;
	}

}
#reread Graph writing flow
open(INFILE, "<$infile");
open(OUTFILE, ">$infile.out");
while ($line = <INFILE>) {
	chomp $line;
	if ($line =~ /^NODE\s+(-*[0-9]+)\s+(-*[0-9]+)/) {
		$line = $line." ".($nodehash{$1}/2);

	}
	elsif ($line =~ /^ARC\s+(-*[0-9]+)\s+(-*[0-9]+)\s+(-*[0-9]+)/) {
		$line = $line." ".($archash{$1.",".$2}/2);
	}
	print OUTFILE "$line\n";

}

close (OUTFILE);
close(INFILE);
exit(0);
