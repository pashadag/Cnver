#! /usr/bin/perl

#use Math::CDF qw(:all);

sub round{
    my($number) = shift;
    return int($number+0.5);
}

sub floor{
    my($number) = shift;
    return int($number);
}

sub ceil{
    my($number) = shift;
    return int($number+.9999);
}

sub evval{
    my($x) = shift;
    my($k) = shift;
    #print "zz $x  $k\n";
    if ($x == 0) {
	return 999999999;
    }
    $z =  -$k * log($x) + $x * $lambda;
    #print STDERR "eval $z\n";
    return $z;
}

$contigcost = 10000000;
$multfactor = 1000;
$baseflow = $ARGV[0]; # some kind of node parameter
$global_read_lambda = $ARGV[1]; # read arrival every this base pairs
$nodecnt = $ARGV[2];
#$arccnt = $ARGV[3];
$maxlines = $ARGV[3];
$freelimit = $ARGV[4];

$allownone = 1;

if (defined $ARGV[5]) {
   $allownone = $ARGV[5];
}

$supersource = int($nodecnt/2)+1;

#print "p min $nodecnt $arccnt\n";

while ($line=<STDIN>) {
    $line =~ /([0-9]+) ([0-9]+) ([0-9.]+) (-*[0-9]+) (.)\s*([0-9]*)/;
    chomp($line);
    @arr = split(/ /, $line);
    $from = $arr[0];
    $to = $arr[1];
    $k = $arr[2];
    $len = $arr[3];
    $type = $arr[4];

    if (defined($arr[5])) {
	$edge_lambda = $arr[5];
	chomp $edge_lambda;
    }
    else {
	$edge_lambda = 0;
    } 

    if ($type eq "k") {
	print "a $from $to 0 99999 0\n";
#	if ($allownone > 0) {
#	    $slope = int($multfactor*($k * $allownone));
#	    print "a $to $from 0 1 $slope\n";
#	}

	next;
    }
    elsif ($type eq "r") {
	if ($len < 0 && $k == 0) {
	    print "a $from $to $baseflow $baseflow 0\n";
	    next;
	}
	if ($len >= 0 && $len <= $freelimit) {
	    print "a $from $to 0 99999 1\n"; #make them pay something for "free" edges to prevent huge loops
		next;
	}
	$global_lambda = $global_read_lambda;
#all the stuff below is actually this case; 
    }
    else {
	print STDERR "INVALID EDGE TYPE $type in $line\n";
	exit(1);
    }



    #$lambda = $len/$global_lambda;
    $lambda=$edge_lambda*$len/$baseflow;	
	#print "E: $edge_lambda $len $k\n";
	#$b=$k/$lambda;
	#print "B: $b\n";
    $up = ceil($k/$lambda+0.0001);
    $down = floor($k/$lambda+0.0001);
    #$yyy=evval($up,$k);
    #$zzz=evval($down,$k); 
   #print "TT ,$yyy,$zzz,$len,$lambda, $k, $global_lambda, $up, $down, $k/$lambda";   
 if (evval($up,$k) > evval($down,$k)) {
	$start = $down;
    }
    else {
	$start = $up;
    }

    if ($forced) {
	$end = $start+1;
	$slope = int($multfactor *(evval($end,$k)  - evval($start,$k))/($end-$start));
	print "a $from $to $start $start $slope\n";

	if ($forced > $start) {
	    $slope = int($multfactor *(evval($forced,$k)  - evval($start,$k))/($forced-$start));
	    $zz = $forced - $start;
	    print "a $from $to $zz $zz $slope\n";
	}
	elsif ($forced < $start) {
		$slope =  int($multfactor * (-1 *  (evval($forced,$k)  - evval($start,$k))/($forced-$start)));
		$zz= $start-$forced;
		print "a $to $from $zz $zz $slope\n";
	}
	next;
    }

    $end = $start + 1;
    $slope =  int($multfactor*(evval($end,$k)  - evval($start,$k))/($end-$start));
    if ($slope < 0) {
	$z = evval($end,$k);
	$y =  evval($start,$k);
	print STDERR "NEGSLOPE $lambda $start $end $k $len ($z $y) $slope\n"

    }
    print "a $from $to $start $end $slope\n";

    $istart = $start;
    $i = 1;
    while ($i < $maxlines) {
	if ($i == $maxlines -1) {
	    $step = 2;
	} else {
	    $step = 1;
	}
	$start = $end; 
	$end = $start + $step;	
	$slope = int($multfactor *(evval($end,$k)  - evval($start,$k))/($end-$start));
	#print "$start vs $end, $k\n";
	print "a $from $to 0 $step $slope\n";
	$i = $i + 1;
    }

    $start = $end; 
    $end = $start + 10;	
    $slope = int($multfactor *(evval($end,$k)  - evval($start,$k))/($end-$start));
    print "a $from $to 0 99 $slope\n";
    $slope = 0;
    $i = 0;
    $start = $istart;
    while ($i < $maxlines) {
	if ($start > 2 && $i == $maxlines -1) {
	    $step = 2;
	} else {
	    $step = 1;
	}
	$end = $start;
	$start = $end-$step;
	if ($start < 1) { $start = 1; last; }
	$slope =  int($multfactor * (-1 *  (evval($end,$k)  - evval($start,$k))/($end-$start)));
	print "a $to $from 0 $step $slope\n";
	$i = $i + 1;
    }
    if ($start > 1) {
	$end = $start;
	$start = $end-10;
	if ($start < 1) { $start = 1; }
	$slope =  int($multfactor * (-1 *  (evval($end,$k)  - evval($start,$k))/($end-$start)));
        $step = $end-$start;
	print "a $to $from 0 $step $slope\n";
    }
    if ($allownone > 0) {
	$end = $start;
	$start = 1/($allownone+1);
        $slope =  int($multfactor * (-1 *  (evval($end,$k)  - evval($start,$k))/($end-$start)));
	$step = $end;
	print "a $to $from 0 $step $slope\n";
    }
}
