#!/bin/perl
use strict;
$main::multi=100;
&main;

sub main{
    my @result;
    @result=&read_table;
    &show_table(@result);
}



sub read_table{
    my $line;
    my @data;
    my @name;
    my @table;
    my $i;

    while($line=<stdin>){
  if($line=~/^#/){
	    print $line;
	    next;
	}
	$line=~s/\s+//g;
	@data=split(/,/, $line, 21);
	if($data[0] eq ""){
	    @name=@data;
	    $i=0;
	}else{
	  $i++;
	  @{$table[$i]}=@data;
	}
    }
    return(\@name, \@table);
}

sub show_table{
    my @names=@{$_[0]};
    my @data=@{$_[1]};
    my @oneLine;
    my $i; my $j;
    my $value;

#    print "# at a scale of ln(2)*$main::multi\n";
    print "# at a scale of real numner\n";

    print " ";
    for($i=0; $i<@names; $i++){
	print $names[$i]."    ";
    }
    print "\n";

    for($i=1; $i<@data; $i++){
	print $data[$i][0]." ";

	for($j=1; $j<@{$data[$i]}; $j++){
	    $value=0.00;
	    if($data[$i][$j] ne ""){
		$value=$data[$i][$j];
	    }elsif($data[$j][$i] ne ""){
		$value=$data[$j][$i];
	    }
	    $value=1-$value/3.5;
	    if($value<0.01){
		$value=0.01;
	    }
#	    $value=int(log($value)/log(2)*$main::multi);
#	    printf("%4d ", $value);
	    printf("%1.4f ", $value);
	}
	print "\n";
    }
}

