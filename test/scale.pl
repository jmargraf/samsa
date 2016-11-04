#!/usr/bin/perl

use strict;
use warnings;
use PDL;

$PDL::BIGPDL = 1; 

my $a; 
my $v;
my $t0;
my $t1;
my $dt;
my $i;


for($i=1000;$i<=100000;$i=$i+1000){
  $t0 = time;
  $a = ones $i,$i;
  $v = zeroes $i;

#  print $a;

  ($a,$v) = eigens_sym $a;

#  print $a;
  $t1=time;
  $dt=$t1-$t0;
  print "$i  $dt\n";
}
