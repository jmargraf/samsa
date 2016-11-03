package CCSD;
use strict;
use warnings;
use PDL;

use Exporter 'import';

our $VERSION = '1.00';
our @EXPORT  = qw(SpinInts ECCSD);

##################################
#  Calc CCSD energy
##################################
sub ECCSD{

}

##################################
#  AO to MO integral transformation
##################################
sub SpinInts{
  my $p;
  my $q;
  my $r;
  my $s;
  my $pr;
  my $qs;
  my $ps;
  my $qr;
  my $prqs;
  my $psqr;
  my $value2;
  my $value1;
  my $quart=0;
  my $half=0;
  my $twoth=0;

  print "    \n";
  print "  Spin-Orbital Integral Transform\n";
  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  Spin-Orbital Integral Transform\n";
  close LOG;


  for($p=0; $p < $::dim*2; $p++){
    for($q=0; $q < $::dim*2; $q++){
      for($r=0; $r < $::dim*2; $r++){
        for($s=0; $s < $::dim*2; $s++) {
          { use integer;
          $pr = Index2e($p/2,$r/2);
          $qs = Index2e($q/2,$s/2);
          }
          $prqs = Index2e($pr,$qs);
          $value1 = $::MOijkl[$prqs] * ($p%2 == $r%2) * ($q%2 == $s%2);
          { use integer;
          $ps = Index2e($p/2,$s/2);
          $qr = Index2e($q/2,$r/2);
          }
          $psqr = Index2e($ps,$qr);
          $value2 = $::MOijkl[$psqr] * ($p%2 == $s%2) * ($q%2 == $r%2);
          $::SpinInts[$p][$q][$r][$s] = $value1 - $value2;
        }
      }
    }
  }

}

sub Index2e{
  my $i = $_[0];
  my $j = $_[1];
  my $ij;

  if($i>$j){
    $ij = $i*($i+1)/2 + $j;
  }else{
    $ij = $j*($j+1)/2 + $i;
  }

  return($ij);
}

