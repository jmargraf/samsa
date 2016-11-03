package CCSD;
use strict;
use warnings;
use PDL;

use Exporter 'import';

our $VERSION = '1.00';
our @EXPORT  = qw(SpinInts ECCSD FockSpin);

##################################
#  Calc CCSD energy
##################################
sub ECCSD{
  my @t1;
  my @t2;

  my $indi;
  my $indj;
  my $inda;
  my $indb;

  my $i;
  my $j;
  my $a;
  my $b;
  my @Dijab;
  my @Dia;

  my $iCC;
  my $maxCC=10;

  my $E_CC;

  my $eVorder = zeroes(2*$::dim) ;
  $eVorder = qsorti $::ES;

  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  running CCSD\n";
  close LOG;
  print     "\n";
  print     "  running CCSD\n";

# initialize t1 amplitudes
  for($indi=0;$indi<2*$::nocc;$indi++){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda++){
      $a = at($eVorder,$inda);
      $t1[$i][$a] = 0.0;
    }
  }

# initialize t2 amplitudes
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<2*$::nocc;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<2*$::dim;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $Dijab[$i][$j][$a][$b] = at($::ES,$i)+at($::ES,$j)-at($::ES,$a)-at($::ES,$b);
          $t2[$i][$j][$a][$b] = ($::SpinInts[$i][$j][$a][$b]-$::SpinInts[$i][$j][$b][$a])
                                /$Dijab[$i][$j][$a][$b];
        }
      }
    }
  }

# MBPT2 vi t2 amplitudes
  $E_CC = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_CC += ($::SpinInts[$i][$j][$a][$b]-$::SpinInts[$i][$j][$b][$a])*$t2[$i][$j][$a][$b];
        }
      }
    }
  }

  $E_CC = $E_CC/4.0;

  print "    E2 = $E_CC \n";

}

##################################
#  SpinOrbit Fock
##################################
sub FockSpin{
  my $p;
  my $q;
  my $m;
  my $ph;
  my $indp;
  my $indq;
  my $indm;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

# Transform to MO basis  
  $::Fock = transpose($::Coeff) x $::Fock x $::Coeff;

  for($indp=0;$indp<2*$::dim;$indp=$indp+2){
    $p = 2*at($eVorder,$indp/2);
    for($indq=0;$indq<2*$::dim;$indq=$indq+2){
      $q = 2*at($eVorder,$indq/2);
      $::f_s[$p][$q] = at($::Fock,$p/2,$q/2); # at($::Hcore,$p/2,$q/2)/2;
#     for($indm=0;$indm<2*$::nocc;$indm=$indm+2){
#       $m = 2*at($eVorder,$indm/2);
#       $::f_s[$p][$q] += $::SpinInts[$p][$m][$q][$m];
#        $::f_s[$p+1][$q+1] = $::f_s[$p][$q];
#     }
      $::f_s[$p+1][$q+1] = $::f_s[$p][$q];
    }
  }

  open(LOG,">>","$::name.out");
  print LOG "\n";

  $::ES= zeroes(2*$::dim);

  for($indp=0;$indp<2*$::dim;$indp=$indp+2){
    $p = 2*at($eVorder,$indp/2);
    $ph = index($::ES,$p);
    $ph .= $::f_s[$p][$p];
    $ph = index($::ES,$p+1);
    $ph .= $::f_s[$p+1][$p+1];
    print     "    $::f_s[$p][$p] ";
    print     "    $::f_s[$p+1][$p+1] \n";
    print LOG "    $::f_s[$p][$p] ";
    print LOG "    $::f_s[$p+1][$p+1] \n";
  }
  close LOG;

}

##################################
#  SpinOrbit integral transformation
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

