package MP2;
use strict;
use warnings;
use PDL;

use Exporter 'import';

our $VERSION = '1.00';
our @EXPORT  = qw(MOints EMP2 dEMP2 MOints2 dEMBPT ENPT2);

##################################
#  Calc MP2 IP correction (SpinOrbs)
##################################
sub dEMBPT{
  my $i;
  my $j;
  my $a;
  my $b;
  my $k;
  my $l;
  my $indi;
  my $indj;
  my $inda;
  my $indb;
  my $indk;
  my $indl;
  my $Dijab;

  my $E2_OS  = 0.0;
  my $E2_SSa = 0.0;
  my $E2_SSb = 0.0;
  my $E2     = 0.0;

  my @dE2;
  my @dEhpp;
  my @dEhpp_OS;
  my @dEhpp_SSa;
  my @dEhpp_SSb;
  my @dEhhp;
  my @dEhhp_OS;
  my @dEhhp_SSa;
  my @dEhhp_SSb;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

  # check bounds!
  for($indi=0; $indi<2*$::nocc; $indi=$indi+2) {
    $i = 2*at($eVorder,$indi/2);
    for($inda=2*$::nocc; $inda < 2*$::dim; $inda=$inda+2) {
      $a = 2*at($eVorder,$inda/2);
      for($indj=0; $indj < 2*$::nocc; $indj=$indj+2) {
        $j = 2*at($eVorder,$indj/2);
        for($indb=2*$::nocc; $indb < 2*$::dim; $indb=$indb+2) {
          $b = 2*at($eVorder,$indb/2);
          $Dijab = (at($::Eps,$i/2) + at($::Eps,$j/2) - at($::Eps,$a/2) - at($::Eps,$b/2));
          $E2_OS  += ($::SpinInts[$i][$j+1][$a][$b+1]*$::SpinInts[$i][$j+1][$a][$b+1])/$Dijab;

          if($indb<$inda and $indj<$indi){
            $E2_SSa += 0.5*$::SpinInts[$i][$j][$a][$b]*($::SpinInts[$i][$j][$a][$b]-$::SpinInts[$i][$j][$b][$a])/($Dijab);

            $E2_SSb += 0.5*$::SpinInts[$i+1][$j+1][$a+1][$b+1]
                     *($::SpinInts[$i+1][$j+1][$a+1][$b+1]-$::SpinInts[$i+1][$j+1][$b+1][$a+1])/($Dijab);
          }
        }
      }
    }
  }

  $E2 = $E2_OS + $E2_SSa + $E2_SSb;
  print "\n";
  print "  MBPT2 in Spin Orbital Basis: \n";
  print "    E2     = $E2 \n";
  print "    E2_OS  = $E2_OS \n";
  print "    E2_SSa = $E2_SSa \n";
  print "    E2_SSb = $E2_SSb \n";
  print "\n";
  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  MBPT2 in Spin Orbital Basis: \n";
  print LOG "    E2     = $E2 \n";
  print LOG "    E2_OS  = $E2_OS \n";
  print LOG "    E2_SSa = $E2_SSa \n";
  print LOG "    E2_SSb = $E2_SSb \n";
  print LOG "\n";
  close LOG;

# 1 hole 2 particles
  for($indi=0; $indi<2*$::dim; $indi=$indi+2) {
    $i = 2*at($eVorder,$indi/2);
    $dEhpp_OS[$indi/2]  = 0.0;
    $dEhpp_SSa[$indi/2] = 0.0;
    for($inda=2*$::nocc; $inda < 2*$::dim; $inda=$inda+2) {
      $a = 2*at($eVorder,$inda/2);
      for($indj=0; $indj < 2*$::nocc; $indj=$indj+2) {
        $j = 2*at($eVorder,$indj/2);
        for($indb=2*$::nocc; $indb < 2*$::dim; $indb=$indb+2) {
          $b = 2*at($eVorder,$indb/2);
          $Dijab = (at($::Eps,$i/2) + at($::Eps,$j/2) - at($::Eps,$a/2) - at($::Eps,$b/2));

          if($a!=$i){
            $dEhpp_OS[$indi/2]  += -($::SpinInts[$i][$j+1][$a][$b+1]
                                    *$::SpinInts[$i][$j+1][$a][$b+1])/$Dijab;
          }
          if(($b!=$i) and ($a!=$i) and $b!=$a and $j!=$i ){
            $dEhpp_SSa[$indi/2] +=  -0.5*$::SpinInts[$i][$j][$a][$b]
                                   *($::SpinInts[$i][$j][$a][$b]-$::SpinInts[$i][$j][$b][$a])/($Dijab);
          }
        }
      }
    }
  }

# 2 holes 1 particle
  for($indi=0; $indi<2*$::dim; $indi=$indi+2) {
    $i = 2*at($eVorder,$indi/2);
    $dEhhp_OS[$indi/2]  = 0.0;
    $dEhhp_SSa[$indi/2] = 0.0;
    for($inda=2*$::nocc; $inda < 2*$::dim; $inda=$inda+2) {
      $a = 2*at($eVorder,$inda/2);
      for($indk=0; $indk < 2*$::nocc; $indk=$indk+2) {
        $k = 2*at($eVorder,$indk/2);
        for($indl=0; $indl < 2*$::nocc; $indl=$indl+2) {
          $l = 2*at($eVorder,$indl/2);
          $Dijab = (at($::Eps,$k/2) + at($::Eps,$l/2) - at($::Eps,$a/2) - at($::Eps,$i/2));

          if($k!=$i){
            $dEhhp_OS[$indi/2]  += ($::SpinInts[$k][$l+1][$i][$a+1]
                                   *$::SpinInts[$k][$l+1][$i][$a+1])/($Dijab);
          }
          if(($l!=$i) and ($k!=$i) and ($l!=$k) and ($a!=$i)){
#            print "$indk $indl $indi $inda\n";
            $dEhhp_SSa[$indi/2] +=  $::SpinInts[$k][$l][$i][$a]
                                  *($::SpinInts[$k][$l][$i][$a]-$::SpinInts[$k][$l][$a][$i])/(8.0*$Dijab);
          }
        }
      }
    }
  }

  open(LOG,">>","$::name.out");
  print     "    Orb       hpp                        hhp                        DE-IP\n";
  print LOG "    Orb       hpp                        hhp                        DE-IP\n";
  for($i=0; $i<$::dim; $i++) {
    $dE2[$i] = $dEhpp_SSa[$i] + $dEhpp_OS[$i] + $dEhhp_SSa[$i] + $dEhhp_OS[$i] ;

    print     "    $i     $dEhpp_SSa[$i] + $dEhpp_OS[$i] + $dEhhp_SSa[$i] + $dEhhp_OS[$i]    $dE2[$i]    \n";
    print LOG "    $i     $dEhpp_SSa[$i] + $dEhpp_OS[$i] + $dEhhp_SSa[$i] + $dEhhp_OS[$i]    $dE2[$i]    \n";
  }
  close LOG;


}

##################################
#  Calc Epstein-Nesbet PT2 Energy
##################################
sub ENPT2{
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

  my $E_EN;

  my $eVorder = zeroes(2*$::dim) ;
  $eVorder = qsorti $::ES;

  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  running ENPT2\n";
  close LOG;
  print     "\n";
  print     "  running ENPT2\n";

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
          $t2[$i][$j][$a][$b] = ($::SpinInts[$i][$j][$a][$b]-$::SpinInts[$i][$j][$b][$a]);
        }
      }
    }
  }

# MBPT2 vi t2 amplitudes
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/$Dijab[$i][$j][$a][$b];
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "    E_EN (MP2)  = $E_EN \n";

# add -<ab||ab> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "    E_EN (-abab) = $E_EN \n";

# add -<ij||ij> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (-ijij) = $E_EN \n";

# add +<ai||ai> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]  
                                         +$::SpinInts[$a][$i][$a][$i]-$::SpinInts[$a][$i][$i][$a]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (+aiai) = $E_EN \n";

# add +<bi||bi> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]
                                         +$::SpinInts[$a][$i][$a][$i]-$::SpinInts[$a][$i][$i][$a]  
                                         +$::SpinInts[$b][$i][$b][$i]-$::SpinInts[$b][$i][$i][$b]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (+bibi) = $E_EN \n";

# add +<aj||aj> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]
                                         +$::SpinInts[$a][$i][$a][$i]-$::SpinInts[$a][$i][$i][$a]  
                                         +$::SpinInts[$b][$i][$b][$i]-$::SpinInts[$b][$i][$i][$b]
                                         +$::SpinInts[$a][$j][$a][$j]-$::SpinInts[$a][$j][$j][$a]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (+ajaj) = $E_EN \n";

# add +<bj||bj> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]
                                         +$::SpinInts[$a][$i][$a][$i]-$::SpinInts[$a][$i][$i][$a]
                                         +$::SpinInts[$b][$i][$b][$i]-$::SpinInts[$b][$i][$i][$b]
                                         +$::SpinInts[$a][$j][$a][$j]-$::SpinInts[$a][$j][$j][$a]  
                                         +$::SpinInts[$b][$j][$b][$j]-$::SpinInts[$b][$j][$j][$b]);
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (+bjbj) = $E_EN \n";

# add +<bj||bj> to denominator
  $E_EN = 0.0;
  for($indi=0;$indi<2*$::nocc;$indi=$indi+1){
    $i = at($eVorder,$indi);
    for($inda=2*$::nocc;$inda<2*$::dim;$inda=$inda+1){
      $a = at($eVorder,$inda);
      for($indj=0;$indj<$indi;$indj=$indj+1){
        $j = at($eVorder,$indj);
        for($indb=2*$::nocc;$indb<$inda;$indb=$indb+1){
          $b = at($eVorder,$indb);
          $E_EN += $t2[$i][$j][$a][$b]*$t2[$i][$j][$a][$b]/
                  ($Dijab[$i][$j][$a][$b]-$::SpinInts[$a][$b][$a][$b]+$::SpinInts[$a][$b][$b][$a]
                                         -$::SpinInts[$i][$j][$i][$j]+$::SpinInts[$i][$j][$j][$i]
                                         +$::SpinInts[$a][$i][$a][$i]-$::SpinInts[$a][$i][$i][$a]
                                         +$::SpinInts[$b][$i][$b][$i]-$::SpinInts[$b][$i][$i][$b]
                                  +0.555*($::SpinInts[$a][$j][$a][$j]-$::SpinInts[$a][$j][$j][$a]
                                         +$::SpinInts[$b][$j][$b][$j]-$::SpinInts[$b][$j][$j][$b]));
        }
      }
    }
  }

  $E_EN = $E_EN/4.0;

  print "         (scale) = $E_EN \n";

}

##################################
#  Calc MP2 IP correction
##################################
sub dEMP2{
  my $i;
  my $j;
  my $a;
  my $b;
  my $k;
  my $l;
  my $ia;
  my $ja;
  my $ib;
  my $jb;
  my $iajb;
  my $ibja;
  my $ik; #= Index2e($i,$k);
  my $al; #= Index2e($a,$l);
  my $il; #= Index2e($i,$l);
  my $ak; #= Index2e($a,$k);
  my $ikal; #= Index2e($ik,$al);
  my $ilak; #= Index2e($il,$ak);
  my $indi;
  my $indj;
  my $inda;
  my $indb;
  my $indk;
  my $indl;
  my $spin;

  my @dEmp2; 
  my @dEhhp;
  my @dEhpp;
  my @dEhhp_SS;
  my @dEhpp_SS;
  my @dEhhp_OS;
  my @dEhpp_OS;
  my @dEmp2_s;

  my $Factor = 1.0;

  print "    \n";
  print "  Calculating MP2 IP corrections \n";
  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  Calculating MP2 IP corrections\n";
  close LOG;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

  for($indi=0; $indi<$::dim; $indi++) {
    $i = at($eVorder,$indi);
#    print "$indi -> $i\n";
    $dEhpp[$indi] = 0.0;
    $dEhpp_OS[$indi] = 0.0;
    $dEhpp_SS[$indi] = 0.0;
    for($inda=$::nocc; $inda < $::dim; $inda++) {
      $a = at($eVorder,$inda);
      $ia = Index2e($i,$a);
      for($indj=0; $indj < $::nocc; $indj++) {
        $j = at($eVorder,$indj);
        $ja = Index2e($j,$a);
        for($indb=$::nocc; $indb < $::dim; $indb++) {
          $b = at($eVorder,$indb);
          $jb = Index2e($j,$b);
          $ib = Index2e($i,$b);
          $iajb = Index2e($ia,$jb);
          $ibja = Index2e($ib,$ja);
          $dEhpp[$indi] += -$Factor*$::MOijkl[$iajb]*(2.0*$::MOijkl[$iajb] - $::MOijkl[$ibja])/
                         (at($::Eps,$i) + at($::Eps,$j) - at($::Eps,$a) - at($::Eps,$b));
          $dEhpp_OS[$indi] += -1.0*$::MOijkl[$iajb]*$::MOijkl[$iajb]/
                         (at($::Eps,$i) + at($::Eps,$j) - at($::Eps,$a) - at($::Eps,$b));
          $dEhpp_SS[$indi] += -1.0*$::MOijkl[$iajb]*($::MOijkl[$iajb] - $::MOijkl[$ibja])/
                         (at($::Eps,$i) + at($::Eps,$j) - at($::Eps,$a) - at($::Eps,$b));
        }
      }
    }
  }

# What if k==i or l==i?

  for($indi=0; $indi < $::dim; $indi++) {
    $i = at($eVorder,$indi);
    $dEhhp[$indi] = 0.0;
    $dEhhp_OS[$indi] = 0.0;
    $dEhhp_SS[$indi] = 0.0;
    for($inda=$::nocc; $inda < $::dim; $inda++) {
      $a = at($eVorder,$inda);
      for($indk=0; $indk < $::nocc; $indk++) {
        $k = at($eVorder,$indk);
        for($indl=0; $indl < $::nocc; $indl++) {
          $l = at($eVorder,$indl);
          $ik = Index2e($i,$k);
          $al = Index2e($a,$l);
          $il = Index2e($i,$l);
          $ak = Index2e($a,$k);
          $ikal = Index2e($ik,$al);
          $ilak = Index2e($il,$ak);
          $dEhhp[$indi] += $Factor*$::MOijkl[$ilak]*(2.0*$::MOijkl[$ilak] - $::MOijkl[$ikal])/
                        (at($::Eps,$k) + at($::Eps,$l) - at($::Eps,$i) - at($::Eps,$a));
          $dEhhp_OS[$indi] += $::MOijkl[$ilak]*$::MOijkl[$ilak]/
                        (at($::Eps,$k) + at($::Eps,$l) - at($::Eps,$i) - at($::Eps,$a));
          $dEhhp_SS[$indi] += $::MOijkl[$ilak]*($::MOijkl[$ilak] - $::MOijkl[$ikal])/
                        (at($::Eps,$k) + at($::Eps,$l) - at($::Eps,$i) - at($::Eps,$a));
        }
      }
    }
  }

  open(LOG,">>","$::name.out");
  print     "    Orb       hpp                        hhp                        DE-IP\n";
  print LOG "    Orb       hpp                        hhp                        DE-IP\n";
  for($i=0; $i<$::dim; $i++) {
    $dEmp2[$i] = $dEhhp[$i] + $dEhpp[$i] ;
    $dEmp2_s[$i] = $dEhhp_SS[$i] + $dEhhp_OS[$i] + $dEhpp_SS[$i] + $dEhpp_OS[$i] ;

    print     "    $i     $dEhpp[$i]    $dEhhp[$i]    $dEmp2[$i]    $dEmp2_s[$i] \n";
    print LOG "    $i     $dEhpp[$i]    $dEhhp[$i]    $dEmp2[$i]    $dEmp2_s[$i]\n";
  }
  close LOG;


}

##################################
#  Calc MP2 energy
##################################
sub EMP2{
  my $i;
  my $j;
  my $a;
  my $b;
  my $ia;
  my $ja;
  my $ib;
  my $jb;
  my $iajb;
  my $ibja;

  my $indi;
  my $indj;
  my $inda;
  my $indb;

  my $E_SS = 0.0;
  my $E_SSx = 0.0;
  my $E_SSc = 0.0;
  my $E_OS = 0.0;
  my $Emp2 = 0.0;

  print "    \n";
  print "  Calculating MP2 correlation energy\n";
  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "  Calculating MP2 correlation energy\n";
  close LOG;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

  for($indi=0; $indi<$::nocc; $indi++) {
    $i = at($eVorder,$indi);
    for($inda=$::nocc; $inda < $::dim; $inda++) {
      $a = at($eVorder,$inda);
      $ia = Index2e($i,$a);
      for($indj=0; $indj < $::nocc; $indj++) {
        $j = at($eVorder,$indj);
        $ja = Index2e($j,$a);
        for($indb=$::nocc; $indb < $::dim; $indb++) {
          $b = at($eVorder,$indb);
          $jb = Index2e($j,$b);
          $ib = Index2e($i,$b);
          $iajb = Index2e($ia,$jb);
          $ibja = Index2e($ib,$ja);
          $E_OS += -($::MOijkl[$iajb] * $::MOijkl[$iajb])/
                    (at($::Eps,$a) + at($::Eps,$b) - at($::Eps,$i) - at($::Eps,$j));
          $E_SS += -($::MOijkl[$iajb] - $::MOijkl[$ibja])*$::MOijkl[$iajb]/
                    (at($::Eps,$a) + at($::Eps,$b) - at($::Eps,$i) - at($::Eps,$j));
          $E_SSx += +($::MOijkl[$ibja] * $::MOijkl[$iajb])/
                    (at($::Eps,$a) + at($::Eps,$b) - at($::Eps,$i) - at($::Eps,$j));
          $E_SSc += -($::MOijkl[$iajb] * $::MOijkl[$iajb])/
                    (at($::Eps,$a) + at($::Eps,$b) - at($::Eps,$i) - at($::Eps,$j));

#          $Emp2 += $::MOijkl[$iajb] * (2.0 * $::MOijkl[$iajb] - $::MOijkl[$ibja])/
#                   (at($::Eps,$i) + at($::Eps,$j) - at($::Eps,$a) - at($::Eps,$b));
        }
      }
    }
  }

  $Emp2 = $E_SS + $E_OS;

  print " \n";
  print "    Emp2  =   $Emp2 \n";
  print "    E_SS  =   $E_SS \n";
  print "    E_OS  =   $E_OS \n";
  print "    E_SSx =   $E_SSx \n";
  print "    E_SSc =   $E_SSc \n";
  open(LOG,">>","$::name.out");
  print LOG "\n";
  print LOG "    Emp2 =   $Emp2\n";
  print LOG "    E_SS =   $E_SS \n";
  print LOG "    E_OS =   $E_OS \n";

  close LOG;

}

##################################
#  AO to MO integral transformation (4N^5)
##################################
sub MOints{
  my $i;
  my $j;
  my $k;
  my $l; #AO indices
  my $p; 
  my $q;
  my $r;
  my $s; #MO indices
  my $ij;
  my $kl;
  my $pq;
  my $rs;
  my $pqrs;
  my $ijkl;
  my @interi;
  my @interj;
  my @interk;
  my @interl;
  my $temp;
  my $quart=0;
  my $half=0;
  my $twoth=0;

  print "    \n";
  print "  Running integral transformation\n";
  open(LOG,">>","$::name.out");
  print LOG "    \n";
  print LOG "  Running integral transformation\n";
  close LOG;

  for($i=0,$ijkl=0; $i < $::dim; $i++) {
    for($j=0; $j <= $i; $j++) {
      for($k=0; $k <= $i; $k++) {
        for($l=0; $l <= ($i==$k ? $j : $k); $l++,$ijkl++) {
          $::MOijkl[$ijkl]=0.0;
        }
      }
    }
  }

  for($l=0; $l < $::dim; $l++) {
    for($p=0; $p < $::dim; $p++) {
      for($q=0; $q < $::dim; $q++) {
        for($r=0; $r < $::dim; $r++) {
          for($s=0; $s < $::dim; $s++) {
            $pq = Index2e($p,$q);
            $rs = Index2e($r,$s);
#            $rl = Index2e($r,$l);
            $pqrs = Index2e($pq,$rs);
#            $pqrl = Index2e($pq,$rl);

            $interl[$p][$q][$r][$l] += at($::Coeff,$l,$s)*at($::Iijkl,$pqrs);          
          }
        }
      }
    }

  }

  print "      25% ...\n";

  for($k=0; $k < $::dim; $k++) {

    for($p=0; $p < $::dim; $p++) {
      for($q=0; $q < $::dim; $q++) {
        for($r=0; $r < $::dim; $r++) {
          for($l=0; $l < $::dim; $l++) {
#            $pq = Index2e($p,$q);
#            $rl = Index2e($r,$l);
#            $kl = Index2e($k,$l);
#            $pqrl = Index2e($pq,$rl);
#            $pqkl = Index2e($pq,$kl);

            $interk[$p][$q][$k][$l] += at($::Coeff,$k,$r)*$interl[$p][$q][$r][$l];
          }
        }
      }
    }

  }

  @interl = 0;
  print "      50% ...\n";

  for($j=0; $j < $::dim; $j++) {

    for($p=0; $p < $::dim; $p++) {
      for($q=0; $q < $::dim; $q++) {
        for($k=0; $k < $::dim; $k++) {
          for($l=0; $l < $::dim; $l++) {
#            $pq = Index2e($p,$q);
#            $pj = Index2e($p,$j);
#            $kl = Index2e($k,$l);
#            $pqkl = Index2e($pq,$ll);
#            $pjkl = Index2e($pj,$kl);

            $interj[$p][$j][$k][$l] += at($::Coeff,$j,$q)*$interk[$p][$q][$k][$l];
          }
        }
      }
    }

  }

  @interk = 0;
  print "      75% ...\n";


  for($i=0; $i < $::dim; $i++) {

    for($p=0; $p < $::dim; $p++) {
      for($j=0; $j < $::dim; $j++) {
        for($k=0; $k < $::dim; $k++) {
          for($l=0; $l < $::dim; $l++) {
#            $pj = Index2e($p,$j);
#            $ij = Index2e($i,$j);
#            $kl = Index2e($k,$l);
#            $pjkl = Index2e($pj,$kl);
#            $ijkl = Index2e($ij,$kl);

            $interi[$i][$j][$k][$l] += at($::Coeff,$i,$p)*$interj[$p][$j][$k][$l];
          }
        }
      }
    }

  }

  @interj = 0;
  print "      100% ...\n";

  for($i=0; $i < $::dim; $i++) {
    for($j=0; $j <= $i; $j++) {
      $ij = Index2e($i,$j);
      for($k=0; $k < $::dim; $k++) {
        for($l=0; $l <= $k; $l++) {
          $kl = Index2e($k,$l);
          $ijkl = Index2e($ij,$kl);
          $::MOijkl[$ijkl] = $interi[$i][$j][$k][$l];
         }
       }
     }
   } 
       
}


##################################
#  AO to MO integral transformation (N^8)
##################################
sub MOints2{
  my $i;
  my $j;
  my $k;
  my $l; #MO indices
  my $p;
  my $q;
  my $r;
  my $s; #AO indices
  my $ijkl;
  my $pq;
  my $rs;
  my $pqrs;
  my $temp;
  my $quart=0;
  my $half=0;
  my $twoth=0;

  print "    Running integral transformation\n";
  open(LOG,">>","$::name.out");
  print LOG "    \n";
  print LOG "    Running integral transformation\n";
  close LOG;

  for($i=0,$ijkl=0; $i < $::dim; $i++) {
    for($j=0; $j <= $i; $j++) {
      for($k=0; $k <= $i; $k++) {
        for($l=0; $l <= ($i==$k ? $j : $k); $l++,$ijkl++) {
          $::MOijkl[$ijkl]=0;
        }
      }
    }
  }



  for($i=0,$ijkl=0; $i < $::dim; $i++) {
    for($j=0; $j <= $i; $j++) {
      for($k=0; $k <= $i; $k++) {
        for($l=0; $l <= ($i==$k ? $j : $k); $l++,$ijkl++) {
        # Check progress
          if($ijkl/$::dim2 >= 0.24 and 
             $ijkl/$::dim2 <= 0.26 and
             $quart==0){
            $quart=1;
            print "      25% ($ijkl)\n";
          }
          if($ijkl/$::dim2 >= 0.49 and 
             $ijkl/$::dim2 <= 0.51 and
             $half==0){
            $half=1;
            print "      50% ($ijkl)\n";
          }
          if($ijkl/$::dim2 >= 0.74 and 
             $ijkl/$::dim2 <= 0.76 and
             $twoth==0){
            $twoth=1;
            print "      75% ($ijkl)\n";
          }
          for($p=0; $p < $::dim; $p++) {
            for($q=0; $q < $::dim; $q++) {
              $pq = Index2e($p,$q);
              for($r=0; $r < $::dim; $r++) {
                for($s=0; $s < $::dim; $s++) {
                  $rs = Index2e($r,$s);
                  $pqrs = Index2e($pq,$rs);
                   
                  $::MOijkl[$ijkl] += at($::Coeff,$i,$p) * 
                                      at($::Coeff,$j,$q) * 
                                      at($::Coeff,$k,$r) * 
                                      at($::Coeff,$l,$s) * 
                                      at($::Iijkl,$pqrs);

                }
              }
            }
          }
#          print "MO($ijkl) = $::MOijkl[$ijkl] \n";
 
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

