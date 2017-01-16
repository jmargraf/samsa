package LinAl;
use strict;
use warnings;
use PDL;
#use PDL::LinearAlgebra::Real;
#use PDL::Ufunc;
#use PDL::Slatec;

use Exporter 'import';

our $VERSION = '1.00';
our @EXPORT  = qw(DiaS BuildS12 CalcDens CalcEnergy CalcFock DiaF FtoMO FtoAO CalcCorP);

##################################
#  Diagonalize Overlap
##################################
sub DiaS{ 
  my $info=null;

  open(LOG,">>","$::name.out");

  print LOG "  Diagonalizing Overlap: \n";

  $::S_vec = zeroes($::dim,$::dim);
  $::S_val = zeroes($::dim);

  $::S_vec .= $::Sij;

#  syev($::S_vec, 1,0, $::S_val,  $info);
#  print  "Info: $info\n";
  ($::S_vec,$::S_val) = eigens_sym $::S_vec ;
  print LOG "    ... done \n";

  if($::Debug){
    print LOG "\n";
    print LOG "Eigenvalues \n";
    print LOG $::S_val;
    print LOG "\n";
    print LOG "Eigenvectors\n";
    print LOG $::S_vec;
    print LOG "\n";
  }

  close LOG;
}

##################################
#  Build Orthogonalization Matrix
##################################
sub BuildS12{
  my $Temp = zeroes($::dim,$::dim);
  my $Smat = zeroes($::dim,$::dim);

  open(LOG,">>","$::name.out");

  print LOG "  Building S12: \n";

  $::S12   = zeroes($::dim,$::dim);

  $::S_val = $::S_val**(-0.5);
  $Smat = stretcher($::S_val);

#  gemm($Smat, 0, 1, $::S_vec, pdl(1.0), pdl(0.0), $Temp);
#  gemm($::S_vec, 0, 0, $Temp, pdl(1.0), pdl(0.0), $::S12);

  $::S12 = $::S_vec x $Smat x transpose($::S_vec);

#  $Temp  = $Smat x transpose($::S_vec);  
#  $::S12 = $::S_vec x $Temp;

  print LOG "    ... done \n";

  if($::Debug){
    print LOG "\n";
    print LOG "Smat: \n";
    print LOG $Smat;
    print LOG "\n";
    print LOG "\n";

    print LOG "Sval: \n";
    print LOG $::S_val;
    print LOG "\n";
    print LOG "\n";
    print LOG "S-12: \n";
    print LOG $::S12;
    print LOG "\n";
    print LOG "\n";
  }

  close LOG;
}


##################################
#  Transform Fockian to MO basis
##################################
sub FtoMO{
  my $i;
  my $fii;

#  open(LOG,">>","$::name.out");
#  print LOG "  Transform Fock to MO:\n";
#  print LOG "\n";
#  close LOG;

# Orthogonalize Fockian
#  $::Fock = transpose($::S12) x $::Fock x $::S12; # $Temp;
  print "Before: \n";
  print $::Fock;

# Transform to MO basis  
  $::Fock = transpose($::Coeff) x $::Fock x $::Coeff;

#  open(LOG,">>","$::name.out");
#  for($i=0;$i<$::dim;$i++){
#    $fii = at($::Fock,$i,$i);
#    print LOG "F' $i $i = $fii \n" 
#  } 
#  print LOG "\n";
#  close LOG;

}

##################################
#  Transform Fockian to AO basis
##################################
sub FtoAO{
  my $i;
  my $fii;

#  open(LOG,">>","$::name.out");
#  print LOG "  Transform Fock to AO:\n";
#  print LOG "\n";
#  close LOG;

# Orthogonalize Fockian
#  $::Fock = transpose($::S12) x $::Fock x $::S12; # $Temp;

# Transform to AO basis  
  $::Fock = $::Sij x $::Coeff x $::Fock x transpose($::Coeff) x $::Sij;

  print "After: \n";
  print $::Fock;

#  open(LOG,">>","$::name.out");
#  for($i=0;$i<$::dim;$i++){
#    $fii = at($::Fock,$i,$i);
#    print LOG "F' $i $i = $fii \n"
#  }
#  print LOG "\n";
#  close LOG;

}

##################################
# Extended Hueckel Initial Guess
##################################
sub EHT {
  my $Temp = zeroes($::dim,$::dim);
  my $i;
  my $j;
  my %H;
  my $basi;
  my $basj;
  my $ph1;
  my $K = 1.75;

  $H{"H 1S"}   = -13.06 ;
  $H{"H 2P"}   = -3.4  ;
  $H{"Li 2S"}  = -5.41 ;
  $H{"Li 2P"}  = -3.61 ;
  $H{"Be 2S"}  = -9.33 ;
  $H{"Be 2P"}  = -5.88 ;
  $H{"B 2S"}   = -14.0 ;
  $H{"B 2P"}   = -8.24 ;
  $H{"C 2S"}   = -19.42;
  $H{"C 2P"}   = -10.7 ;
  $H{"N 2S"}   = -25.58;
  $H{"N 2P"}   = -13.25;
  $H{"O 2S"}   = -32.49;
  $H{"O 2P"}   = -15.88;
  

  for($i=0;$i<$::dim;$i++){
    $basi = "$::atom[$::basis[$i]] $::bastyp[$i]";
    if(
       (($::atom[$::basis[$i]] eq "H") and $::bastyp[$i] ne "1S" ) or
       (($::atom[$::basis[$i]] ne "H") and $::bastyp[$i] eq "1S" ) or 
       ( $::bastyp[$i] eq "Pol P" ) or
       ( $::bastyp[$i] eq "Pol D" ) 
                              ){
      next;
    }
#    print "$basi \n";
    $ph1  = index2d($::Fock,$i,$i);
    $ph1 .= $K*at($::Sij,$i,$i)*
           ($H{$basi}+$H{$basi})/2.0 ;
#    $ph1 = 0;
#    print "$ph1 \n";
  }

#  print $::Fock;

}

##################################
#  Diagonalize Fockian
##################################
sub DiaF {
  my $Temp = zeroes($::dim,$::dim);

  open(LOG,">>","$::name.out");

# Initial Guess
  if($::iSCF<1){
    my $pert = randsym($::dim,$::dim) x 0.01 ;
    print LOG "  Initial Guess:\n";
    $::Fock  = zeroes($::dim,$::dim);
    $::Coeff = zeroes($::dim,$::dim);
    $::Eps   = zeroes($::dim);
    $::Dens  = zeroes($::dim,$::dim);
    if($::UHF==1){
      $::FockB = zeroes($::dim,$::dim) ;
      $::CoeffB = zeroes($::dim,$::dim);
      $::EpsB   = zeroes($::dim);
      $::DensB  = zeroes($::dim,$::dim);
      
    }

    $::Fock  .= $::Hcore;
    if($::EHT eq "true"){
      print LOG "    Hueckel\n";
      EHT();
    }

    if($::UHF==1){
      $::FockB=$::Fock+$pert;
    }#

    print LOG "    ... done\n";
    print LOG "\n";
  }

# Orthogonalize Fockian
#  gemm($::Fock, 0, 0, $::S12, pdl(1.0), pdl(0.0), $Temp);
#  gemm($::S12, 1, 0, $Temp, pdl(1.0), pdl(0.0), $::Coeff);
  $::Coeff = transpose($::S12) x $::Fock x $::S12; # $Temp;
  if($::UHF==1){
    $::CoeffB = transpose($::S12) x $::FockB x $::S12; # $Temp;
  }

# Diagonalize Fock Matrix
#   my $ifail = zeroes($::dim);
#   syevx($::Coeff,1,0,1,0,0,0,0,0.0, $::dim, $::Eps, $::Coeff ,$ifail, $info);
#   syevd($::Coeff, 1,1, $::Eps,  $info);
#   syev($::Coeff, 1,0, $::Eps,  $info);
#   print  "Info: $info\n";
#   ($::Eps,$::Coeff) = eigsys($::Coeff);
  ($::Coeff,$::Eps) = eigens_sym $::Coeff ;
  if($::UHF==1){
    ($::CoeffB,$::EpsB) = eigens_sym $::CoeffB ;
  }


# Transform back into original basis
#  $Temp .= $::Coeff;
#  gemm($::S12, 0, 0, $Temp, pdl(1.0), pdl(0.0), $::Coeff);
  $::Coeff = $::S12 x $::Coeff;
  if($::UHF==1){
    $::CoeffB = $::S12 x $::CoeffB;
  }

  if($::Debug){
    print LOG "  Transformed Fockian:\n";
    print LOG $::Coeff;
    print LOG "\n";

    print LOG "  WF (ort. Basis):\n";
    print LOG $::Coeff;
    print LOG "\n";

    print LOG "  S12:\n";
    print LOG $::S12;
    print LOG "\n";

    print LOG "  WF:\n";
    print LOG $::Coeff;
    print LOG "\n";

    print LOG "  Eigenvalues:\n";
    print LOG $::Eps;
    print LOG "\n";
  }

  close LOG;
}

##################################
#  Calculate Density Matrix
##################################
sub CalcDens {
  my $i;
  my $j;
  my $k;
  my $l;
  my $ph1;
  my $ph2;
  my $ind;
  my $Densold = zeroes($::dim,$::dim);
  my $DensoldB = zeroes($::dim,$::dim);
  
  $::Drmsd = 0.0;
  $Densold = $::Dens;
  $::Dens  = zeroes($::dim,$::dim);
  if($::UHF==1){
    $::DrmsdB = 0.0;
    $DensoldB = $::DensB;
    $::DensB  = zeroes($::dim,$::dim);
  }

  my $eVorder = zeroes($::dim) ;
  my $eVorderB = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;
  if($::UHF==1){
    $eVorderB = qsorti $::EpsB;
  }


  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){
      $ph1  = index2d($::Dens,$i,$j);
      $ph1 .= 0.0;
      for($k=0;$k<$::nocc;$k++){      
        $ind = at($eVorder,$k);
        $ph1 .= $ph1 + index2d($::Coeff,$ind,$i)*index2d($::Coeff,$ind,$j) ;
      }
      $ph2  = index2d($::Dens,$j,$i);
      $ph2 .= $ph1;
      $::Drmsd = $::Drmsd + (at($::Dens,$i,$j)-at($Densold,$i,$j))**2;
    }
  } 
  if($::UHF==1){
    for($i=0;$i<$::dim;$i++){
      for($j=0;$j<=$i;$j++){
        $ph1  = index2d($::DensB,$i,$j);
        $ph1 .= 0.0;
        for($k=0;$k<$::noccB;$k++){
          $ind = at($eVorderB,$k);
          $ph1 .= $ph1 + index2d($::CoeffB,$ind,$i)*index2d($::CoeffB,$ind,$j) ;
        }
        $ph2  = index2d($::DensB,$j,$i);
        $ph2 .= $ph1;
        $::DrmsdB = $::DrmsdB + (at($::DensB,$i,$j)-at($DensoldB,$i,$j))**2;
      }
    } 
  }

  $::Drmsd = sqrt($::Drmsd);
  if($::UHF==1){
    $::DrmsdB = sqrt($::DrmsdB);
  }

  if($::Debug){
    open(LOG,">>","$::name.out");
    print LOG "  Density Matrix:\n";
    print LOG "    ...done \n";
    print LOG $::Dens;
    print LOG "\n";
    print LOG "\n";
    close LOG;
  }
}

##################################
#  Calculate Energy
##################################
sub CalcEnergy {
  my $i;
  my $j;

  $::Eelec  = 0.0;
  $::Etotal = 0.0;

  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<$::dim;$j++){
      $::Eelec = $::Eelec + index2d($::Dens,$i,$j)
                          *(index2d($::Hcore,$i,$j) 
                          + index2d($::Fock,$i,$j));
    }
  }
  if($::UHF==1){
    $::Eelec  = 0.0;
    $::Etotal = 0.0;
    for($i=0;$i<$::dim;$i++){
      for($j=0;$j<$::dim;$j++){
        $::Eelec = $::Eelec + index2d($::Dens,$i,$j)
                            *(index2d($::Hcore,$i,$j)
                            + index2d($::Fock,$i,$j))*0.5
                            + index2d($::DensB,$i,$j)
                            *(index2d($::Hcore,$i,$j)
                            + index2d($::FockB,$i,$j))*0.5;
      }
    }
  }

  $::Etotal =  $::Eelec + $::Enuclear;
  $::dE = $::Etotal - $::Eold;


  open(LOG,">>","$::name.out");
  if($::UHF==0){
    printf LOG "    SCF   %2d   %.8f   %.8f   %.8f   %.8f  \n",$::iSCF,$::Eelec,$::Etotal,$::dE,$::Drmsd ;
    printf     "    SCF   %2d   %.8f   %.8f   %.8f   %.8f  \n",$::iSCF,$::Eelec,$::Etotal,$::dE,$::Drmsd ;
  }elsif($::UHF==1){
    printf LOG "    SCF   %2d   %.8f   %.8f   %.8f   %.8f   %.8f \n",$::iSCF,$::Eelec,$::Etotal,$::dE,$::Drmsd,$::DrmsdB ;
    printf     "    SCF   %2d   %.8f   %.8f   %.8f   %.8f   %.8f \n",$::iSCF,$::Eelec,$::Etotal,$::dE,$::Drmsd,$::DrmsdB ;
  }
  close LOG;
}

##################################
#  Calculate Fock Matrix
##################################
sub CalcFock {
  my $i;
  my $j;
  my $k;
  my $l;
  my $ij;
  my $kl;
  my $ik;
  my $jl;
  my $ijkl;
  my $ikjl;
  my $ph1 = null;
  my $ph2;
  my $ph1B = null;
  my $ph2B;
  my $Jint;
  my $Kint;
  my $Cint;
  my $Dkl;
  my $DklB;
  my $Fockold = zeroes($::dim,$::dim);
  my $FockoldB = zeroes($::dim,$::dim);
  my $betaeff;
  my $weff;

#  open(LOG,">>","$::name.out");
#  print LOG "  Fock Matrix SCF$::iSCF:\n";

  $Fockold .= $::Fock;
  $::Fock   = zeroes($::dim,$::dim);
  $::Fock  .= $::Hcore;
  if($::UHF==1){
    $FockoldB .= $::FockB;
    $::FockB   = zeroes($::dim,$::dim);
    $::FockB  .= $::Hcore;
  }
  $ph1 =0.0;
  $ph2 =0.0;
  if($::UHF==1){
    $ph1B =0.0;
    $ph2B =0.0;
  }


  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<$::dim;$j++){
      $ph1 =0.0;
      $ph2  = 0.0;
      $ph1  = index2d($::Fock,$i,$j);
      if($::UHF==1){
        $ph1B = 0.0;
        $ph2B = 0.0;
        $ph1B = index2d($::FockB,$i,$j);
      }
      
      for($k=0;$k<$::dim;$k++){
        for($l=0;$l<$::dim;$l++){
          $ij   = $i *($i +1)/2 + $j;
          $kl   = $k *($k +1)/2 + $l;
          $ik   = $i *($i +1)/2 + $k;
          $jl   = $j *($j +1)/2 + $l;
          if($j>$i){
            $ij   = $j *($j +1)/2 + $i;
          } 
          if($l>$k){ 
            $kl   = $l *($l +1)/2 + $k;
          }
          if($k>$i){
            $ik   = $k *($k +1)/2 + $i;
          }
          if($l>$j){
            $jl   = $l *($l +1)/2 + $j;
          }
          $ijkl = $ij*($ij+1)/2 +$kl;
          $ikjl = $ik*($ik+1)/2 +$jl;
          if($jl>$ik){
            $ikjl = $jl*($jl+1)/2 +$ik;
          }
          if($kl>$ij){
            $ijkl = $kl*($kl+1)/2 +$ij;
          }         

          $Jint = at($::Iijkl,$ijkl);
          $Kint = at($::Iijkl,$ikjl);
#          $Jint = 2.0;
#         $Kint = 0.0;
          $Dkl  = at($::Dens,$k,$l);
          if($::UHF==1){
            $DklB  = at($::DensB,$k,$l);
          }

#          print "$i $j $k $l $ijkl $Jint $ikjl $Kint  \n";

          $Cint = 1.0;
          if($::IntA==17){
            $Cint = 1.0 + (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0)* 
#                           at($::Dens,$i,$j)*  
                           at($::Sij,$i,$j)**3;
          }
          if($::IntA==19){
            $Cint = 1.0 + (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0)*
#                           at($::Dens,$i,$j)*  
                           at($::Sij,$i,$j);
          }

          if($::UHF==0){
            $ph2 += $Dkl*(2.0*$Jint - $Kint) * $Cint;
          }elsif($::UHF==1){
            $ph2  += ($Dkl*$Jint + $DklB*$Jint - $Dkl*$Kint)  * $Cint;
            $ph2B += ($Dkl*$Jint + $DklB*$Jint - $DklB*$Kint) * $Cint;
          }

#          $ph1 .= $ph1 + $Dkl*(2.0*$Jint - $Kint);
#          print "$ph1\n";

        }
      }
      if($::IntA == 1){
        $Cint = at($::Sij,$i,$j)**3 * ($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0 ;
      }elsif($::IntA == 2){
        $Cint = at($::Sij,$i,$j)**2 * ($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0 ;
      }elsif($::IntA == 3){
        $Cint = at($::Sij,$i,$j)**2 * ($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0
              + at($::Sij,$i,$j)**3 * ($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0;
      }elsif($::IntA == 4){
        $Cint = at($::Sij,$i,$j)**2 * (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0)
              - at($::Sij,$i,$j)**3 * (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0) ;
      }elsif($::IntA == 5){
        $Cint = at($::Sij,$i,$j)**2 * (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0)
              - at($::Sij,$i,$j)**3 * (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0) * 0.5 ;      
      }elsif($::IntA == 6){
        $Cint = at($::Sij,$i,$j)**3 * 
                sqrt((1.0/$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]} * (1.0/$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]});
      }elsif($::IntA == 7){
        $Cint = at($::Sij,$i,$j)**3 *
                sqrt($::Beta{$::bastyp[$i]} * $::Beta{$::bastyp[$j]});
      }elsif($::IntA == 8){
        $Cint = at($::Sij,$i,$j)**3 
                * ((1.0/$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}+(1.0/$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$j]})/2.0 ;
      }elsif($::IntA == 9){
        $Cint = at($::Sij,$i,$j)**3
                * ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
      }elsif($::IntA == 10){
        $Cint = at($::Sij,$i,$j)**2
                * ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
      }elsif($::IntA == 11){
        $betaeff = ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                 +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
        $Cint    = $betaeff*exp(-(at($::Sij,$i,$j)**2) / (2.0*$::width**2) ) ; 
      }elsif($::IntA == 12){
        $betaeff = ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                 +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
        $weff    = ($::w{$::bastyp[$i]} + $::w{$::bastyp[$j]})/2;
        $Cint    = $betaeff*exp(-(at($::Sij,$i,$j)**2) / (2.0*$weff**2) ) ;
      }elsif($::IntA == 13){
        $betaeff = ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                 +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
        $Cint    = $betaeff*(abs(at($::Dens,$i,$j)))**3
      }elsif($::IntA == 14){
        $Cint    = 0.44*abs(at($::Dens,$i,$j))/(7.8+((3.0/(4.0*pi()))*abs(at($::Dens,$i,$j)))**(1/3));
      }elsif($::IntA == 15){
        $betaeff = ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
                 +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
        $Cint    = $betaeff*( at($::Dens,$i,$j)*at($::Sij,$i,$j) )**2 ;
      }else{
        $Cint    = 0.0;
      }

      $ph1 .= $ph1 + $ph2 + $Cint; 
      if($::UHF==1){
        $ph1B .= $ph1B + $ph2B + $Cint;
      }

      if($::IntA==18){
        $Cint = 1.0 - (($::Beta{$::bastyp[$i]}+$::Beta{$::bastyp[$j]})/2.0)*
                       at($::Sij,$i,$j);
        $ph1 .= $ph1*$Cint;
      }

      

#     print "$i $j  $ph2 \n";
    }
  }

  $::Fock .= $::Damp*$::Fock + (1.0-$::Damp)*$Fockold; 
  if($::UHF==1){
    $::FockB .= $::Damp*$::FockB + (1.0-$::Damp)*$FockoldB;
  }

  if($::Debug){
    open(LOG,">>","$::name.out");
    print LOG "  Fock Matrix SCF$::iSCF:\n";

    print LOG "    ...done \n";
    print LOG $::Fock;
    print LOG "\n";
    print LOG "\n";

    close LOG;
  }

}

##################################
#  Calculate Correlation Potential
##################################
sub CalcCorP {
# TODO UHF
  my $i;
  my $j;
  my $k;
  my $l;
  my $Cint;
  my $ph1;
  my $ind;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

  for($i=0;$i<$::nocc;$i++){
#  for($i=0;$i<$::dim;$i++){
    $ind = at($eVorder,$i);
    $ph1 =0.0;
    $ph1  = index2d($::Fock,$ind,$ind);

#    $betaeff = ((5.0-$::Arad[$::basis[$i]])*$::Beta{$::bastyp[$i]}
#             +  (5.0-$::Arad[$::basis[$j]])*$::Beta{$::bastyp[$j]})/2.0 ;
#    $Cint    = $betaeff*( at($::Dens,$i,$j)*at($::Sij,$i,$j) )**2 ;
    $Cint = $::Beta{"1S"}*($ph1) + $::Beta{"2S"}  ;
    $ph1 .= $ph1 + $Cint ;
  }
}


