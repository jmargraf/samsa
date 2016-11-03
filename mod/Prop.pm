package Prop;
use strict;
use warnings;
use PDL;

use Exporter 'import';

our $VERSION = '1.00';
our @EXPORT  = qw(Mulliken PrintEigen);

##################################
#  Population Analysis
##################################
sub Mulliken{
  my $i;
  my $j;
  my $temp;

  for($i=1;$i<=$::natoms;$i++){
    $::PopQ[$i] = 0.0;
  }

  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<$::dim;$j++){
#      print " Basisfunction $i is on atom $::basis[$i]  \n";
      $::PopQ[$::basis[$i]] += -at($::Dens,$i,$j)*at($::Sij,$i,$j) ;
      $::PopQ[$::basis[$j]] += -at($::Dens,$i,$j)*at($::Sij,$i,$j) ;
    }
  }

  open(LOG,">>","$::name.out");
  print LOG "  Mulliken Charges: \n";
  print LOG "\n";

  for($i=1;$i<=$::natoms;$i++){
    $::PopQ[$i] += Mol::coreq($::atom[$i])  ;
    print LOG "    $::atom[$i]$i   $::PopQ[$i]  \n";
  }

#  for($i=0;$i<$::dim;$i++){
#    for($j=0;$j<$::dim;$j++){
#      $temp = at($::Dens,$i,$j)*at($::Sij,$i,$j);
#      print LOG "$temp ";
#    }
#    print LOG "\n";
#  }

  print LOG "\n";
  close LOG;

}

##################################
#  Print Orbital Eigenvalues
##################################
sub PrintEigen{
  my $i;
  my $temp;
  my $temp2;

  my $eVorder = zeroes($::dim) ;
  $eVorder = qsorti $::Eps;

  open(LOG,">>","$::name.out");
  print LOG "  Orbital Energies: \n";
  print LOG "    n      E[Ha]        E[eV] \n";
  print LOG "\n";

  for($i=0;$i<$::dim;$i++){
    $temp  = at($::Eps,at($eVorder,$i));
    $temp2 = $temp*$Constants::Ha2eV;
    if($i==($::nocc-1)){
    printf LOG "    %3s    %4.6f    %4.6f  HOMO \n",$i,$temp,$temp2;
    }else{
    printf LOG "    %3s    %4.6f    %4.6f \n",$i,$temp,$temp2;
    }
  }

  print LOG "\n";
  close LOG;

#  open(LOG,">>","$::name.out");
#  print LOG "  Orbital Energies (Internal Order): \n";
#  print LOG "    n      E[Ha]        E[eV] \n";
#  print LOG "\n";

#  for($i=0;$i<$::dim;$i++){
#    $temp  = at($::Eps,$i);
#    $temp2 = $temp*$Constants::Ha2eV;
#    if($i==($::nocc-1)){
#    printf LOG "    %3s    %4.6f    %4.6f  HOMO \n",$i,$temp,$temp2;
#    }else{
#    printf LOG "    %3s    %4.6f    %4.6f \n",$i,$temp,$temp2;
#    }
#  }

#  print LOG "\n";
#  close LOG;


}
