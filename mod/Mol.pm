package Mol;
use strict;
use warnings;

use Exporter 'import';
use Constants;

our $VERSION = '1.00';
our @EXPORT  = qw(readxyz coreq dimensions);

##################################
#  Read Geometry
##################################
sub readxyz{
  my $scratch;
  my @temp;
  my $m;

  open(IN,"<","$::name.xyz");

# Read number of Atoms to $natoms
  $::natoms=<IN> ;
  chomp $::natoms ;

# Read Charge + Mult
  $scratch=<IN> ;
  chomp $scratch ;
  @temp =split / +/, $scratch ;
  $::charge = $temp[0];
  $::mult   = $temp[1];

# Read  $natoms coordinates to xyz
  for($m=1;$m<=$::natoms;$m++){
    $scratch=<IN> ;
    chomp $scratch ;
    @temp =split / +/, $scratch ;
    $::atom[$m]=  $temp[0] ;
    $::xyz[$m][1]=$temp[1] *$Constants::ang2bohr ;
    $::xyz[$m][2]=$temp[2] *$Constants::ang2bohr ;
    $::xyz[$m][3]=$temp[3] *$Constants::ang2bohr ;
    $::Arad[$m] = &atomrad($::atom[$m]) ;
  }

  close IN;                 

}

##################################
#  Get Core Charge
##################################
sub coreq{
#  print "$_[0] \n";
  $_ = $_[0];

  if($_ eq "H"){
  return(1)}
  elsif($_ eq "He"){
  return(2)}
  elsif($_ eq "Li"){
  return(3)}
  elsif($_ eq "Be"){
  return(4)}
  elsif($_ eq "B"){
  return(5)}
  elsif($_ eq "C"){
  return(6)}
  elsif($_ eq "N"){
  return(7)}
  elsif($_ eq "O"){
  return(8)}
  elsif($_ eq "F"){
  return(9)}
  elsif($_ eq "Ne"){
  return(6)}
}

##################################
#  Get Atomic Radius
##################################
sub atomrad{
  $_ = $_[0];

  if($_ eq "H"){
  return(1.0015547741)}
  elsif($_ eq "He"){
  return(0.585815056)}
  elsif($_ eq "Li"){
  return(3.155842401)}
  elsif($_ eq "Be"){
  return(2.116493107)}
  elsif($_ eq "B"){
  return(1.64406161)}
  elsif($_ eq "C"){
  return(1.266116412)}
  elsif($_ eq "N"){
  return(1.058246554)}
  elsif($_ eq "O"){
  return(0.907068475)}
  elsif($_ eq "F"){
  return(0.793684915)}
  elsif($_ eq "Ne"){
  return(0.718095876)}

}

##################################
#  Get Basis Dimensions and Info
##################################
sub dimensions{
  my $i;
  my $j;
  my $k;
  my $l;
  my $ij;
  my $kl;
  my $ijkl;
  my $core;

# Dimension of 1e Matrices
  $::dim = 0;

  if($::basis_set eq "min"){
    for($i=1;$i<=$::natoms;$i++){
      if($::atom[$i] eq "H" or $::atom[$i] eq "He"){
        $::dim=$::dim+1;
      }else{
        $::dim=$::dim+5;
      }
    }
  }elsif($::basis_set eq "svp"){
    for($i=1;$i<=$::natoms;$i++){
      if($::atom[$i] eq "H" or $::atom[$i] eq "He"){
        $::dim=$::dim+5;
      }else{
        $::dim=$::dim+14;
      }
    }  
  }

# Number of unique 2e Integrals
  $::dim2 = 0;
  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){
      for($k=0;$k<$::dim;$k++){
        for($l=0;$l<=$k;$l++){
          $ij = $i*($i+1)/2 + $j;
          $kl = $k*($k+1)/2 + $l;
          if($ij>=$kl){
            $::dim2++;
            $ijkl = $ij*($ij+1)/2+$kl;
#            print "  $i     $j     $k     $l   $ijkl \n"; 
          }
        } 
      }
    }
  }

# Number of electrons
  $::nel = 0;
  $::nocc = 0;
  for($i=1;$i<=$::natoms;$i++){
    $core = coreq($::atom[$i]);
    $::nel = $::nel + $core;  
  }

  $::nocc = $::nel/2;

  $j = -1;
  open(LOG,">>","$::name.out");

  if($::basis_set eq "min"){
    for($i=1;$i<=$::natoms;$i++){
#   $basis[$i] = atom number ;
      if($::atom[$i] eq "H" or $::atom[$i] eq "He"){
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "1S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";      
      }else{
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "1S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
      }
    }
  }elsif($::basis_set eq "svp"){
    for($i=1;$i<=$::natoms;$i++){
#     $basis[$i] = atom number ;
      if($::atom[$i] eq "H" or $::atom[$i] eq "He"){
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "1S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "1S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
      }else{
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "1S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2S";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "2P";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol D";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol D";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol D";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol D";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
        $j++;
        $::basis[$j] = $i;
        $::bastyp[$j] = "Pol D";
        print LOG "    Basisfunction $j is $::atom[$i] $::bastyp[$j] \n";
      }
    }
  }
  print LOG "\n";
  close LOG;

  open(LOG,">>","$::name.out");
  print LOG "  1e Matrix dimension = $::dim \n";
  print LOG "  No of unique 2e int = $::dim2 \n";
  print LOG "  No of electrons     = $::nel \n";
  print LOG "  No of occ orbitals  = $::nocc \n";
  print LOG "\n";
  close LOG;
}

