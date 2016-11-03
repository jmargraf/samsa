package Int1e;
use strict;
use warnings;
use PDL;
#use PDL::LinearAlgebra::Real;

use Exporter 'import';
#use Mol;

our $VERSION = '1.00';
our @EXPORT  = qw(CoreCore Overlap Kinetic CoreEl CalcHcore RunPsi);

##################################
#  Core-Core repulsion
##################################
sub CoreCore{
  my $i;
  my $j;
  my $core1;
  my $core2;

  $::Enuclear=0.0;

  open(LOG,">>","$::name.out");

  print LOG " \n";
  print LOG "  Distance Matrix: \n";
  print LOG " \n";


  for($i=1;$i<=$::natoms;$i++){
    $core1 = Mol::coreq($::atom[$i]);
#    print "core $atom[$i] = $core1 \n";
    for($j=$i+1;$j<=$::natoms;$j++){
      $core2 = Mol::coreq($::atom[$j]);

      $::rAB[$i][$j]=sqrt(
                   ($::xyz[$i][1]-$::xyz[$j][1])**2+
                   ($::xyz[$i][2]-$::xyz[$j][2])**2+
                   ($::xyz[$i][3]-$::xyz[$j][3])**2
                   );

       $::Enuclear = $::Enuclear + ($core1*$core2)/$::rAB[$i][$j];

      print LOG "    rAB $::atom[$i] $i  $::atom[$j] $j  = $::rAB[$i][$j]  \n";
    }
  }

  print LOG " \n";
  print LOG "  E_nuc = $::Enuclear \n";
  print LOG " \n";

  close LOG;

}

##################################
#  Overlap Matrix
##################################
sub Overlap{
  my $i;
  my $j;
  my $scratch;
  my @temp;
  my $k;

  $::Sij = zeroes($::dim,$::dim);

  open(LOG,">>","$::name.out");

  print LOG "  Overlap Matrix: \n";

  open(IN,"<","ints");

  while(<IN>){
    $scratch=$_;
    chomp $scratch;
  # remove whitespaces
    $scratch =~ s/^\s+|\s+$//g ;
    @temp = split / +/, $scratch ;
    if(defined($temp[0]) and ($temp[0]  eq "OVERLAP")){
      print LOG "    ... found ";
      last;
    }
  }

  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){

    # Read Integrals
      $scratch=<IN> ;
      chomp $scratch ;
      @temp =split / +/, $scratch ;
      $k = index2d($::Sij,$j,$i); 
      $k .= $temp[3];
      $k = index2d($::Sij,$i,$j); 
      $k .= $temp[3];
    }
  }

  close(IN);
  if($::Debug){
  print LOG $::Sij ;
  }
  print LOG " \n";

  close LOG;
}


##################################
#  Kinetic Energy Integral Matrix
##################################
sub Kinetic{
  my $i;
  my $j;
  my $scratch;
  my @temp;

  my $k;

  $::Tij = zeroes($::dim,$::dim);

  open(LOG,">>","$::name.out");
  print LOG "  Kinetic Energy Integral Matrix: \n";

  open(IN,"<","ints");

  while(<IN>){
    $scratch=$_;
    chomp $scratch;
  # remove whitespaces
    $scratch =~ s/^\s+|\s+$//g ;
    @temp = split / +/, $scratch ;
    if(defined($temp[0]) and ($temp[0]  eq "KINETIC")){
      print LOG "    ... found ";
      last;
    }
  }

  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){

    # Read Integrals
      $scratch=<IN> ;
      chomp $scratch ;
      @temp =split / +/, $scratch ;
      $k = index2d($::Tij,$j,$i);
      $k .= $temp[3];
      $k = index2d($::Tij,$i,$j);
      $k .= $temp[3];
    }
  }

  close(IN);
  if($::Debug){
    print LOG $::Tij ;
  }
  print LOG " \n";

  close LOG;
}

##################################
#  Core-Electron Attraction Matrix
##################################
sub CoreEl{
  my $i;
  my $j;
  my $scratch;
  my @temp;

  my $k;

  $::Vij = zeroes($::dim,$::dim);
  open(LOG,">>","$::name.out");
  print LOG "  Core-Electron Attraction Matrix: \n";

  open(IN,"<","ints");
    
  while(<IN>){
    $scratch=$_;
    chomp $scratch;
  # remove whitespaces
    $scratch =~ s/^\s+|\s+$//g ;
    @temp = split / +/, $scratch ;
    if(defined($temp[0]) and ($temp[0]  eq "POTENTIAL")){
      print LOG "    ... found ";
      last;
    }
  }

  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){

    # Read Integrals
      $scratch=<IN> ;
      chomp $scratch ;
      @temp =split / +/, $scratch ;
      $k = index2d($::Vij,$j,$i);
      $k .= $temp[3];
      $k = index2d($::Vij,$i,$j);
      $k .= $temp[3];
    }
  }

  close(IN);

  if($::Debug){
    print LOG $::Vij ;
  }
  print LOG " \n";

  close LOG;
}

##################################
#  Core Hamiltonian
##################################
sub CalcHcore{

  $::Hcore = zeroes($::dim,$::dim);
  $::Hcore = $::Tij + $::Vij;

  open(LOG,">>","$::name.out");
  print LOG "  Hcore Matrix: \n";

  if($::Debug){
    print LOG $::Hcore ;
  }
  print LOG "    ... found ";
  print LOG " \n";

  close LOG;
}

##################################
#  Run Psi4 Integral Package
##################################
sub RunPsi{
  my $i;

  open(LOG,">>","$::name.out");
  print LOG "  Running Psi4 Integral Package \n";
  close LOG;

  open(PSI,">","$::name\_PSI4.inp");
  print PSI "import numpy as np \n";
  print PSI " \n";
  print PSI "memory 250 mb \n";
  print PSI " \n";
  print PSI "set globals {  \n";
  if($::basis_set eq "min"){
    print PSI "  basis sto-3g \n";
  }elsif($::basis_set eq "svp"){
    print PSI "  basis def2-SVP \n";
  }
  print PSI "  units bohr  \n";
  print PSI "}\n";
  print PSI " \n";
  print PSI "molecule $::name {\n";
  for($i=1;$i<=$::natoms;$i++){
    print PSI "  $::atom[$i] $::xyz[$i][1] $::xyz[$i][2] $::xyz[$i][3] \n";
  }
  print PSI "}\n";
  print PSI " \n";
  print PSI "ref_wfn = psi4.new_wavefunction($::name, psi4.get_global_option('BASIS')) \n";
  print PSI "mints = MintsHelper(ref_wfn.basisset()) \n";
  print PSI " \n";
  print PSI "S = mints.ao_overlap() \n";
  print PSI "np_S = np.array(S) \n";
  print PSI "print(\"OVERLAP\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $::dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %f' % (i,j,np_S[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "T = mints.ao_kinetic() \n";
  print PSI "np_T = np.array(T) \n";
  print PSI "print(\"KINETIC\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $::dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %f' % (i,j,np_T[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "V = mints.ao_potential() \n";
  print PSI "np_V = np.array(V) \n";
  print PSI "print(\"POTENTIAL\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "while (i < $::dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    print '  %d  %d    %f' % (i,j,np_V[i][j]) \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI "ERI = mints.ao_eri() \n";
  print PSI "np_ERI = np.array(ERI) \n";
  print PSI "print(\"REPULSION\") \n";
  print PSI " \n";
  print PSI "i=0 \n";
  print PSI "k=0 \n";
  print PSI "count = 0 \n";
  print PSI " \n";
  print PSI "while (i < $::dim): \n";
  print PSI "  j=0 \n";
  print PSI "  while (j <= i): \n";
  print PSI "    k=0 \n";
  print PSI "    while (k < $::dim): \n";
  print PSI "      l=0 \n";
  print PSI "      while( l<=k ): \n";
  print PSI "        ij = i*(i+1)/2+j \n";
  print PSI "        kl = k*(k+1)/2+l \n";
  print PSI "        if ij>=kl: \n";
  print PSI "          print '  %d  %d  %d  %d   %f' % (i,j,k,l,np_ERI[i][j][k][l]) \n";
  print PSI "          count = count + 1 \n";
  print PSI "        l=l+1 \n";
  print PSI "      k=k+1 \n";
  print PSI "    j = j+1 \n";
  print PSI "  i = i+1 \n";
  print PSI " \n";
  print PSI " \n";

  close PSI;

  `Psi4 $::name\_PSI4.inp $::name\_PSI4.out > ints`;

  open(LOG,">>","$::name.out");
  print LOG "    ... done  \n";
  close LOG;
}


