#! /usr/bin/perl
use strict;
use warnings;

use PDL;

use FindBin;                  
use lib "$FindBin::Bin/mod";  

use Mol;
use Constants;
use Int1e;
use Int2e;
use LinAl;
use Prop;
use MP2;
use CCSD;

####### Declaration of Variables
#
# calculation details
  our $name;
  our $natoms;
  our $charge;
  our $mult;
  our @atom;
  our @xyz;
  our @rAB;
  our $dim;
  our $dim2;
  our $nel;
  our $nocc;
  our @basis;
  our @bastyp;
  our @basat;
  our @PopQ;
  our $Debug = 0;
  our $Res =1;
  our @Arad;

# counting variables
  my  $i;
  my  $j;
  our $iSCF=0;
  my $time;
  my $starttime;
  my $endtime;

# Energies
  our $Enuclear;
  our $Eelec;
  our $Etotal;
  our $Eold;
  our $dE=0.0;
  our $Drmsd;

# 1el Integrals
  our $Sij;
  our $Tij;
  our $Vij;
  our $Hcore;

# S eigenvectors and values
  our $S12;
  our $S_vec;
  our $S_val;

# 2el Integrals
  our $Iijkl;
  our @MOijkl;
  our @SpinInts;

# SCF Matrices 
  our @f_s;
  our $Fock;
  our $Coeff;
  our $Dens;
  our $Eps;
  our $ES;

# SCF options
  my $MaxIter = 1000;
  my $Econv   = 1.0e-9;
  my $Dconv   = 1.0e-9;
  our $Damp   = 0.5;
  our $DoDamp = "on";
  our $EHT    = "false";
  our %Beta ;
  our %w;
  our $IntA;
  our $width;

  our $basis_set = "min";
#
####### End Declaration of Variables

##################################
#  Main Program  
##################################
# program starts

  $name = $ARGV[0];
  chomp $name;
  open(LOG,">","$name.out");

  if($DoDamp eq "off"){
    $Damp = 1.0;
  }

  if(@ARGV > 1) {
    $IntA       = $ARGV[1];
    $Beta{"1S"} = $ARGV[2];
    $Beta{"2S"} = $ARGV[3];
    $Beta{"2P"} = $ARGV[4];
    if($IntA==11){
      $width    = $ARGV[5];
    }
    if($IntA==12){
      $w{"1S"}  = $ARGV[5];
      $w{"2S"}  = $ARGV[6];
      $w{"2P"}  = $ARGV[7];
    }
  }else{
    $IntA = 666;
    $Beta{"1S"} = 0.0;
    $Beta{"2S"} = 0.0;
    $Beta{"2P"} = 0.0;
  }

#  $starttime = time;

# Basic Output
  print LOG "                                \n";
  print LOG "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  print LOG "  SAMSA - a Semiempirical and Ab Initio Molecular Simulation Application \n";
  print LOG "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  print LOG "\n";
  print LOG "  (c) 2016 Johannes Margraf\n";
  print LOG "\n";
  print LOG "  SCF calculation for $name  \n";
  print LOG "\n";

  print  LOG "\n";
  printf LOG " b(1s)=%1.4f b(2s)=%1.4f b(2p)=%1.4f Integral Approximation: %s\n",$Beta{"1S"},$Beta{"2S"},$Beta{"2P"},$IntA ;
  print  LOG "\n";


# Read input
  readxyz();

# Print geometry 
  print LOG "  Geometry (Bohr) :\n";
  print LOG "\n";
  for($i=1;$i<=$natoms;$i++){
    printf LOG "    %2s    %2.6f  %2.6f  %2.6f (at rad. = %3.6f) \n",$atom[$i],$xyz[$i][1],$xyz[$i][2],$xyz[$i][3],$Arad[$i];
  }
  print LOG "\n";
  print LOG "  1A = $Constants::ang2bohr Bohr  \n";

  close LOG;

# Calculate Distance Matrix and Core-Core Repulsion
  CoreCore();

# Calculate matrix dimensions 
  dimensions(); 

# Run Psi integral package
  RunPsi();

# Calculate overlap matrix
  Overlap();

# Kinetic energy integrals
  Kinetic();

# Nuclear attraction integrals
  CoreEl();

# Calculate Hcore
  CalcHcore();

# Get repulsion integrals
  Repulsion();

# Diagonalize Overlap Matrix
  DiaS();

# Build S-12 
  BuildS12();

# Run SCF Loop
  SCF();

# Run Population Analysis
  Mulliken(); 

# Print Orbitals
  PrintEigen();

# Run MO integral transformation and calculate MP2 energy
  $starttime = time;
  MOints();
  EMP2();
  $endtime = time;
  $time = $endtime - $starttime;

  print "    MP2 time = $time s\n";
  open(LOG,">>","$name.out");
  print LOG "    MP2 time = $time s\n";
  print LOG "\n";
  close LOG;

# Transform MOs to Spin Integrals
  SpinInts();

# Calculate PT2 IP corrections
  dEMBPT();
#  dEMP2();

  FockSpin();

#  ECCSD();

  ENPT2();

##################################
#  End Main Program   
##################################


##################################
#  Subroutines:
##################################

##################################
#  SCF Routine
##################################
sub SCF{

  $starttime = time;

# Make initial Guess and Diagonalize Fockian
  DiaF();

# Calculate Density Matrix
  CalcDens();

  open(LOG,">>","$name.out");
  printf LOG "    SCF   %2s   %12s   %12s   %10s   %10s  \n","It","E_elec","E_tot","dE","Drms" ;
  printf     "    SCF   %2s   %12s   %12s   %10s   %10s  \n","It","E_elec","E_tot","dE","Drms" ;
  close LOG;

# Calculate Energy
  CalcEnergy();

  $Eold = 0.0;

# Start SCF
  for($iSCF=1;$iSCF<=$MaxIter;$iSCF++){
    
    $Eold = $Etotal;

    CalcFock();
    if($IntA==16){
      FtoMO();
      CalcCorP();
      FtoAO();
    }
    DiaF();
    CalcDens();
    CalcEnergy();

    if(abs($dE)<1.0e-3 and $Damp != 1.00 and $Damp != 0.7 and $iSCF>3 and $DoDamp eq "on"){
      $Damp = 0.7;
      open(LOG,">>","$name.out");
      print LOG "    Damping set to 0.7 \n";
      print     "    Damping set to 0.7 \n";
      close LOG;
    }

    if(abs($dE)<1.0e-5 and $Damp != 1.00 and  $DoDamp eq "on"){
      $Damp = 1.00;
      open(LOG,">>","$name.out");
      print LOG "    Damping turned off \n";
      print     "    Damping turned off \n";
      close LOG;
    }


    if(abs($dE)<$Econv and abs($Drmsd)<$Dconv){
      open(LOG,">>","$name.out");
      print     "\n";
      print     "    Self-consistency achieved! \n";
      print     "\n";
      print LOG "\n";
      print LOG "    Self-consistency achieved! \n";
      print LOG "\n";
      close LOG;
      last;
    }
  }
# End SCF Loop

  open(LOG,">>","$name.out");
  print LOG "    E_final = $Etotal \n";
  print LOG "\n";
  close LOG;

  $endtime = time;
  $time = $endtime - $starttime;
  open(LOG,">>","$name.out");
  print LOG "    SCF time = $time s\n";
  print LOG "\n";
  close LOG;
  print "    SCF time = $time s\n";

}

##################################
#  End Subroutines
##################################
