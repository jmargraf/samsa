package Int2e;
use strict;
use warnings;
use PDL;
#use PDL::LinearAlgebra::Real;

use Exporter 'import';
#use Mol;

our $VERSION = '1.00';
our @EXPORT  = qw(Repulsion);

##################################
#  Electron-Electron repulsion
##################################
sub Repulsion{
  my $i;
  my $j;
  my $k;
  my $l;
  my $ij;
  my $kl;
  my $ijkl;

  my $scratch;
  my @temp;
  my $ph;

  $::Iijkl = zeroes($::dim2);

  open(LOG,">>","$::name.out");

  print LOG "  Repulsion Integrals: \n";

  open(IN,"<","ints");
    
  while(<IN>){
    $scratch=$_;
    chomp $scratch;
  # remove whitespaces
    $scratch =~ s/^\s+|\s+$//g ;
    @temp = split / +/, $scratch ;
    if(defined($temp[0]) and ($temp[0]  eq "REPULSION")){
      print LOG "    ... found \n";
      last;
    }
  }

#  for(my $count=1;$count<=228;$count++){
#    $scratch = <IN>;
#    chomp $scratch;
#    $scratch =~ s/^\s+|\s+$//g ;
#    @temp = split / +/, $scratch ;
#    $i = $temp[0]-1;
#    $j = $temp[1]-1;
#    $k = $temp[2]-1;
#    $l = $temp[3]-1;
#    $ij = $i*($i+1)/2 + $j;
#    $kl = $k*($k+1)/2 + $l;
#    $ijkl = $ij*($ij+1)/2+$kl;
#    $ph = index($::Iijkl,$ijkl);
#    $ph .= $temp[4];
#  }


  for($i=0;$i<$::dim;$i++){
    for($j=0;$j<=$i;$j++){
      for($k=0;$k<$::dim;$k++){
        for($l=0;$l<=$k;$l++){
          $ij = $i*($i+1)/2 + $j;
          $kl = $k*($k+1)/2 + $l;
          if($ij>=$kl){
            $ijkl = $ij*($ij+1)/2+$kl;
            $scratch=<IN> ;
            chomp $scratch ;
            @temp =split / +/, $scratch ;
            $ph = index($::Iijkl,$ijkl);
            $ph .= $temp[5];
          }
        }
      }
    }
  }


  close(IN);

#  print LOG $::Iijkl ;
#  print LOG " \n";
#  print LOG " \n";

  close LOG;
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
}
