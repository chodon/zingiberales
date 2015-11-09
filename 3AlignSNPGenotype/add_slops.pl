#!/usr/bin/perl -w

use warnings;
use strict;


my @lib = ('bdjq', 'cost', 'curc', 'jnub', 'knkv', 'lskk', 'musa', 'tzns', 'uoel');


foreach my $lib (@lib){
my $dir = '/global/scratch/cdspecht/chodon/zingiberales/sep_exons_full_cds/extended/';
# name an out file
my $newbed = $dir . $lib . '.slopped.bed';
open(BED, ">$newbed");
#name infile
my $sloptmp = $dir. $lib . '.bedsslop.tmp';
open(IN, "<$sloptmp");


# loop through temp slop file
my @slop = <IN> ;  
   
         
  # logic is to add or subtract slop from bed if there is sequence available 
   

    foreach my $line (@slop) {
	chomp $line;
          my @tmp = split(/\t/,$line);
          
          my $c1 = $tmp[0];
          my $c2 = $tmp[1];
          my $c3 = $tmp[2];
          my $c4 = $tmp[3];
          my $c5 = $tmp[4];
          my $c6 = $tmp[5];
          my $c7 = $tmp[6];         
    
          my $newstart;
          my $newend;

              if ($c6 < ($c3 - $c1)){
                  $newstart = ($c3 - $c1);
                  }       
             
              else { 
              $newstart = $c6;
	  }

              if ($c7 > ($c4 + $c1)){
                  $newend = ($c4 + $c1);
              }  
              else { 
              $newend = $c7;
	  }
              

 print BED $c2 . "\t" . $newstart . "\t" . $newend . "\n";
           
	      
	      }



   # print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
   # print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
    # process each seq
    # print $seq->id . ' = '.$seq->seq()."\n";

close BED;
close IN;
}


