#!/usr/bin/perl -w

use warnings;
use strict;
use lib "/global/home/users/cdspecht/chodon/programs/local_perl/BioPerl-1.6.1";
use Bio::SeqIO;
use Getopt::Std;


##define libraries, home directory and extension to home directory plus library where bams and vcfs are
  my %opts = (H=>undef, G=>undef); 
  getopts('H:G:', \%opts);

my @lib = qw(CS38 CS02 CS06 CS07 CS09 CS12 CS13 CS14 CS15 CS16 CS17 CS18 CS19 CS20 CS21 CS22 CS23 CS24 CS25 CS27 CS28 CS29 CS30 CS31 CS32 CS33 CS34 CS35 CS36 CS37 CS39 CS40 CS41 CS42 CS43 CS44 CS45 CS46 CS47 CS48 CS49 CS50 CS51 CS52 CS53 CS54 CS55 CS57 CS58 CS69 CS70 CS71 CS72 BAAU BRUD BYPY CAGR CHSE DAS2 DASP DITH HAMA HWUP KIAU NSPR PPQR QGLJ TRPA TZNS JQCX OYLU KNKV JNUB LSKK UOEL BDJQ COST);

#my @lib = qw(CS72 CS30);

  my $home = $opts{H};

  my $vcfs = $opts{G};

my $selectseqs = '/global/home/users/cdspecht/bin/selectseqs.pl'; 

foreach my $lib (@lib){

my $dir = $home . $lib . $vcfs;

# name a reference sequence file
my $refs = $home . $lib . '/' . 'map2' . '/' . $lib . '.clean4.fa';

## this ignores deletions and insertions ##

# name an out file
my $hap0fil = $dir . $lib . '.haplo0O.fasta';
open(HAP0, ">$hap0fil");
my $hap1fil = $dir . $lib . '.haplo1O.fasta';
open(HAP1, ">$hap1fil");
my $indellist = $dir . $lib . 'indel.list';
open(INDEL, ">$indellist");
# name vcf file
my $rawvcf = $dir. $lib . '.filter.vcf';
my $processvcf = $dir. $lib . '.shortO.vcf';
system("grep -v '\#' $rawvcf \| awk \'length(\$4) < 2\' > $processvcf");
my $vcffil = $dir. $lib . '.shortO.vcf';
open(VCF, "<$vcffil");
my @vcf = <VCF>;

# name coverage file

my $coverage =  $home . $lib . '/' . 'map2' . '/' . 'contigs2' . '/' . $lib . '_m2n240_4.cov';    


# handle IO with bioperl
my $seqio = Bio::SeqIO-> new(
                             -file     => $refs,
                             -format => 'FASTA',
                             );

# loop through reference sequences
while( my $seq = $seqio->next_seq ) {
  
    my $name = $seq->id;

    # split into two haplotypes
    my $h0name = $name;
    my $h1name = $name;
    my $h0seq = $seq->seq();
    my $h1seq = $seq->seq();


    # get vcs for contig
    my @refvcs = ();

    my $h0indlen = 0;
    my $h1indlen = 0;

    foreach my $vcline (@vcf)   
    {
          my @tmp = split(/\t/,$vcline);
          my $prea0; my $prea1; my $join; my $newletter;
          my $a0; my $a1;     
          if ($name eq $tmp[0])
          {
              
## for coverage less than 10, turn to Ns unless qual, then leave ref, or homozygous 1.0 then leave alt

              if ($tmp[6] =~ /dp10/){

                 if ($tmp[6] !=~ /altdef|qual/){
                  $a0 = 'N';
                  $a1 = 'N';
                  }
  
	         if ($tmp[6] =~ m/altdef/) {
                  $a0 = $tmp[4];
                  $a1 = $tmp[4];
                  }

	         if ($tmp[6] !=~ m/altdef/){
		      if ($tmp[6] =~ m/qual/) {
                      $a0 = $tmp[3];
                      $a1 = $tmp[3];
	      }
                   }
	     }
## for coverage less than 20, turn to IUPAC, using consensus later, unless qual, then leave ref, or homozygous 1.0, then leave alt
              
              if ($tmp[6] =~ m/dp20/) {
 
                  if ($tmp[6] !=~ m/altdef|qual/){
                  $prea0 = $tmp[3];
                  $prea1 = $tmp[4];
                  $join = "$prea0$prea1";
                    unless (length($join) > 2) {
                    for ($join ){ 
                      if (/TC|CT/) {
		      $a0 = 'Y';
		      $a1 = 'Y';
                      }
		      if (/AG|GA/){
                      $a0 = 'R';
                      $a1 = 'R';
		      }
		      if (/CA|AC/){
                      $a0 = 'M';
                      $a1 = 'M';
		      }
		      if (/TG|GT/){
                      $a0 = 'K';
                      $a1 = 'K';
		      }
		      if (/TA|AT/){
                      $a0 = 'W';
                      $a1 = 'W';
		      }
       		      if (/GC|CG/){
                      $a0 = 'S';
                      $a1 = 'S';
		      }
                    print "$join for $lib at $tmp[0] position $tmp[1] is $a0 or $a1 \n";
                    }
           	  }
		  
	          }
                  if ($tmp[6] =~  m/altdef/) {
                  $a0 = $tmp[4];
                  $a1 = $tmp[4];
                  }
              
                  if ($tmp[6] =~ m/qual/){
                    
                    if ($tmp[6] =~ m/hethigh/) { 
		        $a0 = $tmp[3];
                        $a1 = $tmp[3];
		    }
                    if ($tmp[6] =~ m/hetlow|altpos|homhigh|homlow|altdef/) {
		       $a0 = $tmp[4];
		       $a1 = $tmp[4];
	       	    }
                    else {
		       $a0 = $tmp[3];
		       $a1 = $tmp[3];
                    }
		 }
	      }
 ## for qual, use 50 percent rule
               
	      if ($tmp[6] =~ m/qual/) {
                 
                  if ($tmp[6] !=~ m/dp20|dp10/) {
                  
                      if ($tmp[6] =~ m/hethigh/) {
		      $a0 = $tmp[3];
		      $a1 = $tmp[3];
		      }
                  
                      if ($tmp[6] =~ m/hetlow|altpos|altdef|homhigh|homlow/) {
		      $a0 = $tmp[4];
		      $a1 = $tmp[4];
		      }
                  
                      else {
                      $a0 = $tmp[3];
                      $a1 = $tmp[3];
                       }
                   }              
	      }
## no dp and no qual

	      if ($tmp[6] !=~ m/qual|dp20|dp10/) {

                  if ($tmp[6] =~ m/hethigh/) {
                      $a0 = $tmp[3];
                      $a1 = $tmp[3];
	              }

                  if ($tmp[6] =~ m/hetlow|homhigh|altdef/) {
                      $a0 = $tmp[4];
                      $a1 = $tmp[4];
		      }
		  
                  if ($tmp[6] =~ m/altpos|PASS|homlow/) {
                      $prea0 = $tmp[3];
		      $prea1 = $tmp[4];
		      $join = "$prea0$prea1";
                      unless (length($join) > 2) { 
		      for ($join ){
			  if (/TC|CT/) {
			      $a0 = 'Y';
			      $a1 = 'Y';
			  }
			  if (/AG|GA/){
			      $a0 = 'R';
			      $a1 = 'R';
			  }
			  if (/CA|AC/){
			      $a0 = 'M';
			      $a1 = 'M';
			  }
			  if (/TG|GT/){
			      $a0 = 'K';
			      $a1 = 'K';
			  }
			  if (/TA|AT/){
			      $a0 = 'W';
			      $a1 = 'W';
			  }
			  if (/GC|CG/){
			      $a0 = 'S';
			      $a1 = 'S';
			  }
			  print "$join for $lib at $tmp[0] position $tmp[1] is $a0 or $a1\n";
		      }
		  }
                    
	          }
	     } 

## deal with if two alternate bases

	      if ($tmp[4] =~ m/,/) {
		  
                  $a0 = 'N';
                  $a1 = 'N';
              
              print "$tmp[0] for $lib has more than one alternate base";
                  }  
               
### deal with insertions

	      if (length($tmp[4]) > 1) {
                  $a0 = $tmp[3];
                  $a1 = $tmp[4];
		  print "$tmp[0] for $lib has an indel\n";
	      }

### deal with phasing inconsistent sites !!! note that I am using ambiguity code 'D' so I can count these sites, but the sites are not D, will be converted to N's later
              unless ($tmp[6] =~ m/qual|dp10|dp20/) { 
              for ($tmp[6] =~ m/PASS|homlow|altpos/) {
		  if ($tmp[7] =~ m/PhasingInconsistent/) {
		      $a0 = 'D';
		      $a1 = 'D';
		      print "$tmp[0] for $lib at position $tmp[1] is phasing inconsistent\n";
                  }
	      }
	  }
#this does indel counting


              substr($h0seq, $tmp[1]+$h0indlen-1, 1) = $a0;
              substr($h1seq, $tmp[1]+$h1indlen-1, 1) = $a1;

              $h0indlen = $h0indlen + length($a0)-1;
              $h1indlen = $h1indlen + length($a1)-1;
              print $tmp[0] ."\t" . $tmp[1] . "\t" . $h0indlen ."\t" . $h1indlen ."\n";
	      print INDEL $tmp[0] . "\t" . $tmp[1] . "\t" . $h0indlen . "\n";
          } 
      }



    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
    # process each seq
    # print $seq->id . ' = '.$seq->seq()."\n";
}
close HAP0;
close HAP1;
close INDEL;
close VCF;

##now filter out high coverage sites (max coverage = greater than 4x average - this is too conservative, some nuclear genes will be filtered, should fix this by creating a blast subroutine for finding the chloroplast and mitochondrial genes later)

my $tempdir = $home . $lib . '/map2/contigs2/temp/';
my $highnumber =  qx{awk \'\{sum+=\$3\} END \{print 4*(sum/NR)\}\' $coverage};
print $highnumber;
my $highContigs = qx{awk \'\$3>$highnumber\' $coverage \| awk \'\{print \$1\}\' \| sort \| uniq};
my $highlist = $tempdir . 'high.list';
open (OUTL, ">$highlist");
print OUTL $highContigs;
close OUTL;

## now make under 1x coverage bedfile and mask fasta for under 5x coverage, then cut out under 1x coverage from masked fasta                                                          
my $tempcov0a = $coverage . '0a';
my $tempcov0b = $coverage . '0b';
my $tempcov5a = $coverage . '5a';
my $tempcov5b = $coverage . '5b';
my $over4fasta = $home . $lib . '/map2/contigs2/' . $lib . '.over4.fasta';
my $masked = $home . $lib . '/map2/contigs2/' . $lib . '.covMask.fasta';
my $seqlist = $home . $lib . '/map2/contigs2/' . $lib . '.seqlist';
#make 0 cov bed file
system("awk \'\$3 > 0\' $coverage \| awk '\{print \$1\"\\t\"\(\$2-1\)\"\\t\"\$2\}' > $tempcov0a");
system("mergeBed -i $tempcov0a > $tempcov0b");
#make less than 5 cov bedfile
system("awk \'\$3 < 5\' $coverage \| awk '\{print \$1\"\\t\"\(\$2-1\)\"\\t\"\$2\}' > $tempcov5a");
system("mergeBed -i $tempcov5a > $tempcov5b");
#mask with less than 5x cov
system("maskFastaFromBed -fi $hap0fil -bed $tempcov5b -fo $over4fasta");
#remove 0x cov
#select sequences in 0coverage file - if they are in there, they have some coverage at all
system("awk \'\{print \$1\}\' $tempcov0b | sort | uniq > $seqlist");
system("$selectseqs -f $seqlist $over4fasta > $masked");
system("rm $tempcov0a $tempcov5a");


##count phasing inconsistent sites
#### count ambiguities

my $fastaIn = $masked;
my $ambiguityOut = $tempdir . 'ambCount.txt';
my $phaseInconOut = $tempdir . 'phaseInconCount.txt';
open (OUT1, ">$ambiguityOut");
open (OUT2, ">$phaseInconOut");
my $in = Bio::SeqIO->new(-file => $fastaIn, -format => 'fasta');
while (my $seqobj = $in->next_seq){
    my $id = $seqobj->id;
    my $seq = $seqobj->seq;
    my $length = $seqobj->length;
    my $countA = 0;
    my $countP = 0;
    for (my $i = 0; $i < $length; $i++) {
        my $sub = substr($seq,$i,1);
        if ($sub =~ /Y|M|R|S|W|K/i) {
            $countA++;
	}    
        if ($sub =~ /D/i) {
            $countP++;
        }
    }
    my $ambiguities = sprintf("%.1f",$countA * 100 /$length);
    print OUT1 $id,"\t",$ambiguities,"\n";
    my $phasingInconsistent = sprintf("%.1f",$countP * 100 /$length);
    print OUT2 $id,"\t",$phasingInconsistent,"\n";
}
close OUT1;
close OUT2;
##change D's from phaseInconsistent counting to Ns - note this should work even if contig names contain the letter D
my $maskedB = $home . $lib . '/map2/contigs2/' . $lib . '.covMaskB.fasta';
my $replaceDinSed = '/^[ACTGNYMRSWKactgnymrswk]/';
system("sed \'$replaceDinSed s\/D\/N\/g\' $masked > $maskedB");

## now remove low coverage (approximate average less than 11x calculated with number of reads * 100, although this is not exactly the coverage b/c merged reads and adapter trimming, but close enough for now

my $bamtobed = '/global/home/users/cdspecht/bin/bedtools2-2.20.1/bin/bamToBed';
my $coverageBed = '/global/home/users/cdspecht/bin/bedtools2-2.20.1/bin/coverageBed';
my $bamForCov = $home . $lib . '/map2/contigs2/' . $lib . '_m2n240_4r.bam';
my $lowlist = $tempdir . 'low.list';

system("$bamtobed -i $bamForCov \| $coverageBed -a - -b $tempcov0b \| awk \'\(\(\$4*100\)\/\$5\) < 11\' \| awk \'\{print \$1\":\"\$2\"-\"\$3\}\' > $lowlist");

###/global/home/users/cdspecht/bin/bedtools2-2.20.1/bin/bamToBed -i CS02_m2n240_4r.bam | /global/home/users/cdspecht/bin/bedtools2-2.20.1/bin/coverageBed -a - -b CS02_m2n240_4.cov0b | awk '(($4*100)/$5) < 11' | awk '{print $1":"$2"-"$3}'

##make list to remove
my $removeA = $tempdir . 'remove.tempA';
my $removeB = $tempdir . 'remove.tempB';
system("awk \'\$2 > 1.1\' $ambiguityOut >> $removeA");
system("awk \'\$2 > 0.3\' $phaseInconOut >> $removeA");
#system("cat $highlist $lowlist $removeA | awk \'\{print \$1\}\' | sort | uniq > $removeB"); 
system("cat $highlist $removeA | awk \'\{print \$1\}\' | sort | uniq > $removeB");

#remove them using external perl script
my $firstFilter = $home . $lib . '/map2/contigs2/' . $lib . '.4filts_good.fasta';
my $filtered    = $home . $lib . '/map2/contigs2/' . $lib . '.4filts_bad.fasta';
system("$selectseqs -v -p -f $removeB $maskedB > $firstFilter");
system("$selectseqs -p -f $removeB $maskedB > $filtered");

my $finalFasta = '/global/scratch/cdspecht/zingiberales/nucMapAFfastas/all_150123/' . $lib . '.AF.fa';
my $substitution = '>' . $lib . '_';
my $getRidofNs = '\'{gsub(/^[N]+|[N]+$/,"");print}\'';
my $getRidofEmptys = '\'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}\''; 
my $getnoemptylines = '\'^$\'';
#system("sed \'s\/\>\/$substitution\/\' $firstFilter \| awk $getRidofNs | grep -v $getnoemptylines | awk $getRidofEmptys > $finalFasta");
system("sed \'s\/\>\/$substitution\/\' $maskedB \| awk $getRidofNs | grep -v $getnoemptylines | awk $getRidofEmptys > $finalFasta");
}


