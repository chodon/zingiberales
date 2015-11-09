use warnings;
use strict;
use lib "/global/home/users/cdspecht/chodon/programs/local_perl/BioPerl-1.6.1";
use Bio::SeqIO;

my @lib = ('bdjq');

#my @lib = ('cost', 'curc', 'jnub', 'knkv', 'lskk', 'uoel');

#my @lib = ('tzns');
my $all_covered_exon_list = '/global/scratch/cdspecht/chodon/zingiberales/sep_exons_full_cds/all_covered_150bp_exons_list';

foreach my $lib (@lib){
    my $dir = '/global/scratch/cdspecht/chodon/zingiberales/';
# name an out file
my $hap0fil = $dir . 'sep_exons_full_cds/' . $lib . '.haplo0.fasta';
open(HAP0, ">$hap0fil");
my $hap1fil = $dir. 'sep_exons_full_cds/' . $lib . '.haplo1.fasta';
open(HAP1, ">$hap1fil");

# name vcf file
my $preraw = $dir . 'sep_exons_full_cds/' . $lib . '.varscan.test';
my $rawvcf = $dir . 'sep_exons_full_cds/' . $lib . '.varscan_150bp_only.test';
    system("grep -f $all_covered_exon_list $preraw \> $rawvcf");
my $processvcf = $dir . 'sep_exons_full_cds/' . $lib . '.varscan.process150.vcf';
system("grep -v \'Chrom\' $rawvcf \| sed \'s\/\%\/\/\' \> $processvcf");
my $vcffil = $dir . 'sep_exons_full_cds/' . $lib . '.varscan.process150.vcf';
open(VCF, "<$vcffil");
my @vcf = <VCF>;
    chomp;
# name reference sequence file
my $refs = $dir . 'sep_exons_full_cds/' . 'all_exons.fa';
    
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

    foreach my $vcline (@vcf){

          my @tmp = split(/\t/,$vcline);
          my $a0; my $a1;     
      
          if ($name eq $tmp[0]){
              
             # if (substr($tmp[9],0,1) == 0){
              #    $a0 = $tmp[4];
               #   $a1 = $tmp[3];
             # }

             # if (substr($tmp[9],0,1) == 1){
              #    $a0 = $tmp[3];
               #   $a1 = $tmp[4];
             # }
              
	    #  if ((substr($tmp[7],0,5) eq 'ABHet') && (substr($tmp[7],6,4) >= 0.8) && (substr($tmp[9],0,1) == 0)){
	#	      $a0 = $tmp[3];
	#	      $a1 = $tmp[3];
	 #     }

	  #    if ((substr($tmp[7],0,5) eq 'ABHet') && (substr($tmp[7],6,4) >= 0.8) && (substr($tmp[9],0,1) == 1)){
	#	      $a0 = $tmp[4];
	#	      $a1 = $tmp[4];
	 #     }

         #     if (substr($tmp[9],0,3) eq "1|1"){
	#	  $a0 = $tmp[4];
	#	  $a1 = $tmp[4];
	 #     }  
	      if (($tmp[3] eq $tmp[2]) && ($tmp[6] < 50)) {                                                                                                               $a0 = $tmp[3];                     
                  $a1 = $tmp[3];                                  
              }              

	      if (($tmp[3] ne $tmp[2]) && ($tmp[6] > 50)) {
                      $a0 = substr($tmp[18],0,1);
		      $a1 = $tmp[2];
              }
              
	      if (($tmp[3] ne $tmp[2]) && ($tmp[6] <= 50)) {    
                  $a0 = $tmp[2];
                  $a1 = substr($tmp[18],0,1);
      	      }

	      if (($tmp[3] eq $tmp[2]) && ($tmp[6] >= 50)) {

		  $a0 =$tmp[3];
                  $a1 = $tmp[3];
              }

              if (($tmp[3] eq "N") && ($tmp[4] == 0) && ($tmp[5] == 1)) {
                  
                  if (substr($tmp[18],0,1) eq "-") {
		      $a0 = $tmp[2];
		      my $delete_pre = $tmp[18];
		      my ($delete_post) = substr($delete_pre,1) =~ /([A|C|G|T]+)/;
		      $a1 = $delete_post;
		      my $a0length = length($a1);
		      my $a0toDelete = "Z" x $a0length;
		      $a0 = $tmp[2] . $a0toDelete;
		      print $lib . "\t" . $tmp[0] . "\t" . 'delete' . "\t" . $tmp[1] . "\t" . $a0length . "\n";
                  }

                  elsif (substr($tmp[18],0,1) eq "+") {
		      my $insert_pre = $tmp[18];
		      my ($insert_post) = substr($insert_pre,1) =~ /([A|C|G|T]+)/;
		      $a0 = $tmp[2] . $insert_post;
		      $a1 = $tmp[2];
		      print $lib . "\t" . $tmp[0] . "\t" . 'insert' . "\t" . $tmp[1] . "\t" . length($insert_post) . "\n";
                  }
                  
                  else {

                  $a0 = substr($tmp[18],0,1);
                  $a1 = $tmp[2];

	          }
              }

              if (($tmp[3] eq "N") && ($tmp[4] == 1) && ($tmp[5] == 0)) {
                  $a0 = $tmp[2];
                  $a1 = $tmp[2];
              }

	      if (($tmp[3] eq "N") && ($tmp[6] > 50) && (($tmp[4] + $tmp[5]) > 1 )) {

                  if (substr($tmp[18],0,1) eq "-") {
                      $a0 = $tmp[2];
                      my $delete_pre = $tmp[18];
                      my ($delete_post) = substr($delete_pre,1) =~ /([A|C|G|T]+)/;
                      $a1 = $delete_post;
                      my $a0length = length($a1);
                      my $a0toDelete = "Z" x $a0length;
                      $a0 = $tmp[2] . $a0toDelete;
                      print $lib . "\t" . $tmp[0] . "\t" . 'delete' . "\t" . $tmp[1] . "\t" . $a0length . "\n";
                  }

                  elsif (substr($tmp[18],0,1) eq "+") {
                      my $insert_pre = $tmp[18];
                      my ($insert_post) = substr($insert_pre,1) =~ /([A|C|G|T]+)/;
                      $a0 = $tmp[2] . $insert_post;
                      $a1 = $tmp[2];
                      print $lib . "\t" . $tmp[0] . "\t" . 'insert' . "\t" . $tmp[1] . "\t" . length($insert_post) . "\n";
                  }

                  else {

		      $a0 = substr($tmp[18],0,1);
		      $a1 = $tmp[2];

                  }
              }



              if (substr($tmp[3],0,1) eq "+") {
                  my $insert_pre = $tmp[3];
                  my ($insert_post) = substr($insert_pre,1) =~ /([A|C|G|T]+)/; 
                  $a0 = $tmp[2] . $insert_post;
		  $a1 = $tmp[2];
		  print $lib . "\t" . $tmp[0] . "\t" . 'insert' . "\t" . $tmp[1] . "\t" . length($insert_post) . "\n";
              }
              
              if ((substr($tmp[3],0,1) eq "*") && (substr($tmp[3],2,1) eq "-")) {
                  
		  $a0 = $tmp[2];
                  my $delete_pre = $tmp[3];
                  my ($delete_post) = substr($delete_pre,3) =~ /([A|C|G|T]+)/;
                  $a1 = $delete_post;
                  my $a0length = length($a1);
                  my $a0toDelete = "Z" x $a0length;
                  $a0 = $tmp[2] . $a0toDelete;
                  print $lib . "\t" . $tmp[0] . "\t" . 'delete' . "\t" . $tmp[1] . "\t" . $a0length . "\n";

	      }

	      if ((substr($tmp[3],0,1) eq "*") && (substr($tmp[3],2,1) eq "+") && ($tmp[6] >= 50)) {
		  my $insert_pre = $tmp[3];
                  my ($insert_post) = substr($insert_pre,3) =~ /([A|C|G|T]+)/;
                  $a0 = $tmp[2] . $insert_post;
                  $a1 = $tmp[2];
                  print $lib . "\t" . $tmp[0] . "\t" . 'insert' . "\t" . $tmp[1] . "\t" . length($insert_post) . "\n";
	      }

	      if ((substr($tmp[3],0,1) eq "*") && (substr($tmp[3],2,1) eq "+") && ($tmp[6] < 50)) {
                  $a0 = $tmp[2];
                  $a1 = $tmp[2];
              }

	      if (substr($tmp[3],0,1) eq "-") {

	          $a0 = $tmp[2];
                  my $delete_pre = $tmp[3];
                  my ($delete_post) = substr($delete_pre,1) =~ /([A|C|G|T]+)/; 
                  $a1 = $delete_post;
                  my $a0length = length($a1);
                  my $a0toDelete = "Z" x $a0length;
                  $a0 = $tmp[2] . $a0toDelete;
                  print $lib . "\t" . $tmp[0] . "\t" . 'delete' . "\t" . $tmp[1] . "\t" . $a0length . "\n";
              }

              if (($tmp[4] < 5) && ($tmp[5] < 5 )){
                  $a0 = lc($a0);
                  $a1 = lc($a1);
                  
              }

              substr($h0seq, $tmp[1]+$h0indlen-1, 1) = $a0;
              substr($h1seq, $tmp[1]+$h1indlen-1, 1) = $a1;
	      # substr($h0seq, $tmp[1]-1, 1) = $a0;
	      # substr($h1seq, $tmp[1]-1, 1) = $a1;

              $h0indlen = $h0indlen + length($a0)-1;
              $h1indlen = $h1indlen + length($a1)-1;
 
#print $tmp[0] ."\t" . $tmp[1] . "\t" . $h0indlen ."\t" . $h1indlen ."\n";
          } 
      }



    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
    # process each seq
    # print $seq->id . ' = '.$seq->seq()."\n";
}
close HAP0;
close HAP1;
close VCF;
}
