########################################################################################################
# a script to use run 8_2 aligner with novoalign           #
# written by Chodon Sass chodon at gmail.com August 2012   
## adapted from:                              #
########################################################################################################
########################################################################################################
# a script to prepare shell scripts to use for assembly with multiple de novo assemblers on TACC       #
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 13 Dec 2011                                 #
########################################################################################################

use warnings;
use strict;

#if you want to use these script, your library should be in the first column of your library file

#home directory with all my files; make sure to end with a backslash
#my $fastqdir = '/media/Data/chodon/zingiberales/';
my $home = '/media/Data/chodon/zingiberales/';
#library name
#my @lib = ('COST');
my @lib = ('TZNS', 'KNKV', 'BDJQ', 'JNUB', 'OYLU', 'JQCX', 'LSKK', 'UOEL');
#already run: 'COST'
#, TZNS, KNKV, BDJQ, JNUB, OYLU, JQCX, LSKK, UOEL');
#'Canna_sp-TZNS', 'Curcuma_olena-JQCX-A', 'Curcuma_olena-OYLU-B', 'Heliconia_sp-KNKV', 'Maranta_leuconeura-JNUB', 'Orchidantha_maxillaroides-LSKK', 'Strelitzia_reginae-UOEL', 'Zingiber_officinale-BDJQ-young_leaves', 'Costus_pulverulentus');

my $insertSize = 200;
my $insertStd = 70;
my $maxScorePE = 502;
my $maxScoreU = 248;



###########################
# behold the subroutines! #
# modify at your own risk #
###########################
	
foreach my $lib (@lib) {
	makeNovo($lib);

sub makeNovo {
	my ($lib) = @_;
    my $fastqdir = $home . $lib . '/' . $lib . '/';
    my $resultsDir =  $home . $lib . '/' . $lib . '_full_cds_novo540';
    my $ref = $home  . 'full_cds_sep_exons.fa';
    my $out = $resultsDir . '/' . $lib . '.sorted';
    my $cov = $resultsDir . '/' . $lib . '.cov'; 
    my $cleanzip1 = $fastqdir . $lib . '_1_final.fastq.gz';
    my $cleanzip2 = $fastqdir . $lib . '_2_final.fastq.gz';
	my $reads1 = $fastqdir . $lib . '_1_final.fastq';
	my $reads2 = $fastqdir . $lib . '_2_final.fastq';
    my $cleanzipu = $fastqdir . $lib . '_u_final.fastq.gz';
        my $readsu = $fastqdir . $lib . '_u_final.fastq';
#    my $cleanReads = $fastqdir . $lib;
    my $contigs = $home . 'musa_cds.fna';
   # my $contigs2 = $home . $lib . '/' . $lib . '_assembly.final.fa';
   # my $bed = $home . $lib . '/' . 'Final_BED.txt';
    
    system("mkdir $resultsDir");
    ##make contigs from abyss end in .fa
 #   system("cp $contigs $contigs2");
    ##unzip clean fastqs for novoalign
    system("gunzip $cleanzip1");
    system("gunzip $cleanzip2");
    system("gunzip $cleanzipu");
    
    #run truncated novoalign perl script
my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";	
system ("novoindex $indexed_assemblies_in_target $ref");

system("novoalign -d $indexed_assemblies_in_target -f $reads1  $reads2 -i PE $insertSize\, $insertStd -t $maxScorePE -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outPairedSam1");

system("novoalign -d $indexed_assemblies_in_target -f $readsu -t $maxScoreU -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outSoloSam1");

#    system("perl /media/Data/passerine/chodon/drakolibs/draco_scrubbed/8_2.Alignment_by_novoalignJP.pl -f $ref -r $cleanReads -a $contigs -o $resultsDir -m $insertSize -v $insertStd -t $maxScore -F STDFQ");
    
    #only printout aligned reads
    system("grep -v \'ZS:Z:NM\' outPairedSam1 > target_pair.sam");
    system("grep -v \'ZS:Z:NM\' outSoloSam1 > target_solo.sam");
    
    #run samtools to bam, merge sort and index
    system("samtools view -bS target_pair.sam > target_pair.bam");
    system("samtools view -bS target_solo.sam > target_solo.bam");
    system("samtools merge -f target.bam target_solo.bam target_pair.bam");
    system("samtools sort target.bam $out");
    system("samtools index $out.bam"); 
    
    #generate a genome file, generate coverage with bedtools
    system("samtools idxstats $out.bam > $out.stats");
    system("awk \< $out.stats \'\{print \$1\"\\t\"\$2 \}\' \> $out.genome");
    system("sed \'\$d\' $out.genome \> $out.genome2");
    system("mv $out.genome2 $out.genome");
    system("genomeCoverageBed -ibam $out.bam -d -g $out.genome \> $cov");
    
    #remove temp files to get ready for next library
    system("rm target_pair.sam target_solo.sam target_pair.bam target_solo.bam outSoloSam1 outPairedSam1 target.bam");
    
    #run samtools to generate fasta
	system("samtools mpileup -uf $ref $out.bam | bcftools view -cg - | vcfutils_CS.pl vcf2fq -d 1 > $lib.fullcds.fasta");

    #zip up
	system("gzip $reads1");
	system("gzip $reads2");
	system("gzip $readsu");

    print "\n\nfinished with $lib \n\n\n";
    }
  }

	
