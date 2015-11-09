########################################################################################################
# a script to use make pileups and run varscan to correct reference sequence          #
# written by Chodon Sass chodon at gmail.com May 2014   

########################################################################################################


use warnings;
use strict;
use lib "/global/home/users/cdspecht/chodon/programs/local_perl/BioPerl-1.6.1";
use Bio::SeqIO;
use List::Util qw(min max);


#this defines the library (L) and family specific reference fasta (Z) so it can be specified at the qsub command

#home directory with all my files; make sure to end with a backslash

my $fastqdir = '/global/scratch/cdspecht/zingiberales/Project_Specht/precleaned/clean/';
my $AFoutDir = '/global/scratch/cdspecht/zingiberales/Project_Specht/precleaned/nucMapAF/';
my $home = '/global/scratch/cdspecht/';


#program paths
my $novoalign = '/global/home/users/cdspecht/chodon/programs/novocraft/novoalign';
my $novoindex = '/global/home/users/cdspecht/chodon/programs/novocraft/novoindex';
my $gatk = '/global/home/users/cdspecht/bin/GenomeAnalysisTK.jar';


#library name
#already done = CS72 CS30

my @lib = qw(CS02 CS06 CS07 CS09 CS12 CS13 CS14 CS15 CS16 CS17 CS18 CS19 CS20 CS21 CS22 CS23 CS24 CS25 CS27 CS28 CS29 CS30 CS31 CS32 CS33 CS34 CS35 CS36 CS37 CS38 CS39 CS40 CS41 CS42 CS43 CS44 CS45 CS46 CS47 CS48 CS49 CS50 CS51 CS52 CS53 CS54 CS55 CS57 CS58 CS69 CS70 CS71 CS72 COST TZNS UOEL OYLU LSKK KNKV JQCX JNUB BDJQ BRUD CAGR KIAU PPQR QGLJ TRPA BAAU BYPY CHSE DAS2 DASP DITH HAMA HWUP NSPR);

###########################
# behold the subroutines! #
# modify at your own risk #
###########################

system("date");
	
foreach my $lib (@lib) {
    makegatk($lib);


sub makegatk {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $AFoutDir . $lib . '/' . 'map2' . '/';
	my $contigDir = $resultsDir . 'contigs2' . '/';
	my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean4.fa';       
        my $bam = $contigDir . $lib . '_m2n240_4r.bam';
        my $vcf = $contigDir . $lib . '.filter.vcf';
        
        my $filterQual = q("QD < 2.0 || HaplotypeScore > 13.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0");
        my $filterDp20 = '"DP < 20 && DP > 10"';
        my $filterDp10 = '"DP < 11"';
        my $filterHethigh = '"ABHet > 0.799"';
        my $filterHetlow = '"ABHet < 0.200"';
        my $filterAltpos = '"ABHet < 0.500 && ABHet > 0.199"';
        my $filterAltdef = '"ABHom == 1.00"';
        my $filterHomhigh = '"ABHom > 0.799 && ABHom < 1.00"';
        my $filterHomlow = '"ABHom > 0.499 && ABHom < 0.8"';
        my $exp = '--filterExpression';
        my $name = '--filterName';
       
###### run gatk
         
        system("java -Xmx8g -jar /global/home/users/cdspecht/bin/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $bam -o $contigDir$lib.hc.vcf");
        system("java -Xmx8g -jar /global/home/users/cdspecht/bin/GenomeAnalysisTK.jar -T ReadBackedPhasing -R $ref -I $bam --variant $contigDir$lib.hc.vcf --min_base_quality_score 10 -o $contigDir$lib.raw.vcf");
        system("java -Xmx8g -jar /global/home/users/cdspecht/bin/GenomeAnalysisTK.jar -T VariantAnnotator -A DepthPerAlleleBySample -A AlleleBalance -A FisherStrand -A HaplotypeScore -A HardyWeinberg -R $ref -I $bam --variant $contigDir$lib.raw.vcf -o $contigDir$lib.annotated.vcf");
        system("java -Xmx8g -jar /global/home/users/cdspecht/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -I $bam -R $ref -o $contigDir$lib.ug.vcf -dcov 250 --genotype_likelihoods_model BOTH");
        system("java -jar /global/home/users/cdspecht/bin/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -o $vcf --variant $contigDir$lib.annotated.vcf $exp $filterQual $name \"qual\" $exp $filterDp20 $name \"dp20\" $exp $filterDp10 $name \"dp10\" $exp $filterHethigh $name \"hethigh\" $exp $filterHetlow $name \"hetlow\" $exp $filterAltpos $name \"altpos\" $exp $filterAltdef $name \"altdef\" $exp $filterHomhigh $name \"homhigh\" $exp $filterHomlow $name \"homlow\"");  
}



}

	
