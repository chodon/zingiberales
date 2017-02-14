########################################################################################################
# a script to use make pileups and run varscan to correct reference sequence          #
# written by Chodon Sass chodon at gmail.com May 2014   

########################################################################################################


use warnings;
use strict;
#?# use lib "/global/home/users/cdspecht/chodon/programs/local_perl/BioPerl-1.6.1";
use Bio::SeqIO;
use List::Util qw(min max);
use Getopt::Std;

#this defines the library (L) and family specific reference fasta (Z) so it can be specified at the qsub command

my %opts = (L=>undef, Z=>undef);
getopts('L:Z:', \%opts);

#home directory with all my files; make sure to end with a backslash

#?# my $fastqdir = '/global/scratch/cdspecht/zingiberales/Project_Specht/precleaned/clean/';
#?# my $mapsembleOutDir = '/global/scratch/cdspecht/zingiberales/Project_Specht/precleaned/nucMapAF/';
#?# my $home = '/global/scratch/cdspecht/';


#program paths
#?# my $VarScan = '/global/scratch/cdspecht/chodon/zingiberales/sep_exons_full_cds/VarScan.v2.3.6.jar';
#?# my $novoalign = '/global/home/users/cdspecht/chodon/programs/novocraft/novoalign';
#?# my $novoindex = '/global/home/users/cdspecht/chodon/programs/novocraft/novoindex';
#?# my $gatk = '/global/home/users/cdspecht/bin/GenomeAnalysisTK.jar';
#?# my $mapsembler = '/global/home/users/cdspecht/bin/mapsembler2/tools/mapsembler';
#?# my $trash1 = '/global/scratch/cdspecht/chodon/programs/zing_nuc_map_align/' . $opts{L} . 'temp/trashmeplease*';
#?# my $trash2 = '/global/scratch/cdspecht/chodon/programs/zing_nuc_map_align/' . $opts{L} . 'temp/index_k*';
#?# my $recipblast = '/global/home/users/cdspecht/bin/7recipBlasting.pl';
## note that -b, -f and -r are set per library in this pipeline
my $picard = '/clusterfs/vector/home/groups/software/sl-6.x86_64/modules/picard/2.4.1/picard.jar';
    

#?# my $originalRef = $home . 'zingiberales/Project_Specht/precleaned/nucAlignFix/' . $opts{L} .  '/' . 'map2/' . $opts{L} . '.clean4.fa';  

#arguments for novoalign
my $insertSize = 230;
my $insertStd = 70;
my $maxScore = 240;

#library name

my @lib = ($opts{L});      #cost
#my @lib = ('','');      #musa
#my @lib = ('','');      #heli
#my @lib = ('','');      #mara
#my @lib = ('','');      #orch
#my @lib = ('','');      #zing
#my @lib = ('','');      #cann
#my @lib = ('','');      #stre

###########################
# behold the subroutines! #
# modify at your own risk #
###########################

system("date");
	
foreach my $lib (@lib) {
    mapsembler($lib);
    makeInitialAlignment($lib);
    makepileups($lib);
    mapsembler2($lib);
    makeSecondaryAlignment($lib);
    makepileups2($lib);
    makeThirdAlignment($lib);
    makepileups3($lib);
    makeFourthAlignment($lib);
    makepileups4($lib);    
    makeFinalAlignment($lib);


sub mapsembler {
    my ($lib) = @_;
    
    my $clean1 = $fastqdir . $lib . '_1_final.txt';
    my $clean2 = $fastqdir . $lib . '_2_final.txt';
    my $cleanu = $fastqdir . $lib . '_u_final.txt';
    
    #define directories
    my $resultsDir = $mapsembleOutDir . $lib . '/';
    system("mkdir $resultsDir");
    my $contigDir = $resultsDir . 'contigs' . '/';
    system("mkdir $contigDir");
    my $mapout = $resultsDir . $lib . '_mapsembler';
    my $mergeout = $contigDir . $lib . '_merge.fasta';
    my $k21 = $resultsDir . 'k21.fasta';
    my $k31 = $resultsDir . 'k31.fasta';
    my $k41 = $resultsDir . 'k41.fasta';
    my $k51 = $resultsDir . 'k51.fasta';
    my $k61 = $resultsDir . 'k61.fasta';
    my $k71 = $resultsDir . 'k71.fasta';
    my $k21map = $mapout . '_k_21_q_25_c_2_t_2.fasta';
    my $k31map = $mapout . '_k_31_q_25_c_2_t_2.fasta';
    my $k41map = $mapout . '_k_41_q_25_c_2_t_2.fasta';
    my $k51map = $mapout . '_k_51_q_25_c_2_t_2.fasta';
    my $k61map = $mapout . '_k_61_q_25_c_2_t_2.fasta';
    my $k71map = $mapout . '_k_71_q_25_c_2_t_2.fasta';
    
    #call mapsembler
    
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 21");
    system("rm $trash1 $trash2");
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 31");
    system("rm $trash1 $trash2");
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 41");
    system("rm $trash1 $trash2");
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 51");
    system("rm $trash1 $trash2");
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 61");
    system("rm $trash1 $trash2");
    system("$mapsembler $originalRef $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 71");
    system("rm $trash1 $trash2");
    
    #get ready for and call 7recipBlasting.pl
    
    #renames sequences so can see kmer in name
    my $mapHeader = 'TODO: read and store the fragment comment';
    system("sed \'s\/$mapHeader\/k21\/\' $k21map > $k21");
    system("sed \'s\/$mapHeader\/k31\/\' $k31map > $k31");
    system("sed \'s\/$mapHeader\/k41\/\' $k41map > $k41");
    system("sed \'s\/$mapHeader\/k51\/\' $k51map > $k51");
    system("sed \'s\/$mapHeader\/k61\/\' $k61map > $k61");
    system("sed \'s\/$mapHeader\/k71\/\' $k71map > $k71");
    
    system("cat $k21 $k31 $k41 $k51 $k61 $k71 > $mergeout");
    system("rm $mapout*");
    #calls sonal's script
    system("perl $recipblast -o $contigDir -a $mergeout -b $originalRef");
    
    print "\n\nfinished with first assembly round for $lib \n\n\n";
    
    }    

sub makeInitialAlignment {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $mapsembleOutDir . $lib . '/';
        my $contigDir = $resultsDir . 'contigs' . '/';
        my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $preref = $contigDir . 'finalBaitAssembly.fa';
        my $ref = $resultsDir . $lib . '_m1.fa';
        system("cp $preref $ref");
        
        my $bam = $contigDir . $lib . '_m1n240';
        
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
	    
        #define ref and bam for next step    
        my $tempBam = $contigDir . $lib . '_m1n240r.bam';
    

###### run novoalign on mapsembled reference
        #make out directories and files

        #generate index
		my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";
		system("$novoindex $indexed_assemblies_in_target $ref");

		#run novoalign
		#define temp files
		
		my $outPairedSam1 = $tempDir . 'ops1';
		my $outSoloSam1 = $tempDir . 'oss1';
		my $target_pair_sam = $tempDir . 'tp.sam';
		my $target_solo_sam = $tempDir . 'ts.sam';
		my $target_pair_bam = $tempDir . 'tp.bam';
		my $target_solo_bam = $tempDir . 'ts.bam';
		my $target_bam = $tempDir . 't.bam';
		my $target_sorted = $tempDir . 't.sorted';
		my $target_sorted_bam = $tempDir . 't.sorted.bam';
		my $target_duped_bam = $tempDir . 't.duped.bam';
		my $target_rg_bam = $tempDir . 't.rg.bam';
		my $target_genome = $tempDir . 't.genome';
		my $target_genome2 = $tempDir . 't2.genome';
		
        system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE $insertSize, $insertStd -t 90 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");

        system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 90 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outSoloSam1");

        #only printout aligned reads
		system("grep -v \'ZS:Z:NM\' $outPairedSam1 > $target_pair_sam");
		system("grep -v \'ZS:Z:NM\' $outSoloSam1 > $target_solo_sam");
    
        #run samtools to bam, merge sort and index
		system("samtools view -bS $target_pair_sam > $target_pair_bam");
		system("samtools view -bS $target_solo_sam > $target_solo_bam");
		system("samtools merge -f $target_bam $target_solo_bam $target_pair_bam");
		system("samtools sort $target_bam $target_sorted");
		system("samtools index $target_sorted_bam"); 
  
		#make readgroups and remove duplicates
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard MarkDuplicates INPUT=$target_sorted_bam OUTPUT=$target_duped_bam METRICS_FILE=$bam.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard AddOrReplaceReadGroups INPUT=$target_duped_bam OUTPUT=$target_rg_bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
        system("samtools sort $target_rg_bam $bam");
	system("samtools index $bam.bam");
        system("samtools idxstats $bam.bam > $bam.stats");

        #generate a genome file, generate coverage with bedtools
       	system("awk \< $bam.stats \'\{print \$1\"\\t\"\$2 \}\' \> $target_genome");
	system("sed \'\$d\' $target_genome \> $target_genome2");
        system("mv $target_genome2 $bam.genome");
        system("genomeCoverageBed -ibam $bam.bam -d -g $bam.genome \> $bam.cov");

        ##### run gatk target realigner
        system("samtools faidx $ref");
        my $dict = substr ($ref, 0, -2) . "dict";
        my $intervals = substr ($ref, 0, -2) . "intervals";
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $ref -I $bam.bam -o $intervals");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $ref -I $bam.bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");
       
       #remove temp files, get ready for phasing 
       system("mv $bam.realigned.bam $tempBam");
       system("samtools index $tempBam");
       system("rm -r $tempDir");#target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.bam target.sorted.bam target.sorted.bam.bai target.duped.bam target.rg.bam target.genome");
         
    print "\n\nfinished with initial alignment for $lib \n\n\n";
    
}
     
    
sub makepileups {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir = $mapsembleOutDir . $lib . '/';
        my $contigDir = $resultsDir . 'contigs' . '/';
        
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '_m1.fa';
        my $tempRef = $resultsDir . $lib . '.clean1.fa';
        my $bam = $contigDir . $lib . '_m1n240r.bam';
        my $tempDir = $contigDir . 'temporary/';
        system("mkdir $tempDir");
        
        #define outfiles
        my $pileup = $tempDir  . $lib . '.pileup';
        my $tempPileup = $tempDir  . $lib . '.temp.pileup';
        my $varfile = $tempDir  . $lib . '.varscan'; 
        my $tempVarfile = $tempDir  . $lib . '.temp.varscan'; 
 	
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
    
    #run samtools to generate pileup
    system("samtools mpileup -AB -f $ref $bam > $pileup");

    #run varscan 
    system("java -Djava.io.tmpdir=$tempDir -jar $VarScan pileup2cns $pileup --min-coverage 1 --min-var-freq 0.4 --min-reads2 1 --p-value 0.17 > $varfile");
     
    ####### parse VarScan to new fasta (note - this was adapted from another file so rather than renaming variables, i just 
            # used the definitions already in the script. therefore, some files are referenced by 2 different handles
    
    # name an out file                                                                                                         
	my $hap0fil = $tempDir . $lib . '.haplo0.fasta';
	open(HAP0, ">$hap0fil");
	my $hap1fil = $tempDir . $lib . '.haplo1.fasta';
	open(HAP1, ">$hap1fil");

    # name vcf file                                                                                                            
	
	my $rawvcf = $varfile;
	my $processvcf = $tempDir . $lib . '.varscan.process.vcf';
	system("grep -v \'Chrom\' $rawvcf \| sed \'s\/\%\/\/\' \> $processvcf");
	my $vcffil = $tempDir . $lib . '.varscan.process.vcf';
	open(VCF, "<$vcffil");
	my @vcf = <VCF>;
	chomp;
    
    # handle IO with bioperl                                                                                                   
        my $seqio = Bio::SeqIO-> new(
                             -file     => $ref,
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

      		    if (($tmp[3] eq $tmp[2]) && ($tmp[6] < 50)) {                                                               
	       			    $a0 = $tmp[3];
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

		    $h0indlen = $h0indlen + length($a0)-1;
		    $h1indlen = $h1indlen + length($a1)-1;
		    
		}
	    }

	    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
	    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
	}
	close HAP0;
	close HAP1;
	close VCF;

##### remove Zs to resolve indels, number of Z's plus equal number of following base pairs need to be removed

	my $rawFasta = $tempDir . $lib . '.haplo0.fasta';
        my $rawFasta2 = $tempDir . $lib . '.haploZ.fasta';
        system("sed \'s\/z\/Z\/g\' $rawFasta > $rawFasta2");
        system("mv $rawFasta2 $rawFasta");

        my $Zs = $tempDir . $lib . '.Zs.txt';
        open(FAS, "<$rawFasta");
        open(OUT, ">$Zs");  
        my @zcount;
        while(<FAS>){
	    chomp(my $line = $_);
            while ( $line  =~ /(Z+)/g) {
        push @zcount, length($1);
        print OUT length($1), "\n";
        }
	}
        my $maxZcount = max @zcount;
        print "$maxZcount \t $lib \n";
        close FAS;
        close OUT; 
        
        my $cleanFasta = $tempDir . $lib . '.clean.fa';
        my $precleanFasta = $tempDir . $lib . '.preclean.fa';
        system("cp $rawFasta $precleanFasta");

        for (my $i = $maxZcount; $i > 0; $i--) {
	    system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #      sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta
            system("mv $cleanFasta $precleanFasta");
	}

         ### do it again, if there were Z's in the way that prevented Zs from getting removed

	for (my $i = $maxZcount; $i > 0; $i--) {
            system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #          sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta                                                                                          
            system("mv $cleanFasta $precleanFasta");
        }
         
        system("mv $precleanFasta $cleanFasta");
        system("cp $cleanFasta $tempRef");

    ## cleanup
    system("rm -r $tempDir");

    print "\n\nfinished with first cleanup for $lib \n\n\n";

}

sub mapsembler2 {
    my ($lib) = @_;
    
    my $clean1 = $fastqdir . $lib . '_1_final.txt';
    my $clean2 = $fastqdir . $lib . '_2_final.txt';
    my $cleanu = $fastqdir . $lib . '_u_final.txt';
    
    my $ref = $mapsembleOutDir . $lib . '/' . $lib . '.clean1.fa';
    
    #define directories
    my $resultsDir = $mapsembleOutDir . $lib . '/' . 'map2' . '/';
    system("mkdir $resultsDir");
    my $contigDir = $resultsDir . 'contigs2' . '/';
    system("mkdir $contigDir");
    my $mapout = $resultsDir . $lib . '_mapsembler';
    my $mergeout = $contigDir . $lib . '_merge.fasta';
    my $k21 = $resultsDir . 'k21.fasta';
    my $k31 = $resultsDir . 'k31.fasta';
    my $k41 = $resultsDir . 'k41.fasta';
    my $k51 = $resultsDir . 'k51.fasta';
    my $k61 = $resultsDir . 'k61.fasta';
    my $k71 = $resultsDir . 'k71.fasta';
    my $k21map = $mapout . '_k_21_q_25_c_2_t_2.fasta';
    my $k31map = $mapout . '_k_31_q_25_c_2_t_2.fasta';
    my $k41map = $mapout . '_k_41_q_25_c_2_t_2.fasta';
    my $k51map = $mapout . '_k_51_q_25_c_2_t_2.fasta';
    my $k61map = $mapout . '_k_61_q_25_c_2_t_2.fasta';
    my $k71map = $mapout . '_k_71_q_25_c_2_t_2.fasta';

    #call mapsembler
    
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 21");
    system("rm $trash1 $trash2");
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 31");
    system("rm $trash1 $trash2");
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 41");
    system("rm $trash1 $trash2");
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 51");
    system("rm $trash1 $trash2");
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 61");
    system("rm $trash1 $trash2");
    system("$mapsembler $ref $clean1 $clean2 $cleanu -E -t 2 -o $mapout -k 71");
    system("rm $trash1 $trash2");
    
    #get ready for and call 7recipBlasting.pl                                                                                                                                                                     
    #renames sequences so can see kmer in name                                                                                                    
    my $mapHeader = 'TODO: read and store the fragment comment';
    system("sed \'s\/$mapHeader\/k21\/\' $k21map > $k21");
    system("sed \'s\/$mapHeader\/k31\/\' $k31map > $k31");
    system("sed \'s\/$mapHeader\/k41\/\' $k41map > $k41");
    system("sed \'s\/$mapHeader\/k51\/\' $k51map > $k51");
    system("sed \'s\/$mapHeader\/k61\/\' $k61map > $k61");
    system("sed \'s\/$mapHeader\/k71\/\' $k71map > $k71");

    system("cat $k21 $k31 $k41 $k51 $k61 $k71 > $mergeout");
    system("rm $mapout*");
    #calls sonal's script
    system("perl $recipblast -o $contigDir -a $mergeout -b $ref");
    
    print "\n\nfinished with second assembly round for $lib \n\n\n";
    
    }    

sub makeSecondaryAlignment {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $mapsembleOutDir . $lib . '/' . 'map2' . '/';
        my $contigDir = $resultsDir . 'contigs2' . '/';
	my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $preref = $contigDir . 'finalBaitAssembly.fa';
        my $ref = $resultsDir . $lib . '_m2.fa';
        system("cp $preref $ref");
        
        my $bam = $contigDir . $lib . '_m2n240';
        
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
	    
        #define ref and bam for next step
	    
        my $tempBam = $contigDir . $lib . '_m2n240r.bam';
    

###### run novoalign on mapsembled reference
        #make out directories and files

        #generate index
		my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";
		system("$novoindex $indexed_assemblies_in_target $ref");

		#run novoalign
		#define temp files
		
		my $outPairedSam1 = $tempDir . 'ops1';
		my $outSoloSam1 = $tempDir . 'oss1';
		my $target_pair_sam = $tempDir . 'tp.sam';
		my $target_solo_sam = $tempDir . 'ts.sam';
		my $target_pair_bam = $tempDir . 'tp.bam';
		my $target_solo_bam = $tempDir . 'ts.bam';
		my $target_bam = $tempDir . 't.bam';
		my $target_sorted = $tempDir . 't.sorted';
		my $target_sorted_bam = $tempDir . 't.sorted.bam';
		my $target_duped_bam = $tempDir . 't.duped.bam';
		my $target_rg_bam = $tempDir . 't.rg.bam';
		my $target_genome = $tempDir . 't.genome';
		my $target_genome2 = $tempDir . 't2.genome';
		
        system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE $insertSize, $insertStd -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");

        system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outSoloSam1");

        #only printout aligned reads
		system("grep -v \'ZS:Z:NM\' $outPairedSam1 > $target_pair_sam");
		system("grep -v \'ZS:Z:NM\' $outSoloSam1 > $target_solo_sam");
    
        #run samtools to bam, merge sort and index
		system("samtools view -bS $target_pair_sam > $target_pair_bam");
		system("samtools view -bS $target_solo_sam > $target_solo_bam");
		system("samtools merge -f $target_bam $target_solo_bam $target_pair_bam");
		system("samtools sort $target_bam $target_sorted");
		system("samtools index $target_sorted_bam"); 
  
	#make readgroups and remove duplicates
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard MarkDuplicates INPUT=$target_sorted_bam OUTPUT=$target_duped_bam METRICS_FILE=$bam.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard AddOrReplaceReadGroups INPUT=$target_duped_bam OUTPUT=$target_rg_bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
	system("samtools sort $target_rg_bam $bam");
	system("samtools index $bam.bam");
        system("samtools idxstats $bam.bam > $bam.stats");

        #generate a genome file, generate coverage with bedtools
	system("awk \< $bam.stats \'\{print \$1\"\\t\"\$2 \}\' \> $target_genome");
	system("sed \'\$d\' $target_genome \> $target_genome2");
        system("mv $target_genome2 $bam.genome");
        system("genomeCoverageBed -ibam $bam.bam -d -g $bam.genome \> $bam.cov");

        
        ##### run gatk target realigner
        system("samtools faidx $ref");
        my $dict = substr ($ref, 0, -2) . "dict";
        my $intervals = substr ($ref, 0, -2) . "intervals";
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $ref -I $bam.bam -o $intervals");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $ref -I $bam.bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");


#remove temp files, get ready for phasing
        
       system("mv $bam.realigned.bam $tempBam");
       system("samtools index $tempBam");
       system("rm -r $tempDir");#target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.bam target.sorted.bam target.sorted.bam.bai target.duped.bam target.rg.bam target.genome");
         

    print "\n\nfinished with secondary alignment for $lib \n\n\n";
    
    

}


sub makepileups2 {
	my ($lib) = @_;
        
        #define directories
	my $resultsDir = $mapsembleOutDir . $lib . '/' . 'map2' . '/';
	my $contigDir = $resultsDir . 'contigs2' . '/';
        
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '_m2.fa';
        my $tempRef = $resultsDir . $lib . '.clean2.fa';
        my $bam = $contigDir . $lib . '_m2n240r.bam';
        my $tempDir = $contigDir . 'temporary/';
        system("mkdir $tempDir");
        
        #define outfiles
        my $pileup = $tempDir  . $lib . '.pileup';
        my $tempPileup = $tempDir  . $lib . '.temp.pileup';
        my $varfile = $tempDir  . $lib . '.varscan'; 
        my $tempVarfile = $tempDir  . $lib . '.temp.varscan'; 
 	
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
	my $clean2 = $fastqdir . $lib . '_2_final.txt';
	my $cleanu = $fastqdir . $lib . '_u_final.txt';
    #run samtools to generate pileup

    system("samtools mpileup -AB -f $ref $bam > $pileup");

    #run varscan 
   
    system("java -Djava.io.tmpdir=$tempDir -jar $VarScan pileup2cns $pileup --min-coverage 1 --min-var-freq 0.4 --min-reads2 1 --p-value 0.17 > $varfile");
   
    ####### parse VarScan to new fasta (note - this was adapted from another file so rather than renaming variables, i just 
            # used the definitions already in the script. therefore, some files are referenced by 2 different handles
    
    # name an out file                                                                                                         
	my $hap0fil = $tempDir . $lib . '.haplo0.fasta';
	open(HAP0, ">$hap0fil");
	my $hap1fil = $tempDir . $lib . '.haplo1.fasta';
	open(HAP1, ">$hap1fil");

    # name vcf file                                                                                                            
	
	my $rawvcf = $varfile;
	my $processvcf = $tempDir . $lib . '.varscan.process.vcf';
	system("grep -v \'Chrom\' $rawvcf \| sed \'s\/\%\/\/\' \> $processvcf");
	my $vcffil = $tempDir . $lib . '.varscan.process.vcf';
	open(VCF, "<$vcffil");
	my @vcf = <VCF>;
	chomp;
    
    # handle IO with bioperl                                                                                                   
        my $seqio = Bio::SeqIO-> new(
                             -file     => $ref,
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

      		    if (($tmp[3] eq $tmp[2]) && ($tmp[6] < 50)) {                                                               
	       			    $a0 = $tmp[3];
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

		    $h0indlen = $h0indlen + length($a0)-1;
		    $h1indlen = $h1indlen + length($a1)-1;
		    
		}
	    }

	    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
	    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
	}
	close HAP0;
	close HAP1;
	close VCF;

##### remove Zs to resolve indels, number of Z's plus equal number of following base pairs need to be removed

	my $rawFasta = $tempDir . $lib . '.haplo0.fasta';
        my $rawFasta2 = $tempDir . $lib . '.haploZ.fasta';
        system("sed \'s\/z\/Z\/g\' $rawFasta > $rawFasta2");
        system("mv $rawFasta2 $rawFasta");

        my $Zs = $tempDir . $lib . '.Zs.txt';
        open(FAS, "<$rawFasta");
        open(OUT, ">$Zs");  
        my @zcount;
        while(<FAS>){
	    chomp(my $line = $_);
            while ( $line  =~ /(Z+)/g) {
        push @zcount, length($1);
        print OUT length($1), "\n";
        }
	}
        my $maxZcount = max @zcount;
        print "$maxZcount \t $lib \n";
        close FAS;
        close OUT; 
        
        my $cleanFasta = $tempDir . $lib . '.clean.fa';
        my $precleanFasta = $tempDir . $lib . '.preclean.fa';
        system("cp $rawFasta $precleanFasta");

        for (my $i = $maxZcount; $i > 0; $i--) {
	    system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #      sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta
            system("mv $cleanFasta $precleanFasta");
	}

         ### do it again, if there were Z's in the way that prevented Zs from getting removed

	for (my $i = $maxZcount; $i > 0; $i--) {
            system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #          sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta                                                                                          
            system("mv $cleanFasta $precleanFasta");
        }
         
        system("mv $precleanFasta $cleanFasta");
        system("cp $cleanFasta $tempRef");

    ## cleanup
    system("rm -r $tempDir");

    print "\n\nfinished with second cleanup for $lib \n\n\n";

}

sub makeThirdAlignment {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $mapsembleOutDir . $lib . '/' . 'map2' . '/';
        my $contigDir = $resultsDir . 'contigs2' . '/';
        my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean2.fa';       
        my $bam = $contigDir . $lib . '_m2n240_2';
        
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
	    
	#define ref and bam for next step
	    
        my $tempBam = $contigDir . $lib . '_m2n240_2r.bam';
    

###### run novoalign on mapsembled reference
        #make out directories and files

        #generate index
		my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";
		system("$novoindex $indexed_assemblies_in_target $ref");

		#run novoalign
		#define temp files
		
		my $outPairedSam1 = $tempDir . 'ops1';
		my $outSoloSam1 = $tempDir . 'oss1';
		my $target_pair_sam = $tempDir . 'tp.sam';
		my $target_solo_sam = $tempDir . 'ts.sam';
		my $target_pair_bam = $tempDir . 'tp.bam';
		my $target_solo_bam = $tempDir . 'ts.bam';
		my $target_bam = $tempDir . 't.bam';
		my $target_sorted = $tempDir . 't.sorted';
		my $target_sorted_bam = $tempDir . 't.sorted.bam';
		my $target_duped_bam = $tempDir . 't.duped.bam';
		my $target_rg_bam = $tempDir . 't.rg.bam';
		my $target_genome = $tempDir . 't.genome';
		my $target_genome2 = $tempDir . 't2.genome';
		
        system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE $insertSize, $insertStd -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");

        system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outSoloSam1");

        #only printout aligned reads
		system("grep -v \'ZS:Z:NM\' $outPairedSam1 > $target_pair_sam");
		system("grep -v \'ZS:Z:NM\' $outSoloSam1 > $target_solo_sam");
    
        #run samtools to bam, merge sort and index
		system("samtools view -bS $target_pair_sam > $target_pair_bam");
		system("samtools view -bS $target_solo_sam > $target_solo_bam");
		system("samtools merge -f $target_bam $target_solo_bam $target_pair_bam");
		system("samtools sort $target_bam $target_sorted");
		system("samtools index $target_sorted_bam"); 
  
		#make readgroups and remove duplicates
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard MarkDuplicates INPUT=$target_sorted_bam OUTPUT=$target_duped_bam METRICS_FILE=$bam.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard AddOrReplaceReadGroups INPUT=$target_duped_bam OUTPUT=$target_rg_bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
	    system("samtools sort $target_rg_bam $bam");
	    system("samtools index $bam.bam");
        system("samtools idxstats $bam.bam > $bam.stats");

        #generate a genome file, generate coverage with bedtools
		system("awk \< $bam.stats \'\{print \$1\"\\t\"\$2 \}\' \> $target_genome");
		system("sed \'\$d\' $target_genome \> $target_genome2");
        system("mv $target_genome2 $bam.genome");
        system("genomeCoverageBed -ibam $bam.bam -d -g $bam.genome \> $bam.cov");

        
##### run gatk target realigner
        system("samtools faidx $ref");
        my $dict = substr ($ref, 0, -2) . "dict";
        my $intervals = substr ($ref, 0, -2) . "intervals";
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $ref -I $bam.bam -o $intervals");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $ref -I $bam.bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");


#remove temp files, get ready for phasing
        
       system("mv $bam.realigned.bam $tempBam");
	   system("samtools index $tempBam");

       system("rm -r $tempDir");#target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.bam target.sorted.bam target.sorted.bam.bai target.duped.bam target.rg.bam target.genome");
         

    print "\n\nfinished with third alignment (2nd novo240 on 2nd mapsembler) for $lib \n\n\n";
    
    

}

sub makepileups3 {
	my ($lib) = @_;
        
        #define directories
	my $resultsDir = $mapsembleOutDir . $lib . '/' . 'map2' . '/';
	my $contigDir = $resultsDir . 'contigs2' . '/';
        
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean2.fa';
        my $tempRef = $resultsDir . $lib . '.clean3.fa';
        my $bam = $contigDir . $lib . '_m2n240_2r.bam';
        my $tempDir = $contigDir . 'temporary/';
        system("mkdir $tempDir");
        
        #define outfiles
        my $pileup = $tempDir  . $lib . '.pileup';
        my $tempPileup = $tempDir  . $lib . '.temp.pileup';
        my $varfile = $tempDir  . $lib . '.varscan'; 
        my $tempVarfile = $tempDir  . $lib . '.temp.varscan'; 
 	
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
    #run samtools to generate pileup

    system("samtools mpileup -AB -f $ref $bam > $pileup");

    #run varscan 
   
    system("java -Djava.io.tmpdir=$tempDir -jar $VarScan pileup2cns $pileup --min-coverage 1 --min-var-freq 0.4 --min-reads2 1 --p-value 0.17 > $varfile");
    
    ####### parse VarScan to new fasta (note - this was adapted from another file so rather than renaming variables, i just 
            # used the definitions already in the script. therefore, some files are referenced by 2 different handles
    
    # name an out file                                                                                                         
	my $hap0fil = $tempDir . $lib . '.haplo0.fasta';
	open(HAP0, ">$hap0fil");
	my $hap1fil = $tempDir . $lib . '.haplo1.fasta';
	open(HAP1, ">$hap1fil");

    # name vcf file                                                                                                            
	
	my $rawvcf = $varfile;
	my $processvcf = $tempDir . $lib . '.varscan.process.vcf';
	system("grep -v \'Chrom\' $rawvcf \| sed \'s\/\%\/\/\' \> $processvcf");
	my $vcffil = $tempDir . $lib . '.varscan.process.vcf';
	open(VCF, "<$vcffil");
	my @vcf = <VCF>;
	chomp;
    
    # handle IO with bioperl                                                                                                   
        my $seqio = Bio::SeqIO-> new(
                             -file     => $ref,
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

      		    if (($tmp[3] eq $tmp[2]) && ($tmp[6] < 50)) {                                                               
	       			    $a0 = $tmp[3];
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

		    $h0indlen = $h0indlen + length($a0)-1;
		    $h1indlen = $h1indlen + length($a1)-1;
		    
		}
	    }

	    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
	    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
	}
	close HAP0;
	close HAP1;
	close VCF;

##### remove Zs to resolve indels, number of Z's plus equal number of following base pairs need to be removed

	my $rawFasta = $tempDir . $lib . '.haplo0.fasta';
        my $rawFasta2 = $tempDir . $lib . '.haploZ.fasta';
        system("sed \'s\/z\/Z\/g\' $rawFasta > $rawFasta2");
        system("mv $rawFasta2 $rawFasta");

        my $Zs = $tempDir . $lib . '.Zs.txt';
        open(FAS, "<$rawFasta");
        open(OUT, ">$Zs");  
        my @zcount;
        while(<FAS>){
	    chomp(my $line = $_);
            while ( $line  =~ /(Z+)/g) {
        push @zcount, length($1);
        print OUT length($1), "\n";
        }
	}
        my $maxZcount = max @zcount;
        print "$maxZcount \t $lib \n";
        close FAS;
        close OUT; 
        
        my $cleanFasta = $tempDir . $lib . '.clean.fa';
        my $precleanFasta = $tempDir . $lib . '.preclean.fa';
        system("cp $rawFasta $precleanFasta");

        for (my $i = $maxZcount; $i > 0; $i--) {
	    system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #      sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta
            system("mv $cleanFasta $precleanFasta");
	}

         ### do it again, if there were Z's in the way that prevented Zs from getting removed

	for (my $i = $maxZcount; $i > 0; $i--) {
            system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #          sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta                                                                                          
            system("mv $cleanFasta $precleanFasta");
        }
         
        system("mv $precleanFasta $cleanFasta");
        system("cp $cleanFasta $tempRef");

    ## cleanup
    system("rm -r $tempDir");

    print "\n\nfinished with third cleanup (2nd fix of 2nd mapsembler plus 1 fix of 1st mapsembler)  for $lib \n\n\n";

}

sub makeFourthAlignment {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $mapsembleOutDir . $lib . '/' . 'map2' . '/';
	my $contigDir = $resultsDir . 'contigs2' . '/';
	my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean3.fa';       
        my $bam = $contigDir . $lib . '_m2n240_3';
        
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
        my $clean2 = $fastqdir . $lib . '_2_final.txt';
        my $cleanu = $fastqdir . $lib . '_u_final.txt';
	    
        #define ref and bam for next step    
        my $tempBam = $contigDir . $lib . '_m2n240_3r.bam';
    

###### run novoalign on mapsembled reference
        #make out directories and files

        #generate index
		my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";
		system("$novoindex $indexed_assemblies_in_target $ref");

		#run novoalign
		#define temp files
		
		my $outPairedSam1 = $tempDir . 'ops1';
		my $outSoloSam1 = $tempDir . 'oss1';
		my $target_pair_sam = $tempDir . 'tp.sam';
		my $target_solo_sam = $tempDir . 'ts.sam';
		my $target_pair_bam = $tempDir . 'tp.bam';
		my $target_solo_bam = $tempDir . 'ts.bam';
		my $target_bam = $tempDir . 't.bam';
		my $target_sorted = $tempDir . 't.sorted';
		my $target_sorted_bam = $tempDir . 't.sorted.bam';
		my $target_duped_bam = $tempDir . 't.duped.bam';
		my $target_rg_bam = $tempDir . 't.rg.bam';
		my $target_genome = $tempDir . 't.genome';
		my $target_genome2 = $tempDir . 't2.genome';
		
        system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE $insertSize, $insertStd -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");

        system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 240 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outSoloSam1");

        #only printout aligned reads
	system("grep -v \'ZS:Z:NM\' $outPairedSam1 > $target_pair_sam");
	system("grep -v \'ZS:Z:NM\' $outSoloSam1 > $target_solo_sam");
    
        #run samtools to bam, merge sort and index
	system("samtools view -bS $target_pair_sam > $target_pair_bam");
	system("samtools view -bS $target_solo_sam > $target_solo_bam");
	system("samtools merge -f $target_bam $target_solo_bam $target_pair_bam");
	system("samtools sort $target_bam $target_sorted");
	system("samtools index $target_sorted_bam"); 
  
       	#make readgroups and remove duplicates
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard MarkDuplicates INPUT=$target_sorted_bam OUTPUT=$target_duped_bam METRICS_FILE=$bam.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard AddOrReplaceReadGroups INPUT=$target_duped_bam OUTPUT=$target_rg_bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
	system("samtools sort $target_rg_bam $bam");
	system("samtools index $bam.bam");
        system("samtools idxstats $bam.bam > $bam.stats");

        #generate a genome file, generate coverage with bedtools
	system("awk \< $bam.stats \'\{print \$1\"\\t\"\$2 \}\' \> $target_genome");
	system("sed \'\$d\' $target_genome \> $target_genome2");
        system("mv $target_genome2 $bam.genome");
        system("genomeCoverageBed -ibam $bam.bam -d -g $bam.genome \> $bam.cov");

        
##### run gatk target realigner
        system("samtools faidx $ref");
        my $dict = substr ($ref, 0, -2) . "dict";
        my $intervals = substr ($ref, 0, -2) . "intervals";
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $ref -I $bam.bam -o $intervals");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $ref -I $bam.bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");


#remove temp files, get ready for phasing
        
       system("mv $bam.realigned.bam $tempBam");
       system("samtools index $tempBam");

       system("rm -r $tempDir");#target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.bam target.sorted.bam target.sorted.bam.bai target.duped.bam target.rg.bam target.genome");
         

    print "\n\nfinished with third alignment (2nd novo240 on 2nd mapsembler) for $lib \n\n\n";
    
    

}

sub makepileups4 {
	my ($lib) = @_;
        
        #define directories
	    my $resultsDir = $mapsembleOutDir . $lib . '/' . 'map2' . '/';
	    my $contigDir = $resultsDir . 'contigs2' . '/';
        
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean3.fa';
        my $tempRef = $resultsDir . $lib . '.clean4.fa';
        my $bam = $contigDir . $lib . '_m2n240_3r.bam';
        my $tempDir = $contigDir . 'temporary/';
        system("mkdir $tempDir");
        
        #define outfiles
        my $pileup = $tempDir  . $lib . '.pileup';
        my $tempPileup = $tempDir  . $lib . '.temp.pileup';
        my $varfile = $tempDir  . $lib . '.varscan'; 
        my $tempVarfile = $tempDir  . $lib . '.temp.varscan'; 
 	
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
	    my $clean2 = $fastqdir . $lib . '_2_final.txt';
	    my $cleanu = $fastqdir . $lib . '_u_final.txt';
    #run samtools to generate pileup

    system("samtools mpileup -AB -f $ref $bam > $pileup");

    #run varscan 
   
    system("java -Djava.io.tmpdir=$tempDir -jar $VarScan pileup2cns $pileup --min-coverage 1 --min-var-freq 0.4 --min-reads2 1 --p-value 0.17 > $varfile");
    
    ####### parse VarScan to new fasta (note - this was adapted from another file so rather than renaming variables, i just 
            # used the definitions already in the script. therefore, some files are referenced by 2 different handles
    
    # name an out file                                                                                                         
	my $hap0fil = $tempDir . $lib . '.haplo0.fasta';
	open(HAP0, ">$hap0fil");
	my $hap1fil = $tempDir . $lib . '.haplo1.fasta';
	open(HAP1, ">$hap1fil");

    # name vcf file                                                                                                            
	
	my $rawvcf = $varfile;
	my $processvcf = $tempDir . $lib . '.varscan.process.vcf';
	system("grep -v \'Chrom\' $rawvcf \| sed \'s\/\%\/\/\' \> $processvcf");
	my $vcffil = $tempDir . $lib . '.varscan.process.vcf';
	open(VCF, "<$vcffil");
	my @vcf = <VCF>;
	chomp;
    
    # handle IO with bioperl                                                                                                   
        my $seqio = Bio::SeqIO-> new(
                             -file     => $ref,
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

      		    if (($tmp[3] eq $tmp[2]) && ($tmp[6] < 50)) {                                                               
	       			    $a0 = $tmp[3];
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

		    $h0indlen = $h0indlen + length($a0)-1;
		    $h1indlen = $h1indlen + length($a1)-1;
		    
		}
	    }

	    print HAP0 ">" . $h0name . "\n" . $h0seq . "\n";
	    print HAP1 ">" . $h1name . "\n" . $h1seq . "\n";
	}
	close HAP0;
	close HAP1;
	close VCF;

##### remove Zs to resolve indels, number of Z's plus equal number of following base pairs need to be removed

	my $rawFasta = $tempDir . $lib . '.haplo0.fasta';
        my $rawFasta2 = $tempDir . $lib . '.haploZ.fasta';
        system("sed \'s\/z\/Z\/g\' $rawFasta > $rawFasta2");
        system("mv $rawFasta2 $rawFasta");

        my $Zs = $tempDir . $lib . '.Zs.txt';
        open(FAS, "<$rawFasta");
        open(OUT, ">$Zs");  
        my @zcount;
        while(<FAS>){
	    chomp(my $line = $_);
            while ( $line  =~ /(Z+)/g) {
        push @zcount, length($1);
        print OUT length($1), "\n";
        }
	}
        my $maxZcount = max @zcount;
        print "$maxZcount \t $lib \n";
        close FAS;
        close OUT; 
        
        my $cleanFasta = $tempDir . $lib . '.clean.fa';
        my $precleanFasta = $tempDir . $lib . '.preclean.fa';
        system("cp $rawFasta $precleanFasta");

        for (my $i = $maxZcount; $i > 0; $i--) {
	    system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #      sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta
            system("mv $cleanFasta $precleanFasta");
	}

         ### do it again, if there were Z's in the way that prevented Zs from getting removed

	for (my $i = $maxZcount; $i > 0; $i--) {
            system("sed \'s\/Z\\\{$i\\\}\[A|C|G|T|N\]\\\{$i\\\}\/\/gI\' $precleanFasta > $cleanFasta");
         #          sed  's /Z \ {7  \ } [A|C|G|T|N ] \ {7  \ } / /gI ' uoel.haplo0.fasta                                                                                          
            system("mv $cleanFasta $precleanFasta");
        }
         
        system("mv $precleanFasta $cleanFasta");
        system("cp $cleanFasta $tempRef");

    ## cleanup
    system("rm -r $tempDir");

    print "\n\nfinished with fourth cleanup (3rd fix of 2nd mapsembler plus 1 fix of 1st mapsembler)  for $lib \n\n\n";

}


sub makeFinalAlignment {
	my ($lib) = @_;
        
        #define directories
        my $resultsDir =  $mapsembleOutDir . $lib . '/' . 'map2' . '/';
	    my $contigDir = $resultsDir . 'contigs2' . '/';
	    my $tempDir = $contigDir . 'temp/';
        system("mkdir $tempDir");
       
        #define reference and alignment bams
        my $ref = $resultsDir . $lib . '.clean4.fa';       
        my $bam = $contigDir . $lib . '_m2n240_4';
        
        #define fastqs
        my $clean1 = $fastqdir . $lib . '_1_final.txt';
	    my $clean2 = $fastqdir . $lib . '_2_final.txt';
	    my $cleanu = $fastqdir . $lib . '_u_final.txt';
	    
	    #define ref and bam for next step
	    
	    my $tempBam = $contigDir . $lib . '_m2n240_4r.bam';
    

###### run novoalign on mapsembled reference
        #make out directories and files

        #generate index
		my $indexed_assemblies_in_target =  substr ($ref, 0, -2) . "nix";
		system("$novoindex $indexed_assemblies_in_target $ref");

		#run novoalign
		#define temp files
		
		my $outPairedSam1 = $tempDir . 'ops1';
		my $outSoloSam1 = $tempDir . 'oss1';
		my $target_pair_sam = $tempDir . 'tp.sam';
		my $target_solo_sam = $tempDir . 'ts.sam';
		my $target_pair_bam = $tempDir . 'tp.bam';
		my $target_solo_bam = $tempDir . 'ts.bam';
		my $target_bam = $tempDir . 't.bam';
		my $target_sorted = $tempDir . 't.sorted';
		my $target_sorted_bam = $tempDir . 't.sorted.bam';
		my $target_duped_bam = $tempDir . 't.duped.bam';
		my $target_rg_bam = $tempDir . 't.rg.bam';
		my $target_genome = $tempDir . 't.genome';
		my $target_genome2 = $tempDir . 't2.genome';
		
        system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE $insertSize, $insertStd -t 90 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outPairedSam1");

        system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 90 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $outSoloSam1");

        #only printout aligned reads
		system("grep -v \'ZS:Z:NM\' $outPairedSam1 > $target_pair_sam");
		system("grep -v \'ZS:Z:NM\' $outSoloSam1 > $target_solo_sam");
    
        #run samtools to bam, merge sort and index
		system("samtools view -bS $target_pair_sam > $target_pair_bam");
		system("samtools view -bS $target_solo_sam > $target_solo_bam");
		system("samtools merge -f $target_bam $target_solo_bam $target_pair_bam");
		system("samtools sort $target_bam $target_sorted");
		system("samtools index $target_sorted_bam"); 
  
		#make readgroups and remove duplicates
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard MarkDuplicates INPUT=$target_sorted_bam OUTPUT=$target_duped_bam METRICS_FILE=$bam.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard AddOrReplaceReadGroups INPUT=$target_duped_bam OUTPUT=$target_rg_bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
	system("samtools sort $target_rg_bam $bam");
	system("samtools index $bam.bam");
        system("samtools idxstats $bam.bam > $bam.stats");

        #generate a genome file, generate coverage with bedtools
		system("awk \< $bam.stats \'\{print \$1\"\\t\"\$2 \}\' \> $target_genome");
		system("sed \'\$d\' $target_genome \> $target_genome2");
        system("mv $target_genome2 $bam.genome");
        system("genomeCoverageBed -ibam $bam.bam -d -g $bam.genome \> $bam.cov");

        
##### run gatk target realigner
        system("samtools faidx $ref");
        my $dict = substr ($ref, 0, -2) . "dict";
        my $intervals = substr ($ref, 0, -2) . "intervals";
        system("java -Djava.io.tmpdir=$tempDir -Xmx32g -jar $picard CreateSequenceDictionary REFERENCE=$ref OUTPUT=$dict");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $ref -I $bam.bam -o $intervals");
        system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $ref -I $bam.bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");
 

#remove temp files, get ready for phasing
        
       system("mv $bam.realigned.bam $tempBam");
	   system("samtools index $tempBam");

       system("rm -r $tempDir");#target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.bam target.sorted.bam target.sorted.bam.bai target.duped.bam target.rg.bam target.genome");
         

    print "\n\nfinished with fourth alignment (novo90) for $lib \n\n\n";
	system("date");
    

}


}

	
