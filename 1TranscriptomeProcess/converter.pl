 ################## pertains to FASTQ parsing only ##################
use strict;
use Bio::SeqIO;

my @lib = ('KNKV');
#, 'OYLU', 'KNKV', 'JNUB', 'LSKK', 'UOEL', 'BDJQ');
#'COST');
#, 'JQCX', 'OYLU', 'KNKV', 'JNUB', 'LSKK', 'UOEL', 'BDJQ'); 

foreach my $lib (@lib) {

    my $dir = '/media/Data/chodon/zingiberales/' . $lib . '/' . $lib . '/';
    my @readnum = ('1', '2');

  # specifies the Illumina variant
    foreach my $readnum (@readnum) {
	my $zip = $dir . $lib . '_' . $readnum . '.fastq.bz2';
        
	system("bunzip2 $zip");
	
	my $in = Bio::SeqIO->new(-format    => 'fastq-illumina',
				 -file      => $dir . $lib . '_' . $readnum . '.fastq');

  # simple 'fastq' format defaults to 'sanger' variant
	my $out = Bio::SeqIO->new(-format    => 'fastq',
				  -file      => '>' . $dir . $lib . '_sanger_' . $readnum . '.fastq');

  # $seq is a Bio::Seq::Quality object
	while (my $seq = $in->next_seq) {
	    $out->write_seq($seq);  # convert Illumina 1.3 to Sanger format
	}

	print "\n\nconverted read number $readnum of $lib\n\n";
}
}
