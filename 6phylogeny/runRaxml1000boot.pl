use warnings;
use strict;
use File::Basename;
#chodon at berkeley dot edu
## note that this ignores the first file in the directory...
die (qq/
(first try! Nov 6 2015)
Runs raxml on individual gene fasta files in a directory
Usage: runRaxml.pl <directory with individual gene alignments> 
    \n/) if !@ARGV;

my $orig_directory;

#this just makes sure that the directory ends with a slash, even if the user didn't type it
if ($ARGV[0] =~ m/\/$/ ){
    $orig_directory = $ARGV[0];
}
else {
    $orig_directory = $ARGV[0] . "/";
}

#names your output directories
my $result_dir = $orig_directory . 'raxmlOut/';
my $result2 = $orig_directory . 'raxmlExtraOut/';
my $result3 = $orig_directory . 'bootExtraOut';
mkdir $result_dir unless -e $result_dir;
mkdir $result2 unless -e $result2;
mkdir $result3 unless -e $result3;

#tells to make an array of all the files in your folder that end with a specified ending
opendir(DIR,$orig_directory);
my $ending = '.fasta';
my @orig_files = grep { $_ =~ /$ending/ } readdir(DIR) ;
closedir(DIR);
chomp (@orig_files);
die(qq/\nHmm...did you specify the correct directory? \n\n/) if (scalar (@orig_files) == 0);

#308 nuclear plus one chloroplast gene alignments
my @geneNumber = (0..308);

foreach my $geneNumber (@geneNumber) {
    my $prefile = $orig_files[$geneNumber];
    my $file = basename($prefile) =~ /(\S+).fasta/;
    my $lib_name = $1;
    my $full_name = $1 . ".fasta";
    my $adjusted = $1 . ".lim.fa";
    print "gene working on is $full_name\n";
    my $getnoemptylines = '\'^-+$\'';
    my $getRidofEmptys = '\'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}\'';
    system("grep -Ev $getnoemptylines $full_name | awk $getRidofEmptys > $adjusted");
    my $out_file = $1;
    my $bootout = "RAxML_bootstrap." . $1;
    my $finalML = "RAxML_bestTree." . $1;
    my $newOut = $geneNumber . "/" . "RAxML_bootstrap.final";
    system("raxmlHPC-SSE3 -m GTRGAMMA -N 1000 -s $adjusted -n $out_file -p 12345 -b 12345");
    system("mkdir $geneNumber");
    my $bootstrapOut = 'RAxML_bootstrap.' . $out_file;
    system("mv $bootstrapOut $newOut");
    my $infoOut = 'RAxML_info.' . $out_file;  
    system("mv $infoOut bootExtraOut/");
    system("mv $adjusted raxmlExtraOut/");
    my $reduced = $adjusted . '.reduced';
    system("mv $reduced raxmlExtraOut/");
}

