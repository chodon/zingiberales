use strict;
use Getopt::Std;

#### works! but 
#### need to add if statement so that i goes only till end of file, not all way to 62###
####

# usage: perl parse_gff_to_bed.pl -i gff3_file -o outputFasta -f inputFasta -b outputBed
my %opts = (i=>undef, o=>undef, f=>undef, b=>undef);
getopts('i:o:f:b:', \%opts);

#my $gff_to_bed = 'musa_cds_to_bed_to_excel_part2';
my $cdsFasta = $opts{f};
my $out = $opts{o};
my $outbed = $opts{b};
my $genefile = 'genefile';
my $input=$opts{i};
my $shortGFF = 'shortGFF';
my $interOut = 'tempOut';
my $posGFF = 'posGFF';
my $negGFF = 'negGFF';
#awk '$3 ~ /mRNA/ {print $9}' coffea_canephora.gff3 | awk 'BEGIN{FS=";"}{print $1}' | sed 's/ID=//g'
system("awk \'\$3 ~ /mRNA/ \{print \$9\}\' $input \| awk \'BEGIN\{FS=\";\"\}\{print \$1\}\' \| sed \'s\/ID=\/\/g\' | sort | uniq > $genefile");
#grep 'CDS' coffea_canephora.gff3 | grep -v 'Name' | awk 'BEGIN{FS="\t"}{print $9"\t"$1"\t"$4"\t"$5}' | sort | sed 's/Parent=//g'
system("grep \'CDS\' $input \| grep -v \'Name\' \| awk \'BEGIN\{FS=\"\t\"\}\$7~/+/\{print \$9\"\t\"\$1\"\t\"\$4\"\t\"\$5\}\' | sort -V | sed \'s/Parent=//g\' > $posGFF");
system("grep \'CDS\' $input \| grep -v \'Name\' \| awk \'BEGIN\{FS=\"\t\"\}\$7~/-/\{print \$9\"\t\"\$1\"\t\"\$4\"\t\"\$5\}\' |  sed \'s/Parent=//g\' > $negGFF");
system("cat $negGFF $posGFF > $shortGFF");

open(IN, "<$shortGFF") or die ("no such file!");
open(OUT,">>$interOut");
open(GENE,"<$genefile");

my @gene = <GENE>;
chomp(@gene);


foreach my $gene (@gene) {

    my $temp = 'temporary';
    system("grep \'$gene\' $shortGFF \> $temp ");
    open (GFF,"<$temp");
    my @line = <GFF>;

    my $first = $line[0];
    chomp $first;
    my @parts0 = split(/\t/, $first);
    my $start0 = 0;
    my $g0 = $parts0[0];
    my $s0 = $parts0[1];
    my $o0 = $parts0[2];
    my $e0 = $parts0[3];
    my $ns0 = (($e0 - $o0)+1); 
    my @nsarray = ();
    push @nsarray, $ns0;

    print OUT   $g0,"\t",$start0,"\t",$ns0,"\n";

    foreach my $i (1..100) {
	my $numb = $line[$i];
	chomp $numb;
	my @partsi = split(/\t/, $numb);
	my $gi = $partsi[0];
	my $oi = $partsi[2];
	my $ei = $partsi[3];
	my $nsi = pop(@nsarray);
	my $nsiend = $nsi + ($ei - $oi) + 1;
	push @nsarray, $nsiend;
	print OUT $gi,"\t",$nsi,"\t",$nsiend,"\n"; 
    }
    
}

close OUT;
close GFF;
close GENE;

#this removes all the extra lines that are from line #57
system("grep \'GSC\' $interOut > $outbed");
#grep -f list coffea_canephora.gff3 | awk '$3~/mRNA/{print $9}' | awk 'BEGIN{FS=";"}{print$1"\t"$3}' | sed 's/ID=//' | sed 's/Name=//'

#this get the gene names from the gff and also the GSC id number
my $genelist2 = 'genelist2';
system(" awk \'\$3~/mRNA/\{print\$9\}\' $input | awk \'BEGIN\{FS=\";\"\}\{print \$1\"\t\"\$3\}\' | sed \'s/ID=//\' | sed \'s/Name=//\' > $genelist2");

#this makes the name and id number into a hash and then replaces the gsc id number with the gene name in the bedfile

my $bedfile = $outbed;
my $hashfile = 'genelist2';
    
open(my $hfile, '<', $hashfile) or die "Cannot open $hashfile: $!";
my %names;
while (my $line = <$hfile>) {
    chomp($line);
    my ($oldname, $newname) = split /\t/, $line;
    $names{$oldname} = $newname;
}
close $hfile;

my $regex = join '|', sort { length $b <=> length $a } keys %names;
$regex = qr/$regex/;


open(my $bfile, '<', $bedfile) or die "Cannot open $bedfile: $!";
open(my $out4, '>', 'out3.bed') or die "Cannot open out3.bed: $!";

while (my $line = <$bfile>) {
    chomp($line);
    $line =~ s/($regex)/$names{$1}/g;
    print $out4 $line, "\n";
}

close $out4;
close $bfile;

system("mv out3.bed $outbed");

#sysem("bedtools getfasta -fi $cdsFasta -bed out3.bed -fo shortOutCoffea.fasta");
