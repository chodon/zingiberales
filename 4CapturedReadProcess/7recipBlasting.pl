use warnings;
use strict;
use Getopt::Std;

###written by sonal dot singhal at gmail dot com in 2012, adapted slightly by chodon at gmail dot com
###note -- now needs to be adapted for updated blast i.e. makeblastdb instead of formatdb etc

my %opts = (b=>undef, a=>undef, o=>undef);
getopts('b:a:o:', \%opts);


my $bait = $opts{b};
my $assembly = $opts{a};
#where you want output to go
my $dir = $opts{o};
#number of processors you can use for this exercise
my $np = 4;
#evalue limit; don't want matches unless they are at least this good
my $evalue = 1e-40;

#################################
# and here is the actual script #
#################################

mkdir($dir) unless(-d $dir);
my $call1 = system("formatdb -i $bait -p F");
my $call2 = system("formatdb -i $assembly -p F");

my $out = $dir . "finalBaitAssembly.fa";
my $blast1 = $dir . "blastAssemblyBait.out";
my $blast2 = $dir . "blastBaitAssembly.out";
my $call3 = system("blastall -p blastn -d $bait -i $assembly -a $np -e $evalue -m 8 -o $blast1 -b 10") unless (-f $blast1);
my $call4 = system("blastall -p blastn -d $assembly -i $bait -a $np -e $evalue -m 8 -o $blast2 -b 10") unless (-f $blast2);

my $r1 = parseBlast($blast1);
my $r2 = parseBlast($blast2);
my $seq1 = parseSeq($assembly);
my $seq2 = parseSeq($bait);
my %seq1 = %{$seq1}; my %seq2 = %{$seq2};
my %r1 = %{$r1}; my %r2 = %{$r2};

open(OUT, ">$out");
foreach my $c (keys %seq2) {
    if ($r2{$c}) {
	my $numMatches = scalar(keys %{$r2{$c}});
	#if it has one high match
	if ($numMatches == 1) {
	    my @match = keys %{$r2{$c}};
	    my $match = $match[0];
	    if ($r1{$match}{$c}) {
		#this is a recip blast match!
		if (length($seq2{$c}) > length($seq1{$match})) {
		    my $diff = -(length($seq1{$match}) - length($seq2{$c}));
		    print "$c: Match but original seq $diff bp longer.\n";
		    print OUT ">", $c, "\n", $seq2{$c}, "\n";
		}
		else {
		    my $diff = length($seq1{$match}) - length($seq2{$c});
		    print "$c: Match and assembled seq $diff bp longer.\n";
		    print OUT ">", $c, "\n", $seq1{$match}, "\n";
		}
	    }
	    else {
		print "$c: No match. Sad story.\n";
		print OUT ">", $c, "\n", $seq2{$c}, "\n";
	    }
	}
	#it has multiple high matches
	else {
	    my @match1 = keys %{$r2{$c}};
	    my @match2;
	    foreach (@match1) {
		push(@match2,$_) if $r1{$_}{$c};
	    }
	    @match2 = sort {length($seq1{$a}) <=> length($seq1{$b}) } @match2;
	    my $match = $match2[$#match2];
	    if (length($seq2{$c}) > length($seq1{$match})) {
		my $diff = -(length($seq1{$match}) - length($seq2{$c}));
		print "$c: Match but original seq $diff bp longer.\n";
		print OUT ">", $c, "\n", $seq2{$c}, "\n";
	    }
	    else {
		my $diff = length($seq1{$match}) - length($seq2{$c});
		print "$c: Match and assembled seq $diff bp longer.\n";
		print OUT ">", $c, "\n", $seq1{$match}, "\n";
	    }
	        
	}
    }
    else {
	print "Bummer! $c did not get assembled apparently. Will put given $c sequence into final file, instead.\n";
	print OUT ">", $c, "\n", $seq2{$c}, "\n";
    }
}

sub parseBlast {
    my ($file) = @_;

    open(IN, "<$file");
    my %r; my $score;
    while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	if ($r{$d[0]}) {
	    if ($d[11] == $score) {
		$r{$d[0]}{$d[1]}++;
	    }
	}
	else {
	    $r{$d[0]}{$d[1]}++;
	    $score = $d[11];
	}
    }
    close(IN);
    return(\%r);
}

sub parseSeq {
    my ($file) = @_;

    open(IN, "<$file");
    my %seq; my $c;
    while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
	    $c = $1;
	}
	else {
	    $seq{$c} .= $line;
	}
    }

    return(\%seq);
}
