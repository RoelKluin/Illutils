# GPL v.2.0, 2011-11-11
# (C) Roel Kluin, NKI
package illumina;
require Exporter;
my @ISA = qw(Exporter);

our %primer2v1Ndx = ("ATCACG" => 1, "CGATGT" => 2, "TTAGGC" => 3, "TGACCA" => 4, "ACAGTG" => 5,
	"GCCAAT" => 6, "CAGATC" => 7, "ACTTGA" => 8, "GATCAG" => 9, "TAGCTT" => 10,
	"GGCTAC" => 11, "CTTGTA" => 12);

# 1-12 have an A as 7th: primer v2.
our %Ndx2primer = (1 => "ATCACGA", 2 => "CGATGTA", 3 => "TTAGGCA", 4 => "TGACCAA", 5 => "ACAGTGA",
	6 => "GCCAATA", 7 => "CAGATCA", 8 => "ACTTGAA", 9 => "GATCAGA", 10 => "TAGCTTA",
	11 =>"GGCTACA", 12 => "CTTGTAA",
	13 => "AGTCAAC", 14 => "AGTTCCG", 15 => "ATGTCAG", 16 => "CCGTCCC", 18 => "GTCCGCA",
	19 => "GTGAAAC", 20 => "GTGGCCT", 21 => "GTTTCGG", 22 => "CGTACGT", 23 => "GAGTGGA",
	25 => "ACTGATA", 27 => "ATTCCTT");

#:18,23s/\("[ACTG]\{7\}"\) => \([0-9]\+\)\([,)]\)/\2 => \1\3/g
our %primer2Ndx = ("ATCACGA" => 1, "CGATGTA" => 2, "TTAGGCA" => 3, "TGACCAA" => 4, "ACAGTGA" => 5,
        "GCCAATA" => 6, "CAGATCA" => 7, "ACTTGAA" => 8, "GATCAGA" => 9, "TAGCTTA" => 10,
        "GGCTACA" => 11, "CTTGTAA" => 12,
	"AGTCAAC" => 13, "AGTTCCG" => 14, "ATGTCAG" =>15, "CCGTCCC" => 16, "GTCCGCA" => 18,
	"GTGAAAC" => 19, "GTGGCCT" => 20, "GTTTCGG" => 21, "CGTACGT" => 22, "GAGTGGA" => 23,
	"ACTGATA" => 25, "ATTCCTT" => 27);

#specifies per illumina version the primer length.
our %v2plen = ('1' => 5, '2' => 6);

sub v1primer {
	my $ret = shift;
	$ret = (substr $ret, 0, 6) if length($ret) == 7 and $ret =~ /A$/;
	return exists $primer2v1Ndx{$ret} ? $ret : undef;
}

sub v2primer {
	my $ret = shift;
	$ret .= 'A' if length($ret) == 6;
	return exists $primer2Ndx{$ret} ? $ret : undef;
}

my @EXPORT_OK = qw(%Ndx2primer %primer2Ndx %primer2v1Ndx %v2plen v2primer v1primer vprimer);
1;
