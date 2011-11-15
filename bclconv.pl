#!/usr/bin/perl

# Manually retrieve dna and/or qualities from bcl files
# GPL v.2.0, 2011-11-11
# (C) Roel Kluin, Arno Velds, NKI
# version 0.5

use strict;
use warnings;

use FindBin;
use FileHandle;
use Getopt::Long;
use Pod::Usage;

# Specify default illumina primer version here
my $defpv = undef;

use lib "$FindBin::Bin/modules";
use illumina qw(Ndx2primer primer2Ndx primer2v1Ndx v2plen v2primer v1primer vprimer);

my $tiles = 1306;
my $nreads = 1000000;
my $readlen = 50;
my ($c1, $c2, $version, $lane, $path, $help, $man, $qual, $both, $uniq, $ill);

Getopt::Long::Configure("no_ignore_case", "prefix_pattern=(--|-)");
GetOptions('lane|l:i' => \$lane, 'start-cycle|s:i' => \$c1, 'end-cycle|e:i' => \$c2, 'tiles|t=s' => \$tiles,
	'version|v:i' => \$version, 'path|p=s' => \$path, 'qualities|q' => \$qual, 'nreads|n:i' => \$nreads,
	'both|b' => \$both, 'uniq|u' => \$uniq, 'read-length|L:i' => \$readlen, '50' => sub { $readlen = 50; },
	'75' => sub { $readlen = 75; }, '100' => sub { $readlen = 100; }, 'v1' => sub { $version = 1; },
	'v2' => sub { $version = 2; }, 'help|?|h' => \$help, 'man' => \$man, 'i|illumina' => \$ill);

pod2usage(-verbose => 2) if defined $man;
pod2usage(0) if defined $help;

sub exit_msg {
	my $msg = shift;
	pod2usage(-exitval => "NOEXIT");
	die "$msg\n\n";
}

$path = shift if not defined $path;
if ((not $lane) && $path && (not -d $path)) {
	$lane = $path;
	$path = shift;
}
exit_msg("No path given") if not defined $path;
exit_msg("No such directory: $path") unless -d $path;

$lane = shift unless defined $lane;

$path =~ m/L00(\d)\/?/;
if ($1) {
	exit_msg("specified lane does not match the one in the given path") if ($lane && ($1 != $lane));
	$lane = $1;
} else {
	unless (defined $lane) {
		$lane = shift;
		exit_msg("please specify lane") unless defined $lane;
	}
	exit_msg("lane $lane is not between 1 and 8") unless $lane > 0 && $lane < 9;

	if (-d $path."/L00$lane") {
		$path .= "/L00$lane";
	} else {
		my @list = `ls -d1rt "$path/Data/Intensities/B"*"/L00$lane/" `;
		exit_msg("found no L00$lane/ subdirectory in $path") if not defined $list[0];
		exit_msg(join "\n", ("found multiple L00$lane subdirectories in $path:",
			@list, "specify the one to use with --path")) if @list > 1;
		$path = $list[0];
		chomp $path;
	}
}
exit_msg("Illumina has no version $version") if not exists $illumina::v2plen{$version};

$c1 = $readlen + 1 if not defined $c1;
if (not defined $c2) {
	if (not defined $version) {
		exit_msg("no illumina primer or end cycle specified") if not defined $defpv;
		warn "no illumina primer or end cycle specified".
			" assuming it's a v$defpv illumina primer.\n";
		# FIXME: can we automate choice between v1 and v2 primers? maybe file size?
		$version = $defpv;
	}
	$c2 = $c1 + $illumina::v2plen{$version};
}


my @files;
warn "reading cyles $c1..$c2 for tiles $tiles\n";
foreach my $tile (split /,/, $tiles) {
	exit_msg("Invalid tile: $tile") if $tile !~ /^[0-9]+$/;
	for ($c1..$c2) {
		my $f = "$path/C$_.1/s_$lane"."_$tile.bcl";
		exit_msg("Invalid tile: $f") if	not -f $f;
		push @files, $f;
	}
}

my @bases = qw/A C G T/;

#print $nclust,"\n";
my (@qual, @res, $func);
sub bases_only {
	my ($b, $r) = @_;
	$b = ord(pack("b8", substr($b,0,2)."00000"));
	$res[$r] .= $bases[$b];
}

sub bases_and_qual {
	my ($b, $r) = @_;
	$qual[$r] .= chr(33 + ord(pack("b8", substr($b,2)."00")));
	$b = ord(pack("b8", substr($b,0,2)."00000"));
	$res[$r] .= $bases[$b];
}

sub qual_only {
	my ($b, $r) = @_;
	$qual[$r] .= chr(33 + ord(pack("b8", substr($b,2)."00")));
}

if (defined $both) {
	@res = ("") x $nreads;
	@qual = ("") x $nreads;
	$func = \&bases_and_qual;
} elsif (defined $qual) {
	@qual = ("") x $nreads;
	$func = \&qual_only;
} else {
	@res = ("") x $nreads;
	$func = \&bases_only;
}
foreach my $f (@files) {
	my $fh = FileHandle->new($f, "r");
	binmode $fh or exit_msg("couldn't open $f: $!\n");

	my $buff;
	my $ret = read($fh ,$buff, 4);
	exit_msg("error".(defined $ret ? " ($ret)" : '')." while reading $f: $!\n") if not $ret;

	my $nclust = unpack("L", $buff);
	warn "nclust = $nclust\n";

	seek($fh, $nreads, 0) or exit_msg("couldn't seek to $nreads bytes in $f: $!\n");
	for (1 .. $nreads) {
		read($fh, $buff, 1) or exit_msg("error while reading $f ($_): $!\n");
		my $bits = unpack("b8", $buff);
		$func->($bits, $_);
	}
}
if (defined $ill) {
	my %t;
	print "count\tindex\tDNA_index\n";
	my ($c2pfunc, $p2nh);
	if ($version == 1) {
		($c2pfunc, $p2nh) = (\&illumina::v1primer, \%illumina::primer2v1Ndx);
	} else {
		($c2pfunc, $p2nh) = (\&illumina::v2primer, \%illumina::primer2Ndx);
	}
	for (0..$nreads) {
		my $r = $c2pfunc->($res[$_]);
		die $r if ((defined $r) && (not exists $p2nh->{$r}));
		$t{defined $r ? $p2nh->{$r}."\t".$r : "-\tunmatched"}++;
	}
	print "$t{$_}\t$_\n" foreach (sort {$t{$a} <=> $t{$b}} keys %t);
} elsif (defined $uniq) {
	my %t;
	$t{($both ? $res[$_]."\t".$qual[$_] : ($qual ? $qual[$_] : $res[$_]))."\n"}++ for (0..$nreads);
	print "$t{$_}\t$_" foreach (sort {$t{$a} <=> $t{$b}} keys %t);
} else {
	print "".($both ? $res[$_]."\t".$qual[$_] :
		($qual ? $qual[$_] : $res[$_]))."\n" for (0..$nreads);
}

__END__

=head1 NAME

bclconv.pl - read sequence from bcl files

=head1 SYNOPSIS

bclconv.pl [options] <path>

=head1 OPTIONS

=over 8

=item B<-l|--lane=>

Specify lane

=item B<-p|--path=>

You can specify the path

=item B<-v1|-v2|-v=|--version=>

Specify an Illumina primer version

=item B<-50|-75|-100|-L=|--read-length=>

Specify read length [default 50]

=item B<-s|--start-cycle=>

Specify start cycle [default <read length> +1]

=item B<-e|--end-cycle=>

Specify end cycle [default <start cycle> + primerlength]

=item B<-t|--tiles=>

Specify a comma seperated list of tiles [default 1306]

=item B<-q|--qualities>

Show qualities instead of bases [default: no ]

=item B<-b|--both>

Show qualities alongside bases (overrides -q) [default: no ]

=item B<-u|--uniq>

Show only uniq  [default: no ]

=item B<-h|-?|--help>

Print options message

=item B<--man>

Print manual page

=item B<-n|--reads>

Specify the number of reads

=item B<-i|--illumina>

only count and print illumina primers, nr and dna.

=back

=head1 DESCRIPTION

Read bcl files and read from specified cycles, tiles the bases and or qualities.
This is especially useful to read the primers and their qualities as the current
Casava version doesn't provide this information.

=head1 EXAMPLE

=over 8

=item count unique index primers (v2) of an 50bp run, lane $lane:

./bclconv.pl -v2 -u -p $path_to_rundir -l $lane

=item count uniq indexes of 5000 version 1 primers in 75bp, tile 1306:

./bclconv.pl -v1 -75 -u -l 7 -n 5000 -t 1101 $path_to_rundir2 $lane

=item show index and quality for the first 100 reads:

./bclconv.pl -v1 -100 -b -n=100 $path_to_rundir3/Data/Intensities/BaseCalls/L001

=back

=cut
