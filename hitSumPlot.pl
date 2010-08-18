#!/usr/bin/env perl
#
# original filename: hitSumPlot.pl
#
# Steven Bickerton
#  Dept. of Astrophysical Sciences, Princeton University
#  bick@astro.princeton.edu
#  Created: Sat Nov  8, 2008  19:31:53 EST
#  Host: bender.astro.Princeton.EDU
#  Working Directory: /Users/bick/usr/src/cdiffracSim
#


use strict;
use warnings;
use File::Basename;

use Getopt::Std;

use Local::Constants;
use Local::Astrotools;

my %opt;
getopts("f:r", \%opt);
my $exe = basename($0);
my $usage = "Usage: $exe [-f fresdir] [-r (relative-n)] hitSumFile n\n";

my ($hitfile, $n) = @ARGV;
die $usage unless $n;

my $fresdir = ($opt{f}) ? $opt{f} : "fresnelFiles";

# read in the hitsfile and get values for entry n
my @values;
my $i = 1;
open(HITFILE, "$hitfile");
while(<HITFILE>) {
    next if /^\#/;
    @values = split;
    my $compare = $opt{r} ? $i : $values[0];
    last if $compare == $n;
    $i++;
}
close(HITFILE);

die "No entry $n found in $hitfile\n" unless scalar @values;

# run eventplot with the values
my ($N, $efile, $index, $nk, $rank, $rad, $AU, $ip, $xcor, $cor0, $nxc, $xchi, $chi0, $rchi2, $rchi2n, $dof, $baryc, $fft) = @values;
my $fresfile = sprintf "$fresdir/fres-%05d_%05d", $rad, $AU;
$efile =~ s/\.hits$//;
my $cmd = "eventplot -i $ip $efile $index $fresfile";

system ($cmd);

exit 0;
