#!/usr/bin/env perl
#
# original filename: sigmaTest.pl
#
# Steven Bickerton
#  Dept. of Astrophysical Sciences, Princeton University
#  bick@astro.princeton.edu
#  Created: Wed Jan 13, 2010  18:06:00 EST
#  Host: hammer5.astro.Princeton.EDU
#  Working Directory: /u/bick/usr/src/cdiffracSim64
#


use strict;
use warnings;
use File::Basename;
use Getopt::Std;

use Local::Constants;
use Local::Astrotools;
use Local::Statistics;

my %opt;
getopts("ehf:", \%opt);

my $exe = basename($0);
my $usage = "Usage: $exe fresfile counts\n";

my ($fresfile, $count) = @ARGV;
die $usage unless $count;
die $usage if $opt{'h'};

my $hz = $opt{'f'} ? $opt{'f'} : 40.0;

#inputs: fresnelfile count-level hz
#inputs: aa AU count-level hz

my $tmp_dir = "$ENV{'HOME'}/tmp";
my $fresDir = $tmp_dir."/fresnelFiles";
my $tsfile = $tmp_dir."/ts.fits";
my $tsAddFile = $tsfile.".add.fits";
my $tsNormFile = $tsfile.".add.norm.fits";
my $nAdd = 20;

my $elong = 180;
my $incl = 0;
my $offset = 0;
my $center = 0;
my $random = 1;

# load the fresnel file params
my %fresPars;
open(FRES, $fresfile);
while(<FRES>) {
    next unless /^\#/;
    my ($key, $expl, $val) = $_ =~ /^\#\s+([\w\d]+)\s+([^\s]+):\s+([^\s]+)/;
    $fresPars{$key} = $val;
}
close(FRES);
my $AU = $fresPars{'AU'};
my $lambda = 0.5*($fresPars{'lambLo'} + $fresPars{'lambHi'});
my $fsu = sqrt($AU*$AU_M*$lambda/2.0);

# make a time series
my $v = `elong2v $elong $incl $AU`;
my $eventWidth = 5.0*$fsu/$v; # in seconds

my $P = $nAdd*$eventWidth*50;
my $N = $P*$hz;
my $cmdTS = "ran_poisson -f $count $N $P > $tsfile";
print "$cmdTS\n" if $opt{'e'};
system($cmdTS);

# add the fresfile n times
my $cmdAdd = "addKBO $elong $incl $fresfile $tsfile $nAdd $offset $center $random";
print "$cmdAdd\n" if $opt{'e'};
system($cmdAdd);

# get a list of the added objects
my $nLines = 60 + $nAdd;
my $cmdList = "fold $tsAddFile | head -$nLines | awk '\$1~/HIT/{print \$4}'";
my @indices = `$cmdList`;

# get the sigma for the n entries
(-d $fresDir) or mkdir $fresDir;
unlink("$fresDir/fres-*");
system("cp $fresfile $fresDir/");
my ($cycles, $NperCycle, $corrTh, $chiTh, $Noff, $norm) = (2, 20, 2.0, 1.0, 0, 0);
my $parFile = "$tmp_dir/params.detect";
open(PAR, ">$parFile");
print PAR "$cycles $NperCycle $corrTh $chiTh $fresDir $Noff $norm\n";
close(PAR);

my $normWidth = 10.0*$eventWidth;
my $cmdNorm = "fftSmooth -b $tsAddFile 1:2 $normWidth";
print $cmdNorm if $opt{'e'};
system($cmdNorm);

my $cmdDet = "detect $elong $incl $parFile $tsNormFile 2> /dev/null";
print $cmdDet if $opt{'e'};
system($cmdDet);

# get the hits
my @dat =();
my $hitFile = $tsNormFile.".hits";
open(HITS, "$hitFile");
my $expect;
while(<HITS>) {
    next if /\#/;
    my @line = split;
    my ($rstar, $aa, $AU, $ip, $cor_i, $cormg, $cormg0, $co_n, $chi_i, $chimg, $chimg0, $ch_n, $rX2, $PX2, $rXnull, $dof, $baryc) = @line;
    for my $i (@indices) {
	push @dat, $cormg if ($cor_i >= $i-1 and $cor_i <= $i+1);
    }
    $expect = $cormg0;
}
close(HITS);
$expect = 0.0 unless $expect;

my $mean = mean(\@dat);
my $rms = rms(\@dat);
my $med = median(\@dat);

printf STDOUT "%.2f %.2f (exp=%.2f) (med=%.2f)\n", $mean, $rms, $expect, $med;

exit 0;
