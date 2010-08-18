#!/usr/bin/env perl
#
# original filename: getRateFromStats.pl
#
# Steven Bickerton
#  Dept. of Astrophysical Sciences, Princeton University
#  bick@astro.princeton.edu
#  Created: Mon Oct  6, 2008  16:32:07 DST
#  Host: bender.astro.Princeton.EDU
#  Working Directory: /Users/bick/usr/src/cdiffracSim
#


use strict;
use warnings;
use File::Basename;

use Getopt::Std;

use Local::KBOdensity;
use Local::Constants;
use Local::Astrotools;

my $exe = basename($0);
my $usage = "Usage: $exe [options] statfile q_S\n".
    "  -h This help message\n".
    "  -k r_knee [25km]\n".
    "  -L q_L (slope for large objects) [4.6]\n".
    "  -S q_S (slope for small objects) [2.0]\n".
    "  -w lambda (wavelength - used to get fresnel scale) [5.5e-7m]\n".
    "\n";

my %options = ();
getopts("hL:S:k:w:", \%options);
die $usage if $options{h};

my $rknee = $options{k} ? $options{k} : 25000.0;
my $qL = $options{L} ? $options{L} : 4.6;
my $qS = $options{S} ? $options{S} : 2.0;
my $lambda = $options{w} ? $options{w} : 5.5e-7;

my ($statfile) = @ARGV;
die $usage unless $statfile;


# ===========================================================
# read in the stats file
my %stat;
my $AU = -1;
open(STATFILE, "$statfile") or die "Unable to open $statfile";
while(<STATFILE>) {
    next if /^\#/;
    my ($Rstar, $r_kbo, $AUtmp, $ip, $nrec, $nadd, $nhit, $vret) = split;
    $stat{$r_kbo} = [] unless $stat{$r_kbo};

    if ( $AU < 0 ) {
	$AU = $AUtmp;
    } else {
	die "Multiple AUs defined in this statfile.  Exiting.\n" 
	    unless $AUtmp = $AU;
    }

    $nrec = $nadd if $nrec > $nadd;

    push @{$stat{$r_kbo}}, [$ip, $nrec/(1.0*$nadd), $vret];
}
close(STATFILE);


# ===========================================================
# get an equivalent impact parameter and vret;
my (%ip_equiv, %vret, %ip_max);

foreach my $r_kbo (sort {$a<=>$b} keys %stat) {
    my $ip_last = 0.0;
    $ip_equiv{$r_kbo} = 0.0;
    $ip_max{$r_kbo} = 0.0;
    $vret{$r_kbo} =  ${ $stat{$r_kbo} }[0]->[-1]; # last val in 1st (0) ip entry
    foreach my $ref ( @{$stat{$r_kbo}} ) {
	my ($ip, $rec, $vret) = @$ref;
	$ip_equiv{$r_kbo} += ($ip - $ip_last) * $rec;
	$ip_max{$r_kbo} = $ip if $ip > $ip_max{$r_kbo};
	$ip_last = $ip;
    }
}


# ===========================================================
# get the recovery in the radius+d_radius slice
my ($mu, $n_ms) = (0.0, 0.0);
my @r_kbo = sort {$a<=>$b} keys %stat;
my $nSas_prev = KBOdensity($r_kbo[0], $rknee, $qL, $qS);

if ( $ip_equiv{$r_kbo[0]} > 0 ) {
    die "non-zero impact parameter in first r-bin.  Include smaller bins.\n";
}
if ( $ip_equiv{$r_kbo[0]} > 0.95 * $ip_max{$r_kbo[0]} ) {
    die "ip_equiv ~= ip_max! ... Include larger impact parameters.\n";
}

my $fsu = sqrt($lambda * $AU * $AU_M / 2.0);
my $r_min = 0.0;

foreach my $r_kbo ( @r_kbo ) {
    
    next if ($r_kbo == $r_kbo[0]); # already have the first one

    my $nSas  = KBOdensity($r_kbo, $rknee, $qL, $qS);
    my $dnSas = $nSas_prev - $nSas;
    my $dn_ms = $dnSas * ($ARCSECperRAD**2) / ($AU * $AU_M)**2;

    if ( ($ip_equiv{$r_kbo} > 0.5 * $fsu) &&  # quite arbitrary
	 ( $r_min < 1.0e-3 ) ) {
	$r_min = $r_kbo;
    }
    $mu += 2.0 * $ip_equiv{$r_kbo} * $vret{$r_kbo} * $dn_ms;

    $nSas_prev = $nSas;

}

printf STDOUT "%.1f %.4g\n", $r_min, $mu;

exit 0;
