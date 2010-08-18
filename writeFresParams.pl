#!/usr/bin/env perl
#
# original filename: writeFresParams.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Mon Sep 17, 2007  11:23:10 DST
#  Host: bender.astro.princeton.edu
#  Working Directory: /Users/bick/sandbox/cdiffracSim
#


use strict;
use warnings;
use File::Basename;

use Local::Constants;
use Local::Stardata;

my $exe = basename($0);
my $usage = "Usage: $exe lamblo lambhi V-mag startype AU [off=x] [ar=x] [order=x]\n";

my ($lamblo, $lambhi, $V, $mk, $AU, @args) = @ARGV;
die $usage unless $AU;

my $number = '[+-]?\d+\.?\d*';
my ($offset, $ar, $order) = (0.0, 1.0, 5);
foreach my $arg (@args) {
  
  if ($arg=~/off=$number/) {
    $offset = $arg =~ /off=($number)/;
  }
  if ($arg=~/ar=$number/) {
    $ar = $arg =~ /ar=($number)/;
  }
  if ($arg=~/order=\d+/) {
    $order = $arg =~ /order=(\d+)/;
  }
}

sub round ($$) {
  my ($value, $unit) = @_;

  my $intpart = int($value/$unit);
  my $fracpart = $value/$unit - $intpart;
  $intpart += 1 if $fracpart > 0.5;

  return $intpart*$unit;
}

my $Nlamb = ($lambhi - $lamblo) * 1.0e8;
$Nlamb = 1 if $Nlamb == 0;

my $lamb = ($lamblo+$lambhi)/2.0;
my $fs = sqrt($lamb*$AU*$AU_M/2.0);

my $x00Step = round($fs/20.0, 1);
my $maxX00 = 200.0*$x00Step;

my $starRo = $star{$mk}{'R'} * $Ro_M;
my $Mv = $star{$mk}{'Mv'};
my $star_dist = 10**((($V-$Mv) + 5.0)/5.0) * $PC_M;
my $star_angle = $starRo / $star_dist;
my $RStarProj = int($star_angle * $AU * $AU_M);

my $star_area = $PI * $RStarProj**2 / $fs**2;
my $Nstar = int(100.0 * $star_area);
$Nstar = 1 if $Nstar == 0;

printf STDOUT "# lambLo lambHi Nlamb maxX00 x00step Rstarproj Nptsrc AU offs AR order\n".
  "$lamblo $lambhi $Nlamb $maxX00 $x00Step $RStarProj $Nstar $AU $offset $ar $order\n";

exit 0;
