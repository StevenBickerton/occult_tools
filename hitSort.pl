#!/usr/bin/env perl
#
# original filename: hitSummary.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Mon Feb 11, 2008  10:21:15 EST
#  Host: bender.astro.princeton.edu
#  Working Directory: /Users/bick/sandbox/analysis
#


use strict;
use warnings;
use File::Basename;
use Local::Constants;
use Local::Astrotools;

my $exe = basename($0);
my $usage = "Usage: $exe eventfile\n";

my ($efile, $ilimit) = @ARGV;
die $usage unless defined($efile);
die "File not found\n" unless -s $efile;

$ilimit = 50 unless defined($ilimit);

my %events;

open (EVENTS, "$efile");
  
while (<EVENTS>) { 
  next if /^\#/; 
  
  my $t_limit = 1.0;
  my @record = split;
  my ($RStar, $aa, $AU, $offset, 
      $icor, $mgcor, $cor0, $corN, 
      $ichi, $mgchi, $chi0, $chiN, 
      $rchi2, $Prchi, $rchi2n, $dof, $baryc, $used_fft) = @record;
  
  $mgcor = 1e-9 if defined($mgcor) and $mgcor==0.0;
  $cor0 = 1e-9 if defined($cor0) and $cor0==0.0;
  $mgchi = 1e-9 if defined($mgchi) and $mgchi==0.0;
  $chi0 = 1e-9 if defined($chi0) and $chi0==0.0;
  my $cor_ratio = ($mgcor>$cor0) ? $cor0/$mgcor : $mgcor/$cor0;
  my $chi_ratio = ($mgchi>$chi0) ? $chi0/$mgchi : $mgchi/$chi0;
  #my $rank = $mgcor*$cor_ratio + abs($mgchi*$chi_ratio);
  my $rank = $Prchi;
  $rank = 0.0 if $Prchi =~ /^nan$/;

  #my $time = sprintf "%.3f", $icor*$t_step;
  push @record, ($rank);
  
  # if it's close to one that we already have, put it with that one
  my $t_key = $icor;
  foreach my $existingT (keys %events) {
    $t_key = $existingT if ( abs($icor-$existingT) < $ilimit );
  }
  
  push @{$events{$t_key}}, [@record];	
  
}
close (EVENTS);
die "No events in $efile ... exiting\n" unless (keys %events);

# --- print 'em
my $ebase = basename($efile);
my $elen = length($ebase) - 1;
printf STDOUT "# %2s %${elen}s  %7s %3s %9s ".
  "%4s %5s %4s   %5s %5s %3s  %6s %6s  ".
  "%6s %7s %3s  %5s %3s".
  "\n", "N", "efile", "index", "nk", "rank", 
  "rad", "AU", "ip",    "xcor", "cor0", "nxc", "xchi", "chi0",
  "rchi2", "rchi2n", "dof", "baryc", "fft";
my $Nhit = 1;
foreach my $t_key (sort {$b<=>$a} keys %events) {
  my @record = @{$events{$t_key}};
  my $nkern = sprintf("%3d",scalar @record);
  foreach my $hit ( sort {$b->[-1]<=>$a->[-1]} @{$events{$t_key}} ) {
    printf STDOUT "%-3d $ebase  %7d $nkern %5.3e ".
      "%4.0f %5.0f %4.0f   %5.2f %5.2f %3d  %6.2f %6.2f  ".
	"%6.2f %7.2f %3d  %5.2f %3d".
	"\n", $Nhit++, @$hit[4, -1, 1,2,3, 5,6,7, 9,10,12,14,15,16, 17];
  }
}

exit 0;
