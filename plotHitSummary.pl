#!/usr/bin/env perl
#
# original filename: plotHitSummary.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Tue Feb 19, 2008  17:32:47 EST
#  Host: bender.astro.princeton.edu
#  Working Directory: /Users/bick/sandbox/cdiffracSim
#


use strict;
use warnings;
use File::Basename;

use Local::Constants;
use Local::Astrotools;

my $exe = basename($0);
my $usage = "Usage: $exe hitSummaryfile fresdir xcor_thresh [outform(png|eps)]\n";

my ($hitfile, $fresdir, $xcor_thresh, $outform) = @ARGV;
die $usage unless $xcor_thresh;

$outform = "png" unless $outform;
die "Only eps and png accepted as outforms.\n" unless $outform=~/^(png|eps)$/;

my $ip_tmpfile = "$ENV{HOME}/tmp/plotHitSumm.tmp.dat";


sub loadFresHead($) {
  my ($fresfile) = @_;
  my %head;
  open(FRESFILE, "$fresfile") or die "opening $fresfile\n";
  while(<FRESFILE>) {
    my @line = split;
    if ( $_=~/\#/ ) {
      $head{$line[1]} = $line[3];
    } else {
      last;
    }
  }
  close(FRESFILE);
  return %head;
}

sub resetHeadIP($$$) {
  my ($fresfile, $label, $value) = @_;
  my $tmpfile = "$ENV{HOME}/tmp/plotHitSumm.tmp2";
  open(FRESFILE, "$fresfile") or die "opening $fresfile\n";
  open(TMPFILE, ">$tmpfile") or die "opening $tmpfile\n";
  while(<FRESFILE>) {
    if ($_=~/^\#/) {
      my @line = split;
      if ( $line[1] eq $label) {
	print TMPFILE "# $label  $line[2] $value\n";
      } else {
	print TMPFILE "$_";
      }
    } else {
      print TMPFILE $_;
    }
  }
  close(TMPFILE);
  close(FRESFILE);
  rename $tmpfile, $fresfile;
}

# load the hitfile
my @hits;
open(HITFILE, "$hitfile") or die "Opening $hitfile\n";
while(<HITFILE>) {
  next if /^\#/;
  push @hits, [split];
}
close(HITFILE);


# loop over each hit
foreach my $hit (@hits) {

  my ($ID, $hits, $index, $nk, $rank, $aa, $AU, $ip, $xcor, $cor0, $nxc, $xchi, $chi0, $rchi2, $rchi2n, $dof, $bary, $fft) = @$hit;

  next if $xcor < $xcor_thresh;

  (my $tsfile = $hits) =~ s/.hits$//;

  # -- make the fresfile name
  my $aa_s = sprintf("%05d", $aa);
  my $AU_s = sprintf("%05d", $AU);
  my $fresfile = "${fresdir}/fres-${aa_s}_${AU_s}";

  my %head = loadFresHead($fresfile);
  my $fsu = sqrt( 0.5*($head{lambLo}+$head{lambHi})*$AU*$AU_M/2.0 );
  my $outfile;
  my $IDs = sprintf("%03d",$ID);

  if ($tsfile=~/.*\.fits$/) {
    ($outfile = $tsfile) =~ s/fits$/${IDs}.${outform}/;
  } else {
    $outfile = $tsfile . ".${IDs}.${outform}";
  }

  # -- offset the fresfile
  system("offsetPattern $fresfile $ip 1 | awk \'\$1~/#/{print} \$1!~/#/{print \$2,\$1,\$3}\' > $ip_tmpfile");
  resetHeadIP($ip_tmpfile, "offset", $ip);

  # -- make the plot
  my $elong = 180.0;
  my $incl = 0.0;
  my $v = `elong2v $elong $incl $AU`;
  #my $hwid = (20*$fsu+2.0*$head{aa}+2.0*$head{RStar})/($v*$dt); # a dummy value
  my $hwid;
  my $kwid;
  if ($head{lambLo} > 1.0e-7) {
      $hwid = 50;
      $kwid = 0.25*$fsu/$v;
  } else {
      $hwid = 200;
      $kwid = 0.25*$fsu/$v;
  }
  my $wid = 2.0*$hwid;
  my $rank_s = sprintf("%.2g",$rank);
  my $cmd = "eventplot -i $ip -t \"Nk=$nk  $xcor/$cor0  $xchi/$chi0  (Xv/o=$rchi2/$rchi2n, Bary=$bary)\" -w $wid -k $kwid $tsfile $index $fresfile $outfile";
  system($cmd);
  #print $cmd."\n";
}

exit 0;
