#!/usr/bin/env perl
#
# original filename: plotStars.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Sat Mar 18, 2006  19:07:14 EST
#  Host: kuiper
#  Working Directory: /1/home/bickersj/sandbox/cdiffracSim
#

use strict;
use warnings;
use File::Basename;

use lib "$ENV{HOME}/libperl";
use Local::Constants;

my $exe = basename($0);
my $usage = "Usage: $exe starfile [epsfile]\n";

my ($starfile,$epsfile) = @ARGV;
die $usage unless $starfile;

my ($aa, $AU) = $starfile =~ /fres-(\d\d\d\d\d)_(\d\d\d\d\d).*/;
$AU *= 1.0;

my $RStar = 0;
my $n=0;
open (STAR, "$starfile");
while (<STAR>) {
    if (/^\#/) { $RStar = (split)[3]; next; }
    $n++;
}
close(STAR);


my $terminal = "x11 enhanced font \"Helvetica,20\"\n";
my $setoutput = "";
my $size = 1.0;
my $lw = 1.0;
my $ps = 1.0;
my $pt = 3;
if ($epsfile) {
    $terminal = "postscript eps enhanced color \"Helvetica\" 48";
    $setoutput = "set output \'$epsfile\'\n";
    $size = 2.0;
    $lw = 4.0;
    $ps = 2.0;
    $pt = 3;
}

my $xyextra = 1.3;

my $gp = "#!/usr/bin/env gnuplot\n".
    "set terminal $terminal\n".
    "$setoutput\n".
    "set size square $size,$size\n".
    "set xlabel \'Star x-dimension (m)\'\n".
    "set ylabel \'Star y-dimension (m)\'\n".
    "lw2 = $lw\n".
    "set border lw lw2\n".
    "set label 1 \'$n sources\' at graph 0.95,0.94 right\n".
    "set label 2 \"R* = $RStar m\" at graph 0.05,0.94 left\n".
    "set label 3 \"D = $AU AU\" at graph 0.5,0.94 cen\n".
    #"set polar\n".
    #"set rrange [0:$xyextra*$RStar]\n".
    "set parametric\n".
    "set xrange [-$xyextra*$RStar:$xyextra*$RStar]\n".
    "set yrange [-$xyextra*$RStar:$xyextra*$RStar]\n".
    "set samples 400\n".
    "plot \'$starfile\' u 1:2 t \"\" w p pt $pt ps $ps lc \"black\", $RStar*sin(t),$RStar*cos(t) t \"\" lt 1 lw $lw lc \"black\"\n";




open (GNUPLOT, "| gnuplot -persist");
printf GNUPLOT "$gp\n";
close (GNUPLOT);


exit 0;
