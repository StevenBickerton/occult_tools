#!/usr/bin/env perl
#
# original filename: rip_stat.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Thu Jun  5, 2008  17:02:18 DST
#  Host: bender.astro.princeton.edu
#  Working Directory: /Users/bick/usr/bin/cdiffracSim
#


use strict;
use warnings;
use File::Basename;

use Local::Constants;
use Local::Astrotools;

my $exe = basename($0);
my $usage = "Usage: $exe statfile other-args\n";

my ($statfile, @other) = @ARGV;
die $usage unless $statfile;

my $other = join(" ", @other);
system("plot $statfile [:] [-500:] 3d u2:4:5 w p pt 7 ps 2 palette set:view=map sett:xlabel=\'radius (m)\' sett:ylabel=\'Im. Par. (m)\' $other");

exit 0;
