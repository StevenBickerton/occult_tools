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

# dump the fit
my $results = `rmin_bmax_from_stats $statfile -d junk`;

my $other = join(" ", @other);
system("plot junk [0:] [0:] 3d u1:2:'(\$4>-0.1?\$4:1/0)':'(2.0)' w p pt 7 ps variable palette set:view=map sett:xlabel=\'radius (m)\' sett:ylabel=\'Im. Par. (m)\' unset:key notitle $other");

exit 0;
