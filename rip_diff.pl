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
my $usage = "Usage: $exe statfile snr_ref other-args\n";

my ($statfile, $snr_ref, @other) = @ARGV;
die $usage unless $snr_ref;

my $other = join(" ", @other);
system("plot $statfile [:] [-500:] 3d u2:4:'((\$10)>$snr_ref?(\$9-\$10):1/0)':'(2.0)' u2:4:'(\$10<$snr_ref?(\$9-\$10):1/0)':'(1.5)' w p pt 7 ps variable palette set:view=map sett:xlabel=\'radius (m)\' sett:ylabel=\'Im. Par. (m)\' set:cbrange=[-8:8] unset:key $other");

exit 0;
