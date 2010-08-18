#!/usr/bin/env perl
#
# original filename: animFresPatt.pl
#
# Steven Bickerton
#  Dept. of Physics/Astronomy, McMaster University
#  bick@physics.mcmaster.ca
#  Made with makeScript, Sat Mar 18, 2006  10:39:15 EST
#  Host: kuiper
#  Working Directory: /1/home/bickersj/sandbox/cdiffracSim
#

use strict;
use warnings;
use File::Basename;

use lib "$ENV{HOME}/libperl";
use Local::Constants;

my $exe = basename($0);
my $usage = "Usage: $exe outfile.mpg fresnelFiles ...\n";

my ($outfile, @fresfiles) = @ARGV;
die $usage if $ARGV[0] and $ARGV[0] =~ /help/;
die $usage unless ($outfile and $outfile =~ /.*\.mpg/);
die $usage unless $fresfiles[0] and -s $fresfiles[0];


# read in the last file and get the dimmensions
my %head;
open (FILE1, "$fresfiles[-1]") or die "Unable to open last fresnelfile\n";
while (<FILE1>) {
    my @line = split;

    if ( /^\#/ ) {
	$line[3] *= 1e9 if ($line[1] =~ /^lamb/);
	$head{$line[1]} = $line[3];
	next;
    }
}
close (FILE1);

#die "maxX00: $head{maxX00}\n";
 

# make the plot files
my @tmpbases;
my $tmpDir = "$ENV{HOME}/tmp";
my $count = 1;
my $N = scalar @fresfiles;
foreach my $fresfile (@fresfiles) {

    my $fresbase = basename $fresfile;
    my ($rad, $AU, $RStar) = map {sprintf "%d", $_ } $fresbase =~ /fres-(\d\d\d\d\d)_(\d\d\d\d\d)_(\d\d\d\d\d)/;

    my $scale = 1.0;
    my $xl = $scale*0.75;
    my $yl = $scale*0.85;
    my $dxl = $scale*0.0;
    my $dyl = $scale*(-0.07);
    my $lw = 3.0;

    my $tmpbase = sprintf "$tmpDir/animFresPatt%05d",$count;
    my $tmpfile = "${tmpbase}.eps";
    push @tmpbases, $tmpbase;
    my $gp = "#!/usr/bin/env gnuplot\n".
	"set term postscript eps enhanced color \"Helvetica\" 22\n".
	"set output \'$tmpfile\'\n".
	"set size 1,1\n".
	"set xrange [-$head{maxX00}/1000:$head{maxX00}/1000]\n".
	"set yrange [0:1.5]\n".
	"set xlabel \'Diffraction Shadow Dimension (km)\' 0,0\n".
	"set ylabel \'I/<I>\' 1,0\n".
	#"set xtics rotate by -60\n".
	"set border lw $lw\n".
	#"set tmargin 2\n".
	#"set bmargin 4\n".
	#"set lmargin 7\n".
	#"set rmargin 17\n".
	"set multiplot\n".
	"set size 0.75,1.0\n".
	"set origin 0,0\n".
	"set grid\n".
	"\n".
	
	# add some info about the patterns
	"set label 1 \"R* = $RStar m\" at screen $xl+0.0*$dxl,$yl+0.0*$dyl left\n".
	"set label 2 \"D = $AU AU\" at screen $xl+1.0*$dxl,$yl+1.0*$dyl  left\n".
	"set label 3 \"a_{KBO} = $rad m\" at screen $xl+2.0*$dxl,$yl+2.0*$dyl left\n".
	"set label 4 \"{/Symbol l} = $head{lambLo} - $head{lambHi} nm\" at screen $xl+3.0*$dxl,$yl+3.0*$dyl left\n".

	# plot
	"plot \'$fresfile\' u (\$1/1000):3 t \"\" w l lt 1 lw $lw lc \"black\", \'\' u (-\$1/1000):3 t \"\" w l lt 1 lw $lw lc \"black\"\n";
    
    open (GP, "|gnuplot") or die "Unable to open Gnuplot";
    printf GP "$gp\n";
    close (GP);

    system ("convert $tmpfile ${tmpbase}.mpg");

    printf STDERR "plot $count of %d done\n", $N;

    $count++;
}

# build the animation
unlink "$outfile";
foreach my $tmpbase (@tmpbases) {
    system ("cat ${tmpbase}.mpg >> $outfile");
}

# clean up
foreach my $tmpbase (@tmpbases) { unlink "${tmpbase}.eps", "${tmpbase}.mpg"; }

exit 0;
