#!/usr/bin/env perl
#
#

use strict;
use warnings;

my $usage = "$0 infile epsfile [(dump)arrows]\n";

my ($infile) = shift @ARGV;
die $usage unless $infile and -s $infile;

my ($epsfile, $dumparrows) = ("", 0);
foreach my $arg (@ARGV) {
    if ($arg =~ /.*\.eps$/) {
	$epsfile = $arg;
    }
    if ($arg =~ /dump/) {
	$dumparrows = 1;
    }
}

my $count = 1;
sub box ($$$$) {

    my ($x, $y, $dx, $dy) = @_;

    my $left = $x-$dx;
    my $right = $x+$dx;
    my $top = $y+$dy;
    my $bot = $y-$dy;

    my $str = sprintf "set arrow ${count}0 from %.1f,%.1f to %.1f,%.1f nohead ls 1\n", $left, $bot, $left, $top;
    $str .= sprintf "set arrow ${count}1 from %.1f,%.1f to %.1f,%.1f nohead ls 1\n", $left, $top, $right, $top;
    $str .= sprintf "set arrow ${count}2 from %.1f,%.1f to %.1f,%.1f nohead ls 1\n", $right, $top, $right, $bot;
    $str .= sprintf "set arrow ${count}3 from %.1f,%.1f to %.1f,%.1f nohead ls 1\n", $right, $bot, $left, $bot;

    return $str;
}


system ("echo \"0,0\" > ./zero");

my $xyextra = 1.3;
my $xyarrow = 1.15;
my $xylabel = 1.05;
my ($xmax, $ymax) = (0.0,0.0);
my $boxes;
open (INFILE, "$infile");
while (<INFILE>) {

    my ($x, $y, $dx, $dy) = split /\s+/, $_;

    $boxes .= sprintf "%s", box($x,$y,$dx,$dy);
    
    $count++;

    $xmax = ($x+$dx) if (($x+$dx)>$xmax);
    $ymax = ($y+$dy) if (($y+$dy)>$ymax);

}
close (INFILE);

if ($dumparrows) {
    print $boxes;
    exit 0;
}

my $max = ($xmax>$ymax) ? $xmax : $ymax;
$max *= $xyextra;

my $aspRat = sprintf "%.2f", $xmax/$ymax;
$xmax = sprintf "%d",$xmax;
$ymax = sprintf "%d",$ymax;


my $terminal = "x11 enhanced font \"Helvetica,20\"\n";
my $setoutput = "";
my $size = 1.0;
my $lw = 1.0;
if ($epsfile) {
    $terminal = "postscript eps enhanced color \"Helvetica\" 48";
    $setoutput = "set output \'$epsfile\'\n";
    $size = 2.0;
    $lw = 4.0;
}

$count--;
my $gp = "#!/usr/bin/env gnuplot\n".
    "set terminal $terminal\n".
    "$setoutput\n".
    "set size ratio 1 $size,$size\n".
    "set xlabel \'KBO x-dimension (m)\'\n".
    "set ylabel \'KBO y-dimension (m)\'\n".
    "set label 1 \"$count boxes\" at graph 0.95,0.93 right\n".
    "set label 2 \"Aspect Ratio = $aspRat\" at graph 0.07,0.93 left\n".
    "lw2 = $lw\n".

    "set style line 1 lt 1 lw $lw lc rgb \"black\"\n".
    
    # vertical arrow
    "set arrow 1 from -$xyarrow*$xmax,0 to -$xyarrow*$xmax,$ymax lt 1 lw $lw lc \"black\"\n".
    "set label 3 \"$ymax m\" at -$xylabel*$xmax,0.95*$ymax left\n".

    # horiz arrow
    "set arrow 2 from 0,-$xyarrow*$ymax to $xmax,-$xyarrow*$ymax lt 1 lw $lw lc \"black\"\n".
    "set label 4 \"$xmax m\" at 0.75*$xmax,-$xylabel*$ymax cen\n".
    "set border lw lw2\n".
    "$boxes\n".
    "plot [-$max:$max] [-$max:1.1*$max] \'zero\' t \"\" w dot\n";



open (GNUPLOT, "| gnuplot -persist");
printf GNUPLOT "$gp\n";
close (GNUPLOT);


unlink ("./zero");
