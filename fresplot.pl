#!/usr/bin/env perl
#
#

use strict;
use warnings;
use Math::Trig;
use File::Basename;

my ($fresfile, @args) = @ARGV;
($fresfile and -s $fresfile) or die "usage: $0 fresfile [otherFresfiles] [elong] [outfile] [dump]\n";

my $PI = 3.141592654;
my $AU_M = 149600000000;
my $number = '[+-]?\d+\.?\d*';
my $elong = 180.0;
my $outfile;
my ($eps, $png) = (0,0);
my @fresfile = ("$fresfile");
my $dumpcode = 0;
my $xrange = 0;
foreach my $arg (@args) {
  if ($arg=~/^$number$/) {
    $elong   = $arg;
  } elsif ($arg=~/^.*\.eps$/ ) {
    $eps = 1;
    $outfile = $arg;
  } elsif ($arg=~/^.*\.png$/ ) {
    $png = 1;
    $outfile = $arg;
  } elsif ($arg=~/dump/) {
    $dumpcode = 1;
  } elsif ($arg=~/x=${number}/) {
      ($xrange) = $arg =~ /x=($number)/;
  } else {
    push @fresfile, $arg;
  }
}

# get the yrange of the data
my ($Imin,$Imax) = (1.0,1.0);
my ($xmin,$xmax) = (0.0,0.0);
my $has_negatives = 0;
my %head;
my $k = 1;
foreach my $fresfile2 (@fresfile) {
    open (DATA, "$fresfile2");
    while (<DATA>) {
    
    # load only the data for the 1st fresfile
    if ($_ =~ /^\#/ && $k == 1) {
        my ($param, $definition, $value) = (split)[1,2,3];
        $value *= 1e9 if ($param =~ /^lamb/); # convert to nm
        $head{$param} = $value;
        #printf "$param $definition $value\n";
        next;
    }
    next if /^\#/;

    # get the max's and min's of _all_ fresfiles
    my ($x, $Ihole, $Idisk) = split;
    $has_negatives = 1 if $x < 0.0;

    if ($Idisk<$Imin) {
        $Imin = $Idisk;
        $xmin = $x;
    }
    if ($Idisk>$Imax) {
        $Imax = $Idisk;
        $xmax = $x;
    }
    }
    close (DATA);
    $k++;
}

# a subroutine to convert elongation (radians) to retrograde velocity
sub elong2v ($$) {

    my ($elong, $AU) = @_;  # elong in radians

    my $G = 6.67259e-11;
    my $Me = 5.974e24;
    my $Mo = 1.989e30;

    my $Ve = sqrt($G*$Mo/$AU_M);

    # earth's velocity projected on the line of sight
    my $alpha = $PI - $elong;
    my $V_e_perp = $Ve * cos($alpha);

    # proposed KBO velocity projected on the line of sight
    my $Vk = $Ve / sqrt($AU);
    my $beta = asin ( sin($alpha)/$AU );
    my $V_k_perp = $Vk * cos($beta);

    # the velocity measured as seen from earth (ie. probably retrograd)
    my $vRet = ($V_k_perp - $V_e_perp);   #/100 to convert to metres
        
    return $vRet;  # in metres
}


my ($aa, $AU) = ($head{"aa"}, $head{"AU"} );
(my $filename = $fresfile) =~ s/_/\\_/g;


my $vRet = abs( elong2v($PI*$elong/180.0,$AU) );
my $vRet_s = sprintf "%.1f", $vRet/1000;

my $angle = sprintf "%.0f", $head{angle} || 0.0;

my $lambda = 1e-9*($head{lambHi} + $head{lambLo})/2.0 ;
my $fresnelscale = sqrt($lambda*$AU*$AU_M/2.0);
my $fsu = sprintf "%d", $fresnelscale;
$xrange = 8.3*$fresnelscale + 2.0*$head{aa} + 2.0*$head{RStar} unless $xrange > 0;
my $trange = $xrange/$vRet;

# set range 30% above/below the peaks and troughs
my $IrangeMin = 1.0 - 1.8*(1.0 - $Imin);
my $IrangeMax = 1.0 + 1.5*(1.0 - $Imin);

my $geoShadowMin = ($IrangeMin>0) ? $IrangeMin : 0;

# put a width arrow 10% above the peak
my $IpkArrow = 1.0 + 1.1*($Imax - 1.0);
my $IpkLabel = $IpkArrow + 0.2*($IrangeMax - 1.0);
my $IpkLabel2 = $IpkArrow + 0.1*($IrangeMax - 1.0);
my $width = 2*$xmax;
my $twidth = sprintf "%.3f", $width/$vRet;

# put a height arrow and label it
my $xhArrow = 2.0*$fresnelscale + $head{aa} + $head{RStar}; #$xmax + 0.1*($xrange - $xmax);
#my $xhLabel = $xmax + 0.5*($xrange - $xmax);
my $xhLabel = 2.5*$fresnelscale + $head{aa} + $head{RStar}; #$xhArrow + 500.0;
my $peak = sprintf "+%6.3f", 100 * ($Imax - 1.0);
my $depth = sprintf "-%6.3f", 100 * (1.0 - $Imin);


my $xh = -0.95*$xrange;
my $dyh = 0.11*(1.0-$IrangeMin);
my $yh = $IrangeMin + 0.2*(1.0 - $IrangeMin);

my $IgArr = 1.0 - 0.90*(1.0 - $IrangeMin);
my $IgLab = 1.0 - 0.80*(1.0 - $IrangeMin);
my $IgLab2 = 1.0 - 0.90*(1.0 - $IrangeMin);
my $diam = 2.0*$head{aa};
my $tdiam = sprintf "%.3f", $diam/$vRet;


# get the modification time of the file
my @month = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my $mtime = (stat($fresfile))[9];
my ($ss,$mm,$hh,$mdy,$M,$Y,$wdy,$ydy,$isdst) = localtime($mtime);
$Y += 1900;
#$Y = $Y % 100;
my $infile_mtime = sprintf "$Y $month[$M] $mdy %02d:%02d:%02d", $hh,$mm,$ss;



my $terminal = "x11 enhanced font \"Courier,14\"";
my $setoutput = "\n";
my $basename = basename($filename);
my $title = "file: $basename  ($infile_mtime)";
my $labelSize = 14;
if (defined($outfile) ) {
  if ($eps) {
    $terminal = "postscript eps enhanced color \"Courier\" 14";
    $labelSize = 16;
  }
  if ($png) {
    $terminal = "png";
    $labelSize = 12;
  }
  $setoutput = "set output \'$outfile\'";
  $title = "Diffraction Shadow for a $head{aa}m KBO at $head{AU} AU";
}
my $labelSize2 = $labelSize + 4;
my $titleSize = $labelSize + 4;


my $x2label = "Time (sec) (Assumes ${elong}{/Symbol \260} Solar Elongation)";
$x2label = "Time (sec) (Assumes ${elong}d Solar Elongation)" if $png;
my $lamblabel = "{/Symbol=${labelSize} l}";
$lamblabel = "lambda" if $png;

my $gpt_string = "#!/usr/bin/env gnuplot\n".
    "#\n".
    "set term $terminal\n".
    "$setoutput\n".

    "set style line 2 lt 2 lw 3 lc rgb \"forest-green\"\n".
    "set style line 3 lt 3 lw 3 lc rgb \"blue\"\n".

    "set tmargin 6\n".
    "set bmargin 4\n".
    "set rmargin 9\n".
    "set lmargin 9\n".
    "set title \'$title\' font \"Courier,${titleSize}\"\n".
    "set xlabel \'Distance from Shadow Centre (m)\' font \"Courier,${labelSize2}\"\n".
    "set x2label \'$x2label\' font \"Courier,${labelSize2}\"\n".
    "set ylabel \'I/<I>\' font \"Courier,${labelSize2}\"\n".
    "set xrange [-$xrange:$xrange]\n".
    "set x2range [-$trange:$trange]\n".
    "set xtic nomirror\n".
    "set x2tic nomirror\n".
    "set yrange [$IrangeMin:$IrangeMax]\n".
    "set y2range [$IrangeMin:$IrangeMax]\n".
    "set y2tic\n".

    # fresnel scale
    "yfsarrow=1.0+0.8*($IrangeMax-1.0)\n".
    "yfslabel=1.0+0.65*($IrangeMax-1.0)\n".
    "yfslabel2=1.0+0.55*($IrangeMax-1.0)\n".
    "xfsarrow=0.8*$xrange\n".
    "xfslabel=xfsarrow-$fresnelscale/2.0\n".
    "set arrow 10 from xfsarrow,yfsarrow to xfsarrow-$fresnelscale,yfsarrow lw 2 heads\n".
    "set label 10 \"1 Fsu\" at xfslabel,yfslabel center\n".
    "set label 11 \"($fsu m)\" at xfslabel,yfslabel2 center\n".

    # geometric shadow lines
    "set arrow 1 from -$xrange,1 to -$aa,1 ls 2 nohead\n".
    "set arrow 2 from -$aa,1 to -$aa,$geoShadowMin ls 2 nohead\n".
    "set arrow 3 from $aa,$geoShadowMin to $aa,1 ls 2 nohead\n".
    "set arrow 4 from $aa,1 to $xrange,1 ls 2 nohead\n".
    ($geoShadowMin ? "" : "set arrow 5 from -$aa,$geoShadowMin to $aa,$geoShadowMin ls 2 nohead\n").

    # Star size lines
    "set arrow 11 from -$head{RStar},1 to $head{RStar},1 ls 3 nohead\n".
    "set arrow 12 from -$head{RStar},1+0.2*$depth/100.0 to -$head{RStar},1-0.2*$depth/100.0 ls 3 nohead\n".
    "set arrow 13 from $head{RStar},1+0.2*$depth/100.0 to $head{RStar},1-0.2*$depth/100.0 ls 3 nohead\n".


    # peak width arrow and label
    "set arrow 100 from -$xmax,$IpkArrow to $xmax,$IpkArrow lt 1 lw 2 heads\n".
    "set label 100 \'${width}m\' at 0,$IpkLabel font \"Courier,${labelSize}\" center\n".
    "set label 101 \'(${twidth}s)\' at 0,$IpkLabel2 font \"Courier,${labelSize}\" center\n".

    # peak height arrows and labels
    "set arrow 200 from -$xhArrow,$Imax to -$xhArrow,$Imin lt 1 lw 2 heads\n".
    "set label 202 \'$peak %%\' at -$xhLabel,$Imax font \"Courier,${labelSize}\" right\n".
    "set label 203 \'$depth %%\' at -$xhLabel,($Imin+1.0)/2.0 font \"Courier,${labelSize}\" right\n".

    # geometric width and labels
    "lgarrow=0.2\n".
    "set arrow 300 from -$aa-lgarrow*($xrange-$aa),$IgArr to -$aa,$IgArr lt 1 lw 2 head\n".
    "set arrow 301 from $aa+lgarrow*($xrange-$aa),$IgArr to $aa,$IgArr lt 1 lw 2 head\n".
    "set label 300 \'Geometric Shadow\' at -$aa-lgarrow*($xrange-$aa),$IgLab font \"Courier,${labelSize}\" right\n".
    "set label 301 \'${diam}m (${tdiam}s)\' at -$aa-lgarrow*($xrange-$aa),$IgLab2 font \"Courier,${labelSize}\" right\n".

    # R-side header labels
    "set label 401 \'${lamblabel} = $head{lambLo}-$head{lambHi} nm ($head{nLambda})\' at 0.92*$xrange,$yh+0.0*$dyh font \"Courier,${labelSize}\" right\n".
    "set label 402 \'v_{retro} = $vRet_s m/s\' at 0.92*$xrange,$yh+1.0*$dyh font \"Courier,${labelSize}}\" right\n".
    "set label 403 \'R*= $head{RStar} m ($head{nStars})\' at 0.92*$xrange,$yh+2.0*$dyh font \"Courier,${labelSize}\" right tc ls 3\n".
    "set label 404 \'R_{KBO}= $head{aa} m\' at 0.92*$xrange,$yh+(3.0)*$dyh font \"Courier,${labelSize}}\" right tc ls 2\n".
    "set label 405 \'D = $head{AU} AU\' at 0.92*$xrange,$yh+(4.0)*$dyh font \"Courier,${labelSize}\" right\n".
    "set label 406 \'offset = $head{offset} m\' at 0.92*$xrange,$yh+(5.0)*$dyh font \"Courier,${labelSize}\" right\n".
    "set label 407 \'AR/ang = $head{aspectRatio}/$angle d\' at 0.92*$xrange,$yh+(6.0)*$dyh font \"Courier,${labelSize}\" right\n".

    "set grid\n".
    "plot ";


my $bits = 255;
my $colorIncr = $bits / (scalar @fresfile);
$k = 0;
foreach my $fresfile2 (@fresfile) {
    my $lc = sprintf "lc rgb \"#%x0000\"", ($bits - $k*$colorIncr);
    if ($has_negatives) {
      $gpt_string .= "\'$fresfile2\' using 1:3 t \"\" w lp lt 1 pt 1 $lc\\\n";
    } else {
      $gpt_string .= "\'$fresfile2\' using 1:3 t \"\" w lp lt 1 pt 1 $lc, \'$fresfile2\' using (-\$1):3 t \"\" w lp lt 1 pt 1 $lc,\\\n";
    }
    $k++;
}
$gpt_string =~ s/,\\$/\n/;

open (GNUPLOT, "|gnuplot -persist");
printf GNUPLOT "$gpt_string";
printf STDOUT "$gpt_string" if $dumpcode;
close (GNUPLOT);

#printf "$fresnelscale\n";

exit 0;
