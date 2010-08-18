#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Getopt::Std;

my %opt;
getopts("i:k:l:t:w:y:", \%opt);
my ($infile, $index, $fresfile, $outfile) = @ARGV;

die "usage: $0 infile index fresfile [outfile]\n".
    "\t\t set fresfile = x if not showing fresfile\n".
    "Options: \n".
    "         -i impact param\n".
    "         -l label\n".
    "         -t title\n".
    "         -w width\n".
    "         -y ymin:ymax\n".
    "\n" unless $fresfile;
die "File $infile not found.\n" unless -s $infile;

my $lab = ($opt{l}) ? "$opt{l}" : "";
my $title = ($opt{t}) ? "$opt{t}" : "";
foreach my $arg (@ARGV) {
  if ( (my $tmp) = $arg=~/^label=(.*)$/) {
    $lab = $tmp;
  }
  if ( (my $tmp) = $arg=~/^title=(.*)$/) {
    $title = $tmp;
  }
}
my $ip = ($opt{i}) ? "$opt{i}" : "0";
my $hwid_default = 100;
my $hwid = ($opt{w}) ? "$opt{w}" : $hwid_default;
my ($ymin, $ymax);
if ($opt{y}) {
    my $delim = ":";
    ($ymin, $ymax) = split /$delim/, $opt{y};
}


my $tmpdir = "$ENV{HOME}/tmp";
my $tmpfile1 = "$tmpdir/eventplot.tmp1.dat";
my $tmpfile2 = "$tmpdir/eventplot.tmp2.dat";

my ($elong, $incl) = (180.0, 0.0);


# dump the data to a file to get stats
system("fitsTSdump -c $infile $index $hwid > $tmpfile1");
my $TSstats = `cstats $tmpfile1`;
my @TSstats = split /\n/, $TSstats;

my $t0 = (split /\s+/, $TSstats[10])[1];
my $tf = (split /\s+/, $TSstats[12])[1];
my $dt = ($tf - $t0) / ($hwid - 1.0);


my $rms = (split /\s+/, $TSstats[9])[2];

# sort out the output files and relevant variables
my $eps = ($outfile and $outfile=~/.*.eps$/) ? 1 : 0;
my $png = ($outfile and $outfile=~/.*.png$/) ? 1 : 0;    

my $date = `date`;
chomp $date;
my $term = "x11 enhanced font 24";
my $output = "";
my $lw = 1;
if ($eps) {
  $term = "postscript eps enhanced color font 20";
  $output = "set out \'$outfile\'";
  $lw = 1;
}
if ($png) {
  $term = "png large";
  $output = "set out \'$outfile\'";
  $lw = 1;
}


# deal with the fresnel file, if present
my $showfres = ( -s $fresfile ) ? 1 : 0;

my $fres_dy = 3.0*$rms;
my $fresplot = "";
my $labels = "";
my ($tmin, $tmax) = ("","");
my ($Imin, $Imax) = (1.0, 1.0);
my $shadeStarRad = "";
my $fsu_arrow = "";
my $v;
my %head;
my $hwid2 = $hwid;
my $kwid = 2.0*$dt;
if ($showfres) {

  # get the velocity
  open(FRESFILE, "$fresfile");
  while(<FRESFILE>) {
    my @line = split;
    next unless @line;
    if ($_=~/^\#/) {
      $head{$line[1]} = $line[3];
    } else {
      $Imax = $line[2] if $line[2] > $Imax;
      $Imin = $line[2] if $line[2] < $Imin;
    }
  }
  close(FRESFILE);


  $v = `elong2v $elong $incl $head{AU}`;
  chomp $v;
  my $overdump = 1.05;
  my $fsu = sqrt( 0.5*($head{lambLo}+$head{lambHi})*$head{AU}*1.5e11/2 );

  $kwid = ($opt{k}) ? "$opt{k}" : $fsu / $v;

  if ($hwid==$hwid_default) {
      $tmax = 1.5 * (16.3*$fsu+2.0*$head{aa}+2.0*$head{RStar})/$v;
      $tmin = -$tmax;
  } else {
      $tmax = $hwid * $dt / 2.0;
      $tmin = -$hwid * $dt / 2.0;
  }
  $hwid2 = $overdump * ($tmax-$tmin) / $dt;

  my $xrange = $tmax;
  $fres_dy = 3.0*$rms + (1.0-$Imin);

  my $dt_10 = $dt / 10.0;
  my $kernfile = basename($fresfile);
  $kernfile =~ s/fres/kern/;
  $kernfile = $tmpdir . '/' . $kernfile; 
  system("offsetPattern $fresfile $ip 1 > $kernfile");
  $fresplot = ",\\\n".
      "\'$kernfile\' u (\$2/$v):(\$3+$fres_dy) w l ls 3 t \"Theory\",\\\n".
      "\'\' u (-\$2/$v):(\$3+$fres_dy) w l ls 3 t \"\"\\\n";

  # shaded regions for star and Rkbo
  my $starlayer = "behind";
  my $radlayer = "back";
  if ($head{aa} > $head{RStar}) {
    $starlayer="back";
    $radlayer = "behind";
  }
  my $shadeStar = 
    "set obj 1 rect from -$head{RStar}/$v,1+1.1*$fres_dy to $head{RStar}/$v,$Imax+$fres_dy\n".
      "set obj 1 $starlayer fillstyle solid nobord fc rgb \"#ff5500\"\n";
  my $shadeRad =
    "set obj 2 rect from -$head{aa}/$v,$Imin+$fres_dy to $head{aa}/$v,1+0.9*$fres_dy\n".
      "set obj 2 $radlayer fillstyle solid nobord fc rgb \"#00cc00\"\n";
  $shadeStarRad = "${shadeStar}${shadeRad}";

  # Fsu dimmension arrow
  $fsu_arrow = "set arrow 1 from -0.5*$fsu/$v, 1+$fres_dy to 0.5*$fsu/$v, 1+$fres_dy heads\n";

  # R-side header labels
  my $dyh = 0.07;
  my $yh = 0.15;
  my $labelSize = ($eps) ? 16 : 14;
  my $vRet_s = sprintf "%.1f", $v/1000.0;
  my $angle = ($head{angle}) ? $head{angle} : 0.0;
  my $lambsym = ($png) ? "lambda" : "{/Symbol=${labelSize} l}";
  my $vsym = ($png) ? "v" : "v_{retro}";
  my $Rksym = ($png) ? "Rkbo" : "R_{KBO}";
  my $Rlabels =  "set label 402 \'$vsym = $vRet_s m/s\' at gr 0.96,$yh+0.0*$dyh right\n".
      "set label 401 \'${lambsym} = $head{lambLo}-$head{lambHi} m ($head{nLambda})\' at gr 0.96,$yh+1.0*$dyh right\n".
      "set label 404 \'R*= $head{RStar} m ($head{nStars})\' at gr 0.96,$yh+2.0*$dyh right tc rgb \"#ff3300\"\n";
  
  my $Llabels = "set label 405 \'D = $head{AU} AU\' at gr 0.04,$yh+(0.0)*$dyh left\n".
      "set label 406 \'i.p. = $ip m\' at gr 0.04,$yh+(1.0)*$dyh left\n".
      "set label 403 \'$Rksym= $head{aa} m\' at gr 0.04,$yh+(2.0)*$dyh left tc rgb \"#009900\"\n";


  $labels = $Llabels . $Rlabels;
    

} else {
  $fresplot = "\n";
}


# shaded negatives
my $shadeNegative = 
    "set obj 5 rect from gr 0, 0 to gr 1, first 0\n".
    "set obj 5 behind fillstyle solid nobord fc rgb \"#dddddd\"\n";


system("fitsTSdump -c $infile $index $hwid2 > $tmpfile1");
system("fftSmooth $tmpfile1 1:2 $kwid > $tmpfile2");

my $file_str = $infile;
$title = "$file_str at index $index (${date})" unless length($title) > 1;

($ymin, $ymax) = ( -6.0*$rms, 1.0+$fres_dy+3.0*($Imax-1.0) ) unless ( $opt{y});


#die $title;
#  write the plot script
my $kwid_s = sprintf("%.5g",$kwid);
my $gpscript = "#!/usr/bin/env gnuplot\n".
  "set term $term\n".
  "$output\n".
  "set grid\n".
  "set key bottom center horizontal\n".
  "set title \"$title\" offset 0,0 noenhanced\n".
  "set xlabel \"Time (s)\"\n".
  "set x2label \"Size (km)\"\n".
  "set ylabel \"Intensity\"\n".
  "$labels\n".
  "$shadeStarRad\n".
  "$shadeNegative\n".  
  "$fsu_arrow\n".
  "set label 100 \'$lab\' at gr 0.5,0.95 center\n".
  "set xrange [$tmin:$tmax]\n".
  "set x2range [$tmin*$v:$tmax*$v]\n".
  "set yrange [$ymin:$ymax]\n".
  "set xtic nomirror\n".
  "set x2tic nomirror\n".
  "set ytic 0.2\n".
  "set style line 1 lt 1 lw $lw lc rgb \"red\"\n".
  "set style line 2 lt 1 lw $lw lc rgb \"black\"\n".
  "set style line 3 lt 1 lw $lw lc rgb \"blue\"\n".
  "plot \'$tmpfile2\' u 1:4 w histep ls 1 t \"Raw\",\\\n".
  "\'\' u 1:2 w histep ls 2 t \"$kwid_s s smooth\"\\\n". #,\\\n".
    #"0 w line ls 2 t \"\"\\\n".
  #"\'\' u 1:(\$4-\$2) w histep ls 1 t \"Smooth Resid\"".
  $fresplot;


open(GP, "|gnuplot -persist");
print GP $gpscript;
close (GP);

#print $gpscript;

exit 0;
