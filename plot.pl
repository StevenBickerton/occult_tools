#!/usr/bin/env perl
#
#  plot.pl - quick and dirty plotter using gnuplot.
#            can pipe data from command to get a plot quickly.
#  Steve Bickerton, Sun Sept 25, 2005;
#

use strict;
use warnings;

my $number = '[+-]?\d+\.?\d*[eE]?[+-]?\d*';

# see if there's an input file
( $ARGV[0])                    and 
    ( $ARGV[0] !~ /.*\.eps$/ ) and 
    ( -s $ARGV[0] )            or 
    unshift(@ARGV, "-");

my $infile = shift;

my @string = @ARGV;

my $epsfile = "";
my $pngfile = "";
if (@string) {
    $epsfile = pop @string if ($string[-1] =~ /.*\.eps$/ );
}
if (@string) {
    $pngfile = pop @string if ($string[-1] =~ /.*\.png$/ );
}

my $gnuterm = ($ENV{'GNUTERM'}) ? $ENV{'GNUTERM'} : "x11";
my $term = "set term $gnuterm enhanced\n";
my $output = "\n";
if ($epsfile) {
    $term = "set terminal postscript eps enhanced color \"Helvetica\" 20\n";
    $output = "set output \'$epsfile\'\n";
}
if ($pngfile) {
    $term = "set terminal png enhanced\n";
    $output = "set output \'$pngfile\'\n";
}

my $string = join (" ", @string);

my $usage = "\nusage: $0 [infile] [gnuplot plot string] [file.eps]\n".
    "\n".
    " uX:Y         = sets x and y data columns to X and Y\n".
    " cCOL:LO:HI   = sets conditionals column:low:high (limits inclusive)\n".
    " eCOL:LO:HI   = sets exclusion column:low:high (limits inclusive)\n".
    " hl           = sets header labels ... axis labels taken from header\n".
    " set:par=val  = sets a gnuplot parameter to value\n".
    " sett:par=val = sets a text parameter to \'value\' (note the quotes)\n".
    " lcN          = sets a line colour column: N\n".
    " logx         = sets x axis to logscale (can use x,y,xy)\n".
    " 3d           = 3d plot ... need 3 columns\n".
    " labelN:foo[:vh]  = sets label #N \"foo\" at corner vh\n".
    "                    vh is ul,ur(dflt),ll,lr\n".
    "                    -must quote (\'\') or escape (\\) non-letters\n".
    " vline:x      = puts a vertical line at X\n".
    " hline:y      = puts a horizontal line at Y\n".
    " vshade:x1:x2 = vertically shade between x1 and x2\n".
    " hshade:y1:y2 = horizontally shade between y1 and y2\n".
    "\n".
    "eg.   plot file.dat \"[0:100] u1:2 w p\"\n".
    "      cat file.dat | plot\n".
    "      cat file.dat | plot [:] [10:] u1:3 w l logy\n".
    "      cat file.dat | plot u1:4 out.eps\n".
    "      plot file.dat u1:3 sett:xlabel=\"Seconds\"\n".
    "      cat file.dat | plot u1:2 u1:3 u1:4 w l logxy\n".
    "      cat file.dat | plot u1:2 lc3\n".
    "      cat file.dat | plot u1:2 u1:3 u1:4 lc5 lc6 lc7\n".
    "      cat file.dat | plot 3d u1:2:3\n".
    "      cat file.dat | plot u1:3 hl\n".
    "      cat file.dat | plot u1:3 c2:0:5\n".
    "      plot file.dat u1:2 vline:0 hline:0\n".
    "      plot file.dat u1:2 vshade:0:10 hl logy\n".
    "      plot file.dat u1:3 label1:a=3.142 label2:b='sin(x)'\n".
    "\n";

die $usage if $string=~/help/;



# -------------------------------------------------------------------------
# variables which one might potentially want to change
my $nbins = 8;

my $tmpdir = "$ENV{HOME}/tmp/";
( -d $tmpdir ) or $tmpdir = "";

my $tmpfile = $tmpdir."plottmp";          # this file _WILL_ be deleted
my $scriptfile = $tmpdir."plotscript.gp"; # this file _WON'T_ be delected


# make a header for the script file ... in case you want to keep it.
my $date = `date`;
chomp $date;
my $host = ($ENV{HOST}) ? "$ENV{HOST}": 
    ( ($ENV{HOSTNAME}) ? "$ENV{HOSTNAME}" : "unknown" ) ;
my $cwd = ( $ENV{PWD} ) ? "$ENV{PWD}" : `pwd`;
$cwd = "unknown" if ( ! $cwd );

my $header = "#!/usr/bin/env gnuplot\n".
    "#\n".
    "# Steven Bickerton\n".
    "# Dept. of Physics/Astronomy, McMaster University\n".
    "# bick\@physics.mcmaster.ca\n".
    "# Made with plot.pl $date\n".
    "# Host: $host\n".
    "# Working Directory: $cwd\n".
    "#\n";



# -------------------------------------------------------------------------
# get some relevant info from the arguments
my @range;
$string = "";
my $set = "";
my $splot = 0;
my @plots;
my @columns;
my @lc_col;
my @conditions;
my @exclusions;
my $script_warning = "";
my $set_labels = 0;

# user labels
my $set_user_label = "";
my %nlabel = ( "ul" => 1, "ll" => 1, "ur" => 1, "lr" => 1);

# user lines
my $vlines = "";
my $hlines = "";
my $nvline = 1;
my $nhline = 100;

# shading
my $nvshade = 1;
my $vshade = "";
my $nhshade = 100;
my $hshade = "";
my $nrshade = 200;
my $rshade = "";

# harmonics
my $pharm = "";
my $npharm = 200;

# dump
my $dump_stdout = 0;

my $no_title = 0;

sub scriptWarn ($) {my ($msg) = @_; chomp $msg; 
                    $script_warning .= "# WARN: $msg\n";}

foreach my $part (@string) {

    # snip out the range information so we can put it in before the filename
    # in the printf to gnuplot
    if ($part =~ /^\[.*:.*\]$/ ) {
        push @range, "$part ";

        # make it possible to set parameters
    } elsif ($part =~ /^set:/ ) {
        my ($param, $value) = $part =~ /^set:(\S+)=(\S*)$/;
        die "Error setting parameter: $part\n\n$usage\n" unless defined $value;
        $set .= "set $param $value\n";
        
        # make it possible to unset parameters
    } elsif ($part =~ /^unset:/ ) {
        my ($param) = $part =~ /^unset:(\w+)$/;
        die "Error unsetting parameter: $part\n\n$usage\n" unless $param;
        $set .= "unset $param\n";
        
        # make a shortform for logs
    } elsif ($part =~ /^log/ ) {
        my ($value) = $part =~ /^log(\w\w?\w?)$/;
        die "Error setting log scale: $part\n\n$usage\n" 
            unless $value =~ /(x|y|z|xy|xz|yz|xyz)/;
        $set .= "set logscale $value\n";
        
        # make it possible to set text parameters
    } elsif ($part =~ /^sett:/ ) {
        my ($param, $value) = $part =~ /^sett:(\S+)=(.*)$/;
        die "Error setting text parameter: $part\n\n$usage\n" unless $param;
        $set .= "set $param \'$value\'\n";
        
        # set splot
    } elsif ($part =~ /^3d$/ ) {
        $splot = 1;
        
        # set a line colour column
    } elsif ($part =~ /^lc/ ) {
        my ($lc_col) = $part =~ /^lc(\d+)$/;
        die "Error in lc: $part\n\n$usage\n" unless $lc_col;
        push @lc_col, $lc_col;
        
        # set header labels
    } elsif ($part =~ /^(headlab|hl)$/ ) {
        $set_labels = 1;
        
        # set conditionals
    } elsif ($part =~ /^c\d+/ ) {
        my ($col,$lo,$hi) = $part =~ /^c(\d+):(\S+):(\S+)$/;
        die "Error in conditional: $part\n\n$usage\n" unless $col;
        scriptWarn("Conditional $part is not handled by this script\n");
        push @conditions, [$col,$lo,$hi];
        
        # set exclusions
    } elsif ($part =~ /^e\d+/ ) {
        my ($col,$lo,$hi) = $part =~ /^e(\d+):(\S+):(\S+)$/;
        die "Error in exclusion: $part\n\n$usage\n" unless $col;
        scriptWarn("Exclusion $part is not handled by this script\n");
        push @exclusions, [$col,$lo,$hi];
        
        # extract the 'using' statements
    } elsif ($part =~ /^u\S+:/) {
        (my $part2 = $part) =~ s/u/using /;
        push @columns, $part2;
        
        # extract user labels
    } elsif ($part =~ /^label\d+:/) {
        my ($num,$label) = $part =~ /label(\d+):(.+)$/;
        die "Error in label: $part\n\n$usage\n" unless $label;
        my ($ver,$hor) = $label =~ /\w+:([ul])([lr])$/;
        $label =~ s/:[ul][rl]$//;
        
        $ver = 'u' unless $ver;
        $hor = 'r' unless $hor;
        my $nlab = $nlabel{$ver.$hor}++;
        my $v_incr = 0.05;
        my ($xpos,$xjust) = ($hor eq 'l') ? (0.05,"l") : (0.95,"r");
        my ($ypos) = ($ver eq 'u') ? (1.0 - $nlab*$v_incr) : ($nlab*$v_incr);
        $set_user_label .= "set label $num \"$label\" at graph $xpos,$ypos $xjust\n";
        
        # extract vertical lines
    } elsif ($part =~ /^vline:\S+/) {
        my ($xpos) = $part =~ /^vline:(\S+)$/;
        die "Error in vline: $part\n\n$usage\n" 
            unless ($xpos=~/^$number$/  || $xpos=~/^${number}e${number}$/ );
        
        $vlines .= "set arrow $nvline from $xpos,graph 0 to $xpos,graph 1.0 nohead\n";
        $nvline++;
        
        # extract horizontal lines
    } elsif ($part =~ /^hline:\S+/) {
        my ($ypos) = $part =~ /^hline:(\S+)$/;
        die "Error in hline: $part\n\n$usage\n"
            unless ($ypos=~/^$number$/  || $ypos=~/^${number}e${number}$/ );
        
        $vlines .= "set arrow $nhline from graph 0,first $ypos to graph 1.0,first $ypos nohead\n";
        $nhline++;
        
        
        # extract shading 
    } elsif ($part =~ /^vshade:/) {
        my ($xlo,$xhi) = $part =~ /^vshade:(\S+):(\S+)$/;
        die "Error in vshade: $part\n\n$usage\n" 
            unless ( ($xlo=~/^$number$/  || $xlo=~/^${number}e${number}$/ ) &&
                     ($xhi=~/^$number$/  || $xhi=~/^${number}e${number}$/ ) );
        
        $vshade .= "set obj $nvshade rect from $xlo,gr 0 to $xhi,gr 1.0\n".
            "set obj $nvshade behind ".
            "fillstyle solid nobord fc rgb \"#bbbbbb\"\n";
        $nvshade++;
        
    } elsif ($part =~ /^hshade:/) {
        my ($ylo,$yhi) = $part =~ /^hshade:(\S+):(\S+)$/;
        die "Error in hshade: $part\n\n$usage\n"
            unless ( ($ylo=~/^$number$/  || $ylo=~/^${number}e${number}$/ ) &&
                     ($yhi=~/^$number$/  || $yhi=~/^${number}e${number}$/ ) );
        
        $hshade .= "set obj $nhshade rect ".
            "from gr 0,fir $ylo to gr 1.0,fir $yhi\n".
            "set obj $nhshade behind fillstyle solid noborder ".
            "fc rgb \"#bbbbbb\"\n";
        $nhshade++;
        
    } elsif ($part =~ /^rshade:/) {
        my ($xlo,$ylo,$xhi,$yhi) = $part =~ /^rshade:(\S+):(\S+):(\S+):(\S+)$/;
        die "Error in rshade: $part\n\n$usage\n"
            unless ( ($xlo=~/^$number$/  || $xlo=~/^${number}e${number}$/ ) &&
                     ($xhi=~/^$number$/  || $xhi=~/^${number}e${number}$/ ) &&
                     ($ylo=~/^$number$/  || $ylo=~/^${number}e${number}$/ ) &&
                     ($yhi=~/^$number$/  || $yhi=~/^${number}e${number}$/ ) );
        
        $rshade .= "set obj $nrshade rect ".
            "from $xlo,$ylo to $xhi,$yhi\n".
            "set obj $nrshade behind fillstyle solid noborder ".
            "fc rgb \"#bbbbbb\"\n";
        $nrshade++;
        
        
        # show period harmonics
    } elsif ($part =~ /^pharm:/) {
        my ($period, $nharm) = $part =~ /^pharm:($number):?(\d+)?$/;
        die "Error in pharm: $part\n\n$usage\n"
            unless ($period);
        
        foreach my $harmonic ( map {$period/$_} (1 .. $nharm) ) {
            
            my $harmonicS = sprintf "%.2g", $harmonic;
            $pharm .= "set arrow $npharm from $harmonic, gr 1.00 to $harmonic, gr 0.95\n";
            $pharm .= "set label $npharm \"$harmonicS\" at $harmonic, gr 1.02 left rotate by 90\n";
            $npharm++;
        }
        
        # show frequency harmonics
    } elsif ($part =~ /^fharm:/) {
        my ($period, $nharm) = $part =~ /^fharm:($number):?(\d+)?$/;
        die "Error in fharm: $part\n\n$usage\n"
            unless ($period);
        
        foreach my $harmonic ( map {$period*$_} (1 .. $nharm) ) {
            
            my $harmonicS = sprintf "%.2g", $harmonic;
            $pharm .= "set arrow $npharm from $harmonic, gr 1.00 to $harmonic, gr 0.95\n";
            $pharm .= "set label $npharm \"$harmonicS\" at $harmonic, gr 1.02 left rotate by 90\n";
            $npharm++;
        }
        
        # set variable to dump gpstring to STDOUT
    } elsif ($part =~ /^dump$/ ) {
        $dump_stdout = 1;

    } elsif ($part =~ /^notitle$/) {
        $no_title = 1;
        
        # what's left is the plot string
    } else {
        $string .= "$part ";
    }
}

if (! $columns[0]) {
    $columns[0] = ($splot) ? "using 1:2:3" : "using 1:2";
}

foreach my $j (0 .. $#columns) {
    $lc_col[$j] = $lc_col[0] if ($columns[$j] && ! $lc_col[$j]);
}
my $lc_col = scalar @lc_col;



# -------------------------------------------------------------------------
# get the column info for the first plot to make axis labels
my ($xcol, $ycol, $zcol) = $columns[0] =~ /^using (\S+):(\S+)?:?(\S+)?$/;
$xcol = "1" unless $xcol;
$ycol = "2" unless $ycol;
$zcol = "3" unless $zcol;

my $xcol_int = ($xcol =~ /^\d+$/) ? $xcol : 0;
my $ycol_int = ($ycol =~ /^\d+$/) ? $ycol : 0;
my $zcol_int = ($zcol =~ /^\d+$/) ? $zcol : 0;


# -------------------------------------------------------------------------
# we'll put the data in a tmp file ... this is necessary if we're reading
#   the data from STDIN ... or could use '-' as plotfile instead.
my @lc_vals;
my @lc_bin = (0)     x $lc_col;
my @lc_min = (1e99)  x $lc_col;
my @lc_max = (-1e99) x $lc_col;
my ($xlabel, $ylabel, $zlabel);
# dump the data to a tmp file
open (INFILE, "$infile");
open (TMP, ">$tmpfile");
READIN: while (<INFILE>) { 

    # grab column labels if they exist (last for comments split values)
    if ($_ =~ /^\#/) {
        next READIN unless $set_labels;
        my @line = split;
        my $xlab = "\'$line[$xcol]\'" if $xcol_int && $line[$xcol_int];
        my $ylab = "\'$line[$ycol]\'" if $ycol_int && $line[$ycol_int];
        my $zlab = "\'$line[$zcol]\'" if $zcol_int && $line[$zcol_int];
        ($xlabel,$ylabel) = ($xlab,$ylab) if ($xlab && $ylab);
        $zlabel = "set zlabel ".$zlab if ($xlab && $ylab && $zlab);
        next READIN;
    }
    
    # keep the blank lines as they mark data blocks
    if ($_ =~ /^\s+$/) { printf TMP "$_"; next READIN;}

    my @line = split; 

    # collect info on the values in the line colour column
    foreach my $j (0 .. $#lc_col) {

        if ($lc_col[$j]) {
            
            # if the number of different values is less than nbins
            #   add a new value to the list
            if ( scalar keys %{$lc_vals[$j]} <= $nbins ) {
                ${$lc_vals[$j]}{$line[$lc_col[$j]-1]} = 1;
                
                # otherwise ... set the flag to bin data in this column
            } else { 
                $lc_bin[$j] = 1; 
            }
            
            my $lc = $line[$lc_col[$j]-1];  # the value of this data point
            $lc_min[$j] = $lc if ($lc<$lc_min[$j]);
            $lc_max[$j] = $lc if ($lc>$lc_max[$j]);
        }
        
    }
    
    
    # filter data according to the conditional
    foreach my $condition (@conditions) {
        my ($col, $lo, $hi) = @$condition;
        next READIN 
            unless $line[$col-1] >= $lo and $line[$col-1] <= $hi;
    }
    
    # filter data according to the exclusions
    foreach my $exclusion (@exclusions) {
        my ($col, $lo, $hi) = @$exclusion;
        next READIN 
            unless $line[$col-1] <= $lo or $line[$col-1] >= $hi;
    }
    
    printf TMP "$_"; 
}
close (TMP);
close (INFILE);


# set default labels if they couldn't be found in the header.
$xlabel = "\'X-column(s): $xcol\'" unless ($xlabel and $ylabel);
$ylabel = "\'Y-column(s): $ycol\'" unless ($xlabel and $ylabel);
$zlabel = "set zlabel \'Z-column(s): $zcol\'" 
    unless ($xlabel and $ylabel and $zlabel);



# -------------------------------------------------------------------------
# make it possible to splot
my $plot = $splot ? "splot" : "plot";





# -------------------------------------------------------------------------
# build the set range statements
my $range = "";
my @axis = ("x", "y", "z");
foreach my $q (0 ..2) {
    if ($range[$q]) {
        my @values = $range[$q] =~ /^\[(.*):(.*)\] $/;
        
        my $ax = $axis[$q];
        my ($min, $max) = @values;
        
        my $minl = ($min =~ /^$number$/) ? "${ax}min" : "";
        my $maxl = ($max =~ /^$number$/) ? "${ax}max" : "";
        
        $range .= sprintf "%s${ax}min = $min\n", ($min =~ /^$number$/) ? 
            "" : "# ";
        $range .= sprintf "%s${ax}max = $max\n", ($max =~ /^$number$/) ? 
            "" : "# ";
        
        $range .= "set ${ax}range [$minl:$maxl]\n";
        
    }
}
foreach my $q (0 ..2) {
}




# -------------------------------------------------------------------------
# change line color for different columns

my $bits = 255;
my @hex = ("#%x0000", "#00%x00", "#0000%x");
my $mod = 3;

foreach my $j (0 .. $#columns) {
    
    #my $colorIncr = $bits / ( ( scalar keys %{$lc_vals[$j]} ) - 1);
    my $colorIncr = ($lc_vals[$j]) ? 
        $bits / ( scalar keys %{$lc_vals[$j]} ) : $bits;
    
    # grab the string for a single plot with these columns
    #   if lc_col is set, we'll overwrite this string with the first
    #      line in the colour sequence
    
    $columns[$j] =~ s/:$//;  # if it's just a single column, lose the :
    
    if ($lc_col[$j] && $lc_bin[$j]) {
        push @{$plots[$j]}, sprintf 
            ",\\\n\'$tmpfile\' $columns[$j] $string lc rgb \"$hex[$j % $mod]\"", $bits;
    } else {
        push @{$plots[$j]}, sprintf ",\\\n\'$tmpfile\' $columns[$j] $string";
    }
    
    # if we're binning the data for line colours
    if ($lc_col[$j] && $lc_bin[$j]) {
        
        my ($xcolumn, $ycolumn) = $columns[$j] =~ /^using (\S+):(\S+)$/;
        $xcolumn = "\$".$xcolumn if $xcolumn =~ /^\d+$/;
        my $lc_col = $lc_col[$j];
        my $range = $lc_max[$j] - $lc_min[$j];
        my $step = $range / ($nbins-1.0);
        
        # wipe out the original
        shift @{$plots[$j]};
        
        # loop over the number of bins and push the plot string onto the array
        foreach my $q (0 .. $nbins-1) {
            my $lo = sprintf "%.3g", $q*$step;
            my $hi = sprintf "%.3g", ($q+1)*$step;
            my $rgb = sprintf "\"$hex[$j % $mod]\"", ($bits - $q*$colorIncr);
            my $col_string = "((\$$lc_col>=$lo && \$$lc_col<$hi) ? $xcolumn : (1/0)):$ycolumn";
            push @{$plots[$j]}, ",\\\n\'\' using $col_string $string t \"$lo < \$$lc_col < $hi  $xcolumn:$ycolumn\" lc rgb $rgb";
        }
        
        
        # if we're not binning the data
        #   ie. there are only a few values and we're just giving each value
        #   it's own colour ... colour does not depend on the value itself.
    } elsif ($lc_col[$j] && !$lc_bin[$j]) {
        
        my ($xcolumn, $ycolumn) = $columns[$j] =~ /^using (\S+):(\S+)$/;
        $xcolumn = "\$".$xcolumn if $xcolumn =~ /^\d+$/;
        my $lc_col = $lc_col[$j];  
        my $k=0;
        
        # wipe out the original
        shift @{$plots[$j]};
        
        # loop over each different value and push the plot string to the array
        foreach my $lc_val ( sort {$a<=>$b} keys %{$lc_vals[$j]} ) {
            
            my $rgb = sprintf "\"$hex[$j % $mod]\"", ($bits - $k*$colorIncr);
            my $col_string = "((\$$lc_col==$lc_val) ? $xcolumn :(1/0)):$ycolumn";
            push @{$plots[$j]}, ",\\\n\'\' using $col_string $string t \"\$$lc_col=$lc_val $xcolumn:$ycolumn\" lc rgb $rgb";
            $k++;
        }
    }
}



# -------------------------------------------------------------------------
# write the plot script
my @allplots;
my $nplot = 0;
foreach my $plotref (@plots) { 
    $nplot += scalar @$plotref; 
    if ($splot) {
        foreach my $plot (@$plotref) {
            $plot =~ s/\s+lc.*//;
        }
    }
    push @allplots, join("",@$plotref);     
}

# replace the leading comma and blank filename
for ($allplots[0]) {
    s/^,//;
    s/^\\\n\'\'/\\\n\'$tmpfile\'/;
}

# unset the key if it's only a single plot
my $setkey = (scalar keys %{$lc_vals[0]} or $columns[1]) ? "" : "unset key\n";
$setkey = "set key below\n" if ( $nplot >= 6 && $nplot <=12 );
$setkey = "unset key\n" if ( $nplot > 12);

# handle underscores in the filename to print it in the title
(my $unenhanced_infile = $infile) =~ s/_/\\_/g;

# get the date/time and add it to the title
my @month = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my $mtime = ($infile eq '-') ? time() : (stat($infile))[9];
my ($ss,$mm,$hh,$mdy,$M,$Y,$wdy,$ydy,$isdst) = localtime($mtime);
$Y += 1900;
#$Y = $Y % 100;
my $infile_mtime = sprintf "$Y $month[$M] $mdy %02d:%02d:%02d", $hh,$mm,$ss;

my $title = ($infile eq '-') ? 
    "<STDIN>  ($infile_mtime)" : "$unenhanced_infile   ($infile_mtime)";

# get host and user to put in the title
my $title_info = "Host: $host User: $ENV{USER}";
my $title_line = "";
if ($no_title > 0) {
    $title_line = "";
    $title = "";
} else {
    $title_line = "set title sprintf(\"%s\\n\\n%s\\n\",\'$title\', \"$title_info\")";
}


my $gp_string = "$header".
    "$script_warning".
    "#\n".
    "$term".
    "$output".
    "\n".
    "set border lw 2\n".
    "set grid lw 2\n".
    "\n".
    "$title_line\n".
    "set xlabel $xlabel\n".
    "set ylabel $ylabel\n".
    "$zlabel\n".
    "\n".
    "set style line 1 lt 1 lw 2 lc rgb \"black\"\n".
    "set style line 2 lt 2 lw 2 lc rgb \"black\"\n".
    "set style line 3 lt 3 lw 2 lc rgb \"black\"\n".
    "set style line 4 lt 4 lw 2 lc rgb \"black\"\n".
    "\n".
    "$set_user_label\n".
    "\n".
    "$vlines".
    "$hlines".
    "\n".
    "$vshade".
    "$hshade".
    "$rshade".
    "\n".
    "$pharm".
    "\n".
    "$setkey".
    "\n".
    "$set".
    "\n".
    "$range".
    "\n".
    "$plot".
    join("",@allplots);

$gp_string .= "\n";



# -------------------------------------------------------------------------
# pipe the plot command to gnuplot
open (GNUPLOT, "|gnuplot -persist") or 
    die "gnuplot pipe failed\n";
print GNUPLOT "$gp_string";
close (GNUPLOT);


$gp_string =~ s/$tmpfile/$infile/g if $infile!~/^-$/;

open (SCRIPT, ">$scriptfile") or 
    die "gnuplot open scriptfile: $scriptfile failed\n";
print SCRIPT "$gp_string";
close(SCRIPT);

print STDOUT "$gp_string" if $dump_stdout;

unlink ("$tmpfile");
