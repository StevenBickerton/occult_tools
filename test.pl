#!/usr/bin/env perl
#
#

use warnings;
use strict;

my ($program2test) = @ARGV;
my $all = 1 unless $program2test;

my %program;
foreach my $program2test (@ARGV) {   $program{$program2test} = 1;  }

my $buildDir = `pwd`;
chomp $buildDir;

my $log = "$buildDir/test_commands.log";
system ("rm -f $log") if ( -f $log);

sub mysystem($;$) {
    my ($cmd,$doLog) = @_;
    if (defined $doLog) {
        open (LOG, ">>$log");
        #printf LOG "echo \"$cmd\"\n";
        printf LOG "$cmd\n";
        close(LOG);
    }
    return system("$cmd");
}

(-d "test") or mkdir("test");
mysystem("cp ts test/");
mysystem("cp ts.fits test/");
chdir ("test/");

my $fresDir = "fresnelFiles";
my $fresDirT = "fresnelFilesT";
my $fresDirTA = "fresnelFilesTA";
my $fresPref = "fres-";
my $fresfile = "fres-01000_00320";
my $tsfile = "ts";
my $tsfits = "ts.fits";
my $tsfileTA = "tsTA";
my $ts_add_file = "ts.add";
my $ts_add_fits = "ts.fits.add.fits";
my ($command, $command1, $command2);
my $paramfile;
my $elong = 180.0;
my $incl = 0.0;
my $offset = 500;
my $dt = 0.025;
my $random = 1;
my $norandom = 0;
my $center = 0;
my $maxOffset = 10000;
my $Noffset = 10;
my $corrThresh = 5;
my $chiThresh = -2;



##################################################
#  test fresnelBox
#    make a pattern for: g' filter, 100m star, 40 AU
#
##################################################

my $program = "fresnelBox";
if ( $program{"$program"} or $all ) {
    
    printf STDERR "Testing $program ... ";
    
    
    mysystem ("rm -rf $fresDir") if ( -d $fresDir );
    
    mkdir("$fresDir");
    chdir("$fresDir");
    
    mysystem ("cp $buildDir/$program ./");
    
    $paramfile = "params.$program";
    open (PARAMS, ">${paramfile}1");
    printf PARAMS "# lambLo lambHi Nlambda   maxX00 x00Step    RStarProj Nptsrc  AU    offset  aspectRatio order\n";
    printf PARAMS "4.0e-7 5.5e-7 10    20000 50   100.0 25 40.0   0.0 1.0 5\n";
    close (PARAMS);

    open (PARAMS, ">${paramfile}2");
    printf PARAMS "# lambLo lambHi Nlambda   maxX00 x00Step    RStarProj Nptsrc  AU    offset  aspectRatio order\n";
    printf PARAMS "4.0e-7 5.5e-7 10    20000 50   100.0 25 320.0   0.0 1.0 5\n";
    close (PARAMS);
    
    # Usage: ./fresnelBox aa paramFile [dump]
    $command = "./$program 300 ${paramfile}1 1 2> error.${program};".
	"./$program 1000 ${paramfile}2 2>> error.${program};";
 
    my $out1 = "fres-00300_00040";
    my $out2 = "fres-01000_00320";
    my $out3 = "${out1}.stars";
    my $out4 = "${out1}.boxes";
    my $ran = (mysystem($command, 1)==0) ? 1:0;
    my $wrote = (-s $out1 and -s $out2 and -s $out3 and -s $out4);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

    chdir("../");
}




##################################################
#  test fresnelT
#    make a pattern for: g' filter, 100m star, 40 AU
#
##################################################

$program = "fresnelT";
if ( $program{"$program"} or $all ) {
    
    printf STDERR "Testing $program ... ";
    
    
    mysystem ("rm -rf $fresDirT") if ( -d $fresDirT );
    
    mkdir("$fresDirT/");
    chdir("$fresDirT");
    
    mysystem ("cp $buildDir/$program ./");
    
    $paramfile = "params.$program";
    open (PARAMS, ">$paramfile");
    printf PARAMS "# lambLo lambHi Nlambda   maxX00 x00Step    RStarProj Nptsrc  AU    offset  aspectRatio order\n";
    printf PARAMS "4.0e-7 5.5e-7 10    20000 50   100.0 25 40.0   0.0 1.0 5\n";
    close (PARAMS);
    
    # Usage: ./fresnelBox aa paramFile [dump]
    $command = "./$program 300 5770 $paramfile 1 2> error.${program};".
	"./$program 400 5770 $paramfile 2>> error.${program};";
 
    my $out1 = "fres-00300_00040";
    my $out2 = "fres-00400_00040";
    my $out3 = "${out1}.stars";
    my $out4 = "${out1}.boxes";
    my $ran = (mysystem($command, 1)==0) ? 1:0;
    my $wrote = (-s $out1 and -s $out2 and -s $out3 and -s $out4);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

    chdir("../");
}



##################################################
#  test fresnelTA
#    make a pattern for: g' filter, 100m star, 40 AU, TNO rotated 30deg, ip=1km
#
##################################################

$program = "fresnelTA";
if ( $program{"$program"} or $all ) {
    
    printf STDERR "Testing $program ... ";
    
    
    mysystem ("rm -rf $fresDirTA") if ( -d $fresDirTA );
    
    mkdir("$fresDirTA/");
    chdir("$fresDirTA");
    
    mysystem ("cp $buildDir/$program ./");
    
    $paramfile = "params.$program";
    open (PARAMS, ">$paramfile");
    printf PARAMS "# lambLo lambHi Nlambda   maxX00 x00Step    RStarProj Nptsrc  AU    offset  aspectRatio order angle\n";
    printf PARAMS "4.0e-7 5.5e-7 8    20000 50   100.0 16 40.0   1000.0 4.0 4 30.0\n";
    close (PARAMS);
    
    # Usage: ./fresnelBox aa paramFile [dump]
    $command = "./$program 300 5770 $paramfile 1 2> error.${program};".
      "./$program 400 5770 $paramfile 2>> error.${program};";
 
    my $out1 = "fresTA-00300_00040";
    my $out2 = "fresTA-00400_00040";
    my $out3 = "${out1}.stars";
    my $out4 = "${out1}.boxes";
    my $ran = (mysystem($command, 1)==0) ? 1:0;
    my $wrote = (-s $out1 and -s $out2 and -s $out3 and -s $out4);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";
    
    chdir("../");
}


######################################################
#  test addKBO
#   
#######################################################

$program = "addKBO";
if ($program{"$program"} or $all ) {
    printf "Testing $program ... ";
    # Usage: ./addKBO elong incl fresnelfile timeseries N offset center[1|0] random[1|0]
    $command = "../$program $elong $incl ${fresDir}/$fresfile $tsfile 10 $offset 0 1 2> error.${program}";
    
    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s $ts_add_file);
    printf "%s\n", ($ran and $wrote) ? "success" : "failure";
    
}

#########################################################
#  test detect
#
#########################################################

$program = "detect";
if ($program{"$program"} or $all ) {

    printf "Testing $program ... ";
    
    
    $paramfile = "params.${program}";
    open (PARAMS, ">$paramfile");
    printf PARAMS "# cycles NperCycle  corrThresh  chiThresh  fresDir      Noffset  dv\n";
    printf PARAMS "2        20         5.0         2.0        fresnelFiles 10   0";
    close (PARAMS);
    
    # Usage: ./detect elongation paramFile timeseries
    $command = "../$program $elong $incl params.detect ts.add 2> error.$program";
    
    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s "${ts_add_file}.hits");
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}


######################################################
#  test addKBO (fits)
#   
#######################################################

$program = "addKBO";
if ($program{"$program"} or $all ) {
    printf "Testing $program (fits) ... ";
    # Usage: ./addKBO elong incl fresnelfile timeseries N offset center[1|0] random[1|0]
    $command = "../$program $elong $incl ${fresDir}/$fresfile $tsfits 10 $offset 0 1 2> error.${program}";
    
    ( -s $ts_add_fits ) && unlink $ts_add_fits;
    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s $ts_add_fits);
    printf "%s\n", ($ran and $wrote) ? "success" : "failure";
    
}

#########################################################
#  test detect
#
#########################################################

$program = "detect";
if ($program{"$program"} or $all ) {

    printf "Testing $program (fits) ... ";
    
    
    $paramfile = "params.${program}";
    open (PARAMS, ">$paramfile");
    printf PARAMS "# cycles NperCycle  corrThresh  chiThresh  fresDir      Noffset  dv\n";
    printf PARAMS "2        20         5.0         2.0        fresnelFiles 10   1";
    close (PARAMS);
    
    # Usage: ./detect elongation paramFile timeseries
    $command = "../$program $elong $incl params.detect ts.fits.add.fits 2> error.$program";
    
    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s "${ts_add_fits}.hits");
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}



########################################################
#  test komplete
#
########################################################

$program = "komplete";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    # Usage: ./komplete elongation paramFile timeseries
    $command = "../$program $elong $incl params.detect ts.add  2> error.$program";
    

    my $out1 = "${ts_add_file}.hits";
    my $out2 = "${ts_add_file}.stats";

    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s $out1 and -s $out1);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}





########################################################
#  test basis
#
########################################################

$program = "basis";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    
    #Usage: basis elongation paramFile timeseries offset
    $command = "../$program $elong $incl params.detect $ts_add_file $offset 2> error.$program";
    
    my $out1 = "${ts_add_file}.Bhits";
    my $out2 = "${ts_add_file}.Bstats";
 
    my $ran = (mysystem($command,1)==0) ? 1:0;
    my $wrote = (-s $out1 and -s $out2);
    printf "%s\n", ($ran and $wrote) ? 	"success" : "failed";

}



########################################################
#  test makeKernel
#
########################################################

$program = "makeKernel";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    my $out1 = "kernel_ran_$dt";
    my $out2 = "kernel_noran_$dt";
    #Usage: ./makeKernel fresfile elong dt offset random[0|1] centering[0|1]
    $command1 = "../$program $fresDir/$fresfile $elong $incl $dt $offset $random $center > $out1 2> error.$program";
    $command2 = "../$program $fresDir/$fresfile $elong $incl $dt $offset $norandom $center > $out2 2>> error.$program";
    
    my $ran1 = (mysystem($command1,1)==0) ? 1:0;
    my $wrote1 = (-s $out1);
    my $ran2 = (mysystem($command2,1)==0) ? 1:0;
    my $wrote2 = (-s $out2);

    printf "%s\n", ($ran1 && $ran2 && $wrote1 && $wrote2) ? "success":"failed";

}



########################################################
#  test offsetPattern
#
########################################################

$program = "offsetPattern";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    my $out = "${fresfile}_off";
    #Usage: ./offsetPattern fresfile maxOffset nOffset
    $command = "../$program $fresDir/$fresfile $maxOffset $Noffset > $out  2> error.${program}";
    
    my $ran = (mysystem($command,1) == 0) ? 1:0;
    my $wrote = (-s $out);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}



########################################################
#  test xcorrelate
#
########################################################

$program = "xcorrelate";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";

    #Usage: xcorrelate elong fresnelfile timeseries offset corrThresh
    $command = "../$program $elong $incl $fresDir/$fresfile $ts_add_file $offset $corrThresh 2> error.${program}";

    my $out = "${ts_add_file}.xcor";    
    my $ran = (mysystem($command,1) == 0) ? 1:0;
    my $wrote = (-s $out);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}



########################################################
#  test xchi
#
########################################################

$program = "xchi";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";

    #Usage: xchi elong fresnelfile timeseries offset chiThresh
    $command = "../$program $elong $incl $fresDir/$fresfile $ts_add_file $offset $chiThresh 2> error.${program}";


    my $out = "${ts_add_file}.xchi";    
    my $ran = (mysystem($command,1) == 0) ? 1:0;
    my $wrote = (-s $out);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}


########################################################
#  test hideKBO
#
########################################################

$program = "hideKBO";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    
    #Usage: ../hideKBO elongation paramFile timeseries
    $command = "../$program $elong $incl $fresDir $tsfile 2> error.${program}";


    my $out = "${tsfile}.hide";    
    my $ran = (mysystem($command,1) == 0) ? 1:0;
    my $wrote = (-s $out);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}

$program = "hideKBO_TA";
if ($program{"$program"} or $all ) {
    
    printf "Testing $program ... ";
    
    mysystem("cp $tsfile $tsfileTA");
    
    #Usage: ../hideKBO_TA elongation paramFile timeseries
    $command = "../$program $elong $incl $fresDirTA $tsfileTA 2> error.${program}";


    my $out = "${tsfileTA}.hide";    
    my $ran = (mysystem($command,1) == 0) ? 1:0;
    my $wrote = (-s $out);
    printf "%s\n", ($ran and $wrote) ? "success" : "failed";

}


