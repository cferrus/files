#!/usr/bin/env perl
package EccReduce;  # Start the EccReduce namespace
use strict;
use warnings FATAL=>'all';
use List::Util qw/max/;
use Cwd;
use File::Basename;
use File::Copy qw/move/;
use File::Path;

use lib dirname(Cwd::abs_path(__FILE__));
use Utils;
use Machines;
use SpEC;
use PBandJ;

# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(__FILE__) . "/../bin";

# Parse the output files from OmegaDotEccRemoval.py
# If the eccentricity is large (>0.01), then the spin-spin terms
# should be negligible. For these cases, we do not use the spin-spin
# fits because they will overfit d(Omega)/dt.
sub ParseParams {
  my ($dir,$ignore_spin_spin_for_high_ecc) = @_;
  my $include_target_ecc = not $ignore_spin_spin_for_high_ecc;
  my @Params = ParseParamsDat("$dir/Params_F2cos2_SS.dat",$include_target_ecc);
  if ($ignore_spin_spin_for_high_ecc && $Params[0] > 0.01) {
    @Params = ParseParamsDat("$dir/Params_F2cos2.dat",$include_target_ecc);
  }
  return @Params;
}

# Parse a Params_*.dat file from OmegaDotEccRemoval.py
sub ParseParamsDat {
  my ($file,$include_target_ecc) = @_;
  -f $file || die("Could not find $file");
  my @text = Utils::ReadFile($file);
  my @NewParams = split(/\s/,$text[-1]);
  my @OldParams = split(/\s/,$text[-2]);
  my $Ecc = $OldParams[3];
  my $Omega0 = $NewParams[0];
  my $adot0 = $NewParams[1]*1e-4;
  if($include_target_ecc) {
    my $EccT = $NewParams[3];
    return ($Ecc,$EccT,$Omega0,$adot0);
  }
  return ($Ecc,0.0,$Omega0,$adot0);
}

# Replace value of perl variable in the form "$var = bla;"
sub ReplacePerlVarValue {
  my ($file, $var, $value) = @_;
  my $text = Utils::ReadFile($file);
  my $regex = qr/([\$\@]$var\s*=\s*)[^;]*;/;
  die "Could not find regex $regex in file $file\n" if (not $text =~ m/$regex/);
  $text =~ s/$regex/$1$value;/;
  Utils::OverwriteFile($file,$text);
}

# Update next iteration MakeSubmit.input to increment Jobname.
# Replaces Ecc<N-1> with Ecc<N>, does nothing if no match.
sub ReplaceJobname {
  my ($dir, $N) = @_;
  my $Jobname = Utils::SystemOutput("$RealBin/MakeSubmit.py " .
                                    "-d $dir query Jobname");
  chomp($Jobname);
  $Jobname =~ s/Ecc\d+/Ecc$N/;
  Utils::System("$RealBin/MakeSubmit.py -d $dir update --Jobname $Jobname");
}

# This script estimates orbital eccentricty to get better initial parameters.
# It then creates new initial data with the better parameters.
# Returns (NewWork, Resubmit, IsError).
sub Main {
  my ($work, $scratch) = @_;

  # Determine the RunHome (directory with Ev/ and ID/), Ecc, and Lev
  my $PathRegexp = qr|^(.*)/Ecc(\d+)/(Ev[^/]*)/Lev([-\d]+)_|;
  $work =~ m|$PathRegexp| ||
      die "Path $work does not match expected structure ${PathRegexp}.\n";
  my $EccBase = $1 . "/Ecc";
  my $EccNum  = $2;
  my $RunHome = $EccBase . $EccNum;
  my $EvDir = "$RunHome/$3";
  my $Lev = $4;

  # sanity checks
  die "$RunHome does not exist" if (not -d $RunHome);
  die "$EvDir does not exist" if (not -d $EvDir);
  die "Could not parse $work to find level" if (not defined $Lev);
  print STDERR "RunHome=$RunHome\nLev=$Lev\nEccBase=$EccBase\nEccNum=$EccNum\n";

  # Parse the target ecc (in sci notation) from the termination string
  my $text = Utils::ReadFile("$scratch/EccentricityReduction.output");
  (my ($TargetEcc) = $text =~ /TargetEcc\s*=\s*([0-9eE\.\+\-]+);/) ||
    die "Could not parse TargetEcc from output file";
  (my ($MaxIters) = $text =~ /MaxIts\s*=\s*(\d+);/) ||
    die "Could not parse MaxIts from output file";
  (my ($Continue) = $text =~ /ContinueAfterTargetEccReached\s*=\s*(1|0);/) ||
    die "Could not parse ContinueAfterTargetEccReached from output file";
  (my ($RoughEccReduction) = $text =~ /RoughEccReduction\s*=\s*(1|0);/) ||
    die "Could not parse RoughEccReduction from output file";
  (my ($PBandJ) = $text =~ /PBandJ\s*=\s*(1|0);/) ||
    die "Could not parse PBandJ from output file";

  # Stop here if this is NOT the highest resolution present
  if (max(SpEC::GetLevs($EvDir)) > $Lev) {
    warn("Lev$Lev is not the highest level, so it will not be used for eccentricity reduction.");
    return ($RunHome, 0, 0);
  }

  # Use CombineSegments to combine OmegaMag.dat files for the current level
  my $CombineDir = "$EvDir/JoinedForEcc";
  chdir($EvDir) || die("Could not chdir $EvDir");
  my $sing_cmd = new Machines()->GetSingularityCmdIfNeeded();
  if ($sing_cmd ne "") {
    $sing_cmd .= " "; # Add space between 'singularity exec' and command.
  }
  $sing_cmd .= "$RealBin/CombineSegments.py";
  Utils::System("${sing_cmd} -f Horizons.h5 "
                ."-e h5 -o $CombineDir -L $Lev");
  chdir($CombineDir) || die("Could not chdir $CombineDir");

  # Here we need to determine whether we are doing eccentricity reduction
  # (i.e. we want eccentricity to go to zero), or eccentricity control
  # (i.e. we have a target eccentricity).  To do this, determine if
  # we have a TargetParams.input with non-empty values.
  my $ecc_reduction = 1;
  my $TPFile = "$RunHome/TargetParams.input";
  if(-e $TPFile) {
    my $text = Utils::ReadFile($TPFile);
    # Matches '$MassRatio = bla;' with non-whitespace 'bla'.
    if($text =~ m/\$MassRatio\s*=\s*[^;]+;/) {
      $ecc_reduction = 0;
    }
  }

  if($ecc_reduction) {
    my $IDFile = "$RunHome/ID/EvID/ID_Params.perl";
    Utils::System("$RealBin/OmegaDotEccRemoval.py -t bbh ".
                  "-d=$CombineDir/ApparentHorizons --idperl=$IDFile ".
                  "--improved_Omega0_update --no_check");
  } else {
    my $IDFile = "$RunHome/Params.input";
    Utils::System("$RealBin/InitialDataAdjustment.py -p $IDFile -t $TPFile ".
                  "-d $CombineDir/ApparentHorizons");
  }

  # Read the new parameters and eccentricity from the F2cos2 fit
  my ($Ecc,$EccT,$Omega0,$adot0) = ParseParams($CombineDir,$ecc_reduction);

  my $NewParamsFile;
  if($ecc_reduction) {
    # Create the new Params.input file and store it in the JoinedForEcc dir.
    # Do this even if we will be terminating, in case it's needed later.
    $NewParamsFile = "$CombineDir/NewParams.input";
    Utils::MyCopy("$RunHome/Params.input", $NewParamsFile);
    ReplacePerlVarValue($NewParamsFile, "Omega0", $Omega0);
    ReplacePerlVarValue($NewParamsFile, "adot0", $adot0);
  } else {
    $NewParamsFile = "$CombineDir/Params.input";
  }

  # Successful termination conditions:
  #  * Eccentricity < Target, and using high resolution
  if ( abs($Ecc-$EccT) < $TargetEcc && !$RoughEccReduction) {
    print STDERR "abs($Ecc-$EccT) is below the target $TargetEcc. ".
                 "Successfully completed eccentricity removal!\n";
    if ($Continue) {
      # When eccentricity reduction successfully completes, start the other
      # resolution jobs branching off at EccReduction time
      PBandJ::StartOtherLevs($work, $PBandJ);
      # By passing the original work directory and requesting a resubmit,
      # we trigger the new segment logic in BatchJobTermination.
      return ($work, 1, 0);
    } else {
      return ($RunHome, 0, 0);
    }
  }

  if ($RoughEccReduction && abs($Ecc-$EccT) < $TargetEcc) {
    # Rough eccentricity reduction is complete, switch to final eccentricity
    # reduction at higher resolution for next iteration
    print STDERR "abs($Ecc-$EccT) is below the initial target $TargetEcc. ".
                 "Successfully completed rough eccentricity removal!\n";
    # Change DoMultipleRuns.input to RoughEccReduction=0, so that the next
    # iteration is done at high resolution
    ReplacePerlVarValue("$EvDir/DoMultipleRuns.input", "RoughEccReduction", 0);
  } elsif ($EccNum > 0 && $ecc_reduction) {
    # Failure termination conditions
    #  * Eccentricity > 0.9*[previous ecc]  (ecc reduction only)
    # Don't check this if rough ecc reduction has succeeded.
    my $PrevEccNum = $EccNum - 1;
    my $OldEccDir = "$EccBase$PrevEccNum/Ev/JoinedForEcc";
    print STDERR "Getting previous ecc from $OldEccDir\n";
    my ($OldEcc,$OldEccT,$OldOmega0,$Oldadot0) =
        ParseParams($OldEccDir,$ecc_reduction);
    if (abs($Ecc-$EccT) > abs($OldEcc-$OldEccT)) {
      # Ecc appears to be converging too slowly.
      # However, if we have just changed Lev
      # (e.g. this iteration is a high-res ecc reduction using Lev3 but the
      #  previous iteration is a low-res ecc reduction using Lev1),
      # then we are ok.
      my $OldLev = max(SpEC::GetLevs("$EccBase$PrevEccNum/Ev"));
      if($OldLev == $Lev) {
        print STDERR "Eccentricity converging too slowly. ".
                     "Current ecc error is\n",
                     abs($Ecc-$EccT), " for ecc $Ecc, target $EccT\n".
                     "Old ecc error is ", abs($OldEcc-$OldEccT)," for".
                     " ecc $OldEcc, target $OldEccT.\n";
        if ($RoughEccReduction) {
          # If Ecc is converging too slowly in RoughEccReduction we just switch
          # to high resolution for the next iteration
          print STDERR "Stopping RoughEccReduction and switching to higher ".
                "resolution for next eccentricity reduction iteration.\n";
          ReplacePerlVarValue("$EvDir/DoMultipleRuns.input",
              "RoughEccReduction", 0);
        } else {
          # If Ecc is converging too slowly for high resolution ecc reduction,
          # we exit.
          return ($RunHome, 0, 1);
        }
      }
    }
  }

  if ($EccNum >= $MaxIters) {
    # Another failure condition.  Check no matter what.
    print STDERR "The number of eccentricity-reduction iterations exceeded ".
                 "the allowed maximum, $MaxIters.\n";
    if ($RoughEccReduction) {
      # If MaxIters is reached in RoughEccReduction we just switch
      # to high resolution for the next iteration
      print STDERR "Stopping RoughEccReduction and switching to higher ".
            "resolution for next eccentricity reduction iteration.\n";
      ReplacePerlVarValue("$EvDir/DoMultipleRuns.input",
          "RoughEccReduction", 0);
    } else {
      # If MaxIters is reached for high resolution ecc reduction, we exit.
      return ($RunHome, 0, 1);
    }
  }

  # If we get here, we are starting the next round of ecc reduction.
  print STDERR "abs($Ecc-$EccT) > $TargetEcc, ".
      "starting next round of ecc reduction.\n";

  # Create a new EccN directory to start the next iteration
  # and copy the ID job files there (i.e. don't edit existing files)
  my $NextEccNum = $EccNum + 1;
  my $NextRunHome = "${EccBase}${NextEccNum}";
  my $NextEvDir = "$NextRunHome/Ev";
  mkdir($NextRunHome) || die "Could not make next ecc dir '$NextRunHome'";
  mkdir($NextEvDir) || die "Could not make next Ev dir '$NextEvDir'";

  SpEC::CopyInputFiles($RunHome,$NextRunHome);
  SpEC::CopyInputFiles($EvDir,$NextEvDir);
  unlink "$NextRunHome/SpEC.stderr", "$NextRunHome/SpEC.stdout"; #if they exist
  if (-d "$RunHome/newbin") {
    # This is in case someone ran BFI's RestartRun --update-code.
    # BFI will make a newbin directory in $RunHome in that case,
    # so that it will be used for future ecc runs.
    Utils::LinkToOriginal("$RunHome/newbin", "$NextRunHome/bin");
  } else {
    Utils::LinkToOriginal("$RunHome/bin", "$NextRunHome/bin");
  }
  # We link from $work/bin instead of from $EvDir/bin in case
  # someone has changed the bin directory during the run,
  # e.g. by BFI's RestartRun --update-code.
  Utils::LinkToOriginal("$work/bin", "$NextEvDir/bin");

  unless($ecc_reduction) {
    # Remove the "StartJob.sh" script from the new RunHome because
    # we don't want the user to be calling it!  If the user calls it,
    # it will reset Params.input to its initial value in Ecc0.
    unlink("$NextRunHome/StartJob.sh");
  }

  # Copy the updated Params.input to the next Ecc directory
  Utils::MyCopy($NewParamsFile, "$NextRunHome/Params.input");

  # Change Jobname in MakeSubmit.input to reflect eccentricity iteration
  ReplaceJobname($NextRunHome, $NextEccNum);

  # Resubmit=2 submits an ID job
  return ($NextRunHome, 2, 0);
}

#--------------------------------------------------------------------
1; # This is how to end sources
