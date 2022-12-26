#!/usr/bin/env perl
package PBandJ;  # Start the PBandJ namespace
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
use BatchJobTermination;

my $ME = basename(__FILE__);
# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(Cwd::abs_path(__FILE__)) . "/../bin";


# Replace value of perl variable in the form "$var = bla;"
# In contrast to ReplacePerlVarValue() this does not fail
# if the variable does not exist
sub ReplacePerlVarValueIfExists {
  my ($file, $var, $value) = @_;
  my $text = Utils::ReadFile($file);
  my $regex = qr/([\$\@]$var\s*=\s*)[^;]*;/;
  $text =~ s/$regex/$1$value;/;
  Utils::OverwriteFile($file,$text);
}

# Lev*_AA segments are moved to a backup directory for all Levs != $BaseLev.
# And the segments up to branching time are symlinked over from $BaseLev
# to all other Levs.
#
# When the next segments (for Lev != $BaseLev) are made in
# SubmitOtherLevsIfPBandJ(), the AmrTolerances.input files will be updated
# according to the Lev. To get the AmrTolerances.input for Lev != $BaseLev,
# we still need DoMultipleRuns to be run to generate AmrTolerances.input in
# Lev*_AA. This would already have been done in StartOtherLevs(). In this
# function, we rename Lev*_AA to Lev*_AA_PBandJBackup for Lev != $BaseLev,
# and symlink all segments of $BaseLev to the other Levs. This ensures
# that the data up to the branching time is the same for all Levs.
sub MakeLinksForPBandJ {
  my ($EvDir, $BaseLev, $BaseSegmentFullPath) = @_;

  # Get all Lev numbers, and sort them in increasing order
  my @Levs = sort(Utils::uniq(SpEC::GetLevs($EvDir)));

  foreach my $ThisLev (@Levs) {
    # Ignore $BaseLev
    unless($ThisLev == $BaseLev) {
      # Move Lev*_AA out of the way for all other Levs.
      # We want these segments so we can get the AmrTolerances.input files in
      # SubmitOtherLevsIfPBandJ(). Therefore, we save them as a backup
      # directory.
      my $Lev_AA_dir = "$EvDir/Lev$ThisLev"."_AA";
      my $PBandJBackupDir = "$Lev_AA_dir"."_PBandJBackup";
      if (-d $PBandJBackupDir || -l $PBandJBackupDir) {
        die "$PBandJBackupDir already exists, remove and try again.";
      }

      # If $Lev_AA_dir is a symlink at this stage, we skip $ThisLev as this
      # must mean that PBandJ links have already been made for this Lev and we
      # are presumably starting a different Lev.
      if (-l $Lev_AA_dir) {
          warn "Skipping PBandJ restart for Lev$ThisLev as $Lev_AA_dir ".
            "already exists as a symlink.";
          next;
      }
      rename($Lev_AA_dir, $PBandJBackupDir) || die "Cannot rename $Lev_AA_dir";

      # Make symlinks for all Lev$BaseLev_* segments to Lev$ThisLev_*
      # This way, the segments before branching are the same for all
      # Levs.
      foreach my $BaseLevSegment (glob "$EvDir/Lev$BaseLev"."_*[A-Z]") {
        (my $ThisLevSegment = $BaseLevSegment) =~ s/(Lev)${BaseLev}(_[A-Z]+)$/${1}${ThisLev}${2}/;
        if (-d $ThisLevSegment || -l $ThisLevSegment) {
          die "$ThisLevSegment already exists, cannot create symlink";
        }
        # Make a symlink to a relative path, in case the entire
        # directory tree gets moved.
        Utils::MySymlink(basename($BaseLevSegment), $ThisLevSegment);
        # Don't link segments past the BaseLevSegment.
        last if($BaseLevSegment eq $BaseSegmentFullPath);
      }
    }
  }
}

# For PBandJ (Perform Branching after Junk), when starting the other
# Levs, we start them at the branching time so that they have the same
# post-junk initial data.
# So, if PBandJ=1, we do the following:
#   1. Run DoMultipleRuns but don't start the other resolution jobs. This
#   creates the Lev*_AA segments for these Levs.
#   2. Move these out of the way and rename them Lev*_AA_PBandJBackup; we will
#   need these to get the AmrTolerances.input files for these Levs.
#   3. Create symlinks for all $BaseLev segments to the other Levs
sub StartOtherLevs {
  my ($workdir, $PBandJ) = @_;

  my ($EvDir, $BaseLev, $LastSegmentSuffix)
      = SpEC::SplitWorkDirPath($workdir);

  # Change work^ DoMultipleRuns.input to PBandJ=0
  ReplacePerlVarValueIfExists("$EvDir/DoMultipleRuns.input", "PBandJ", 0);
  # Change work^ DoMultipleRuns.input to EccRedRun=0
  ReplacePerlVarValueIfExists("$EvDir/DoMultipleRuns.input", "EccRedRun", 0);


  # Run the StartJob.sh script to submit the other Lev jobs, but use the
  # "-S" option to DoMultipleRuns and use GetResubCmd to make sure all
  # subtleties (i.e. submitting from compute nodes) are handled properly.
  my $text = Utils::ReadFile("$EvDir/StartJob.sh");

  # We want DoMultipleRuns to generate Lev*_AA for all Levs except $BaseLev
  my $exclude = "\'^Lev$BaseLev" . "_[A-Z]+\$\'";
  if ($PBandJ) {
    # If doing PBandJ, we check whether some PBandJ runs are already going.
    # if so, we tell DoMultipleRuns to exclude them also.
    my @LevsToExclude=();
    foreach my $Lev (sort(Utils::uniq(SpEC::GetLevs($EvDir)))) {
      next if $Lev == $BaseLev;
      if (-d "$EvDir/Lev${Lev}_AA") {
        if(-l "$EvDir/Lev${Lev}_AA") {
          push @LevsToExclude,$Lev;
        } else {
          die "We found Lev${Lev}_AA which is not a symlink.\n" .
              "However, this is a PBandJ run.  Something went wrong.\n";
        }
      }
    }
    if (scalar(@LevsToExclude) > 0) {
      my $excluded = "($BaseLev";
      foreach my $Lev (@LevsToExclude) {
        $excluded .= "|$Lev";
      }
      $excluded .= ")";
      $exclude = "\'^Lev$excluded" . "_[A-Z]+\$\'";
    }
  }
  $text =~ s/(DoMultipleRuns)/$1 -s $exclude/ || die "Could not match regexp";

  # This file is used to safeguard against starting Levs != BaseLev from t=0
  # for PBandJ. DoMultipleRuns would try to do this if you run StartJob.sh.
  # So, StartJob.sh looks for $prohibit_startjob_file and exits with a message
  # if it exists.
  my $prohibit_startjob_file = "ProhibitStartJobReruns.txt";

  if ($PBandJ) {
    # If not doing PBandJ, the new Levs will automatically get started.
    # For PBandJ, we run DoMultipleRuns with the -n options so that Lev*_AA is
    # generated for the new Levs but the jobs don't start. Below, in
    # MakeLinksForPBandJ() we will move these out of the way and symlink the
    # segments up to branching time from $BaseLev to the new Levs.
    $text =~ s/(DoMultipleRuns)/$1 -n/ || die "Could not match regexp";

    # If $prohibit_startjob_file already exists, PBandJ_StartJob.sh will
    # fail along with StartJob.sh. However, $prohibit_startjob_file would
    # already exist if, for example, Lev1 and Lev2 have already been started
    # with PBandJ, and now you want to start Lev4. We remove
    # $prohibit_startjob_file here so that PBandJ_StartJob.sh goes
    # through, but $prohibit_startjob_file will get rewritten below.
    unlink "$EvDir/$prohibit_startjob_file";
  }

  my $SubmitFile = "PBandJ_StartJob.sh";
  Utils::OverwriteFile("$EvDir/$SubmitFile", $text);
  chmod 0755, "$EvDir/$SubmitFile" || die "Could not chmod $EvDir/$SubmitFile";
  my $cmd = new Machines()->GetResubCmd(WorkDir => $EvDir,
                                        SubCmd  => "./$SubmitFile");
  Utils::System($cmd);

  # Setup the segments up to branching time for the new Levs
  if ($PBandJ) {
    MakeLinksForPBandJ($EvDir, $BaseLev, $workdir);

    # Create $prohibit_startjob_file. This is to safeguard against someone
    # running Levs other than $BaseLev from t=0. StartJob.sh will fail if
    # this file exists.
    my $text_prohibit_startjob =
      "This is a PBandJ (Perform Branching after Junk) run, which means that
branching into multiple Levs is done at a later time, e.g. after
eccentricity reduction. You are seeing this because you are presumably trying
to run a new Lev from Startjob.sh.

To run a new Lev, do the following instead:
Run MakeNextSegment from the Lev$BaseLev segment that terminates with
termination condition = EccentricityReduction or condition = PBandJTime.
For this run, that would be:
./bin/MakeNextSegment -d Lev$BaseLev\_$LastSegmentSuffix
This will start the other Levs (Lev != $BaseLev) as well.";

    Utils::OverwriteFile(
      "$EvDir/$prohibit_startjob_file",  $text_prohibit_startjob);
  }
}



# For PBandJ (Perform Branching after Junk), in StartOtherLevs(), symlinks
# would already have been created from segments of $BaseLev to the other
# Levs. In addition, a directory called Lev*_AA_PBandJBackup would have been
# created for these Levs. However, StartOtherLevs() does not start the jobs for
# these Levs when PBandJ=1, so we do that in SubmitOtherLevsIfPBandJ(). If this
# is not a PBandJ run, SubmitOtherLevsIfPBandJ() does nothing.
#
# $BaseLevNewSegment is the newest segment of $BaseLev (e.g. Lev3_AF).
# This should be the segment right after the branching time has been reached.
# In SubmitOtherLevsIfPBandJ(), we make a copy of $BaseLevNewSegment to other
# Levs, for example Lev3_AF gets copied to Lec2_AF and Lev1_AF. Then for these
# segments (for Lev != $BaseLev) we copy the AmrTolerances.input file from
# Lev*_AA_PBandJBackup to Lev*_AF. Finally, we delete the Lev*_AA_PBandJBackup
# dirs and start Lev*_AF for all Levs except $BaseLev. Note that
# Lev$BaseLev_AF will by started by EvolutionWrapper.
sub SubmitOtherLevsIfPBandJ {
  my ($scratch, $BaseLevNewSegment, $DoEverythingExceptFinalSubmission) = @_;

  my ($EvDir, $BaseLev, $LastSegmentSuffix)
    = SpEC::SplitWorkDirPath($BaseLevNewSegment);

  # Get all Lev numbers, and sort them in increasing order
  my @Levs = sort Utils::uniq(SpEC::GetLevs($EvDir));

  ## Loop over all Levs, copy $BaseLevNewSegment over to other Levs,
  ## update the AmrTolerances.input from Lev*_AA_PBandJBackup, delete
  ## Lev*_AA_PBandJBackup, and submit the new segment for Lev != $BaseLev.
  foreach my $ThisLev (@Levs) {
    # Ignore $BaseLev as this will be submitted by EvolutionWrapper anyway
    unless($ThisLev == $BaseLev) {
      # The name of the next segment for $ThisLev
      my $ThisLevNewSegment = "$EvDir/Lev$ThisLev"."_$LastSegmentSuffix";
      # Backup dir containing AmrTolerances.input for $ThisLev
      my $PBandJBackupDir ="$EvDir/Lev$ThisLev"."_AA_PBandJBackup";

      # If $PBandJBackupDir does not exist, that means that $PBandJ=0 for
      # this run, and this function does nothing.
      if (-d $PBandJBackupDir) {
        # Get input files from the old segment
        my $InputFiles
        = BatchJobTermination::ReadInputFiles($BaseLevNewSegment);

        # Make the new segment
        mkdir $ThisLevNewSegment || die "$ME: Cannot mkdir $ThisLevNewSegment";
        BatchJobTermination::LinkSpec($BaseLevNewSegment,$ThisLevNewSegment);

        # Update AmrTolerances.input
        $InputFiles->{"AmrTolerances.input"}
        = Utils::ReadFile("$PBandJBackupDir/AmrTolerances.input");

        # Update Evolution.input
        SpEC::ChangeEvolutionInputToRestartFromLastStep
            ($scratch, 1, $InputFiles, "Evolution.input");
        {
          # And now flag all subdomains as changed.
          my $to_add = "FlagAllSubdomainsAsChanged = yes; " .
                       "ResetFilterFunctionsForAMR = yes; " .
                       "StartAmrWithMinimumExtents = yes; ";
          $InputFiles->{"Evolution.input"} =~
              s|(GlobalVarsCheckpoint\s*=\s*Interpolated\()|${1}$to_add|
              or die "regexp not matching";
        }

        my $domaininput = "GrDomain.input";
        unless(exists($InputFiles->{$domaininput})) {
          $domaininput = "Domain.input";
          unless(exists($InputFiles->{$domaininput})) {
            die "Cannot find Domain.input or GrDomain.input ".
                "in $BaseLevNewSegment\n";
          }
        }

        # Update Lev in Domain.input and AmrDriver.input
        foreach my $file ($domaininput,"AmrDriver.input") {
          $InputFiles->{$file}
          =~ s|^(\s*Level\s*=\s*)${BaseLev};$|${1}${ThisLev};|m
              || die "regexp not matching in $file";
        }

        BatchJobTermination::WriteInputFiles($ThisLevNewSegment,$InputFiles);

        # We don't need PBandJBackupDir anymore
        rmtree("$PBandJBackupDir") || die "Could not rmtree $PBandJBackupDir";

        # Update the Jobname in MakeSubmit.py from $BaseLev to $ThisLev
        my $Jobname = Utils::SystemOutput(
                        "$RealBin/MakeSubmit.py -d $BaseLevNewSegment query " .
                        "Jobname");
        chomp($Jobname);
        $Jobname
            =~ s/_Ev.$BaseLev$/_Ev.$ThisLev/ || die "regexp not matching";
        Utils::System(
            "$RealBin/MakeSubmit.py -d $ThisLevNewSegment update " .
            "--Jobname $Jobname");

        # Submit the new segment for $ThisLev
        unless ($DoEverythingExceptFinalSubmission) {
          my $cmd = new Machines()->GetResubCmd(WorkDir => $ThisLevNewSegment);
          Utils::System($cmd);
        }
      }
    }
  }
}
