#!/usr/bin/env perl
package BatchJobTermination;  # Start the BatchJobTermination namespace
use Exporter;                 # load Exporter module
@ISA=qw(Exporter);            # Inherit from Exporter
@EXPORT=qw(SetupForNextRun GetTerminationCondition);

require 5;
use warnings FATAL => 'all';
use strict;
use Cwd;
use File::Basename;
use File::Path qw(mkpath rmtree);
use POSIX qw(floor);
use List::Util qw(max);
use re 'eval';
# Load SpEC perl modules
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;
use SpEC;
use EccReduce;
use PBandJ;
use Machines;

my $ME = basename(__FILE__);
# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(Cwd::abs_path(__FILE__)) . "/../bin";
my $opt_v = 0;  # default verbosity

# Used by SpECBatchScript after the evolution executable has exited
# OUTPUT ARGS:
#   ContinueJob = bool  # True if continuing in the same batch job
#   WorkDir = string    # location of the new work directory
#   Resubmit = int      # 0 -> restart in the same job
#                       # 1 -> resubmit job
#                       # 2 -> resubmit job, but runs MakeSubmit.py itself
#   DoneSegment = bool  # True if the current segment is complete
sub SetupForNextRun {
  my ($workdir, $scratch, $opt_a, $opt_f, $WalltimeMins,
      $StartTime, $opt_e, $opt_b, $in_v, $JobID) = @_;
  $opt_v = $in_v if (defined $in_v);

  Utils::MakePathsAbsolute(\$workdir,\$scratch);

  print STDERR "workdir='",$workdir,"'\nscratch='",$scratch,"'\n" if $opt_v;
  chdir $workdir || die "$ME: Cannot chdir $workdir";
  my $OrigWorkdir = $workdir;   # Save for later

  my $resubmit = 0;
  my $allow_change_cores = 0;
  my $cleanup_checkpoints = 0;

  # Interpolation option passed to ChangeEvolutionInputToRestartFromLastStep.
  # If it retains this default value, then will not change the restart string.
  # Otherwise, this value must be an integer.
  my $change_restart = "no";

  # Determine if this job is over
  my($ContinueThisRun, $reason, $IsError) = GetTerminationCondition($scratch);

  if($ContinueThisRun) {
    print STDERR "$ME: This run should be continued: $reason\n" if ($opt_v);

    # Decide how to continue this job:
    {
      my $allowedseconds = 60*$WalltimeMins;
      my $curtime        = time;
      my $elapsedtime    = $curtime-$StartTime;
      if($allowedseconds-$elapsedtime < ($opt_e+2*$opt_b)*60) {
        print STDERR "$ME: Not enough wallclock time left; ",
                     "Must submit a new batch job.\n" if ($opt_v);
        $resubmit = 1;
      }
    }

    my $InputFiles = ReadInputFiles($workdir);

    # Check the reason for termination
    if($reason =~ /Only\S+OnSlice/ or $reason =~ /Only\S+OnSphericalBdry/) {
      ($workdir, $change_restart)
          = DropInnerShells($workdir,$scratch,$InputFiles);
    } elsif ($reason =~ /SplitOrMergeSubdomains/) {
      # Make new workdir, replace files in it.
      ($workdir, $change_restart)
          = SplitOrMergeSubdomains($workdir,$scratch,$InputFiles);
      $allow_change_cores = 1;
    } elsif ($reason =~ /NsNsCommonHorizon/) {
      # Always resubmit, since number of nodes may change.
      $resubmit = 2;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      $workdir = StartRingdown($workdir,$scratch,
                               "--hydro --inertialframe SpectralInertial",
                               "AHfinder",1,48);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /^CommonHorizon/) {
      # Always resubmit, since number of nodes may change.
      $resubmit = 2;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      $workdir = StartRingdown($workdir,$scratch,"--autorminfac",
                               "ForContinuation",0,48);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /DividedGridBoxesTouching/) {
      # Always resubmit, since number of nodes may change.
      # Use number of nodes in Submit.input rather than overwriting them on the
      # command line.
      $resubmit = 2;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      $workdir = StartNsNsMergedBoxes($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /SingleStarHorizon/) {
      # Always resubmit, since number of nodes may change.
      $resubmit = 2;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      $workdir = StartRingdown($workdir,$scratch,
                               "--hydro --inertialframe SpectralInertial ".
                               "--nopremergermap",
                               "AHfinder",1,undef);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /PreserveRelativeDeltaR[AB]/) {
      # Make new workdir, replace files in it.
      ($workdir, $change_restart) =
          PreserveRelativeDeltaRAB($workdir,$scratch,$InputFiles);
      $allow_change_cores = 1;
    } elsif ($reason =~ /ChangeCutSphereGridPlaneCoordX/) {
      # Make new workdir, replace files in it.
      ($workdir, $change_restart) =
        ChangeCutSphereGridPlaneCoordX($workdir,$scratch,$InputFiles);
    } elsif ($reason =~ /BhNsMergerBegins/) {
      # Always resubmit, since number of nodes may change.
      $resubmit = 1;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      # Make new workdir, replace files in it.
      $workdir=BhNsMerger($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /BhNsPlungeBegins/) {
      # Always resubmit, since number of nodes may change.
      $resubmit = 1;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      # Make new workdir, replace files in it.
      $workdir=BhNsPlungeBegins($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /BhNsSettleDisk/) {
      # Make new workdir, replace files in it.
      $workdir=BhNsSettleDisk($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /BhNsHydroAMRRestart/ or $reason =~ /NsNsHydroAMRRestart/) {
      UpdateHyDomainInput($scratch,$InputFiles);
      $change_restart = 3;
      $cleanup_checkpoints = 1;
    } elsif ($reason =~ /BhNsExpandGrid/ or $reason =~ /NsNsExpandGrid/) {
      ExpandHydroGrid($scratch,$workdir,$InputFiles);
      $change_restart = 3;
      $cleanup_checkpoints = 1;
    } elsif($reason =~ /NsNsDHGauge/){
      $resubmit = 1;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      $workdir=NsNsDHGauge($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    }elsif ($reason =~ /NsNsMergerBegins/) {
      $resubmit = 1;
      SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
      # Make new workdir, replace files in it.
      $workdir=NsNsMerger($workdir,$scratch);
      $InputFiles={}; # No further input file changes needed.
    } elsif ($reason =~ /FoshRhsIsNotBalanced/ or $reason =~ /MemoryIsNotBalanced/) {
      AddOwnerPolicyOutputToDomainInput($scratch,$InputFiles);
      $change_restart = 0;
    } elsif ($reason =~ /FileExists/) {
      AddOwnerPolicyOutputToDomainInput($scratch,$InputFiles);
    } elsif ($reason =~ /^EccentricityReduction$/) {
      $ContinueThisRun = 0; # special termination
      ($workdir, $resubmit, $IsError) = EccReduce::Main($workdir, $scratch);
      unless ($resubmit and ($OrigWorkdir eq $workdir)) {
        # If (resubmit and origwordir eq workdir), then later the
        # code will try to do a wallclock restart, so we need to keep
        # the input files.  Otherwise, EccReduce::Main has already
        # copied input files around and we need to remove the $InputFiles
        # hashref here so that input files are not overwritten.
        $InputFiles = {};
      }
      SpEC::JobNotify("Stop EccRed", $OrigWorkdir,
                      $scratch, $opt_a, $opt_f, $JobID) if (not $resubmit);
    } elsif ($reason =~ /^PBandJTime$/) {
      $IsError = 0;
      $ContinueThisRun = 0; # special termination
      $resubmit = 1;
      PBandJ::StartOtherLevs($workdir, 1);
      unless ($resubmit and ($OrigWorkdir eq $workdir)) {
        # If (resubmit and origwordir eq workdir), then later the
        # code will try to do a wallclock restart, so we need to keep
        # the input files.
        $InputFiles = {};
      }
    } elsif ($reason =~ /EvolveGeodesicsWallClock/) {
      $workdir =EvolveGeodesicsWallClock($workdir,$InputFiles);
      $resubmit = 1;
    } elsif ($reason =~ /EvolveGeodesicsNewLev/) {
      $workdir =EvolveGeodesicsNewLev($workdir,$InputFiles);
      $resubmit = 1;
    } elsif ($reason =~ /^IngoingCharFieldOnSphericalBdry$/ ||
             $reason =~ /^ProportionalIntegral::MaxDt$/ || 
             $reason =~ /^ProportionalIntegral::MinDt$/) {
      $workdir  = 
          RestartWithBetterControlSystemInitialization($workdir,
                                                       $scratch,$InputFiles);
    } elsif ($reason =~ /^WallClock$/ or $reason =~ /^Preemption$/) {
      $workdir = WallClockContinuation($workdir,$scratch,$InputFiles);
      $resubmit = 1;
      $allow_change_cores = 1;
      $change_restart = 0;
    } elsif ($reason =~ /DropJunkShellTime/) {
      ($workdir, $change_restart)
          = DropJunkShell($workdir, $scratch, $InputFiles);
    } else {
      die "$ME: Don't know about reason $reason";
    }

    # If requested (i.e. for a number of termination criteria
    # only used for hydro with FMR), cleanup extraneous checkpoints
    if ($cleanup_checkpoints) {
      CleanupHydroCheckpointDirectories($OrigWorkdir);
    }

    # If expecting a resubmission, make sure the Work directory increments.
    # This can occur when job runs out of time during a non-WallClock restart.
    if ($resubmit and ($OrigWorkdir eq $workdir)) {
      my $dum = undef;
      $workdir = WallClockContinuation($workdir,$scratch,$InputFiles);
      if ($change_restart eq "no") {
        $change_restart = 0;
      }
    }

    # Change restart string after all other logic that can affect the restart
    # string (i.e. work directory, interpolation option) has concluded.
    if ($change_restart ne "no") {
      SpEC::ChangeEvolutionInputToRestartFromLastStep(
        $scratch, $change_restart, $InputFiles, "Evolution.input");
    }

    WriteInputFiles($workdir,$InputFiles);

  } else {
    print STDERR "$ME: This run should NOT be continued: $reason\n" if($opt_v);

    # Send an email message about termination
    SpEC::JobNotify($reason, $workdir, $scratch, $opt_a, $opt_f, $JobID);
  }

  # The current method to change the number of cores on the fly does not
  # work for hydro-FMR; do not allow changes in core numbers for now.
  # Note that SpEC::IsHydro requires input files, which are not there
  # if $workdir is a higher-level Ecc directory, so we check that case first.
  if ((not $workdir =~ m|/Ecc\d+$|) and SpEC::IsHydro($workdir))
  {
      $allow_change_cores = 0;
  }
  # Potentially change MakeSubmit.input and $resubmit
  $resubmit = UpdateSubmit($resubmit, $allow_change_cores, $workdir);

  # For PBandJ (Perform Branching after Junk) the restart does not
  # submit the jobs for the new Levs, so we do that here. If this is
  # not a PBandJ run, PBandJ::SubmitOtherLevsIfPBandJ() does nothing.
  if (($reason =~ /^EccentricityReduction$/) or ($reason =~ /^PBandJTime$/)) {
    # This is a hack for some specific cases such as
    # Support/Tests/TestAutoEccRed, for which $workdir here is the Ecc dir.
    # PBandJ::SubmitOtherLevsIfPBandJ() only works when $workdir is the new
    # segment created after the branching time has been reached.
    unless ($workdir =~ m|/Ecc\d+$|) {
      PBandJ::SubmitOtherLevsIfPBandJ($scratch, $workdir);
    }
  }

  # Determine if the run is continuing in the same batch job
  my $ContinueJob = ($ContinueThisRun and not $resubmit) ? 1 : 0;

  # Determine if we're done running in this segment
  my $DoneSegment = 0;
  unless ($ContinueThisRun and ($OrigWorkdir eq $workdir)) {
    $DoneSegment = 1;
    print STDERR "$ME: Segment in $OrigWorkdir is complete.\n" if ($opt_v);
  }

  return ($ContinueJob, $workdir, $resubmit, $DoneSegment, $IsError);
}

#==============================================================================
# Private Subroutines
#==============================================================================

# * In some cases, we allow the number of cores to change, which affects
#   all resubmission parameters.
# * If restarting (in mode 1), construct the new Submit.sh
#   Mode 2 restarts are supposed to handle their own resubmission options
sub UpdateSubmit {
  my ($resubmit, $allow_change_cores, $workdir) = @_;

  my $QueryInWorkdir = sub {
    my ($field) = @_;
    my $out = Utils::SystemOutput("$RealBin/MakeSubmit.py " .
                                  "-d $workdir query $field");
    chomp($out);
    return $out;
  };
  my $UpdateInWorkdir = sub {
    my ($opts) = @_;
    Utils::System("$RealBin/MakeSubmit.py -d $workdir update $opts");
  };

  my $ForceCPN   = $QueryInWorkdir->("ForceCoresPerNode") eq "True" ? 1 : 0;
  my $ForceCores = $QueryInWorkdir->("ForceCores") eq "True" ? 1 : 0;
  my $Cores      = $QueryInWorkdir->("Cores");
  my $CPN        = $QueryInWorkdir->("CoresPerNode");

  if (($allow_change_cores and not $ForceCores) || $resubmit==1) {
    my ($Nmin, $Nmax, $N);
    if ($allow_change_cores and not $ForceCores) {
      my $NSubdomains = SpEC::TotalNSubdomains($workdir);
      $N    = floor( (48/62)*$NSubdomains );
      $Nmin = max(floor( (40/62)*$NSubdomains ), 1);
      $Nmax = $NSubdomains;
    } else {
      $Nmin = $Nmax = $N = $Cores;
    }

    # TODO: should also use this if Nprocs in Bbh/DoMultipleRuns.input is undef
    my ($NewN, $NewCPN, $NewResubmit, $queue) = Machines::GetNodeSetup_Resub(
      Nmin       => $Nmin,
      Nmax       => $Nmax,
      N          => $N,
      AllowedCPN => [ new Machines()->GetAllowedCPN() ],
      Resubmit   => $resubmit,
      N0         => $Cores,
      CPN0       => $CPN,
      ForceCPN   => $ForceCPN,
    );
    $resubmit = $NewResubmit;

    # NOTE: we take care above to not choose invalid combinations of
    # CoresPerNode and Cores if either ForceCoresPerNode or ForceCores.
    print STDERR "$ME: (N,CPN) = ($Cores,$CPN)->($NewN,$NewCPN)\n" if ($opt_v);
    if ($queue ne "SpECNoSuchQueue") {
        $UpdateInWorkdir->("--CoresPerNode $NewCPN " .
                           "--Cores $NewN --Queue $queue");
    } else {
        $UpdateInWorkdir->("--CoresPerNode $NewCPN --Cores $NewN");
    }
  }

  return $resubmit;
}

# Returns list of ($continue,$reason)
sub GetTerminationCondition {
  my($scratch)=@_;
  local $_;

  my $continue = undef;
  my $reason = "Unknown";
  my $error = undef;
  {
    my $file = "$scratch/TerminationReason.txt";
    unless(-f $file) {
      warn "$ME: Cannot open $file";
      return ($continue,$reason,$error);
    }
    my @text = Utils::ReadFile($file);
    foreach (@text) {
      if(/^Termination condition\s*(.*)/)    {$reason = $1;}
      if(/^To be continued\.\.\.$/)          {$continue = 1;}
      if(/^This termination is an error\.$/) {$error = 1;}
    }
  }
  return ($continue,$reason,$error);
}

#====================
# TerminationCriteria
#====================

sub RestartWithBetterControlSystemInitializationWork {
  # $ringdown is a bool, 1 if ringdown, 0 otherwise.
  my($oldworkdir,$oldscratchdir,$ringdown,$InputFiles)=@_;

  # This termination criterion will be attempted max_iter times.
  my $max_iter  = 3;

  # Move old work directory out of the way; will be replaced.
  my $workdir     = $oldworkdir;
  my $prevworkdir = $oldworkdir;
  $prevworkdir    =~ s|/Lev(\d_\S\S)$|/ControlSysInit0_Lev$1|;
  my $curr_iter   = 0;
  while(-d $prevworkdir) {
    unless ($prevworkdir =~ m|ControlSysInit(\d)_Lev|) {
      die
          "$ME: Prevworkdir $prevworkdir " .
          "doesn't match ControlSysInit(\\d)_Lev\n";
    }
    $curr_iter = $1+1;
    $prevworkdir =~ s|ControlSysInit(\d)_Lev|ControlSysInit${curr_iter}_Lev|;
    die "$ME: Too many iterations, $curr_iter\n" if($curr_iter > $max_iter);
  }
  rename($workdir,$prevworkdir) ||
      die "$ME: Cannot rename $workdir $prevworkdir\n";

  # Copy from old work directory to new work directory.
  if(-d $workdir) {
    die "$ME : Cannot submit to $workdir because it already exists";
  }
  mkdir $workdir    || die "$ME: Cannot mkdir $workdir";

  LinkSpec($prevworkdir,$workdir);

  # Move shape-initialization files into new work directory.
  my $ShapeInit = "$workdir/ControlSystemInit";
  {
    my $ShapeInitOld  = "$prevworkdir/ControlSystemInit$$";
    rename($ShapeInitOld,$ShapeInit) ||
        die "$ME: Cannot rename $ShapeInitOld $ShapeInit\n";
  }

  # Turn off TimeThresholdForError TerminationCriterion (to avoid
  # infinite loops).
  if($curr_iter>=$max_iter) {
    my $file = "Evolution.input";
    my $text = $InputFiles->{$file};
    # Note the /g modifier to catch all instances of TimeThresholdForError.
    if($text =~ s/(TimeThresholdForError\s*=\s*[^;]*;)/#$1/g) {
      $InputFiles->{$file} = $text;
    } else {
      die "$ME: Cannot find TimeThresholdForError TerminationCriterion\n";
    }
  }

  # Change source of initial shape in SpatialCoordMap.input
  {
    my $file = "SpatialCoordMap.input";
    my $text = $InputFiles->{$file};
    if($ringdown) {
      # Need to do all Init* files separately.
      $text =~ s|(\(\s*File\s*=\s*)[^;]+(/Init-[^;]*);|${1}${ShapeInit}${2};|g;
      $InputFiles->{$file} = $text;
    } else {
      # There is only one MatchDir, for all Init* files.
      if($text =~ s/^([^\#]*MatchDir\s*=\s*)[^;]*;\s*$/${1}${ShapeInit};/m) {
        $InputFiles->{$file} = $text;
      } else {
        die "$ME: Cannot find MatchDir in $file\n";
      }
    }
  }

  # Change SmoothAhRadiusFile in GrDataBoxItems.input
  {
    my $file = "GrDataBoxItems.input";
    my $text = $InputFiles->{$file};
    my $a    = "$ShapeInit" . "/";
    # The /g modifier below catches multiple cases.
    # And the /s modifier handles the case of a newline inbetween.
    if($text =~
       s|(SmoothAhRadiusFile\s*=\s*)[^;]*/([^;]+;)|$1$a$2|sg) {
      $InputFiles->{$file} = $text;
    } else {
      die "$ME: Cannot find SmoothAhRadiusFile in $file\n";
    }
  }

  return $workdir;
}

sub RestartWithBetterControlSystemInitializationRingdown {
  my($oldworkdir,$oldscratchdir,$InputFiles)=@_;

  # Create new shape-initialization files, and put them
  # into $oldworkdir.  Temporarily.  They will be moved later.
  my $ShapeInit0 = "$oldworkdir/ControlSystemInit$$";
  mkdir $ShapeInit0 || die "$ME: Cannot mkdir $ShapeInit0\n";

  # Get rminfac and rmaxfac from Continuation/Init-GridParams.txt
  my $rminfac = undef;
  my $rmaxfac = undef;
  my $time = undef;
  {
    my $file = "$oldworkdir/../Continuation/Init-GridParams.txt";
    my $text = Utils::ReadFile($file);
    $text =~ m|rAH\s*=\s*([^;]+);|s;
    my $rah  = $1;
    $text =~ m|Rmin\s*=\s*([^;]+);|s;
    my $rmin = $1;
    $text =~ m|Rtrans\s*=\s*([^;]+);|s;
    my $rmax = $1;
    $text =~ m|Time\s*=\s*([^;]+);|s;
    $time    = $1;
    $rminfac = $rmin/$rah;
    $rmaxfac = $rmax/$rah;
  }

  # Make new Init-* files.
  {
    my $opts    = "-p -rminfac $rminfac -minorder 0 " .
                  "-rmaxfac $rmaxfac -l=-1 -T $time";
    Utils::System("$RealBin/SurfaceToSpatialCoordMapFiles $opts " .
                  "$oldscratchdir/ApparentHorizons/AhCCoefs.dat");
    rename("Init-SmoothAhRadius.txt",
           "$ShapeInit0/Init-SmoothAhRadius.txt")
        || die "Cannot write $ShapeInit0/Init-SmoothAhRadius.txt\n";
    rename("Init-FuncLambdaFactor.txt",
           "$ShapeInit0/Init-FuncLambdaFactor.txt")
        || die "Cannot write $ShapeInit0/Init-FuncLambdaFactor.txt\n";
    rename("Init-FuncLambdaFactor0.txt",
           "$ShapeInit0/Init-FuncLambdaFactor0.txt")
        || die "Cannot write $ShapeInit0/Init-FuncLambdaFactor0.txt\n";
    rename("Init-FuncLambdaFactor1.txt",
           "$ShapeInit0/Init-FuncLambdaFactor1.txt")
        || die "Cannot write $ShapeInit0/Init-FuncLambdaFactor1.txt\n";
    rename("Init-ShapeTimeScales.txt",
           "$ShapeInit0/Init-ShapeTimeScales.txt")
        || die "Cannot write $ShapeInit0/Init-ShapeTimeScales.txt\n";
    rename("Init-GridParams.txt",
           "$ShapeInit0/Init-GridParams.txt")
        || die "Cannot write $ShapeInit0/Init-GridParams.txt\n";
  }

  # Make new Init-FuncTrans.txt
  {
    # First read AhCInertial.dat to get the center
    my $file = "$oldscratchdir/ApparentHorizons/AhCInertial.dat";
    my $text = Utils::ReadFile($file);
    #   Read header
    $text =~ m|^#\s*\[(\d+)\]\s*=\s*Center_x|m;
    my $center_x_col = $1-1;
    $text =~ m|^#\s*\[(\d+)\]\s*=\s*Center_y|m;
    my $center_y_col = $1-1;
    $text =~ m|^#\s*\[(\d+)\]\s*=\s*Center_z|m;
    my $center_z_col = $1-1;
    #   Write file
    my $outfile = "$ShapeInit0/AhC_Center.dat";
    my $outtext = "# [1] = time\n";
    $outtext .= "# [2] = Center_x\n";
    $outtext .= "# [3] = Center_y\n";
    $outtext .= "# [4] = Center_z\n";
    while($text =~ m/(.+)$/mg) {
      next if $1 =~ m/^#/; # Skip comments
      my @a = split(' ',$1);
      $outtext .= $a[0] . " " . $a[$center_x_col] . " " .
          $a[$center_y_col] . " " . $a[$center_z_col] . "\n";
    }
    Utils::OverwriteFile($outfile,$outtext);
    Utils::System("$RealBin/TrajectoryToSpatialCoordMapFiles " .
                  "-N -T $time $outfile");
    rename("Init-FuncTrans.txt",
           "$ShapeInit0/Init-FuncTransSingleHole.txt")
        || die "Cannot write $ShapeInit0/Init-FuncTransSingleHole.txt\n";
  }

  return RestartWithBetterControlSystemInitializationWork($oldworkdir,
                                                          $oldscratchdir,1,
                                                          $InputFiles);
};

sub RestartWithBetterControlSystemInitialization {
  my($oldworkdir,$oldscratchdir,$InputFiles)=@_;

  if ($oldworkdir =~ m|_Ringdown/Lev\S+_AA$|) {
    return RestartWithBetterControlSystemInitializationRingdown($oldworkdir,
                                                                $oldscratchdir,
                                                                $InputFiles);
  };

  # Sanity check 1: We are in Lev*_AA.
  # This criterion is meant only for the first few timesteps of a run.
  unless ($oldworkdir =~ m/_AA$/) {
    die "$ME: RestartWithBetterControlSystemInitialization on segment>AA\n";
  }

  # Sanity check 2: There is an ID/EvID directory.
  my $EvID = $oldworkdir;
  $EvID    =~ s|/Ev/Lev\d_[A-Z][A-Z]$|/ID/EvID|;
  unless (-d $EvID) {
    die "$ME: Cannot find $EvID\n";
  }

  # Create new shape-initialization files, and put them
  # into $oldworkdir.  Temporarily.  They will be moved later.
  my $ShapeInit0 = "$oldworkdir/ControlSystemInit$$";
  mkdir $ShapeInit0 || die "$ME: Cannot mkdir $ShapeInit0\n";

  #  #  Get rminfac from ID.
  my $rminfacA = undef;
  my $rminfacB = undef;
  {
    Utils::EvalPerlCodeInFile("$EvID/ID_Params.perl");
    if (not defined($Utils::ID_rExcA)) {
      die "ERROR: $EvID/ID_Params.perl does not define 'ID_rExc{A,B}'.\n".
          "See Support/LegacyFileConversion/UpdateIDParams.pl for help.";
    }
    $rminfacA = sqrt($Utils::ID_rExcA/$Utils::ID_rA);
    $rminfacB = sqrt($Utils::ID_rExcB/$Utils::ID_rB);
  }
  foreach my $AB ("A","B") {
    my $var     = "\$rminfac${AB}";
    my $rminfac = eval $var;
    my $opts    = "-p -rminfac $rminfac -minorder 0 -rmaxfac 10 -l=-1 -1 -T 0";
    # What dir is this executed in?
    Utils::System("$RealBin/SurfaceToSpatialCoordMapFiles $opts " .
                  "$oldscratchdir/ApparentHorizons/Ah${AB}Coefs.dat");
    rename("Init-SmoothAhRadius.txt",
           "$ShapeInit0/ID_Init_SmoothAhRadius${AB}.txt")
        || die "Cannot write $ShapeInit0/ID_Init_SmoothAhRadius${AB}.txt\n";
    rename("Init-FuncLambdaFactor.txt",
           "$ShapeInit0/ID_Init_FuncLambdaFactor${AB}.txt")
        || die "Cannot write $ShapeInit0/ID_Init_FuncLambdaFactor${AB}.txt\n";
    rename("Init-FuncLambdaFactor0.txt",
           "$ShapeInit0/ID_Init_FuncLambdaFactor${AB}0.txt")
        || die "Cannot write $ShapeInit0/ID_Init_FuncLambdaFactor${AB}0.txt\n";
  }
  # Also copy ID_Init_FuncTrans.txt which is needed by SpatialCoordMap.input
  {
    foreach my $f ("ID_Init_FuncTrans.txt") {
      Utils::MyCopy("$EvID/$f","$ShapeInit0/$f");
    }
  }

  return RestartWithBetterControlSystemInitializationWork($oldworkdir,
                                                          $oldscratchdir,0,
                                                          $InputFiles);
}

sub ChangeCutSphereGridPlaneCoordX {
  my($oldworkdir,$oldscratchdir,$InputFiles)=@_;

  my $workdir = GetNewSubmissionDirectoryAndLinkSpEC($oldworkdir);

  # read in new value for CutSphereGridPlaneCoordX and for
  # MinimumAllowedSmallerCutRadius
  # the format is
  # CutSphereGridPlaneCoordX=    -10.4403512033531864489;
  # MinimumAllowedSmallerCutRadius=    -10.4403512033531864489;
  my $CutXvalue = Utils::ReadFile("$oldscratchdir/ChangesToCutSphereGeometry.output");
  my $CutRvalue = $CutXvalue;
  # $CutXvalue =~ s/CutSphereGridPlaneCoordX.*=//;
  $CutXvalue =~ s/MinimumAllowedSmallerCutRadius.*;//;
  $CutXvalue =~ s/CutSphereGridPlaneCoordX.*=//;
  $CutXvalue =~ s/;.*//;
  $CutXvalue =~ s/\n//;
  chomp $CutXvalue;
  $CutRvalue =~ s/CutSphereGridPlaneCoordX.*;//;
  $CutRvalue =~ s/MinimumAllowedSmallerCutRadius.*=//;
  $CutRvalue =~ s/;.*//;
  $CutRvalue =~ s/\n//;
  chomp $CutRvalue;

  # replace the current cut sphere setting with the new one
  # in each of the following input files:
  # Domain.input, SpatialCoordMap.input, StateChangers.input
  foreach my $inputfile
  ("GrDomain.input","SpatialCoordMap.input", "GrStateChangers.input", "Evolution.input") {
    my $text = $InputFiles->{$inputfile};
    # the value of this option can either be <<DEFAULT>> or a double
    $text =~ s/(CutSphereGridPlaneCoordX\s*=)\s*<<DEFAULT>>\s*;/${1}$CutXvalue;/g;
    $text =~ s/(CutSphereGridPlaneCoordX\s*=)\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)?/${1}$CutXvalue/g;
    $text =~ s/(MinimumAllowedSmallerCutRadius\s*=)\s*<<DEFAULT>>\s*;/${1}$CutRvalue;/g;
    $text =~ s/(MinimumAllowedSmallerCutRadius\s*=)\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)?/${1}$CutRvalue/g;
    $InputFiles->{$inputfile} = $text;
  }

  # Change Evolution.input: restart with mapped interpolation,
  # as a changed cutX plane implies changed maps.
  my $change_restart = 2;

  # Make sure that the initial number of initial amr iterations in
  # Evolution.input is nonzero.  Otherwise, Domain.input will start from
  # the checkpoint file and will not notice the changes made here.
  AddOneAmrIterationIfCurrentlyZero($InputFiles, "Evolution.input");

  # Adjust OwnerPolicies.
  # RemoveOwnerPolicyFromDomainInput($workdir);
  AddOwnerPolicyOutputToDomainInput($oldscratchdir,$InputFiles);

  return ($workdir, $change_restart);
}


#---------------
sub PreserveRelativeDeltaRAB {
  my($oldworkdir,$oldscratchdir,$InputFiles)=@_;

  my ($workdir) = GetNewSubmissionDirectoryAndLinkSpEC($oldworkdir);

  # Read the set of Radii that should remain (i.e., not dropped).
  my %Radii;
  {
    my $file = "$oldscratchdir/SpEC.out";
    my @text = Utils::ReadFile($file);
    foreach (@text) {
      foreach my $ab ("A","B") { # Do for both A and B.
        if(/^PreserveRelativeDeltaR${ab}: spherical shell Radii\s*=(.*)$/){
          $Radii{$ab} = $1; chomp $Radii{$ab};
        }
      }
    }
    die "$ME: Cannot find Radii" unless (defined $Radii{"A"} or
                                         defined $Radii{"B"});
  }

  # Change Domain.input:
  # 1) remove NShells, RMin, RMax, RelativeDeltaR, Radii
  # 2) insert the Radii that we read from SpEC.out
  {
    my $file = basename(SpEC::GetGrDomainInputPath($oldscratchdir));
    my $text = $InputFiles->{$file};
    foreach my $ab ("A","B") { # Do for both A and B.
      next unless(defined $Radii{$ab});

      ### $bal matches balanced parens with possible nested parens inside.
      our $bal;  $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;

      $text      =~ m|(Sphere$ab\s*=\s*)($bal)|;
      my $head   = $1;
      my $inside = $2;
      my $pre    = $`;
      my $post   = $'; #'

      # Remove leading/trailing parens from $inside.
      $inside =~ s/^\(//;
      $inside =~ s/\)$//;

      # Split $inside at semicolons with a newline
      $inside =~ s|;([\t ]*[\S \t]+)|;\n$1|g;

      # We will use the 'Radii' option to CutSpheres::Sphere[AB], so
      # we comment out other options that conflict with 'Radii'.
      foreach my $i ("RMin","RMax","NShells","RelativeDeltaR",
                     "ScalingOfRadii") {
        $inside=~s/^(\s*)($i\s*=\s*[\S \t]+)$/$1#$2/m;
      }

      # Replace or add Radii option.
      unless($inside =~ s/^(\s*Radii\s*=\s*)[\S \t]+$/${1}$Radii{$ab};\n/m) {
        if(substr($inside,0,1) eq "\n") {
          $inside = "Radii=" . $Radii{$ab} . ";" . $inside;
        } else {
          $inside = "Radii=" . $Radii{$ab} . ";\n" . $inside;
        }
      }

      # If there is only one radius, this means that there are no more
      # spheres.  However, the Sphere[AB] option must still be specified
      # because it indicates the inner boundary for the wedges.
      # But in this case, we need to comment out many options that
      # are undefined, so we do that here.
      unless($Radii{$ab} =~ m/,/) { # Check only one radius by absence of comma
        foreach my $i ("L","RadialExtents","RadialMap",
                       "ReplaceSpheresWithWedges",
                       "ComputeWedgeGridSizesBasedOnL",
                       "RotationAboutZ", "RotationAboutRotatedY",
                       "RotationAboutRotatedZ") {
          $inside=~s/^(\s*)($i\s*=\s*[\S \t]+)$/$1#$2/m;
        }
      }

      # Reconstruct the contents of the input file
      $text = $pre . $head . "(" . $inside . ")" . $post;
    }
    $InputFiles->{$file} = $text;
  }

  # Change RefinementOptionHistory.input:
  # remove any references to subdomains beginning with
  #       Sphere$ab$ShellNumber
  #       Cylinder$ab$ShellNumber
  #       CylinderS$ab$ShellNumber
  #       FilledCylinder$ab$ShellNumber
  #       Wedge$ab$ShellNumber
  # where $ShellNumber is >= the number of remaining shells.
  {
    my $file = "RefinementOptionHistory.input";
    my $text = $InputFiles->{$file};
    my $pat  = "(?:Sphere|Cylinder|CylinderS|FilledCylinder|Wedge)";
    foreach my $ab ("A","B") { # Do for both A and B.
      next unless(defined $Radii{$ab});
      my $nshells = $Radii{$ab} =~ tr/,//; # Number of commas in $Radii{$ab}
      $text =~ s/^($pat$ab(\d+)\S+[,;]$)/($2>=$nshells ? "#" : "") .$1/mge;
    }
    # Comma at end becomes semicolon.
    $text =~ s/,$/;/;
    $text =~ s/,((\n#[^\n]*)+)$/;$1/;  # ..even if last line is commented...
    $InputFiles->{$file} = $text;
  }

  # Change Evolution.input: restart with mapped interpolation
  # Note that PreserveRelativeDeltaR can change the outer
  # radius of the shells as a means of preserving the
  # distance from these to the cutting plane.  This will
  # prevent ChangeNumberOfSphericalShells from triggering
  # for most of the run.
  my $change_restart = 2;

  # Make sure that the initial number of initial amr iterations in
  # Evolution.input is nonzero.  Otherwise, Domain.input will start from
  # the checkpoint file and will not notice the changes made here.
  AddOneAmrIterationIfCurrentlyZero($InputFiles, "Evolution.input");

  # Adjust OwnerPolicies.
  # RemoveOwnerPolicyFromDomainInput($workdir);
  AddOwnerPolicyOutputToDomainInput($oldscratchdir,$InputFiles);

  return ($workdir, $change_restart);
}

#---------------
sub BhNsPlungeBegins {
  my($oldworkdir,$oldscratch)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my ($workdir,$lev) = MakeSubdirWithSuffix($oldworkdir,"Plunge");

  #------------------------------------------------------------
  # Get parent of old workdir
  # and files copied from parent (i.e. before DoMultipleRuns).
  #------------------------------------------------------------
  my $oldworkdirparent="$oldworkdir/..";
  Utils::MakePathsAbsolute(\$oldworkdirparent);
  my @fromparent=("GaugeItems.input","Constraints.input",
                  "HySetupAndEvolution.input","HyDataBoxItems.input",
                  "HyStateChangers.input","HyObservers.input",
		  "NeutrinoLeakageItems.input");
  #---------------------
  # Copy to new workdir
  #---------------------
  my %filestocopy;
  foreach my $file (glob "$oldworkdir/*.input") {
    my $ignore=undef;
    foreach my $ff(@fromparent) {
      if($file =~ /$ff$/) {$ignore=1;last;}
    }
    next if($ignore);
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}=$file;
  }
  foreach my $file (@fromparent) {
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}="$oldworkdirparent/$file";
    die "Cannot find ",$filestocopy{$file2} unless (-f $filestocopy{$file2});
  }
  foreach my $file (keys %filestocopy) {
    my $file2 = $file;
    if($file2 =~ s/Plunge$// || !exists($filestocopy{"${file}Plunge"})) {
      Utils::MyCopy($filestocopy{$file},"$workdir/${file2}.input");
    }
  }
  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  #-----------------------------------------
  # Read time of checkpoint from oldscratch
  #-----------------------------------------
  my $tcheckpoint=FindLatestCheckpointTimeFromInputDir($oldscratch);

  #-----------------------------------------
  # CD to workdir and modify input files
  #-----------------------------------------
  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  # Modify input files
  {
    my $file = "DoMultipleRuns.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|(^\s*\$INSPIRALDIR\s*=\s*)[^;]*;|${1}"${oldscratch}";|m;
      $text    =~ s|(^\s*\$TDISRUPT\s*=\s*)[^;]*;|${1}${tcheckpoint};|m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

    $file = "DoMultipleRunsSettleDisk.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|(^\s*\$INSPIRALDIR\s*=\s*)[^;]*;|${1}"${oldscratch}";|m;
      $text    =~ s|(^\s*\$TDISRUPT\s*=\s*)[^;]*;|${1}${tcheckpoint};|m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

    $file = "Evolution.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s/UseInterpolatedRestartAtTime\s*=\s*\d+.\d+;/UseInterpolatedRestartAtTime=${tcheckpoint};/;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

    $file = "GrObservers.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|SphereA|JuggleA|m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

    $file = "GrDataBoxItems.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|#__FORPLUNGE__||m;
      $text    =~ s|__TDISRUPT__|${tcheckpoint}|m;
      $text    =~ s|#__FORPLUNGE__||m;
      $text    =~ s|__TDISRUPT__|${tcheckpoint}|m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

    $file = "WaveExtraction.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      my $re   = "DeltaT = __DTOBS__;";
      $text =~ s/^([^\#]*)(ObservationsPerOrbit\s*=\s*\d+;)\s*$/${1}$re/m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }
  }

  # Run DoMultipleRuns and switch to a new workdir.
  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub BhNsSettleDisk {
  my($oldworkdir,$oldscratch)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my ($workdir,$lev) = MakeSubdirWithSuffix($oldworkdir,"SettleDisk");

  #---------------------
  # Copy to new workdir
  #---------------------
  my %filestocopy;
  foreach my $file (glob "$oldworkdir/*.input") {
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}=$file;
  }
  foreach my $file (keys %filestocopy) {
    my $file2 = $file;
    if($file2 =~ s/SettleDisk$// || !exists($filestocopy{"${file}SettleDisk"})) {
      Utils::MyCopy($filestocopy{$file},"$workdir/${file2}.input");
    }
  }
  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  #-----------------------------------------
  # Read time of checkpoint from oldscratch
  #-----------------------------------------
  my $tcheckpoint=FindLatestCheckpointTimeFromInputDir($oldscratch);

  #-----------------------------------------
  # CD to workdir and modify input files
  #-----------------------------------------
  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  # Modify input files
  {
    my $file = "DoMultipleRuns.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|(^\s*\$PLUNGEDIR\s*=\s*)[^;]*;|${1}"${oldscratch}";|m;
      $text    =~ s|(^\s*\$TSETTLEDISK\s*=\s*)[^;]*;|${1}${tcheckpoint};|m;
      Utils::OverwriteFile($file,$text);
    }

    $file = "Evolution.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s/UseInterpolatedRestartAtTime\s*=\s*\d+.\d+;/UseInterpolatedRestartAtTime=${tcheckpoint};/;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }
  }

  # Run DoMultipleRuns and switch to a new workdir.
  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub NsNsDHGauge {
  my($oldworkdir,$oldscratch)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my ($workdir,$lev) = MakeSubdirWithSuffix($oldworkdir,"DHGauge");

  #------------------------------------------------------------
  # Get parent of old workdir
  # and files copied from parent (i.e. before DoMultipleRuns).
  #------------------------------------------------------------
  my $oldworkdirparent="$oldworkdir/..";
  Utils::MakePathsAbsolute(\$oldworkdirparent);
  my @fromparent=("AmrDriver.input");

  #---------------------
  # Copy to new workdir
  #---------------------
  my %filestocopy;
  foreach my $file (glob "$oldworkdir/*.input") {
    my $ignore=undef;
    foreach my $ff(@fromparent) {
      if($file =~ /$ff$/) {$ignore=1;last;}
    }
    next if($ignore);
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}=$file;
  }
  foreach my $file (@fromparent) {
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}="$oldworkdirparent/$file";
    die "Cannot find ",$filestocopy{$file2} unless (-f $filestocopy{$file2});
  }
  foreach my $file (keys %filestocopy) {
    my $file2 = $file;
    if($file2 =~ s/DHGauge$// || !exists($filestocopy{"${file}DHGauge"})) {
      Utils::MyCopy($filestocopy{$file},"$workdir/${file2}.input");
    }
  }
  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  #-----------------------------------------
  # Read time of checkpoint from oldscratch
  #-----------------------------------------
  my $tcheckpoint=FindLatestCheckpointTimeFromInputDir($oldscratch);

  #-----------------------------------------
  # CD to workdir and modify input files
  #-----------------------------------------
  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  # Modify input files
  {
    my $file = "DoMultipleRuns.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|(^\s*\$TGaugeChange\s*=\s*)[^;]*;|${1}${tcheckpoint};|m;
      $text    =~ s|(^\s*\$LastWorkDir\s*=\s*)[^;]*;|${1}"${oldscratch}";|m;
      Utils::OverwriteFile($file,$text);
    } else {
      die "$ME: Cannot find $file";
    }

  }

  # Run DoMultipleRuns and switch to a new workdir.
  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub NsNsMerger {
  my($oldworkdir,$oldscratch)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my ($workdir,$lev) = MakeSubdirWithSuffix($oldworkdir,"Merger");

  #---------------------
  # Copy to new workdir
  #---------------------
  my %filestocopy;
  foreach my $file (glob "$oldworkdir/*.input") {
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}=$file;
  }
  foreach my $file (keys %filestocopy) {
    my $file2 = $file;
    if($file2 =~ s/Merger$// || !exists($filestocopy{"${file}Merger"})) {
      Utils::MyCopy($filestocopy{$file},"$workdir/${file2}.input");
    }
  }
  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  #-----------------------------------------
  # Read time of checkpoint from oldscratch
  #-----------------------------------------
  my $tcheckpoint=FindLatestCheckpointTimeFromInputDir($oldscratch);

  #-----------------------------------------
  # CD to workdir and modify input files
  #-----------------------------------------
  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  # Modify DoMultipleRuns.input
  {
    my $file = "DoMultipleRuns.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      $text    =~ s|(^\s*\$LastWorkDir\s*=\s*)[^;]*;|${1}"${oldscratch}";|m;
      $text    =~ s|(^\s*\$TMergerBegins\s*=\s*)[^;]*;|${1}${tcheckpoint};|m;
      Utils::OverwriteFile($file,$text);
    }
  }

  # Run DoMultipleRuns and switch to a new workdir.
  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub BhNsMerger {
  my($oldworkdir,$oldscratch)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my $lev;
  my $workdir;
  my $oldworkdirbase;

  if($oldworkdir =~ m|(.*/)([^/]+)$|) {
    $oldworkdirbase = $1; # Part before and including last slash
    my $post        = $2; # Part after last slash
    if($oldworkdirbase =~ m|(.*/)([^/]+)/$|) {
      my $prepre      = $1;
      my $workdirbase = $2;

      if($workdirbase eq "Inspiral") {
        $workdir = "Merger";
      } else {
        $workdir = "$workdirbase" . "_" . "Merger";
      }
    } else {
      die "$ME: Please tell Mark about this error";
    }
    if($post =~ m|(.*)_([A-z]+)$|) {
      $lev      = $1;
      $lev      =~ s/^Lev(\d+)/$1/ || die "$ME: Cannot find lev in $lev";
    } else {
      $lev      = $post;
      $lev      =~ s/^Lev(\d+)/$1/ || die "$ME: Cannot find lev in $lev";
    }
  } else {
    die "$ME : Error with workdir '",$oldworkdir,"'";
  }
  if(-d $workdir) {
    die "$ME : Cannot submit to $workdir because it already exists";
  }
  mkdir $workdir || die "$ME: Cannot mkdir $workdir";

  unless(-d "$oldworkdirbase/MergerInputFiles") {
    die "$ME: Cannot find $oldworkdirbase/MergerInputFiles";
  }

  #---------------------
  # Copy to new workdir
  #---------------------
  foreach my $file (glob "$oldworkdirbase/MergerInputFiles/*.input") {
    my $file2=$file; $file2 =~ s|.*/||;
    Utils::MyCopy($file,"$workdir/$file2");
  }

  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  #------------------------------------------
  # Modify new input files
  #------------------------------------------

  # Insert $lev and initial data into DoMultipleRuns
  {
    my $file = "DoMultipleRuns.input";
    my $text = Utils::ReadFile($file);

    # Change Lev
    $text    =~ s|^\s*(\$LEV\s*=\s*)\S+;|${1}${lev};|m;

    # Change initial data location
    my $tstart = 0.0;
    {
      my $cpdir      = $oldscratch . "/" . "Checkpoints";
      my $biggestint = -1;
      {
        opendir(DIR,"$cpdir") || die "$ME: Cannot opendir $cpdir";
        my @dum = readdir(DIR);
        closedir(DIR)    || die "$ME: Cannot closedir $cpdir";
        foreach my $dir (@dum) {
          next unless ($dir =~ m/^\d+$/);
          if(-d $cpdir/$dir) {
            if($dir > $biggestint) {
              $biggestint = $dir;
            }
          }
        }
      }
      $tstart = FindLatestCheckpointTimeFromInputDir("$cpdir/$biggestint");
      $cpdir .= "/" . $biggestint . "/";
      $text   =~ s|^\s*(\$INSPIRALCPDIR\s*=\s*)\S+;|${1}"${cpdir}";|m;
    }

    # Directory with history files
    $text =~ s|^\s*(\$INSPIRALRUNDIR\s*=\s*)\S+;|${1}"${oldscratch}";|m;

    # Tstart
    $text =~ s|^\s*(\$TSTART\s*=\s*)\S+;|${1}"${tstart}";|m;

    # Excision boundary
    my $rmin = GetMinRadiusOfAh("$oldscratch/ApparentHorizon/GridAhA.dat");
    $text    =~ s|^\s*(\$RADIUSAH\s*=\s*)\S+;|${1}"${rmin}";|m;

    Utils::OverwriteFile($file,$text);
  }

  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}


#---------------
sub StartNsNsMergedBoxes {
  my($oldworkdir,$oldscratch)=@_;

  my $workdir = "$oldworkdir-Merger";
  if(-d $workdir) {
    die "$ME : Cannot submit to $workdir because it already exists";
  }

  my $NPROCS=undef;
  {
    my $text = Utils::ReadFile("$oldscratch/SpEC.out");
    while($text =~ m/MPI starting on (\d+) processors/g) {
      $NPROCS=$1;
    }
  }

  my $Checkpointdir=undef;
  {
    opendir(DIR,"$oldscratch/Checkpoints")
        || die "Cannot opendir $oldscratch/Checkpoints\n";
    my(@entry)=sort {$b <=> $a} (grep {$_ =~ m/^[0-9]+$/} readdir(DIR));
    closedir(DIR);
    $Checkpointdir="$oldscratch/Checkpoints/$entry[0]";
  }

  my $cwd=cwd();

  #### This was transcribed from old setup_merger.sh
  {
    mkdir($workdir) || die "Cannot mkdir $workdir";
    LinkSpec($oldworkdir,$workdir);
    Utils::MakePathsAbsolute(\$workdir);
    # Copy input files
    foreach my $file (glob "$oldworkdir/../*.input") {
      my $file2=$file; $file2 =~ s|.*/||;
      Utils::MyCopy($file,"$workdir/$file2");
    }
    # Swap inspiral and merger input files
    chdir($workdir) || die "$ME: Cannot chdir $workdir";
    foreach my $file (glob "*Merger.input") {
      my $file2=$file; $file2 =~ s|Merger\.input$||;
      if(-r "${file2}.input") {
        rename("${file2}.input","${file2}Inspiral.input") || die;
      }
      rename($file,"${file2}.input") || die;
    }
    Utils::System("$RealBin/DoMultipleRuns -n -B " .
                  "-f $RealBin/GetNsNsInspiralSubstitutions",
                  $opt_v);
    Utils::System("cd $Checkpointdir && $RealBin/RecoverOldCheckpoints",
                  $opt_v);
    Utils::System("$RealBin/GetNsNsMergedBoxFiles" .
                  " $oldscratch $Checkpointdir $workdir",
                  $opt_v);
  }

  $workdir = MergerCommon($workdir);

  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub StartRingdown {
  my($oldworkdir,$oldscratch,$extraopts,$ContinuationDirName,$usehydro,
     $NSubdomainsGR)=@_;

  #------------------------------------------------------
  # Create new workdir. Get the lev while we are at it.
  #------------------------------------------------------
  my ($workdir,$lev) = MakeSubdirWithSuffix($oldworkdir,"Ringdown");

  #---------------------
  # Copy to new workdir
  #---------------------
  my %filestocopy;
  foreach my $file (glob "$oldworkdir/*.input") {
    my $file2=$file; $file2 =~ s|.*/||; $file2 =~ s|\.input$||;
    $filestocopy{$file2}=$file;
  }
  foreach my $file (keys %filestocopy) {
    my $file2 = $file;
    if($file2 =~ s/Ringdown$// || !exists($filestocopy{"${file}Ringdown"})) {
      Utils::MyCopy($filestocopy{$file},"$workdir/${file2}.input");
    }
  }
  LinkSpec($oldworkdir,$workdir);
  Utils::MakePathsAbsolute(\$workdir);

  #---------------------------------------------------
  # Get glob pattern for last several segments.
  #---------------------------------------------------
  my $pat = GlobForDirsWithAhC($oldworkdir,$oldscratch,$ContinuationDirName);

  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  #------------------------------------------------------------------
  # Modify WaveExtraction.input by removing the ObservationsPerOrbit
  # and replacing it with output of every 0.1 (recall the total mass
  # is always one in our runs).
  #------------------------------------------------------------------
  {
    my $file = "GrWaveExtraction.input";
    if (-e $file) {
      my $text = Utils::ReadFile($file);
      my $re   = "DeltaT = 0.1;";
      if($text =~ s/^([^\#]*)(ObservationsPerOrbit\s*=\s*\d+;)\s*$/${1}$re/m) {
        Utils::OverwriteFile($file,$text);
      }
    } else {
      die "$ME: Cannot find $file";
    }
  }

  #---------------------------------------
  # Decide whether to use AMR in ringdown
  #---------------------------------------
  my $RingdownUsesAmr=undef;
  {
    my $text = Utils::ReadFile("Evolution.input");
    $text    =~ s/#.*$//mg; # Remove comments.
    if($text =~ m/AmrDriver\s*=/m) {
      $RingdownUsesAmr=1;
    }
  }

  #---------------------------------
  # Update the submission input file
  #---------------------------------
  {
    # Determine maximum number of cores
    my $MaxCores = undef;
    my $MinCores = undef;
    my $TargetCores = undef;
    my $NSubdomainsHY = 0;
    if ($usehydro) {
      if (not defined $NSubdomainsGR) {
        # This is the SingleStarHorizon case
	my $Old_NSubdomainsHY = SpEC::NSubdomains("$oldscratch/HyDomain.input");
        $NSubdomainsGR = SpEC::NSubdomains("$oldscratch/GrDomain.input");
        $NSubdomainsHY = $Old_NSubdomainsHY;
        $MaxCores = $NSubdomainsGR + $NSubdomainsHY;
	$MinCores = floor($NSubdomainsGR/2) + $NSubdomainsHY;
	$TargetCores = $MaxCores;
      } else {
        # This is the NsNsHorizon case
	# Currently, number of Sd post-merger is hard-coded to 960
	# i.e. 2x8^3-4^3 for two 8x8x8 nested grids.
        $NSubdomainsHY = 960;
	# Estimate number of gr procs
        $NSubdomainsGR = 70;
	$MaxCores = $NSubdomainsGR + $NSubdomainsHY;
	$MinCores = floor($NSubdomainsGR/2) + floor($NSubdomainsHY/6);
	$TargetCores = $NSubdomainsGR + floor($NSubdomainsHY/4);
      }
    } else {
      $MaxCores = $NSubdomainsGR;
      $MinCores = floor($MaxCores/1.5);
      $TargetCores = $MaxCores;
    }

    # Determine number of cores
    # Request all cores on a node if $MaxCores > $cpn
    my ($Cores, $cpn, $queue) = Machines::GetNodeSetup(
      Nmin       => $MinCores,
      Nmax       => $MaxCores,
      N          => $TargetCores,
    );
    if ($usehydro) {
      my $GrCores = $NSubdomainsGR;
      if (not defined $NSubdomainsGR) {
        $GrCores = $Cores - $NSubdomainsHY;
      }
      $extraopts .= " --ngrprocs=$GrCores --ntotalprocs=$Cores --autoncyl=0";
    }

    # Determine Jobname
    my $Jobname = Utils::SystemOutput("$RealBin/MakeSubmit.py query Jobname");
    chomp($Jobname);
    $Jobname = "RD_$Jobname";

    # Update file
    if ($queue ne "SpECNoSuchQueue") {
        Utils::System("$RealBin/MakeSubmit.py update -f --Jobname $Jobname ".
                      "--CoresPerNode $cpn --Cores $Cores --Queue $queue");
    } else {
        Utils::System("$RealBin/MakeSubmit.py update -f --Jobname $Jobname ".
                      "--CoresPerNode $cpn --Cores $Cores");
    }
  }

  #-----------------------------------
  # Run ConstructPostMergerDataFromAhC
  #-----------------------------------
  my $command   =
      "$RealBin/ConstructPostMergerDataFromAhC " .
      "$extraopts " .
      "--contdirname \"$ContinuationDirName\" " .
      "--inputdirs  \"$pat\" " .
      "--inputdir   $oldscratch " .
      "--vout 0.01 --autovout " .
      "--v --rminfac 0.94 " .
      "--stencil 3 " .
      "--rwidfac 1.5 --rwidth 20.0 --derivorderfvt 2 " .
      "--expdecay 1 --rotdecay 20 " .
      ">DoContinuation.out 2>&1";
  Utils::OverwriteFile("DoContinuation",$command);
  Utils::System($command, $opt_v);

  #---------------------------------------------------------
  # Copy output of ConstructPostMergerDataFromAhC to workdir
  #---------------------------------------------------------
  foreach my $file (glob "Continuation/*.input") {
    my $file2 = $file; $file2 =~ s|.*/||;
    if($usehydro) {
      $file2 =~ s/^StateChangers.input/GrStateChangers.input/;
      $file2 =~ s/^Domain.input/GrDomain.input/;
    }
    Utils::MyCopy("$file","$file2");
  }

  #---------------------------------------------------------
  # The DoMultipleRuns.input here is the one from Continuation
  #---------------------------------------------------------
  {
    my $file = "DoMultipleRuns.input";
    my $text = Utils::ReadFile($file);

    # Change the levs to the same as inspiral
    $text =~ s/(\$M(in|ax)Lev\s*=\s*)[^;]*;/$1$lev;/gm;

    # Put ringdown resolution into new run
    if($RingdownUsesAmr) {
      # Get truncation error from the Plunge.
      if(-f "$oldscratch/AmrDriver.input") {
        my $read = Utils::ReadFile("$oldscratch/AmrDriver.input");
        $read    =~ s/#.*$//mg; # Remove comments
        # Remove 'RestrictionsOnExtents=' if it exists.
        # This is because it may include a 'TruncationErrorMax'.
        our $bal; $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;
        $read    =~ s/RestrictionsOnExtents\s*=\s*([^(]+${bal}[,;])+//mg;
        if($read =~ m|^\s*TruncationErrorMax\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?);|m) {
          # Here TruncationErrorMax is specified as just a number.
          my $terrmax = $1;
          $text =~ s/^([^#]*\$TruncationErrorMax\s*=\s*).*$/${1}$terrmax;/m;
        } elsif($read =~
                m|^\s*TruncationErrorMax\s*=\s*((?:Add${bal}[,;]\s*)+)|ms){
          # Here it is a list of 'Add(...)'s.
          my $text2   = $1;
          my $terrmax = undef;
          # Find the 'Add(...)' with a center of zero.
          while($text2 =~ m|Add($bal)|mg) { # Loop over each 'Add...'
            my $text3 = $1;
            if($text3 =~ m|Center\s*=\s*(\S+);|m) {
              # Check if the center is zero
              my @dum = split(',',$1);
              my $isnonzero = undef;
              foreach my $dum (@dum) {
                if($dum != 0.0) {$isnonzero=1;}
              }
              # If center is zero, set the truncation error, and don't
              # keep looking.
              if(!$isnonzero) {
                if($text3 =~ m|Value\s*=\s*(\S+);|m) {
                  $terrmax = $1;
                  # The value may not be a number if we are using
                  # AmrTolerances.input.  So check for this.
                  if($terrmax =~ m|ReadFromAmrTolerancesInputFile\((\S+)\)|) {
                    # Now we need to get a key from AmrTolerances.input.
                    my $key = $1;
                    my $amr_tol_text
                        = Utils::ReadFile("$oldscratch/AmrTolerances.input");
                    if($amr_tol_text =~ m|$key\s*=\s*(\S+);|m) {
                      $terrmax = $1;
                    } else {
                      die "$ME: Cannot find $key in "
                          . "$oldscratch/AmrTolerances.input\n";
                    }
                  }
                  last;
                }
              }
            }
          }
          die "$ME: Cannot find TruncationErrorMax in "
              . "$oldscratch/AmrDriver.input\n" unless defined $terrmax;
          $text    =~
              s/^([^#]*\$TruncationErrorMax\s*=\s*).*$/${1}$terrmax;/m;
        } else {
          die "$ME: Cannot find TruncationErrorMax in "
              . "$oldscratch/AmrDriver.input\n";
        }
      } else {
        # OK, merger didn't use AMR. So set truncation error to something
        # 'sensible'.
        my $terrmax = 0.001;
        warn "$ME: Cannot find TruncationErrorMax in "
            . "$oldscratch/AmrDriver.input. Setting it to $terrmax\n";
        $text =~ s/^([^#]*\$TruncationErrorMax\s*=\s*).*$/${1}$terrmax;/m;
      }
    }

    Utils::OverwriteFile($file,$text);
  }

  $workdir = MergerCommon($workdir);

  #---------------------------------------------------------
  # Clean up and exit
  #---------------------------------------------------------
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return $workdir;
}

#---------------
sub SplitOrMergeSubdomains {
  my($oldworkdir,$oldscratch,$InputFiles)=@_;

  my $workdir = GetNewSubmissionDirectoryAndLinkSpEC($oldworkdir);

  # Cd to workdir, saving cwd.
  my $cwd=cwd();
  chdir($workdir) || die "$ME: Cannot chdir $workdir";

  # Check that Domain.input has an entry
  # ReadFromFiles=RefinementOptionHistory.input;
  {
    my $file = basename(SpEC::GetGrDomainInputPath($oldscratch));
    my $text = $InputFiles->{$file};
    $text =~ s|#.*$||gm; # Remove comments.
    unless($text =~ m|ReadFromFiles[^;]*RefinementOptionHistory\.input|ms) {
      die "File $file does not contain appropriate ReadFromFiles option";
    }
    # We don't care to change Domain.input here, so don't rewrite
    # $InputFiles->{$file};
  }

  # Append AppendToRefinementOptionHistory.output;
  if(-f "$oldscratch/AppendToRefinementOptionHistory.output") {
    my $text2=Utils::ReadFile("$oldscratch/AppendToRefinementOptionHistory.output");
    $text2 = "\n$text2"; # Want to start on a new line.
    AddCommaSeparatedEntriesToFile("RefinementOptionHistory.input",$InputFiles,
                                   $text2);
  }

  # Change refinement options for JuggleBalls
  if(-f "$oldscratch/NewAdaptiveJuggleBallsTree.input"){
    $InputFiles->{"AdaptiveJuggleBallsTree.input"} =
        Utils::ReadFile("$oldscratch/NewAdaptiveJuggleBallsTree.input");
  }

  my $change_restart = 1;

  # RemoveOwnerPolicyFromDomainInput($workdir);
  AddOwnerPolicyOutputToDomainInput($oldscratch,$InputFiles);

  # Cd back.
  chdir($cwd) || die "$ME: Cannot chdir $cwd";

  return ($workdir, $change_restart, $InputFiles);
}

#---------------
# We do not set $change_restart here to avoid overwriting a prior value
sub WallClockContinuation {
  my($oldwork,$oldscratch,$InputFiles)=@_;

  my $work = GetNewSubmissionDirectoryAndLinkSpEC($oldwork);

  AddOwnerPolicyOutputToDomainInput($oldscratch,$InputFiles);

  return $work;
}

#---------------
sub DropInnerShells {
  my($oldwork,$oldscratch,$InputFiles)=@_;

  my $work = GetNewSubmissionDirectoryAndLinkSpEC($oldwork);

  # Cd to workdir, saving cwd.
  my $cwd=cwd();
  chdir($work) || die "$ME: Cannot chdir $work";

  my @sdnames = SubdomainNamesToDrop("$oldscratch/SpEC.out");
  RemoveSubdomainsFromDomainInput($InputFiles,@sdnames);

  # Comment out any per-subdomain resolution changes.
  foreach my $file (keys %$InputFiles) {
    if($file eq "Domain.input") {
      CommentOutLinesBeginningWithSubdomainNames($file,
                                                 $InputFiles,@sdnames);
      AddCommaSeparatedEntriesToFile($file,$InputFiles,undef);
    }
    if ($file eq "DataBoxItems.input" or $file =~ /..DataBoxItems.input/) {
      CommentOutLinesBeginningWithSubdomainNames($file,$InputFiles,@sdnames);
    }
  }

  my $change_restart = 0;
  AddOwnerPolicyOutputToDomainInput($oldscratch,$InputFiles);

  # Cd back.
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return ($work, $change_restart);
}

#---------------
sub DropJunkShell {
  my($oldwork,$oldscratch,$InputFiles)=@_;

  my $work = GetNewSubmissionDirectoryAndLinkSpEC($oldwork);

  # Cd to workdir, saving cwd.
  my $cwd=cwd();
  chdir($work) || die "$ME: Cannot chdir $work";
  my $change_restart = 0;
  AddOwnerPolicyOutputToDomainInput($oldscratch,$InputFiles);

  # Find $RmaxAfterDrop
  my $DomainFile = "SpatialCoordMap.input";
  my $text = $InputFiles->{$DomainFile};
  my $RmaxAfterDrop;
  if ($text =~ /OuterBdryRadius *= *(\d+\.*\d*);/){
    $RmaxAfterDrop = $1;
  }
  else {
    die "SpatialCoordMap.input was improperly formatted for DropJunkShell.";
  }

  # Find the number of shells to drop, and the list of original subdomain radii
  $DomainFile = "GrDomain.input";
  $text = $InputFiles->{$DomainFile};
  my $NumRadiiToDrop = 0;
  my @OldShellRadii;
  if ($text =~ /SphereC *= *\(Radii *= *([\d+\.*\d*\,*]+)/) {
    @OldShellRadii = split(',', $1);
    for my $R (@OldShellRadii) {
      # We set a pretty generous tolerance here because we are reading in and
      # comparing values from different input files, which might be goverened
      # differently regarding rounding. This is fine because the  the SphereC
      # shells all have the same width, which is FAR larger than 1e-5.
      if ($R >= $RmaxAfterDrop + 1.0e-5) {
        $NumRadiiToDrop += 1;
      }
    }
  }
  else {
    die "GrDomain.input was improperly formatted for DropJunkShell.";
  }
  # Get the list of new subdomain radii
  my @NewShellRadii = @OldShellRadii;
  @NewShellRadii = splice @NewShellRadii, 0, -$NumRadiiToDrop;
  # Data for rewriting GrDomain.input
  my $OldShellRadiiText = join ',', @OldShellRadii;
  my $NewShellRadiiText = join ',', @NewShellRadii;
  my $NumRadii = scalar @NewShellRadii;
  $NumRadii -= 1;
  # Rewrite GrDomain.input
  for (my $i=$NumRadii; $i<=$NumRadii+$NumRadiiToDrop; $i++) {
    $text  =~ s|SphereC$i|# SphereC$i|m;
  }
  $text    =~ s|$OldShellRadiiText|$NewShellRadiiText|m;
  $InputFiles->{$DomainFile} = $text;

  # Rewrite Evolution.input
  $DomainFile = "Evolution.input";
  $text = $InputFiles->{$DomainFile};
  $text =~ s|DropJunkShellTime\(|# DropJunkShellTime\(|m;
  $InputFiles->{$DomainFile} = $text;

  # Rewrite GrObservers.input
  $DomainFile = "GrObservers.input";
  $text = $InputFiles->{$DomainFile};
  $text =~ s|Rmax *= *$OldShellRadii[-1]|Rmax = $RmaxAfterDrop|m;
  $InputFiles->{$DomainFile} = $text;

  # Cd back.
  chdir($cwd) || die "$ME: Cannot chdir $cwd";
  return ($work, $change_restart);
}

#---------------
sub EvolveGeodesicsNewLev {
  my($OldWorkDir,$InputFiles)=@_;
  (my $OldLevDir = $OldWorkDir) =~ s|/+[^/]+/*$||;

  my $NewLevDir = SpEC::GetNewSubmissionDirectory($OldLevDir);
  my $NewWorkDir = "$NewLevDir/Run_AA";

  mkpath($NewWorkDir) || die "$ME: Cannot mkdir $NewWorkDir";

#Copy the restart information
  Utils::MyCopy("$OldWorkDir/Run/Restart.output", "$NewWorkDir/Restart.input");
  delete($InputFiles->{"Restart.input"});
  mkpath("$NewWorkDir/RestartInput") ||
    die "$ME: Cannot mkdir $NewWorkDir/RestartInput";
  foreach my $file (glob "$OldWorkDir/Run/RestartOutput/*") {
    my $file2=$file;  $file2 =~ s|.*/||;
    Utils::MyCopy($file,"$NewWorkDir/RestartInput/$file2");
  }

  LinkSpec($OldWorkDir, $NewWorkDir);

  return ($NewWorkDir);
}

#---------------
sub EvolveGeodesicsWallClock {
  my($OldWorkDir,$InputFiles)=@_;
  my $NewWorkDir = GetNewSubmissionDirectoryAndLinkSpEC($OldWorkDir);

  $InputFiles->{"Restart.input"} =
      Utils::ReadFile("$OldWorkDir/Run/Restart.output");
  mkpath("$NewWorkDir/RestartInput") ||
    die "$ME: Cannot mkdir $NewWorkDir/RestartInput";
  foreach my $file (glob "$OldWorkDir/Run/RestartOutput/*") {
    my $file2=$file;  $file2 =~ s|.*/||;
    Utils::MyCopy($file,"$NewWorkDir/RestartInput/$file2");
  }

  return $NewWorkDir;
}

#=========================
# Medium-level Subroutines
#=========================

# Common subfunction for all merger-like routines
sub MergerCommon {
  my ($workdir) = @_;

  Utils::System("$RealBin/DoMultipleRuns -n -L -b ./bin", $opt_v);
  {
    # Change workdir to new Lev* run
    my @levs = (glob "Lev*");
    die "$ME: What to use for Lev? @levs" if(scalar(@levs)!=1);
    $workdir .= "/" . $levs[0];
  }
  return $workdir;
}

#---------------
# Given a file of the form
#   Bla =
#     entry,
#     entry,
#     ....
#     lastentry;
#
# and text of the form
#     newentry,
#     newentry,
#     ...
#     lastnewentry;
#
# Concatenates the new entries to the end of the old entries,
# and overwrites the file with the list of concatenated entries,
# making sure that each line but the last ends in a comma, and the
# last line ends in a semicolon.  If there are commented lines,
# ignore them for the purpose of semicolon/comma.
#
# If 'text' is undefined, simply cleans up the given file by
# making sure that every noncommented line except the last ends in
# a comma, and the last noncommented line ends in a semicolon.
sub AddCommaSeparatedEntriesToFile {
  my($file,$InputFiles,$text)=@_;

  # Read the file
  my $text1 = $InputFiles->{$file};
  # Remove blank lines
  $text1 =~ s/^\s*$//m;
  # Remove the semicolon if there are no options
  $text1 =~ s/(\S+\s*=\s*);/$1/;
  # Replace trailing semicolon AT END OF FILE by comma.
  $text1 =~ s/;\s*$/,/;
  $text1 =~ s/;((\n#[^\n]*)+)$/,$1/; # and on last non-commented line
  if(defined $text) {
    # I have something to add.
    # Replace trailing comma at end of extra text by semicolon.
    $text  =~ s/,\s*$/;/;
    $InputFiles->{$file} = $text1.$text."\n";
  } else {
    # I'm adding nothing, just cleaning commas/semicolons.
    # Replace trailing comma AT END OF last non-commented line by semicolon
    $text1 =~ s/,((\n#[^\n]*)+)$/;$1/;
    # Replace trailing comma AT END OF FILE by semicolon.
    $text1 =~ s/,\s*$/;/;
    $InputFiles->{$file} = $text1 . "\n";
  }
}

#---------------
# Returns a list of lines that were changed.
sub CommentOutLinesBeginningWithSubdomainNames {
  my($file,$InputFiles,@names)=@_;

  my $output;
  my @lines = split(/\n/,$InputFiles->{$file});
  foreach my $line (@lines) {
    foreach my $name(@names) {
      if($line =~ /^\s*$name([xyzXYZ][mp].+)?\s*\(/) {
        $line = "#" . $line;
        last; # Don't need to check other names for this line
      }
    }
    $output .= $line;
  }
  $InputFiles->{$file} = $output;
}

#---------------
# Reads $specoutfile, looks for 'only *going Charspeeds' string, and
# returns names of all subdomains that should be dropped.
sub SubdomainNamesToDrop {
  my($specoutfile)=@_;

  my %names;
  my @text = Utils::ReadFile($specoutfile);
  foreach (@text) {
    if(/^Slice\s+(SliceUFF\.\S+)\s+has only ingoing CharSpeeds\.$/) {
      my $name =  $1;
      $name    =~ s/SliceUFF\.//;#Remove 'Slice' stuff at the beginning.
      # The name that triggers termination is the same one dropped.
      if($name =~ m/(\d+)([XYZ][mp].+)?$/) {
        $name  = $` . $1;
      }
      $names{$name}=1; # Use hash instead of array for uniqueness.
    }
    if(/^Slice\s+(SliceLFF\.\S+)\s+has only outgoing CharSpeeds\.$/) {
      my $name =  $1;
      $name    =~ s/SliceLFF\.//;#Remove 'Slice' stuff at the beginning.
      # The name that triggers termination is > than the one dropped.
      if($name =~ m/(\d+)([XYZ][mp].+)?$/) {
        my $pre=$`;
        my $num=$1;
        --$num;
        $name = $pre . $num;
      }
      $names{$name}=1; # Use hash instead of array for uniqueness.
    }
  }

  return keys(%names);
}

#---------------
sub RemoveSubdomainsFromRefOptHist {
  my($file,$InputFiles,@names)=@_;

  my $text = $InputFiles->{$file};
  # $text =~ s/;\s*$//; # Remove ';' at end.

  # Remove the semicolon if there are no options
  $text =~ s/(RefinementOptionHistory\s*=\s*);/$1/;
  # Replace trailing semicolon AT END OF FILE by comma.
  $text =~ s/;\s*$/,/;
  $text =~ s/;((\n#[^\n]*)+)$/,$1/; # and on last non-commented line

  foreach my $name(@names) {
    $text .= "\n" . $name . "(Remove=yes),";
  }
  # Replace trailing comma AT END OF FILE by semicolon.
  $text =~ s/,\s*$/;/;

  $InputFiles->{$file} = $text;
}


#---------------
sub RemoveSubdomainsFromDomainInput {
  my($InputFiles,@names)=@_;

  # Read old Domain.input
  my $changed     = undef;
  my $filesfound  = 0;
  my $refophist   = undef;
  my $domainfile  = undef;
  foreach my $file ("Domain.input","GrDomain.input") {
    next unless exists($InputFiles->{$file});
    ++$filesfound;
    $domainfile = $file;
    my $olddomaininput = $InputFiles->{$file};

    # Is there a RefinementOptionHistory.input?
    {
      $olddomaininput    =~ s|\s+||g; # Remove whitespace, including newlines.
      if($olddomaininput =~ m|ReadFromFiles=([^;]+);|) {
        my @files=split(/,/,$1);
        foreach my $f (@files) {
          if($f =~ m|RefinementOptionHistory|) {
            $refophist=$f;
          }
        }
      }
    }

    my $newdomaininput=
        ChangeMaskStringInDomainInputWork('$sdlabel',\@names,
                                          split(/\n/,$olddomaininput));
    unless(defined $newdomaininput) {
      $newdomaininput=
          ChangeMaskStringInDomainInputWork('$sdbasename',\@names,
                                            split(/\n/,$olddomaininput));
    }
    unless(defined $newdomaininput) {
      $newdomaininput=
          ChangeMaskStringInDomainInputWork('$sdname',\@names,
                                            split(/\n/,$olddomaininput));
    }

    if(defined $newdomaininput) {
      $changed = 1;
      $InputFiles->{$file} = $newdomaininput;
    }
  }

  die ("$ME: Should find either Domain.input or GrDomain.input")
      unless (1==$filesfound);

  unless($changed) {
    unless(defined $refophist) {
      # There's no RefinementOptionHistory, so make one.
      $refophist =  $domainfile;
      $refophist =~ s/Domain/RefinementOptionHistory/;
      my $text   = $InputFiles->{$domainfile};
      $text     .= "\nReadFromFiles=$refophist;\n";
      $InputFiles->{$domainfile} = $text;
      $InputFiles->{$refophist}  = "RefinementOptionHistory=;";
    }
    RemoveSubdomainsFromRefOptHist($refophist,$InputFiles,@names);
  }
}

#---------------
# Returns string to print into the new Domain.input file.
# Returns undef if it does not find a mask string to change.
sub ChangeMaskStringInDomainInputWork {
  my($suffix,$spherenames,@lines)=@_;
  local $_;

  my $outputstring="";
  my $changed=undef;
  foreach (@lines) {
    foreach my $sdname (@$spherenames) {
      my $sdbasename=$sdname;
      $sdbasename   =~s/\d.*//;# Take off the first number and beyond.
      my $sdlabel   = substr($sdbasename,-1);

      # The following line is a perl trick to evaluate Mask$suffix TWICE
      # and put it into $pat.
      my $pat; eval "\$pat = \"Mask$suffix\"";
      if(/^(\s*${pat}\s*=\s*)(.*);/) {
        my $var  = $1;
        my $list = $2;
        my @list = split(/\s*,\s*/,$list);
        my $outval=0;
        my $comma="";
        foreach my $val (@list) {
          $var  .= $comma . $outval;
          $comma = ",";
          $outval=$val; # This increases the number of zeros in the list.
        }
        print STDERR "For $suffix, $sdname , $sdlabel , $sdbasename, $pat:\n"
            . "Changing $_ to " if($opt_v);
        $_ = $var . ";\n";
        print STDERR $_ if($opt_v);
        $changed=1;
        last; # Don't need to check other spherenames for this line
      }
    }
    $outputstring .= $_;
  }

  if($changed) {
    return $outputstring;
  } else {
    return undef;
  }
}

#---------------
sub RemoveOwnerPolicyFromDomainInput {
  my($newworkdir)=@_;
  my $domainfile = SpEC::GetGrDomainInputPath($newworkdir);
  my $text       = Utils::ReadFile($domainfile);
  ### $bal matches balanced parens with possible nested parens inside.
  our $bal;  $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;
  $text=~s|OwnerPolicy\s*=\s*WeightBalanced\s*${bal};||g;
  $text=~s|OwnerPolicy\s*=\s*LoadBalanced\s*${bal};||g;
  $text=~s|OwnerPolicy\s*=\s*BalancedAndDistributed\s*${bal};||g;

  Utils::OverwriteFile($domainfile,$text);
}

#---------------
# If an $oldscratchdir/NewOwnerPolicy.output exists,
# then incorporate it into Domain.input in $newworkdir.
sub AddOwnerPolicyOutputToDomainInput {
  my($oldscratchdir,$InputFiles)=@_;

  my $file = "$oldscratchdir/NewOwnerPolicy.output";
  return 0 unless(-f $file);

  ### $bal matches balanced parens with possible nested parens inside.
  our $bal;  $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;
  my $reg = qr/^\s*OwnerPolicy\s*=\s*([^();]+$bal?\s*;)\s*$/m;

  my $otext  = Utils::ReadFile($file);
  if($otext =~ m/$reg/m) {
    my $replacement = $1;

    my $dfile = basename(SpEC::GetGrDomainInputPath($oldscratchdir));
    my $dtext = $InputFiles->{$dfile};

    if($dtext =~ m/^\s*OwnerPolicy\s*=\s*/m) {
      my $re = qr/(^\s*OwnerPolicy\s*=\s*)[^();]+$bal?\s*;\s*$/m;
      unless($dtext =~ s/$re/${1}$replacement/m) {
        die "$ME: Cannot set OwnerPolicy in $file from $dfile";
      }
    } else {
      $dtext = "OwnerPolicy=$replacement\n" . $dtext;
    }
    $InputFiles->{$dfile} = $dtext;
  }
}

#---------------
# If an $oldscratchdir/NextHyDomain.output exists,
# then use it to replace HyDomain.input
sub UpdateHyDomainInput {
  my($oldscratchdir,$InputFiles)=@_;

  my $file = "$oldscratchdir/NextHyDomain.output";
  return 0 unless(-f $file);

  # LastHyDomain.input is the source domain for
  # the interpolated restart
  $InputFiles->{"LastHyDomain.input"} = $InputFiles->{"HyDomain.input"};
  my $nexttext  = Utils::ReadFile($file);
  $InputFiles->{"HyDomain.input"} = $nexttext;
}

#---------------
# Automated rescaling of the hydro grid
# for FMR runs
sub ExpandHydroGrid {
  my($oldscratchdir,$newworkdir,$InputFiles)=@_;

  # LastHyDomain.input is the source domain for
  # the interpolated restart
  $InputFiles->{"LastHyDomain.input"} = $InputFiles->{"HyDomain.input"};

  my $cptime;
  my $cpdir;
  ($cptime,$cpdir) = SpEC::FindLatestCheckpointTimeFromInputDir($oldscratchdir);

  my $dmfile = $InputFiles->{"Ev_Params.input"};
  my $evfile = $InputFiles->{"Evolution.input"};
  my $fmrfile = $InputFiles->{"FMRItems.input"};
  return 0 unless ($dmfile =~ m/FmrRegridScale\s*=\s*(\d+.\d+)/);
  my $FmrScale = $1;
  return 0 unless ($evfile =~ m/MonitorTargetToExpandGrid\s*=\s*(\d+.\d+);/);
  my $NewLimit = $1*$FmrScale;
  # This gives us an updated version of Evolution.input
  $evfile =~ s/MonitorTargetToExpandGrid\s*=\s*\d+.\d+;/MonitorTargetToExpandGrid=${NewLimit};/;
  return 0 unless ($evfile =~ m/UseInterpolatedRestartAtTime\s*=\s*\d+.\d+;/);
  $evfile =~ s/UseInterpolatedRestartAtTime\s*=\s*\d+.\d+;/UseInterpolatedRestartAtTime=${cptime};/;
  return 0 unless ($fmrfile =~ m/HalfLengthFinestLevel\s*=\s*(\d+.*\d*),(\d+.*\d*),(\d+.*\d*);/);
  my $newLx = $1/$FmrScale;
  my $newLy = $2/$FmrScale;
  my $newLz = $3/$FmrScale;
  # This gives an updated version of FMRItems.input
  $fmrfile =~ s/HalfLengthFinestLevel\s*=\s*(\d+.*\d*),(\d+.*\d*),(\d+.*\d*);/HalfLengthFinestLevel=$newLx,$newLy,$newLz;/;

  # Now, we need to copy the new files in the new work directory
  $InputFiles->{"Evolution.input"} = $evfile;
  $InputFiles->{"FMRItems.input"} = $fmrfile;

  # Finally, we get the new HyDomain.input
  my $file = "$oldscratchdir/NextHyDomain.output";
  return 0 unless(-f $file);
  $InputFiles->{"HyDomain.input"} = Utils::ReadFile($file);
}

#---------------
# Make sure that the initial number of initial amr iterations in
# Evolution.input is nonzero.  Otherwise, Domain.input will start from
# the checkpoint file and will not notice the changes made here.
sub AddOneAmrIterationIfCurrentlyZero {
  my ($InputFiles, $file) = @_;

  my $text    = $InputFiles->{$file};
  my $replace = 2; # >0 if we replace the text, >1 if we write the text.
  if($text =~ s/^(\s*AmrDriverMaxInitIterations\s*=\s*)(\d+)\s*;/${1}1;/ms){
    $replace=1;
    $replace=0 if($2 > 0);
  }
  if($replace == 2) {
    $text = "AmrDriverMaxInitIterations = 1;\n" . $text;
  }
  if($replace > 0) {
    $InputFiles->{$file} = $text;
  }
}

#---------------
# Return the number of processors used in a Hydro run
sub GetHydroProcs {
  my ($oldscratch) = @_;

  # Figure out how many hydro processors were used in the pre-ringdown run.
  my $PlungeProcessors=undef;
  my $PlungeHyProcessors=undef;
  my $PlungeGrProcessors=undef;
  my $file = "$oldscratch/SpEC.out";
  if (-e $file) {
    my $text = Utils::ReadFile($file);
    if($text =~ m/^MPI starting on (\d+) processors at/m) {
      $PlungeProcessors=$1;
    }
    if($text =~ m|Gr:\s*ranks\s*\(([^)]+)\)\s*Hy:\s*ranks\s*\(([^)]+)\)|s) {
      my $grproc=$1;
      my $hyproc=$2;
      if($grproc =~ m|(\d+)-(\d+)|) {
        $PlungeGrProcessors = $2 - $1 + 1;
      } elsif($grproc =~ /,/) {
        my @a=split(/,/,$grproc);
        $PlungeGrProcessors = scalar(@a);
      }
      if($hyproc =~ m|(\d+)-(\d+)|) {
        $PlungeHyProcessors = $2 - $1 + 1;
      } elsif($hyproc =~ /,/) {
        my @a=split(/,/,$hyproc);
        $PlungeHyProcessors = scalar(@a);
      }
    }
  }
  die "$ME: Cannot find num processors" unless(defined $PlungeProcessors);
  return ($PlungeProcessors,$PlungeHyProcessors);
}

#---------------
sub GlobForDirsWithAhC {
  my ($oldworkdir,$oldscratch,$ContinuationDirName) = @_;

  # Choose $oldworkdir/Run instead of $oldscratch because the latter may
  # be on some scratch directory and may not be as easy to find a glob
  # pattern for.  As long as $oldworkdir/Run is symlinked to $oldscratch
  # this should be fine.
  my $pat = "$oldworkdir/Run";
  unless(Cwd::abs_path($pat) eq Cwd::abs_path($oldscratch)) {
    die "$ME: These do not have the same absolute path: $oldscratch $pat";
  }
  unless($pat =~ s|/Run$||) {
    die "$ME: Cannot find '/Run' at the end of $pat";
  }
  if($pat =~ s|^(.*_[A-z]*)[A-Z][A-Z]$|${1}*|) {
    my $matchstring = quotemeta($1);    # used in regexp
    # Restrict the glob: find the first dir from the merger that contains
    # ForContinuation/AhC.dat and include all dirs from that to the end.
    # This will change a glob like Lev2_*/Run into e.g. Lev2_{BB,BC,BD}/Run
    # This is to prevent copying large amounts of unnecessary Hist- files
    # when running ConstructPostMergerDataFromAhC
    my @paths = glob("$pat/Run");
    while (@paths>0 and not -e "$paths[0]/$ContinuationDirName/AhC.dat") {
      shift @paths;
    }
    unless (@paths>0) {
      die "$ME: didnt find any $ContinuationDirName/AhC.dat";
    }
    my $tags = "";
    foreach my $path (@paths) {
      if($path =~ m|$matchstring(.*)/Run|) {
        $tags .= $1.",";
      } else {
        die "$ME: path $path does not match expected pattern";
      }
    }
    $tags =~ s/,$//;
    if(@paths > 1) {
      $pat =~ s/\*/\{$tags\}/;
    } else {
      $pat =~ s/\*/$tags/;
    }
  }
  $pat .= "/Run";
  return $pat;
}

#======================
# Low-level Subroutines
#======================

sub LinkSpec {
  my($src,$dest)=@_;

  # Link the evolution executable
  Utils::LinkToOriginal("$src/SpEC", "$dest/SpEC");
  print STDERR "$ME: SpEC linked to $src/SpEC\n" if($opt_v);

  # Link the bin directory, if it exists
  if(-d "$src/bin") {
    Utils::LinkToOriginal("$src/bin", "$dest/bin");
    print STDERR "$ME: bin directory linked to $src/bin\n" if($opt_v);
  }
}

#---------------
sub WriteInputFiles {
  my($work,$InputFiles)=@_;

  foreach my $file (keys %$InputFiles) {
    Utils::OverwriteFile("$work/$file",$InputFiles->{$file});
  }
}

#---------------
sub ReadInputFiles {
  my($oldworkdir)=@_;
  my $InputFiles;
  foreach my $file (glob "$oldworkdir/*.input") {
    $InputFiles->{basename($file)} = Utils::ReadFile($file);
  }
  return $InputFiles;
}

#---------------
sub GetNewSubmissionDirectoryAndLinkSpEC {
  my($oldworkdir)=@_;
  my $workdir = SpEC::GetNewSubmissionDirectory($oldworkdir);

  if(-d $workdir) {
    die "$ME : Cannot submit to $workdir because it already exists";
  }
  mkdir $workdir || die "$ME: Cannot mkdir $workdir";

  LinkSpec($oldworkdir,$workdir);

  Utils::MakePathsAbsolute(\$workdir);
  return $workdir;
}

#---------------
# Make a subdirectory for major evolution transitions,
# e.g. Ringdown, Plunge, etc.
sub MakeSubdirWithSuffix {
  my ($oldworkdir,$suffix) = @_;
  my $lev;
  my $workdir;
  if($oldworkdir =~ m|(.*/)([^/]+)$|) {
    $workdir = $1;
    my $foo  = $2;
    if($foo =~ m|(.*)_([A-z]+)$|) {
      $lev      = $1;
      $workdir .= $lev . "_$suffix";
      $lev      =~ s/^Lev(-?\d+)/$1/ || die "$ME: Cannot find lev in $lev";
    } else {
      $lev      = $foo;
      $workdir .= $foo . "_$suffix";
      $lev      =~ s/^Lev(-?\d+)/$1/ || die "$ME: Cannot find lev in $lev";
    }
  } else {
    die "$ME : Error parsing old workdir '$oldworkdir'";
  }
  if(-d $workdir) {
    die "$ME : Cannot submit to $workdir because it already exists";
  }
  mkdir($workdir) || die "$ME: Cannot mkdir $workdir";
  return ($workdir,$lev);
}

#---------------
sub GetMinRadiusOfAh {
  my($file)=@_;

  my $rmin = undef;
  my $col  = undef;
  my @text = Utils::ReadFile($file);
  foreach (@text) {
    if(/^\s*#\s*\[(\d+)\]\s*=\s*min\(r\)/) {
      # Find column with min(r)
      $col = $1;
    } elsif(/^\s*#/) {
      # Skip other comments
    } elsif(defined $col) {
      s/^\s*//; # Remove leading whitespace
      my @a  = split;
      $rmin = $a[$col-1];
    }
  }

  die "Cannot find rmin in $file" unless(defined $rmin);
  return $rmin;
}

# Finds filepattern recursively starting in dir.
# Returns full paths (starting with "dir") matching filepattern.
# filepattern should not contain "/" characters.
sub Find {
  my($filepattern,$dir)=@_;

  die "Filepattern should not contain '/' characters" if($filepattern =~m|/|);

  opendir(DIR,"$dir") || die "Cannot opendir $dir\n";
  my(@entry)=readdir(DIR);
  closedir(DIR);

  my @result;
  foreach my $f (@entry) {
    next if($f eq '.' || $f eq '..');
    my $fullpath = "$dir/$f";
    if(-d $fullpath) { # Skip links
      push(@result,Find($filepattern,$fullpath));
    } elsif(-f $fullpath && $f =~ m|^${filepattern}$|) {
      push(@result,$fullpath);
    }
  }
  return @result;
}

#---------------
sub FindLatestCheckpointTimeFromInputDir {
  my($inputdir)=@_;
  my(@files) = Find("Cp-EvolutionLoopControl\.txt",$inputdir);
  unless(scalar(@files)) {
    die "$ME: Cannot find Cp-EvolutionLoopControl.txt under $inputdir";
  }
  my $a= GetLatestCheckpointTime(@files);
  die "Cannot find checkpoint time for $inputdir, "
      . "perhaps no checkpoint file exists?"
      unless defined($a);
  return DecimalNotation($a);
}

#---------------
sub GetLatestCheckpointTime {
  my(@files)=@_;

  my $time = undef;
  foreach my $file (@files) {
    my @text = Utils::ReadFile($file);
    foreach (@text) {
      if(m|^t\s*=\s*([\d.+-e]+);|m) {
        $time=$1 if(!defined($time) || $1 > $time);
      }
    }
  }
  return $time;
}

#---------------
# Reads file, evals commands on ENTIRE FILE read into $_.
# If file not changed prepends 'prepend' to data, overwrites file.
sub ApplyCommandsToFile {
  my($file,$prepend,@commands)=@_;
  local $_;

  $_ = Utils::ReadFile($file);
  my $oldtext=$_;

  foreach my $code(@commands) {
    eval $code;
    if ($@) {die "$ME: Eval error, $@"};
  }

  if($oldtext eq $_) {
    if(defined $prepend) {
      $_ = $prepend . $_;
      Utils::OverwriteFile($file,$_);
    }
  } else {
    Utils::OverwriteFile($file,$_);
  }
}

#---------------
# Returns true if file exists and it contains a non-whitespace option,
# false otherwise
sub FileHasOption {
  my($file)=@_;

  return 0 unless (-f $file);
  my $text = Utils::ReadFile($file);
  $text    =~ s|#.*$||m; # remove comments from each line
  $text    =~ s|\s+||g;  # remove all whitespace
  if($text =~ m|^[^=]+=\S;|) { # Does it contain text=moretext;?
    return 1;
  } else {
    return 0;
  }
}

#---------------
sub DecimalNotation {
  my($num)=@_;

  $num = sprintf("%20.19g",$num);
  $num =~ s|^\s+||; # Remove whitespace at beginning
  return $num;
}

#---------------
# Remove all checkpoints except the latest
# for all segments at the same Lev in the same root directory

sub CleanupHydroCheckpointDirectories {
    my($oldworkdir)=@_;

    my $topdir = -1;
    my $topdirprefix = "NA";
    # 1) Get list of directories at the same level LEV within the directory Prefix
    return 0 unless ($oldworkdir =~ m/(.*)\/Lev(\d)_[A-Z]{2}$/);
    my $Prefix = $1;
    my $LEV = $2;
    print("Hydro clean-up checking Lev $LEV in directory $Prefix\n");
    opendir(ROOTDIR,"$Prefix") or die "Couldn't open root directory $Prefix";
    my @workdirlist = grep(/^Lev${LEV}_[A-Z]{2}$/,readdir(ROOTDIR));
    closedir(ROOTDIR);
    foreach my $workdir (@workdirlist)
    {
	# 2) Look through Run/Checkpoints and Run/ChainingCheckpoints
	my @dirlist = ("$Prefix/$workdir/Run/Checkpoints", "$Prefix/$workdir/Run/ChainingCheckpoints");
	foreach my $directory (@dirlist){
	    if (-d $directory)
	    {
		opendir (DIR, $directory) or die "Couldn't open directory $directory during hydro clean-up.";
		while (my $dirnb = readdir(DIR)) {
		    # Ignore ./
		    next if ($dirnb =~ m/^\./);
		    # New directory is the latest checkpoint found
		    if ($dirnb > $topdir)
		    {
			if($topdir == -1)
			{
			    # First directory under consideration
			    # Just tag it as the most recent checkpoint
			    $topdir = $dirnb;
			    $topdirprefix = $directory;
			} else
			{
			    # $dirnb is a more recent checkpoint. Remove checkpoint
			    # previously tagged as the most recent.
			    CleanupSingleHydroDirectory("$topdirprefix/$topdir");
			    # Tag new checkpoint as the most recent
			    $topdir = $dirnb;
			    $topdirprefix = $directory;
			}
		    } else{
			# New directory is older than the currently tagged latest checkpoint
			# Remove it.
			CleanupSingleHydroDirectory("$directory/$dirnb");
		    }
		}
		closedir(DIR);
	    }
	}
	# 3) Now cleanup source terms in Saved checkpoints
	# This saves most of the memory cost of hydro checkpoints
	# (but we sacrifice exact reproducibility on restarts if we
        #  need to go back to one of these checkpoints).
	my $SaveDir = "$Prefix/$workdir/Run/SavedCheckpoints/";
	next unless (-d $SaveDir);
	opendir(DIR,$SaveDir) or die "Could not open directory $SaveDir";
	while (my $savecpdir = readdir(DIR))
	{
	    # Ignore ./
	    next if ($savecpdir =~ m/^\./);
	    # Find all files of the form Source*.h5, and remove them.
	    opendir(DIR2,"$SaveDir/$savecpdir") or die $!;
	    while (my $savefilename = readdir(DIR2))
	    {
		if ($savefilename =~ m/^Source.*\.h5$/)
		{
		    unlink("$SaveDir/$savecpdir/$savefilename");
		}
	    }
	    closedir(DIR2);
	}
	closedir(DIR);
    }
}

#---------------
# Remove all files within a given directory
# Will ignore sub-directories

sub CleanupSingleHydroDirectory {
    my($dirtorm)=@_;
    print("Hydro clean-up: Removing directory $dirtorm\n");
    # Remove all files within a directory.
    opendir (DIR2, $dirtorm) or die "Couldn't open to-be-deleted directory $dirtorm.";
    while (my $file = readdir(DIR2)) {
	# Makes sure this is a file and not a directory
	next unless (-f "$dirtorm/$file");
	my $filetorm = $dirtorm."/".$file;
	unlink($filetorm);
    }
    closedir (DIR2);
    # This should only work if directory is empty (i.e. there were no
    # sub-directory).
    rmdir $dirtorm;
}

#-------------------------------------------------
1; # This is how to end sources
