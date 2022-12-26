#!/usr/bin/env perl
package SpEC;  # Start the SpEC namespace

require 5;
use warnings FATAL => 'all';
use strict;
use Cwd qw(cwd abs_path);
use File::Basename;
use Math::Trig;
use List::Util qw(min max);
use POSIX ();
# Load SpEC modules
use lib dirname(abs_path(__FILE__));
use Utils;
use Machines;

my $ME    = __FILE__;
# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(Cwd::abs_path(__FILE__)) . "/../bin";

#-------------------------------------------------
# Make bin directory and copy support items
#   src  = path to source bin directory
#   dest = path to destination bin directory
#   exec = main evolution executable (OPTIONAL)
#   LinkSpEC = bool to decide if bin/../SpEC link should be created (OPTIONAL)
#   LinkBin  = bool if true, link items to bin instead of copying (OPTIONAL)
#   Hydro    = bool if true, link hydro items (OPTIONAL)
sub MakeBinDir {
  my($src,$dest,$exec,$LinkSpEC,$LinkBin,$Hydro) = @_;
  Utils::MakePathsAbsolute(\$src,\$dest,\$exec);
  $LinkSpEC = 1 if (not defined $LinkSpEC);
  $LinkBin  = 0 if (not defined $LinkBin);
  $Hydro  = 0 if (not defined $Hydro);

  print "Making bin directory from $src\n";
  die "Cannot make bin directory, dest $dest already exists!\n" if (-d $dest);
  Utils::System("mkdir -p $dest");
  die "Failed 'mkdir -p $dest'!" if (not -d $dest);

  my @FilesForBinDir = (
    # Globs if src is a repo bin dir
    <$src/../Perl/*.pm>,
    <$src/../../MakefileRules/this_machine*.env>,
    # Globs if src is a Run bin dir
    <$src/*.pm>,
    <$src/this_machine*.env>,
    # Explicit files
    "$src/ApplyObservers",
    "$src/BBH_ID.py",
    "$src/BBH_ID_WriteInputFiles.py",
    "$src/BbhDiagnosticsImpl.py",
    "$src/CombineSegments.py",
    "$src/ComputeIDVout",
    "$src/ConstructPostMergerDataFromAhC",
    "$src/DoMultipleRuns",
    "$src/EvolutionWrapper",
    "$src/EvolveGeodesicsWrapper",
    "$src/GetMpiCmd",
    "$src/GetNsNsInspiralSubstitutions",
    "$src/GetNsNsMergedBoxFiles",
    "$src/GetSingularityCmd",
    "$src/JoinH5",
    "$src/JoinDatFiles",
    "$src/MailSpec",
    "$src/MakeBinDirectory",
    "$src/MakeNextSegment",
    "$src/MakeSubmit.py",
    "$src/ManageRuns.py",
    "$src/Nsubdomains",
    "$src/OmegaDotEccRemoval.py",
    "$src/InitialDataAdjustment.py",
    "$src/BHNSGenerateInitialParamsInput.py",
    "$src/GenerateInitialParamsInput.py",
    "$src/ParseIdParams.py",
    "$src/QueryMachines",
    "$src/RecoverOldCheckpoints",
    "$src/SurfaceToSpatialCoordMapFiles",
    "$src/TarDirectory",
    "$src/TrajectoryToSpatialCoordMapFiles",
    "$src/Utils.py",
    "$src/UpdateH5DataVersion.py",
    "$src/UpdateH5DataVersionUtils.py",
    "$src/ZeroEccParamsFromPN",
  );
  push(@FilesForBinDir,"$src/dmrgen_lessfluff") if $Hydro;
  push(@FilesForBinDir, $exec) if $exec;
  foreach my $path (@FilesForBinDir) {
    if (! -f "$path") {
      Utils::Die("Path '$path' does not exist, but is needed for the bin\n".
                 "directory. You may need to either make execs or use\n".
                 "the '-B' flag to skip bin directory creation.\n");
    }
    my $file = basename($path);
    if ($LinkBin) {
      Utils::MySymlink($path,"$dest/$file");
    } else {
      Utils::MyCopy($path,"$dest/$file");
    }
  }

  # Link $dest/SpEC to the $exec in the bin directory
  if($exec && $LinkSpEC) {
    my $ExecName = basename($exec);
    my $DestBinName = basename($dest);
    Utils::MySymlink("$DestBinName/$ExecName","$dest/../SpEC");
  }

  # Save some system info for future use
  Utils::System("printenv > $dest/env.log");
  system("module list > $dest/module.log 2>&1");
}

#-------------------------------------------------
# Sets the environment variable from a SpEC bin directory
sub LoadEnvFromBin {
  my ($BinDir, $opt_v) = @_;
  unless(defined($opt_v)) { $opt_v=0; }   # Default verbosity
  Utils::MakePathsAbsolute(\$BinDir);

  my %NewENV = %ENV;    # default is current environment
  # Source this_machine.env (use MakefileRules version if
  # $BinDir points to repo)
  my $EnvFile = "$BinDir/this_machine.env";
  $EnvFile = "$BinDir/../../MakefileRules/this_machine.env" if (! -e $EnvFile);
  if (-e $EnvFile) {
    %NewENV = Utils::GetEnvFromSubshell(". $EnvFile", $opt_v);
  } else {
    print "No this_machine.env found in $BinDir.\n";
  }
  # Put the bin dir at the top of the path
  $NewENV{'PATH'} = "$BinDir:$NewENV{'PATH'}";

  return %NewENV;
}

#-------------------------------------------------
# Edit Evolution.input files to restart FromLastStep of a checkpoint,
# potentially interpolating the cp data.
#   $CpDir = directory containing Checkpoint dirs (e.g. old scratch dir)
#   $Interp = flag representing interpolation options
#       0: None
#       1: Interpolated(...)
#       2: Remapped(...)
#       3: Hydro remapped restart (use remap only to use InterpolatorScheduler)
#   $InputFiles = arrayref containing all of the input files.
#   @Files = array of files to edit
sub ChangeEvolutionInputToRestartFromLastStep {
  my($CpDir,$Interp,$InputFiles,@Files)=@_;

  my $fnameprefix = "FilenamePrefix=";
  if(-d "$CpDir/ChainingCheckpoints") {
    $fnameprefix .= "$CpDir/ChainingCheckpoints/;";
  } elsif(-d "$CpDir/Checkpoints") {
    $fnameprefix .= "$CpDir/Checkpoints/;";
  } else {
    $fnameprefix .= "$CpDir/;";
  }

  my $GetPrefix = sub {
    GetGrDomainInputPath($CpDir) =~ m|/([^/]*)Domain.input$|;
    return $1;
  };

  my $intrpstring = "";
  if($Interp==1) {
    my $Prefix = $GetPrefix->();
    $intrpstring=
      "${Prefix}GlobalVarsCheckpoint = Interpolated(" .
      "ResolutionChanger=Spectral; " .
      "Interpolator =ParallelAdaptive(TopologicalInterpolator=CardinalInterpolator;); " .
      "DomainFile=${Prefix}Domain.input; DomainDir=$CpDir/);";
  } elsif($Interp==2) {
    my $Prefix = $GetPrefix->();
    $intrpstring=
      "${Prefix}GlobalVarsCheckpoint = Remapped(" .
      "DomainDir=$CpDir/; DomainFile=${Prefix}Domain.input; ".
      "MapPrefixGridToInertial=GridToInertial; " .
      "MapPrefixSrcGridToInertial=GridToInertial; " .
      "Interpolator=ParallelAdaptive(TopologicalInterpolator=CardinalInterpolator;););";
  } elsif($Interp==3){
      my $Prefix = "Hy";
      $intrpstring=
          "${Prefix}GlobalVarsCheckpoint = Remapped(" .
          "DomainDir=./; DomainFile=Last${Prefix}Domain.input; ".
          "MapPrefixGridToInertial=; " .
          "MapPrefixSrcGridToInertial=;".
          "InterpolatorScheduler = NonBlocking(".
          "  TopologicalInterpolator=CappedPolynomial();".
          "  DistributePoints = OnSource();".
          "  ChooseSubdomain  = FinestLevel(MinimumDistance=3);".
          "  );".
          "Extrapolator=SetToConstant(Const=1.e-16);".
          ");";
  } else {
    # Only other allowed value of Interp is 0 (i.e. standard restart)
    Utils::Die("Unexpected value Interp=$Interp") unless ($Interp==0);
  }

  my $replace = "FromLastStep($fnameprefix$intrpstring);";

  foreach my $file(@Files) {
    my $text = $InputFiles->{$file};
    if($text =~ m/^\s*Restart\s*=\s*/m) {# The 'm' means ^ and $ match any line
      # We found a restart line

      # In the next line, the suffix
      # /m changes '^' and '$' to match the start/end of any line in string.
      # Adding /s causes '.' to match newline, which is not what we want now.

      ### OLD VERSION, dies if the Restart is on a multiple line.
      ###  $text =~ s/(^\s*Restart\s*=\s*).*$/${1}$replace/m;

      ### NEW VERSION: Works if 'Restart' spans multiple lines.
      ### The complicated perl regexp matches nested parentheses.
      ### The following should work in perl5.9:
      #my $re=qr/(^\s*Restart\s*=\s*)([^()]++(\(((?:[^()]++|(?2))*)\))?)\s*;$/;
      ### But we are still using perl 5.8 some places, so instead use this:
      ### $bal matches balanced parens with possible nested parens inside.
      our $bal;  $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;
      my $re = qr/(^\s*Restart\s*=\s*)[^();]+$bal?\s*;\s*$/m;

      # Note that the qr//m behavior changed from perl 5.8.x to 5.10.0.
      # The following matches in 5.8.x, doesn't match in 5.10.0:
      # $re = qr/^bar/; "foo\nbar" =~ /$re/m;
      # So we need to be careful about using qr with /m.

      unless($text =~ s/$re/${1}$replace/m) {# 'm' means ^ and $ match any line
        Utils::Die("Cannot restart from last step in $file");
      }
    } else {
      # This is in case we don't find 'Restart'
      $text = "Restart=$replace\n" . $text;
    }
    $InputFiles->{$file} = $text;
  }
}

#-------------------------------------------------
# Given an old submission directory, return the new submission directory
# for the next segment by incrementing the counter suffix, e.g. "_AA"->"_AB"
sub GetNewSubmissionDirectory {
  my($OldSubmissionDir)=@_;

  my $NewSubmissionDir;
  if($OldSubmissionDir =~ m|(.*/)([^/]+)$|) {
    $NewSubmissionDir = $1;
    my $Dir    = $2;
    if($Dir =~ m|(.*_)([A-z]+)$|) {
      my $DirPrefix = $1;
      my $Counter = $2;
      ++$Counter;
      $NewSubmissionDir .= $DirPrefix . $Counter;
    } else {
      $NewSubmissionDir .= $Dir . "_a";
    }
  } else {
    Utils::Die("Error with submission directory '".$OldSubmissionDir."'");
  }
  return $NewSubmissionDir;
}

#-------------------------------------------------
# Notify the user about the termination or merger of a batch job
#   reason  = termination reason
#   work    = path to work directory
#   scratch = path to scratch directory
#   opt_a   = comma-delimited list of e-mail addresses to notify
#   opt_f   = file to append termination info
#   JobID   = job number assigned by the system job manager
sub JobNotify {
  my($reason, $work, $scratch, $opt_a, $opt_f, $JobID) = @_;

  my $Machine;
  eval { $Machine = new Machines(); };

  # Set default values of input arguments
  unless(defined $reason)  { $reason = "Unknown"; }
  unless(defined $work)    { $work = cwd(); }
  unless(defined($scratch)) {   # Must come after $work
    eval { $scratch = $Machine->GetScratchDir(); };
    if ($@) { $scratch = "$work/Run"; }
  }
  unless(defined $opt_a)   { $opt_a = ""; }
  unless(defined $opt_f)   { $opt_f = ""; }
  unless(defined $JobID)   { $JobID = "Unknown"; }

  # Remove all whitespace and newlines from $opt_a and $opt_f
  $opt_a =~ s/(\s+|\n)//g;
  $opt_f =~ s/(\s+|\n)//g;

  if (($opt_a eq "" || $opt_a eq "__EmailAddresses__") &&
      ($opt_f eq "" || $opt_f eq "__TerminationInfoFile__")) {
    print STDERR "JobNotify has nothing to do! Returning...\n";
    return;
  }

  # gather some additional information about the failed run
  my $AdditionalInfo = "";

  my $StartTime = `grep StartTime $scratch/Evolution.input 2>/dev/null`;
  if (defined($StartTime) && $StartTime =~ /StartTime\s*=\s*([\d.]+);/) {
    $StartTime = sprintf("%.3f",$1);
  } else {
    $StartTime = "Not found";
  }

  my $TermTime = `tail -n1 $scratch/TStepperDiag.dat 2>/dev/null \\
    | awk '{print \$1}'`;
  $TermTime = $TermTime ? sprintf("%.3f",$TermTime) : "Not found";

  my $DirBase = $work;
  $DirBase =~ s|/Run$||;
  $DirBase =~ s|[A-Za-z]{1,2}$||;

  my $ProperSep = `find ${DirBase}*/Run/App* -name HorizonSepMeasures.dat \\
    -follow 2>/dev/null | xargs cat | tail -n1 | awk '{print \$2}'`;
  $ProperSep = $ProperSep ? sprintf("%.3f",$ProperSep) : "Not found";

  my $GhCeLinf = `find ${DirBase}*/Run/Constraint* -name GhCe_Norms.dat \\
    -follow 2>/dev/null | xargs cat | tail -n1 | awk '{print \$4}'`;
  $GhCeLinf  = $GhCeLinf  ? sprintf("%.3g",$GhCeLinf)  : "Not found";


  $AdditionalInfo .=
      sprintf("StartTime:       %9s   Proper Sep: %s\n",$StartTime,$ProperSep);
  $AdditionalInfo .=
      sprintf("TerminationTime: %9s   GhCe Linf:  %s\n",$TermTime,$GhCeLinf);

  if ($reason =~ m/^(SpEC failed|Unknown|)$/) {
    # try to find a reason for termination
    my $text = Utils::ReadFile("$scratch/SpEC.out");
    if ($text =~ m/points cannot be interpolated/) {
      $reason = "Non-interpolated points";
    } elsif ($text =~ m/Floating point exception/) {
      $reason = "Floating point exception";
    } elsif ($text =~ m/New ran out of memory/) {
      $reason = "Out of memory";
    } elsif ($text =~ m/Segmentation fault/) {
      $reason = "Segmentation fault";
    } elsif ($text =~ m/REQUIRE\s+FAILED/) {
      $reason = "Require failed";
    } else {
      $reason = "Terminated (no reason found)";
    }
  }

  # Determine jobname
  my $job = Utils::SystemOutput("$RealBin/MakeSubmit.py " .
                                "-d $work query Jobname");
  chomp($job);

  # Create subject and body of message
  my $Host = Utils::Hostname();
  my $Subject = "SpEC: $job -- $reason";
  my $Body =
       "Reason:     $reason\n"
      ."HostName:   $Host\n"
      ."JobID:      $JobID\n"
      ."WorkDir:    $work\n"
      ."ScratchDir: $scratch\n\n"
      .$AdditionalInfo."\n";

  # Include output from Error.txt files if found
  my @files = <$scratch/Error*.txt>;
  if (@files) {
    $Body .= "Contents of $files[0] are below:\n"
      ."----------------------------------------------------------------\n"
      .`cat $files[0]`
      ."----------------------------------------------------------------\n\n";
  } else {
    $Body .= "No Error*.txt files found.\n\n";
  }

  # Include output from SpEC.stderr if found
  my $file = "$work/SpEC.stderr";
  if (-e $file) {
    $Body .= "The last 10 lines of SpEC.stderr are below:\n"
      ."----------------------------------------------------------------\n"
      .`tail -n10 $file`
      ."----------------------------------------------------------------\n\n";
  } else {
    $Body .= "No SpEC.stderr found.\n\n";
  }

  # Include output from SpEC.out if found
  $file = "$scratch/SpEC.out";
  if (-e $file) {
    $Body .= "The last 30 lines of SpEC.out are below:\n"
      ."----------------------------------------------------------------\n"
      .`tail -n30 $file`
      ."----------------------------------------------------------------\n\n";
  } else {
    $Body .= "No SpEC.out found.\n\n";
  }

  # send email notification
  if($opt_a ne "" && $opt_a ne "__EmailAddresses__") {
    my $mailer = ($Machine ? sub {$Machine->SendMail(@_);} : \&Utils::SendMail);
    eval { $mailer->(subject=>$Subject,body=>$Body,to=>$opt_a,v=>1); };
    if($@) { #Catch errors in mailing and warn about them, but don't die
      warn "Non-fatal error caught using Utils::SendMail: $@\n";
    }
  }

  # append to notification file
  if($opt_f ne "" && $opt_f ne "__TerminationInfoFile__") {
    # construct a shorter message than the emailed one
    my $date = `date`; chomp($date);
    my $ShortBody =
        "Job terminated at $date\n"
        ."Reason:  $reason\n"
        ."WorkDir: $work\n"
        .$AdditionalInfo;

    if ($reason eq "" || $reason eq "Require failed") {
      # Check for Error.txt files.  If found, include output from one of them
      my @files = <$scratch/Error*.txt>;
      if (@files) {
        $ShortBody .= "Beginning of $files[0]:\n"
            .`head -n5 $files[0]`;
      } else {
        $ShortBody .= "No Error*.txt files found.\n";
      }
    }
    $ShortBody .= "========================================================\n";

    print STDERR "Appending message about termination to file: $opt_f\n";
    eval { Utils::AppendToFile($opt_f, $ShortBody); };
    if($@) { #Catch errors and warn about them, but don't die
      warn "Non-fatal error caught using Utils::AppendToFile: $@\n";
    }
  }
}

#-------------------------------------------------
# Block comment out code that matches a regexp and replace it with text.
# Do it nicely so that code blocks still line up.
#   text    = Original text
#   regexp  = Regexp to match (must be first text on its line)
#   replace = Text to replace match with
#   count   = Exact number of times regexp should be matched (optional)
sub BlockReplace {
  my($text,$regexp,$replace,$count) = @_;

  my @matches = ($text =~ m/\n\s*$regexp/g);
  if (defined $count) {
    my $N = scalar(@matches);  # number of matches
    die "Matched $regexp $N times, expected $count" if ($count != $N);
  }

  foreach my $match (@matches) {
    # Grab and comment out the old block of code
    $match =~ m/\n(\s*)($regexp)/;
    my $Indent = $1;
    my $Block = "\n$2";
    $Block =~ s/\n/\n#%%%/g;

    $match = quotemeta($match);
    if ($replace eq "") {
      # Replacement is empty, so don't add anything
      $text =~ s/$match/$Block/;
    } else {
      # Replacement is NOT empty, so add it
      $text =~ s/$match/$Block\n$Indent$replace/;
    }
  }
  return $text;
}

#-------------------------------------------------
# Could be called CopyContentsToHere. Copies directory contents here.
# If a file already exists in $dest, it will not be copied.
sub CopyInputFiles {
  my ($source,$dest) = @_;
  Utils::Die("$source is not a directory!") if (not -d $source);
  for my $srcfile (glob("$source/*")) {
    next if (-d $srcfile); # skip directories
    my $destfile = "$dest/" . basename($srcfile);
    Utils::MyCopy($srcfile, $destfile) unless (-e $destfile);
  }
}

#-------------------------------------------------
# This is a workaround for a filesystem problem on Bridges.
# When accessing a file in quick succession
# on Bridges' /pylon2 filesystem, the file becomes corrupted.
# I (Mark S.) have tested that both of the sleep()s and the
# second OverwriteFile are necessary to work around this
# problem, but I do not understand why.
#
sub OverwriteFileWithPossibleSleep {
  my($file,$text) = @_;
  Utils::OverwriteFile($file,$text);

  my $Machine;
  eval { $Machine = new Machines(); };
  if($Machine and $Machine->GetHost() eq "bridges") {
    sleep(10);
    Utils::OverwriteFile($file,$text);
    sleep(10);
  }
}


#---------------
# Returns the absolute path of the relevant Domain.input file.
#   dir = path to look for file (default: cwd)
sub GetGrDomainInputPath {
  my ($dir) = @_;
  $dir = cwd() if (not $dir);

  my $DomainInputRelPath;
  if (-e "$dir/GrDomain.input") {
    if(-e "$dir/Domain.input") {
      die "Found GrDomain.input and Domain.input in $dir.\n".
          "Is this a hydro run or not?  Don't know what to do.\n";
    }
    $DomainInputRelPath = "$dir/GrDomain.input";
  } elsif(-e "$dir/Domain.input") {
    if(-e "$dir/HyDomain.input") {
      die "Found HyDomain.input and Domain.input in $dir.\n".
          "Is this a hydro run or not?  Don't know what to do.\n";
    }
    $DomainInputRelPath = "$dir/Domain.input";
  } else {
    Utils::Die("Could not find GrDomain.input or Domain.input in $dir.\n");
  }
  return abs_path($DomainInputRelPath);
}

#---------------
# Return time and dir of the latest checkpoint within $inputdir
sub FindLatestCheckpointTimeFromInputDir {
  my($inputdir)=@_;
  my(@files) = Utils::Find('Cp-EvolutionLoopControl.txt',$inputdir);
  unless(@files) {
    Utils::Die("Cannot find Cp-EvolutionLoopControl.txt under $inputdir");
  }
  my ($a,$b) = GetLatestCheckpointTime(@files);
  Utils::Die("Cannot find checkpoint time for $inputdir, "
    . "perhaps no checkpoint file exists?") unless defined($a);
  return (DecimalNotation($a),$b);
}

# Helper function
sub DecimalNotation {
  my($num)=@_;
  $num = sprintf("%20.19g",$num);
  $num =~ s|^\s+||; # Remove whitespace at beginning
  return $num;
}

# Helper function
sub GetLatestCheckpointTime {
  my(@files)=@_;

  my $time   = undef;
  my $myfile = undef;
  foreach my $file (@files) {
    my @text = Utils::ReadFile($file);
    foreach (@text) {
      if(m|^t\s*=\s*([\d.+e]+);|m) {
        if(!defined($time) || $1 > $time){
          $time   = $1;
          $myfile = $file;
          $myfile =~ s|(.*)/[^/]+$|$1|;
        }
      }
    }
  }
  return ($time,$myfile);
}

#---------------
# Return the checkpoint directory in $inputdir matching $time to within $eps
sub GetCheckpointDirectory {
  my($inputdir,$time,$eps)=@_;

  my @files    = Utils::Find('Cp-EvolutionLoopControl.txt',$inputdir);
  my $thisfile = undef;
  foreach my $file (@files) {
    my @text = Utils::ReadFile($file);
    foreach (@text) {
      if(m|^t\s*=\s*([\d.+e]+);|m) {
        if(abs($1-$time)<$eps*(1.0+abs($1)+abs($time))) {$thisfile=$file;}
      }
    }
  }
  die "Cannot find checkpoint file with $time within $eps"
      unless defined($thisfile);
  $thisfile =~ s|(.*)/[^/]+$|$1|;
  return $thisfile;
}

#---------------
# Reads a raw domain file conaining a 'SphericalShells3D SphereC' entry.
# Outputs the maximum spherical radius in grid coordinates (units of M).
sub GridOuterBoundary {
  my($dir)=@_;

  my $file = "$dir/D.input";
  $file = "$dir/Dnew.input" unless (-f $file);
  $file = "$dir/Domain.input" unless (-f $file);
  $file = "$dir/GrDomain.input" unless (-f $file);
  die "Cant find {D,Dnew,Domain,GrDomain}.input in $dir" unless (-f $file);

  my $text = Utils::ReadFile($file);
  $text =~ s/#.*$//mg; # Remove comments ('m' means $ matches end of each line)
  $text =~ s/\n//g;    # Remove newlines
  $text =~ s/\s+//g;   # Remove whitespace

  my $Rmax0=undef;
  ### $bal matches balanced parens with possible nested parens inside.
  our $bal;  $bal = qr/\((?:(?>[^()]+)|(??{$bal}))*\)/;
  while($text =~ m/SphericalShells3D\s*(${bal})/g) {
    my $dum = $1;
    if($dum =~ m/SphereC/) {
      if($dum =~ m/Bounds\s*=\s*[^;]*,([^;]+);/) {
        $Rmax0 = eval $1;
      }
    }
  }
  die "Cannot read Rmax0 from $file" unless ($Rmax0);
  return $Rmax0;
}

# Given a directory, determines if this directory contains hydro input
# files.  Currently does this by detecting if a HyDomain.input file exists.
sub IsHydro {
  my($dir)=@_;
  if($dir) {
    $dir .= "/" unless($dir=~m|/$|);
  } else {
    $dir = ""; # In case $dir was undef.
  }
  # Check if this is a lensing or EH find run. They are run in directories with
  # a Geodesic.input file so check for this file.
  if (-f "${dir}Geodesic.input") {
    return 0;
  }
  my $h = "${dir}HyDomain.input";
  my $g = "${dir}GrDomain.input";
  my $d = "${dir}Domain.input";
  unless (-f "$g" or -f "$d") {
    # Demand that there is a Domain.input or GrDomain.input. Because
    # otherwise we are probably not pointing at input files.
    die "$ME: Cannot find $g or $d. Does $dir point to input files?";
  }
  if(-f "$h") {
      unless (-f "$g") {
        # Demand that there is a GrDomain.input. Because otherwise
        # we are probably not pointing at input files.
        die "$ME: Cannot find $g. Does $dir point to input files?";
      }
    return 1; # 1 means true in perl.
  }
  return 0;
}

#---------------
# Returns number of subdomains in a Domain.input file.
sub NSubdomains {
  my($file)=@_;

  my $Machine;
  eval { $Machine = new Machines(); };

  my $SingularityCmd = $Machine ? $Machine->GetSingularityCmdIfNeeded() : "";
  if($SingularityCmd ne "") {
    $SingularityCmd .= " "; # Add space between 'singularity exec' and command.
  }

  my $exec = $SingularityCmd . "$RealBin/Nsubdomains";
  my $N = Utils::SystemOutput("$exec -d $file");
  die "Could not find integer Nsubdomains in '$N'" if (not $N =~ m/^\d+$/);
  chomp($N);
  return $N;
}

# Returns the total number of subdomains in any *Domain.input file in $dir.
sub TotalNSubdomains {
  my ($dir) = @_;
  my $N = 0;
  foreach my $file (glob "$dir/*Domain.input") {
    $N += NSubdomains($file);
  }
  die "Could not find any *Domain.input files in $dir to extract the number ".
    "of subdomains." unless $N > 0;
  return $N;
}

#---------------
# Determines the radii of SphereC subdomain boundaries and takes into account
# possible junkyard shells that will be dropped later in the run.
sub SetSphereCRadii {
  my ($DropJunkShell,$Rmax,$Rmin,$NumSphereCs,$JunkWidth) = @_;

  # If $DropJunkShell = 0 (i.e. false), then $DropJunkShellTime
  # will remain set to 0. The DropJunkShellTime termination
  # criterion interprets this as "Do nothing", so no termination
  # or shell dropping will occur.
  my $DropJunkShellTime = 0;
  my $RmaxAfterDrop = $Rmax;
  # Radial width of a SphereC subdomain
  my $Rwidth = ($Rmax - $Rmin) / $NumSphereCs;
  if ($DropJunkShell) {
    # Increase $Rmax to account for the junkyard shells. The new $Rmax will
    # be set to the first multiple of $Rwidth that is greater than or equal
    # to $JunkWidth. We need to ensure that the junkyard shells have the same
    # width or junk radiation will start to reflect back into the domain at
    # inner boundary of the first junkyard shell.
    $Rmax += $Rwidth * POSIX::ceil($JunkWidth/$Rwidth);
  }

  # Construct an array of the radii of the SphereC subdomain boundaries.
  # This will be used by GrDomain.input to construct the subdomains.
  my @RadiiArray = ();
  for (my $R=$Rmin; $R <= $Rmax+1e-12; $R+=$Rwidth) {
    push(@RadiiArray, $R);
  }
  # Make sure our value of $Rmax is EXACTLY the radius of the outer boundary
  $Rmax = $RadiiArray[-1];
  my $Radii = join ',', @RadiiArray;
  if ($DropJunkShell) {
    # Trigger the run to stop at one light-crossing time
    $DropJunkShellTime = $Rmax;
  }

  if ($DropJunkShell) {
    # informational output
    print "----------------------------------------\n";
    print "Junkyard shells will be used to prevent reflected junk radiation.\n";
    printf("Rmax will be extended from %d to %d until t=%f.\n",
           $RmaxAfterDrop, $Rmax, $DropJunkShellTime);
    print "----------------------------------------\n";
  }

  return $Radii, $Rmax, $DropJunkShellTime;
}

#---------------
# Determine WaveExtraction Rmin based on SphereC and (bound) eccentric Newtonian orbit
sub AutoWaveExtractionRmin {
  my ($RminSphereC,$Omega,$Mtotal,$ExtendExtractionRegion) = @_;
  # We want WaveExtractionRmin to be outside of the inner radius of SphereC0.
  # (Otherwise, the extraction radius crosses many subdomain boundaries and
  #  this creates noise in the extracted quantities).
  # So we choose it 5% larger, i.e. 1.05*RminSphereC
  # We also want WaveExtractionRmin to be an int, so we
  # choose the smallest int that is greater than this value.
  my $WaveExtractionRmin1 = int(1.05*$RminSphereC+1.0);
  # We also want WaveExtractionRmin to be greater than 1 GW wavelength.
  # (Near-field effects fall off as powers of (lambda/(2 pi r)), so 1 wavelength
  # puts us in a regime where these effects are down by a factor of 1/(2pi) ).
  # The GW wavelength is pi/Omega_orbital; pick next largest int.
  # If Omega is zero, then ignore this condition.
  my $WaveExtractionRmin2 = undef;
  if($ExtendExtractionRegion==0){
    $WaveExtractionRmin2 = $Omega==0.0 ? 0.0 : int(pi/$Omega+1.0);
  }
  elsif($ExtendExtractionRegion==1){
    # If we are extracting all the Weyl scalars then we need to extend Rmin
    # inwards as close as possible to reduce noise due to backscattered junk
    # radiation in the lower index Weyl scalars. Here, WaveExtractionRmin is set
    # to the next largest int above two reduced GW wavelengths.
    # If Omega is zero, then ignore this condition.
    $WaveExtractionRmin2 = $Omega==0.0 ? 0.0 : int(1/$Omega+1.0);
  }
  else {
    die "ERROR: ExtendExtractionRegion should be set to either 0 or 1.";
  }
  my $WaveExtractionRmin  = max($WaveExtractionRmin1, $WaveExtractionRmin2);

  # informational output
  print "----------------------------------------\n";
  printf("For Omega0=%f, RminSphereC0=%f, Wave extraction Rmin is %d\n",
         $Omega, $RminSphereC,$WaveExtractionRmin);
  print "----------------------------------------\n";

  if($ExtendExtractionRegion==0){
    # For a fixed tMerger, Omega0 increases as we go towards large q and the
    # above method can result in a very small WaveExtractionRmin (~85M for q=30
    # for 20 orbits). For best results extrapolating h and Psi4, we want all
    # extraction spheres to be far away from the near-zone effects which can have
    # lambda/r as well as Mtotal/r dependence. Therefore, we limit how small
    # WaveExtractionRmin can be.
    return max($WaveExtractionRmin, int(100*$Mtotal));
  }
  else{
    # If we are extracting all the Weyl scalars then we need to go closer than
    # 100*$Mtotal. The extrapolation procedure can handle this for h and Psi4
    # but higher extrapolation orders will be needed. More care should be taken
    # in post-processing to make sure extrapolation converges.
    return $WaveExtractionRmin;
  }
}

# Gets all the Levs in a directory.
sub GetLevs {
  my($dirr) = @_;
  opendir(my $dir, $dirr) || die "can't opendir $dirr: $!";
  my @Levs = map { /^Lev([-\d]+)/ ? $1 : () } readdir($dir);
  closedir $dir;
  return @Levs;
}

# Determine the $EvDir, $Lev and $LastSegmentSuffix from $workdir.
# $workdir is the path to some segment.  $EvDir is the path to Ev,
# $Lev is the Lev, and $LastSegmentSuffix is the part after the
# underscore in the name of the segment.
# For example, for $workdir = <blah>/Ecc2/Ev/Lev3_AF,
#   $EvDir = <blah>/Ecc2/Ev
#   $Lev = 3
#   $LastSegmentSuffix=AF
sub SplitWorkDirPath {
  my ($workdir) = @_;

  my $PathRegexp = qr|^(.*)/(Ev[^/]*)/Lev([-\d]+)_([A-Z]+)|;
  $workdir =~ m|$PathRegexp| ||
    die "$workdir does not match expected structure $PathRegexp";
  my $RunHome = $1;
  my $EvDir = "$RunHome/$2";
  my $Lev = $3;
  my $LastSegmentSuffix = $4;

  return ($EvDir, $Lev, $LastSegmentSuffix);
}

#---------------
# Determine Rmax based on (bound) eccentric Newtonian orbit
sub AutoRmax {
  my ($MA, $MB, $Separation, $Omega, $aDot, $MinRadiusFrac,
      $MinimumRmax) = @_;

  $MinimumRmax = 600 unless defined($MinimumRmax);

  my $Mtotal = $MA + $MB;

  # total energy, normalized by reduced mass
  my $SpecificEnergy = 0.5*$Separation**2*($Omega**2+$aDot**2)
                       -$Mtotal/$Separation;
  if ($SpecificEnergy >= 0) {
    # hyperbolic orbit, die
    die "Newtonian specific energy is positive, indicating a hyperbolic\n"
        ."encounter. If you want to proceed, set HyperbolicOrbit to 1.\n";
  }

  # angular momentum, normalized by $Mtotal*$Mreduced
  my $SpecificAngularMom = ($Mtotal*$Omega)*($Separation/$Mtotal)**2;

  # semi-major axis with units of length
  my $SemiMajorAxis = -$Mtotal/(2*$SpecificEnergy);

  # orbital period in units of length
  my $OrbitalPeriod = 2*pi*$Mtotal*($SemiMajorAxis/$Mtotal)**1.5;

  # compute Newtonian eccentricity; only used for diagnostic output
  my $Ecc = sqrt(1+2*$SpecificEnergy*$SpecificAngularMom**2);

  # (1) Rmax based on GW frequency, estimated from orbital period
  #     Factor '2.5' is a fudge-factor.  Based on two circular
  #     data-points sampled by Harald (q=6, chiA=-0.96, q=8, chi=0),
  #     OrbitalPeriod tends to be ~20% smaller than 2pi/Omega_0,
  #     so the '2.5' should ensure roughly the same Rmax for circular
  #     binaries as before.
  my $RmaxGW = 2.5*$OrbitalPeriod;

  # (2) Rmax based on motion of outer boundary.
  #     By time-of-merger, outer boundary should have moved in
  #     to at most $MinRadiusFrac=0.9 of its initial value

  #   (2a) 1PN merger time for **circular** orbits
  #        Peters, 1964, Eq. (5.10).  Of course, also in
  #        Blanchet, arXiv:1310.1528) Eqs. 230,315,316.
  my $Tc = 5/256*($SemiMajorAxis)**4/($MA*$MB*$Mtotal);

  #   (2b) include safety factor so Tc is greater than the actual Tc,
  #        by at least 5k (relevant for short runs).  For very long runs,
  #        the factor 2 is based on observations by Harald that for two
  #        of his circular runs, SemiMajorAxis was smaller by factor ~1.14
  #        than D0.  1.14^4=1.7, so factor 2 approximately restores
  #        earlier behavior.
  my $TcSafety = max(2*$Tc, $Tc+5000);

  #   (2c) Eccentric runs merge faster, so the value is just yet more
  #        conservative. [To improve, use Peters, 1964].

  #   (2d) Choose Rmax sufficiently large, such that given the
  #        estimated shift at Rmax, the outer boundary moves in
  #        to at most $MinRadiusFrac*$Rbdry
  my $ShiftSafety = 2;  # extra safety factor
  my $RmaxBdry = (4*$ShiftSafety*$TcSafety/(1-$MinRadiusFrac))**(1/3)
                    /$MinRadiusFrac;

  # (3)  Combine all considerations for Rmax into final Rmax
  my $Rmax = max($RmaxGW, $RmaxBdry);
  if ($Rmax < $MinimumRmax) {
    print "WARNING: Rmax=$Rmax is very small.  Resetting Rmax=$MinimumRmax\n";
    $Rmax = $MinimumRmax;
  }
  if ($Rmax > 6000*$Mtotal) {
    print "WARNING: Rmax=$Rmax is very large.  Resetting Rmax=6000\n";
    $Rmax = 6000*$Mtotal;
  }
  # For large q, the above method will result in a very small Rmax (~250M for
  # q=30). This can be a problem for the outer boundary condition. We want the
  # outer boundary to be far away from all near-zone effects which can have
  # lambda/r as well as Mtotal/r dependence. Therefore, we limit how small Rmax
  # can be.
  if ($Rmax < 600*$Mtotal) {
    print "WARNING: Rmax=$Rmax is very small.  Resetting Rmax=600\n";
    $Rmax = 600*$Mtotal;
  }

  # informational output
  print "----------------------------------------\n";
  printf("Newtonian estimates based on D=%f, Omega0=%f, adot0=%f\n",
         $Separation, $Omega, $aDot);
  printf("  eccentricty e=%f, semi-major axis %f, orbital period P=%f\n",
         $Ecc, $SemiMajorAxis, $OrbitalPeriod);
  printf("  0-PN coalescence time Tc=%f, Rmax=%f.\n", $Tc, $Rmax);
  print "----------------------------------------\n";

  return $Rmax, $SemiMajorAxis;
}

# Verify that given Rmax results in sufficiently slow outer boundary speed
sub CheckOuterBdrySpeed {
  my ($Rmax, $OuterBdryDriftSpeed) = @_;

  my $OuterBdrySpeed = abs($OuterBdryDriftSpeed*$Rmax);
  if ($OuterBdrySpeed >= 1) {
    die "Outer boundary motion would be superluminal: ".
        "drift=$OuterBdryDriftSpeed, Rmax=$Rmax";
  } elsif ($OuterBdrySpeed > 0.1) {
    print "WARNING: Outer boundary speed '$OuterBdrySpeed' is large\n";
  }
}

# Choose ideal NSphereC based on given Rmax and SemiMajorAxis
sub AutoNSphereC {
  my ($Rmax, $SemiMajorAxis, $Mtotal) = @_;

  # Determine number of SphereC's.
  #   (1) width <= 2*SemiMajorAxis
  #       [generalized to eccentric orbits from previous input files,
  #        which had 2*$Separation.]
  #   (2) width <= 40*$Mtotal
  #       [to have a width comparable to waveform features of late inspiral]
  my $MaxSphereCWidth = 40*$Mtotal;
  if (defined $SemiMajorAxis) {
    # SemiMajorAxis not defined for hyperbolic orbits
    $MaxSphereCWidth = min($MaxSphereCWidth, 2*$SemiMajorAxis);
  }

  # Assumes interior to SphereC's is about 1 SphereC width
  my $NumSphereCs = int($Rmax/$MaxSphereCWidth)-1;

  # sanity
  my $MaxSphereCs = 100;
  if ($NumSphereCs > $MaxSphereCs) {
    print "WARNING: Would want to use '$NumSphereCs' sphere Cs.\n";
    print "Using $MaxSphereCs instead.\n";
    $NumSphereCs = $MaxSphereCs;
  }

  return $NumSphereCs;
}

#-------------------------------------------------
1; # This is how to end sources
