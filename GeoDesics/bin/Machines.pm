#!/usr/bin/env perl
package Machines;  # Start the Machines namespace

require 5;
use warnings FATAL => 'all';
use strict;
use POSIX qw(ceil);
use List::Util qw(min max);
use Cwd;
use File::Basename;
# Load Utils module
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;
use DetermineQueue;
use DotSpecOptions;

# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(__FILE__) . "/../bin";

# Constructor for class holding info about the current machine
#
# Possible fields are:
#   AllowedCPN      : list of cores per node that can be requested
#   Backup          : executes a backup of the Run directory
#   BatchEnv        : defines the type of submission environment
#   BatchHeader     : submission script header
#   DefaultQueue    : default submission queue
#   DefaultAccount  : default submission account
#                     Can be overridden by the same field in ~/.SpEC
#   AltMail         : bool specifying if MailSpEC workaround should be used.
#                     Set to 1 (True) if the mail command is non-functional,
#                     and then see MailSpEC for black-holes.org instructions.
#   JobID           : batch ID number for an active job
#   Jobname         : batch name for an active job
#   MaxRuntimeHours : wallclock limit for this machine (and queue)
#   MpiCmd          : mpi program and options used to invoke exec
#   MpiThreadingCmd : mpi program and options used in hybrid MPI + OpenMP runs
#                     (useful for the ray-tracing code)
#   Nprocs          : number of procs requested by an active job
#   MasterNodeName  : Function that returns the name of the master node,
#                     i.e. the node that is running the batch script before
#                     'mpirun' is invoked.  Should be called only in a running
#                     batch job; otherwise returns undef.
#   NodeNameList    : Function that returns the names of all the nodes a job
#                     is running on. Should be called only in a running
#                     batch job; otherwise returns undef. 
#   RegExp          : regexp matching head and compute node hostnames.
#                     Please try to make these as specific as possible!
#                     A regexp that matches different machines is a problem.
#   RestrictByQueue : filters possible (cpn,core) configs based on the queue
#   TmpTensorYlmDbPath : If this is undef, then TensorYlmDb should not be
#                        copied to /tmp.  Otherwise, this is the path that
#                        TensorYlmDb should be copied to before running
#                        an evolution, e.g. /tmp/TensorYlmDB.
#   CanonicalTensorYlmDb : Should point to a canonical TensorYlmDb directory
#                          that will be used instead of the user's TensorYlmDb
#                          directory if the user's TensorYlmDb directory is
#                          nonexistent or empty.
#                          New users will always have a nonexistent or
#                          empty TensorYlmDb the first time they run SpEC.
#   ScratchDir      : returns the path to the run's scratch dir
#   SshName         : ssh to this hostname when resubmitting
#                     If empty string is given, do not ssh to resubmit
#   SharedNodeQueue : If defined, returns the name of the queue to use for
#                     shared-node jobs, where the job uses fewer cores than
#                     the total cores on one node.  The other cores on the
#                     node are left free to be used by other users' jobs.
#   SharedNodeMaxCores: If SharedNodeQueue is defined, specifies the max
#                       number of cores that can be requested in that queue.
#   SubmissionCmd   : command to invoke a batch job (e.g. "qsub")
#   SingularityCmd  : command and options to invoke Singularity.
#
# If the field has a default value or a special behavior, this is
# identified in the "PUBLIC METHODS" for the Getter functions below.
sub new {
  # Complete database of all machines
  my %FULLDATABASE = (
  "bluewaters" => {
    # Each BlueWaters XE node has 32 integer cores and 16 floating-point units.
    # By default, a node uses all 32 of its "processors" (i.e. integer cores),
    # with each FP unit shared between 2 MPI processes.
    # For SpEC BBH simulations, we instead prefer to use 16 "procs" per node,
    # so that each MPI process has a dedicated FP unit. This reduces memory use.
    # To return to the BlueWaters default, change CPN to 32, and change Aprun's
    # flags to re-enable both integer cores per FP unit ("-j 2").
    "AllowedCPN"        => [16],
    "BatchEnv"          => "PBS",
    "BatchHeader"       => \&HeaderPBS_BlueWaters,
    "DefaultAccount"    => "balf",
    "DefaultQueue"      => "normal",
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Aprun,
    "MpiThreadingCmd"   => \&Aprun_threads,
    "RegExp"            => '(^h2ologin|^nid[0-9]+)',
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "briaree" => {
    "AllowedCPN"        => [12],
    "BatchEnv"          => "PBS",
    "MaxRuntimeHours"   => "168",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^briaree|^node-[a-g][0-9]-[0-9][0-9]",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "briaree1",
  },
  "bridges" => {
    "AllowedCPN"        => [28],
    "BatchEnv"          => "SBATCH",
    "DefaultQueue"      => "RM",
    "MaxRuntimeHours"   => 48,
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => 'bridges\.psc\.edu',
    "ScratchDir"        => \&LinkedScratchBridges,
    "SshName"           => "",
  },
  "bridges2" => {
    "AllowedCPN"        => [128],
    "BatchEnv"          => "SBATCH",
    "DefaultQueue"      => "RM",
    "SharedNodeQueue"   => "RM-shared",
    "SharedNodeMaxCores"=> 64,
    "MaxRuntimeHours"   => 48,
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => 'bridges2',
    "CanonicalTensorYlmDb" =>
        "/ocean/projects/phy160003p/mscheel/TensorYlmDBFull",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "caltechhpc" => {
    "AllowedCPN"        => [32,56],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_CaltechHPC,
    "DefaultAccount"    => "sxs",
    "DefaultQueue"      => "any",
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => '^(login\d|hpc-\d\d-\d\d)\.cm\.cluster$',
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "comet" => {
    "AllowedCPN"        => [24],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Comet,
    "DefaultQueue"      => "normal",
    "AltMail"           => 1,
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^comet",
    "ScratchDir"        => \&SimpleScratch,
  },
  "datura" => {
    "AllowedCPN"        => [12],
    "BatchEnv"          => "SGE",
    "BatchHeader"       => \&HeaderSGE_Datura,
    "DefaultQueue"      => "daturamon.q",
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => '(^login-damiana|submit|^kop\d+)(\.damiana\.admin)?',
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "login-damiana",
  },
  "deac" => {
    "AllowedCPN"        => [8],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH,
    "DefaultQueue"      => "small",
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "deac",
    "ScratchDir"        => \&SimpleScratch,
  },
  "edison" => {
      "AllowedCPN"        => [24],
      "BatchEnv"          => "SBATCH",
      "DefaultQueue"      => "regular",
      "MaxRuntimeHours"   => "36",
      "MpiCmd"            => \&Srun,
      "RegExp"            => '(^edison)',
      "ScratchDir"        => \&SimpleScratch,
      "SshName"           => "edison",
  },
  "frontera" => {
    "AllowedCPN"        => [56],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Frontera,
    "MasterNodeName"    => \&MasterNodeName_Slurm,
    "NodeNameList"      => \&NodeNameList_Slurm,
    "DefaultAccount"    => "PHY20018",
    "DefaultQueue"      => "small",
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Ibrun,
    "TmpTensorYlmDbPath" => "/tmp/TensorYlmDB",
    "CanonicalTensorYlmDb" => "/scratch3/projects/sxs/ux450022/TensorYlmDB",
    "RegExp"            => '(^frontera|\S+\.frontera\.)',
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "login1.frontera.tacc.utexas.edu",
  },
  "expanse" => {
    "AllowedCPN"        => [128],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Expanse,
    "DefaultAccount"    => "csd275",
    "DefaultQueue"      => "compute",
    "SharedNodeQueue"   => "shared",
    "SharedNodeMaxCores"=> 127,
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Mpirun,
    "CanonicalTensorYlmDb" =>
        "/expanse/lustre/projects/csd275/ux450022/TensorYlmDBFull",
    "RegExp"            => "expanse",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "gpc" => {
    "AllowedCPN"        => [8],
    "BatchEnv"          => "PBS",
    "AltMail"           => 1,
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^(scinet|gpc-f1[0123]|gpc-f14[0123456])",
    "ScratchDir"        => \&LinkedScratchGPC,
    "SshName"           => "gpc01",
   },
  "graham" => {
    "AllowedCPN"        => [16],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Graham,
    "DefaultAccount"    => "rrg-pfeiffer-ac",
    "MaxRuntimeHours"   => 24,
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^gra.*",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "gravity" => {
    "AllowedCPN"        => [12],
    "BatchEnv"          => "PBS",
    "BatchHeader"       => \&HeaderPBS_ScinetGravity,
    "DefaultQueue"      => "gravity",
    "MaxRuntimeHours"   => "12",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^(gravity|gpc-f14[78])",
    "ScratchDir"        => \&LinkedScratchGPC,
    "SshName"           => "gravity01",
  },
  "guillimin" => {
    "AllowedCPN"        => [12],
    "BatchEnv"          => "PBS",
    "DefaultAccount"    => "jsh-251-ae",
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^(lg|sw|hb|lm|gm)-[0-9]r[0-9][0-9]-n[0-9][0-9]",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },
  "hydra" => {
    "AllowedCPN"        => [16,20],
    "BatchEnv"          => "LoadLeveler",
    "DefaultQueue"      => undef,
    "MaxRuntimeHours"   => "24", # not sure if this is really the max
    # SuperMUC's doc recommend this for Intel MPI
    "MpiCmd"            => \&Mpirun,
    # Hydra's doc only show examples for IBM MPI and poe,
    # SuperMUC's doc show this for IBM MPI but also mpiexec
    #"MpiCmd"            => \&Poe,
    "RegExp"            => '^(hydra\d\d|hy\d\d\d\d)$',
    "RestrictByQueue"   => \&RestrictByQueue_Hydra,
    "ScratchDir"        => \&LinkedScratch_Hydra,
    "SshName"           => "",
  },
  "jureca" => {
    "AllowedCPN"        => [24],
    "BatchEnv"          => "SBATCH",
    "DefaultQueue"      => "batch",
    "AltMail"           => 1,
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Srun,
    "RegExp"            => "^jr.*[.]jureca",
    "ScratchDir"        => \&LinkedScratch_Jureca,
  },
  "kepler" => {
    "AllowedCPN"        => [40],
    "BatchEnv"          => "PBS",
    "MaxRuntimeHours"   => "12",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^kepler",
    "ScratchDir"        => \&SimpleScratch,
  },
  "mamouth" => {
    "AllowedCPN"        => [24],
    "BatchEnv"          => "PBS",
    "MaxRuntimeHours"   => "120",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^ip|^cp",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "pfeiffer-mp2.rqchp.ca",
  },
  "maple" => {
    "AllowedCPN"        => [36],
    "BatchEnv"          => "PBS",
    "BatchHeader"       => \&HeaderPBS_Maple,
    "MaxRuntimeHours"   => "9999",
    "MpiCmd"            => \&Mpirun_Maple,
    "RegExp"            => "^maple|^cn[0-9][0-9]",
    "ScratchDir"        => \&SimpleScratch,
  },
  "minerva" => {
    "AllowedCPN"        => [16],
    "BatchEnv"          => "SBATCH",
    "DefaultQueue"      => "nr",
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Mpirun_Minerva,
    # when using OpenMPI
    #"MpiCmd"            => \&Mpirun,
    "RegExp"            => '^login0[12][.]cluster|^node\d\d\d[.]cluster',
    "ScratchDir"        => \&LinkedScratch_Minerva,
    "SshName"           => "",
  },
  "niagara" => {
    "AllowedCPN"        => [40],
    "BatchEnv"          => "SBATCH",
    "MaxRuntimeHours"   => 24,
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^nia.*",
    "ScratchDir"        => \&LinkedScratchNiagara,
    "SshName"           => "nia-login06",
  },
  "ocean" => {
    "AllowedCPN"        => [12],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Ocean,
    "DefaultQueue"      => "orca-0",
    "MaxRuntimeHours"   => "12",
    "MpiCmd"            => \&Mpirun,
    "RegExp"            => "^o(rca|cean)",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "ocean"
  },
  "pleiades" => {
     "AllowedCPN"        => [16],
     "BatchEnv"          => "PBS",
     "BatchHeader"       => \&HeaderPBS_Pleiades,
     "DefaultAccount"    => "s1996",
     "DefaultQueue"      => "long",
     "MaxRuntimeHours"   => "120",
     "MpiCmd"            => \&Mpirun_Pleiades,
     "RegExp"            => '^pfe.*|^r\d\d\di\dn\d.*',
     "ScratchDir"        => \&SimpleScratch,
     "SshName"           => "",
  },
  "stampede2" => {
    "AllowedCPN"        => [48],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Stampede2,
    "DefaultAccount"    => "TG-PHY990007N",
    "DefaultQueue"      => "skx-normal",
    "MaxRuntimeHours"   => "48",
    "MpiCmd"            => \&Ibrun,
    "RegExp"            => '(^stampede2|\S+\.stampede2\.)',
    "CanonicalTensorYlmDb" => "/work2/00207/ux450022/stampede2/TensorYlmDBFull",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "login1.stampede2.tacc.utexas.edu",
  },
  "trillian" => {
    "AllowedCPN"        => [32],
    "BatchEnv"          => "PBS",
    "BatchHeader"       => \&HeaderPBS_Trillian,
    "MaxRuntimeHours"   => "24",
    "MpiCmd"            => \&Aprun_Trillian,
    "RegExp"            => '(^trillian)',
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
  },

"rusty" => {
     "AllowedCPN"        => [128],
     "BatchEnv"          => "SBATCH",
     "BatchHeader"	=> \&HeaderSBATCH_Rusty,
     "DefaultAccount"    => "cca",
     "DefaultQueue"      => "cca",
     "MaxRuntimeHours"   => "24",
     "MpiCmd"            => \&Mpirun,
     "MpiThreadingCmd"   => \&Mpirun_threads,
     "RegExp"            => 'rusty.*|worker\d{4}',
     "ScratchDir"        => \&SimpleScratch,
     "SshName"           => "",
   },

  "wheeler" => {
    "AllowedCPN"        => [24],
    "BatchEnv"          => "SBATCH",
    "BatchHeader"       => \&HeaderSBATCH_Wheeler,
    "DefaultAccount"    => "sxs",
    "MaxRuntimeHours"   => 24,
    "MpiCmd"            => \&Mpirun,
    "MpiThreadingCmd"   => \&Mpirun_threads,
    "RegExp"            => "^wheeler",
    "ScratchDir"        => \&SimpleScratch,
    "SshName"           => "",
    # Note that the .sif file defined on the next line is to be
    # used by everyone on this machine.
    "SingularityCmd"    => "singularity exec " .
        "-B /usr/lib64:/wheeler/usr/lib64 -B /etc/libibverbs.d " .
        "/home/scheel/Singularity/specwheeler.sif",
  },

  # Default machine used if hostname doesn't match anything above.
  "local" => {
    "AllowedCPN"        => \&CoresPerNodeLocal,
    "BatchEnv"          => "none",
    "MaxRuntimeHours"   => "120",
    "MpiCmd"            => \&MpirunLocal,
    "MpiThreadingCmd"   => \&Mpirun_threads,
    "SshName"           => "",
    "ScratchDir"        => \&SimpleScratch,
  },
  );

  # Extract the database for this machine
  my $host = Hostname(\%FULLDATABASE);
  my %Database = %{$FULLDATABASE{$host}};
  $Database{"Host"} = $host;

  # Default sets of fields for batch environments.
  # The values can be changed by overriding these fields in FULLDATABASE.
  my %AllBatchEnvDefaults = (
  "PBS" => {
    "BatchHeader"       => \&HeaderPBS,
    "JobID"             => \&JobidPBS,
    "Jobname"           => \&JobnamePBS,
    "Nprocs"            => \&NprocsPBS,
    "SubmissionCmd"     => "qsub",
  },
  "SGE" => {
    "BatchHeader"       => \&HeaderSGE,
    "DefaultAccount"    => "", # special logic applies, let SGE decide
    "JobID"             => \&JobidSGE,
    "Jobname"           => \&JobnameSGE,
    "Nprocs"            => \&NprocsSGE,
    "SubmissionCmd"     => "qsub",
  },
  "SBATCH" => {
    "BatchHeader"       => \&HeaderSBATCH,
    "DefaultAccount"    => "", # special logic applies, let SBATCH decide
    "JobID"             => \&JobidSBATCH,
    "Jobname"           => \&JobnameSBATCH,
    "Nprocs"            => \&NprocsSBATCH,
    "SubmissionCmd"     => "sbatch",
  },
  "LoadLeveler" => {
    "BatchHeader"       => \&HeaderLoadLeveler,
    "DefaultAccount"    => "", # there are no accounts on Hydra yet
    "JobID"             => \&JobidLoadLeveler,
    "Jobname"           => \&JobnameLoadLeveler,
    "Nprocs"            => \&NprocsLoadLeveler,
    "SubmissionCmd"     => "llsubmit",
  },
  "none" => {
    "BatchHeader"       => \&HeaderLocal,
    "JobID"             => \&JobidDate,
    "Jobname"           => \&JobnameSGE,
    "Nprocs"            => \&NprocsSGE,
    "SubmissionCmd"     => "bash",
  },
  );

  # Apply recessive default fields for the given batch environment
  my %Defaults = %{$AllBatchEnvDefaults{$Database{"BatchEnv"}}};
  foreach my $key (keys %Defaults) {
    $Database{$key} = $Defaults{$key} if (not exists $Database{$key});
  }

  # Create the class
  my $self = \%Database;
  my $class = shift;
  bless($self,$class);
  return $self;
}

#==============================================================================
# PUBLIC METHODS
#==============================================================================

# Special function to return the standard resubmission command.
# INPUT ARGS:
#   WorkDir = string    # Path to run the submission command. [REQUIRED]
#   SubCmd = string     # Supplants the default submission command.
sub GetResubCmd {
  my ($self, %args) = @_;

  # Verify arguments and set defaults
  Utils::VerifyHashArgs(\%args,"SubCmd","WorkDir");
  die "Must specify 'WorkDir'" if (not exists $args{WorkDir});

  my $SubmitCmd = (exists $args{SubCmd}) ?
    $args{SubCmd} : $self->GetSubmissionCmd() . " ./Submit.sh";
  die "Submission command must not be empty" if (not $SubmitCmd);

  my $workdir = $args{WorkDir};
  my $SshName = $self->GetSshName();
  if (not $SshName) {
    # This handles machines whose compute nodes can't ssh to head nodes.
    return "cd $workdir && $SubmitCmd";
  } else {
    return "ssh $SshName \"cd $workdir && bash -l -c '$SubmitCmd'\"";
  }
}

sub GetHost {
  my ($self) = @_;
  return $self->{"Host"};
}
sub GetBatchHeader {
  my ($self, $cpn, $preempt) = @_;
  my $header = "#!/bin/bash -\n" . $self->{"BatchHeader"}->($cpn, $preempt);
  return $header;
}

# Returns a list of allowed CoresPerNode for the machine
sub GetAllowedCPN {
  my ($self) = @_;
  my $field = $self->{AllowedCPN};
  if (ref($field) eq "CODE") {
    return $field->();
  } else {
    return @{$field};
  }
}

sub GetDefaultQueue {
  my ($self) = @_;
  return (exists $self->{"DefaultQueue"})
      ? $self->{"DefaultQueue"} : undef;
}

sub GetSharedNodeQueue {
  my ($self) = @_;
  return (exists $self->{"SharedNodeQueue"})
      ? $self->{"SharedNodeQueue"} : "SpECNoSuchQueue";
}

sub GetSharedNodeMaxCores {
  my ($self) = @_;
  return (exists $self->{"SharedNodeMaxCores"})
      ? $self->{"SharedNodeMaxCores"} : undef;
}

sub GetTmpTensorYlmDbPath {
  my ($self) = @_;
  return (exists $self->{"TmpTensorYlmDbPath"})
      ? $self->{"TmpTensorYlmDbPath"} : undef;
}

sub GetCanonicalTensorYlmDb {
  my ($self) = @_;
  return (exists $self->{"CanonicalTensorYlmDb"})
      ? $self->{"CanonicalTensorYlmDb"} : undef;
}

# Returns the DefaultAccount (~/.SpEC has priority over Machines.pm)
sub GetDefaultAccount {
  my ($self) = @_;
  my $DotSpecAccount = DotSpecOptions::GetValue("DefaultAccount");
  return $DotSpecAccount ? $DotSpecAccount : $self->{DefaultAccount};
}
sub GetMaxRuntimeHours {
  my ($self) = @_;
  return $self->{"MaxRuntimeHours"};
}
sub GetMpiCmd {
  my ($self,$Nprocs) = @_;
  my $cmd = $self->{"MpiCmd"}->($Nprocs);
  # If necessary, we add 'singularity exec ...' to the end of the cmd.
  my $sing_cmd = GetSingularityCmdIfNeeded($self);
  if ($sing_cmd ne "") {
    $cmd .= " " . $sing_cmd;
  }
  return $cmd;
}
sub GetMpiThreadingCmd {
  my ($self,$Nnodes,$NthreadsPerNode) = @_;
  return $self->{"MpiThreadingCmd"}->($Nnodes,$NthreadsPerNode);
}
sub GetSubmissionCmd {
  my ($self) = @_;
  my $cmd = $self->{"SubmissionCmd"};
}

sub GetSingularityCmd {
  my($self) = @_;
  return (exists $self->{"SingularityCmd"}) ? $self->{"SingularityCmd"} : ""; 
}

# Returns the appropriate singularity command, generally
# "singularity exec some/file/name.sif".
# Returns the command only if singularity is requested by
# the use of "SPEC_USE_SINGULARITY_CONTAINER", and if we
# are not already in the container.  Otherwise returns an empty string.
sub GetSingularityCmdIfNeeded {
  my ($self) = @_;
  # Here we use an environment variable "SPEC_USE_SINGULARITY_CONTAINER"
  # that should be defined in this_machine.env if we want all executables
  # to be run inside of a singularity container when the executables are
  # invoked outside the container.
  if(defined $ENV{"SPEC_USE_SINGULARITY_CONTAINER"}) {
    # The environment variable "SINGULARITY_CONTAINER" is defined by
    # Singularity if we are inside a singularity container.
    if(not defined $ENV{"SINGULARITY_CONTAINER"}) {
      # We are outside the container,
      # so we want to invoke singularity.
      if(exists $self->{"SingularityCmd"}) {
        return $self->{"SingularityCmd"};
      } else {
        die "We chose to use the singularity container on this machine\n".
            "by setting SPEC_USE_SINGULARITY_CONTAINER in this_machine.env,\n".
            "yet there is no singularity command defined for this machine.";
      }
    }
  }
  # In all other cases (SPEC_USE_SINGULARITY_CONTAINER is undefined,
  # or we are inside the container already), we just return
  # an empty string.
  return "";
}
# default: Host (the hash key for this machine)
sub GetSshName {
  my ($self) = @_;
  return (exists $self->{"SshName"}) ? $self->{"SshName"} : $self->GetHost();
}
sub GetJobID {
  my ($self) = @_;
  return $self->{"JobID"}->();
}
sub GetJobname {
  my ($self) = @_;
  return $self->{"Jobname"}->();
}
sub GetScratchDir {
  my ($self,@opts) = @_;
  my $scratch = $self->{"ScratchDir"}->(@opts);
  Utils::MakePathsAbsolute(\$scratch);
  return $scratch;
}
sub GetNprocs {
  my ($self,@opts) = @_;
  return $self->{"Nprocs"}->(@opts);
}
# default: no-op
sub Backup {
  my ($self,@opts) = @_;
  $self->{"Backup"}->(@opts) if (exists $self->{"Backup"});
}

sub RestrictByQueue {
  my ($self, %data) = @_;
  if (exists $self->{"RestrictByQueue"}) {
    return $self->{"RestrictByQueue"}->(%data);
  } else {
    return %data;
  }
}

sub GetMasterNodeName {
  my ($self)=@_;
  return $self->{"MasterNodeName"}->();
}

sub GetNodeNameList {
  my ($self)=@_;
  return $self->{"NodeNameList"}->();
}

# Wrapper for Utils::SendMail, with variant if the machine doesn't have "mail"
sub SendMail {
  my ($self,%args) = @_;
  if (exists $self->{AltMail} && $self->{AltMail}==1) {
    $args{mailer} = "$RealBin/MailSpec";
    $args{host} = $self->GetSshName();
  } elsif (system("which mail >/dev/null 2>&1") != 0) {
    print "WARNING: mail is not installed, so SendMail may fail!"
  }
  Utils::SendMail(%args);
}

#==============================================================================
# PUBLIC FUNCTIONS
#==============================================================================

# Function to get a new node setup during resubmission
# Arguments to GetNodeSetup_Resub:
#   Nmin       = int        # minimum allowed cores
#   Nmax       = int        # maximum allowed cores
#   N          = int        # suggested number of cores
#   AllowedCPN = array ref  # list of allowed cores per node
#   N0         = int        # current number of cores
#   CPN0       = int        # current cores per node
#   Resubmit   = bool       # current resubmission status
#                           # (will try to avoid resubmission if False)
# Returns the following:
#   [0]        = int        # new number of cores
#   [1]        = int        # new number of cores per node
#   [2]        = bool       # new resubmission status
#   [3]        = string     # Queue to be used to override the submission
#                           # queue, or "SpECNoSuchQueue" for no override.
sub GetNodeSetup_Resub {
  my (%args) = @_;
  Utils::VerifyHashArgs(\%args, "Nmin", "Nmax", "N", "AllowedCPN",
                                "Resubmit", "N0", "CPN0", "ForceCPN");
  die "CPN0 $args{CPN0} must be in AllowedCPN @{$args{AllowedCPN}}"
    unless (grep {$_ == $args{CPN0}} @{$args{AllowedCPN}});

  # If not planning to resubmit, see if there are valid setups that don't
  # change the number of nodes. If not, will need to resubmit.
  if (not $args{Resubmit}) {
    my $Nodes0 = ceil($args{N0}/$args{CPN0});
    my $Nmin_NoResub = ($Nodes0-1)*$args{CPN0}+1;
    my $Nmax_NoResub = $Nodes0*$args{CPN0};
    my ($NewN, $NewCPN, $Queue) = GetNodeSetup(
      Nmin => max( $args{Nmin}, min($args{Nmax}, $Nmin_NoResub) ),
      Nmax => min( $args{Nmax}, max($args{Nmin}, $Nmax_NoResub) ),
      N    => $args{N},
      AllowedCPN => [$args{CPN0}],
    );
    if ($Queue eq "SpECNoSuchQueue") {
      # We are not in a shared queue situation, so we don't need to
      # resubmit if the new number of cores is in range.
      if (($NewN >= $Nmin_NoResub) and ($NewN <= $Nmax_NoResub)) {
        return ($NewN, $NewCPN, $args{Resubmit}, $Queue);
      }
    } elsif (($NewN == $args{N0}) and ($NewCPN == $args{CPN0})) {
      # We are in a shared queue situation, so we don't need to resubmit
      # if the number of cores hasn't changed.
      return ($NewN, $NewCPN, $args{Resubmit}, $Queue);
    } # else we must resubmit
  }

  # We will have to resubmit anyway to satisfy the (Nmin,Nmax)
  # criteria, or because we have changed the number of cores when in a
  # shared queue, so since we need to resubmit anyway we now can
  # change more properties of the setup.
  my @AllowedCPN = $args{ForceCPN} ? ($args{CPN0}) : @{$args{AllowedCPN}};
  my ($NewN, $NewCPN, $Queue) = GetNodeSetup(
    Nmin => $args{Nmin},
    Nmax => $args{Nmax},
    N    => $args{N},
    AllowedCPN => \@AllowedCPN,
  );

  return ( $NewN, $NewCPN, 2, $Queue );
}

sub GetOptimumCoresAndCpn {
  my($N,$Nmin,$Nmax,@AllowedCPN)=@_;
  my @Cores_base = ($Nmin..$Nmax);
  return (-1,-1) if (not @Cores_base);
  
  # Determine the number of unused cores for each (cpn,cores) pair.
  # Then find the minimum unused cores over all (cpn,cores) pairs.
  my %Unused = ();
  my $min_unused;
  {
    my $CountUnused = sub {
      my ($n, $cpn) = @_;
      return ceil($n/$cpn)*$cpn - $n;
    };
    foreach my $cpn (@AllowedCPN) {
      $Unused{$cpn} = { map { $_ => $CountUnused->($_,$cpn) } @Cores_base };
    }
    # Minimum value in a hash of hash refs
    $min_unused = min( map { min(values(%$_)) } values(%Unused) );
  }
  
  # Select the (cpn,cores) pairs that share the minimum unused cores.
  my %ValidCores = ();
  foreach my $cpn (keys %Unused) {
    my @valid = grep { $Unused{$cpn}{$_}==$min_unused } @Cores_base;
    $ValidCores{$cpn} = \@valid if (@valid);
  }

  # Restrict ValidCores based on information from the machine queue
  %ValidCores = new Machines->RestrictByQueue(%ValidCores);

  # Determine the distance from suggested "N" for valid (cpn,cores) pairs.
  # Then find the minimum distance over valid (cpn,cores) pairs.
  my %Distance = ();
  my $min_distance;
  {
    foreach my $cpn (keys %ValidCores) {
      next if (not @{$ValidCores{$cpn}});
      $Distance{$cpn} = { map { $_ => abs($_-$N) } @{$ValidCores{$cpn}} };
    }
    $min_distance = min( map { min(values(%$_)) } values(%Distance) );
  }

  # Select the (cpn,cores) pairs that have the minimum unused cores,
  # preferentially choosing such a pair that has larger N, then larger CPN.
  my $BestN = -1;
  my $BestCPN = -1;
  foreach my $cpn (keys %Distance) {
    foreach my $cores (keys %{$Distance{$cpn}}) {
      if ($Distance{$cpn}{$cores} == $min_distance) {
        if ($cores > $BestN || ($cores == $BestN and $cpn > $BestCPN)) {
          $BestN = $cores;
          $BestCPN = $cpn;
        }
      }
    }
  }

  # Sanity check
  die "Chose $BestN cores, outside range [$Nmin,$Nmax]"
    if ($BestN < $Nmin || $BestN > $Nmax);

  return ($BestN,$BestCPN);
}

# Function to get a new node setup
# Arguments to GetNodeSetup:
#   Nmin       = int        # minimum allowed cores
#   Nmax       = int        # maximum allowed cores
#   N          = int        # suggested number of cores
#   AllowedCPN = array ref  # list of allowed cores per node (optional)
# Returns the following:
#   [0]        = int        # new number of cores
#   [1]        = int        # new number of cores per node
#   [2]        = string     # If GetSharedNodeQueue is defined and the number
#                           # of cores fits on a single node, then the string
#                           # is the result of GetSharedNodeQueue.  Otherwise
#                           # it is the string "SpECNoSuchQueue".
sub GetNodeSetup {
  my (%args) = @_;
  Utils::VerifyHashArgs(\%args, "Nmin", "Nmax", "N", "AllowedCPN");
  die "Nmin must be an int>=0, not $args{Nmin}" if (not $args{Nmin} =~ /^\d+$/);
  die "Nmax must be an int>=0, not $args{Nmax}" if (not $args{Nmax} =~ /^\d+$/);
  die "N must be an int>0, not $args{N}" if (not $args{N} =~ /^\d+$/);
  die "ERR: Nmin=$args{Nmin} > Nmax=$args{Nmax}" if ($args{Nmin} > $args{Nmax});

  # Default AllowedCPN is the machine default
  my @AllowedCPN;
  if (exists $args{AllowedCPN}) {
    if (not ref($args{AllowedCPN}) eq 'ARRAY') {
      die "AllowedCPN must be array ref";
    }
    @AllowedCPN = @{$args{AllowedCPN}};
  } else {
    @AllowedCPN = new Machines()->GetAllowedCPN();
  }

  # If we have the option of running on a shared node, and if
  # suggested N is less than SharedNodeMaxCores (so that all the cores fit on
  # the node), then the algorithm is easy: just return the suggested N
  # and the correct AllowedCPN.
  my $SharedNodeQueue = new Machines()->GetSharedNodeQueue();
  if($SharedNodeQueue ne "SpECNoSuchQueue") {
    my $SharedNodeMaxCores = new Machines()->GetSharedNodeMaxCores();
    foreach my $cpn (@AllowedCPN) {
      if($args{N} < $cpn and $args{N} <= $SharedNodeMaxCores) {
        return ($args{N},$cpn,$SharedNodeQueue);
      }
    }
  }

  # Compute the best N and the best CPN, not counting any shared queue.
  my ($BestN, $BestCPN) = GetOptimumCoresAndCpn($args{N},$args{Nmin},
                                                $args{Nmax},@AllowedCPN);
  
  if($SharedNodeQueue ne "SpECNoSuchQueue") {
    # Repeat GetOptimumCoresAndCpn, but treat SharedNodeMaxCores as
    # an additional AllowedCPN to see if that produces a better match.  
    my $SharedNodeMaxCores = new Machines()->GetSharedNodeMaxCores();
    my ($optBestN, $optBestCPN) = GetOptimumCoresAndCpn($args{N},$args{Nmin},
                                                        $args{Nmax},
                                                        @AllowedCPN,
                                                        $SharedNodeMaxCores);
    # Use the shared queue if it has been selected, and if the number of cores
    # is less than the max number of cores for that queue.
    if($optBestCPN == $SharedNodeMaxCores and
       $optBestN <= $SharedNodeMaxCores) {
      # Use the shared queue, with $optBestN cores.
      # But CPN should be set to the actual CPN, not $optBestCPN which is
      # not really a valid CPN.  Choose the smallest CPN that is greater than
      # $SharedNodeMaxCores.  
      my $UseThisCPN = 100000000;
      foreach my $cpn (@AllowedCPN) {
        if($cpn < $UseThisCPN and $cpn > $SharedNodeMaxCores) {
          $UseThisCPN = $cpn;
        }
      }
      return ($optBestN, $UseThisCPN, $SharedNodeQueue);
    }
  }
  
  return ($BestN, $BestCPN, "SpECNoSuchQueue");
}

#==============================================================================
# PRIVATE FUNCTIONS
#==============================================================================

# Return the name of the host in the database that matches the real hostname
sub Hostname {
  my ($Database,$opt_v) = @_;
  unless(defined($opt_v)) { $opt_v=0; }   # Default verbosity

  my $FullHostname = Utils::Hostname();
  foreach my $host (keys %{$Database}) {
    next if ($host =~ "local");
    print STDERR "Looking at database entry: $host\n" if ($opt_v>0);
    my $regexp = $Database->{$host}{"RegExp"};
    if($FullHostname =~ m/$regexp/) {
      print STDERR "Found host $host ($FullHostname)\n" if ($opt_v>0);
      return $host;
    }
  }
  # Die if the host wasn't found
  print STDERR "Cannot find host for $FullHostname in database.\n".
               "Assuming host is a 'local' machine.\n" if ($opt_v>0);
  return "local";
}

#--------------------------------------------------------------------

sub HeaderLocal {
  my $header = <<END;
export NSLOTS="__cores__"   # Number of cores
export JOB_NAME="__SpEC__"  # Job Name
exec 1>>`pwd`/SpEC.stdout    # Output file name
exec 2>>`pwd`/SpEC.stderr    # Error file name
END
  return $header;
}

sub HeaderLoadLeveler {
  # see http://www-01.ibm.com/support/knowledgecenter/SSFJTW_5.1.0/com.ibm.cluster.loadl.v5r1.load100.doc/am2ug_jobkey.htm?lang=en
  my ($ppn) = @_;
  my $header = <<END;
# @ shell=/bin/bash                                                                                     #
# @ error = SpEC.stderr
# @ output = SpEC.stdout
# @ job_type = mpich # use parallel for poe
# @ job_name = __SpEC__
# @ node_usage= not_shared
# @ node = __nodes__
# @ tasks_per_node = ${ppn}
# @ resources = ConsumableCpus(1) # 1 thread per MPI rank
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit = __hours__:__minutes__:00
# @ notification = always
# @ queue
END
  # llsubmit does not like comments on the option lines so we remove them.
  $header = join("\n",map {s/(#[^#]*[^#[:space:]])\s*#.*/$1/;$_} split /\n/,$header);
  return $header;
}

sub HeaderSBATCH {
  my ($ppn) = @_;
  my $header = <<END;
#SBATCH -J __SpEC__                   # Job Name
#SBATCH -o SpEC.stdout                # Output file name
#SBATCH -e SpEC.stderr                # Error file name
#SBATCH -n __cores__                  # Number of cores
#SBATCH --ntasks-per-node $ppn        # number of MPI ranks per node
#SBATCH -p __Queue__                  # Queue name
#SBATCH -t __hours__:__minutes__:00   # Run time
#SBATCH -A __Account__                # Account name
#SBATCH --no-requeue
END
  return $header;
}

sub HeaderSBATCH_Comet {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --partition=compute\n";
  $header .= "#SBATCH --ntasks-per-node $ppn\n";
  return $header;
}

sub HeaderSBATCH_CaltechHPC {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  if($ppn == 32) {
      $header .= "#SBATCH --constraint=skylake\n";
  } elsif($ppn == 56) {
      $header .= "#SBATCH --constraint=cascadelake\n";
  } else {
      die "Unexpected ppn $ppn\n";
  }
  return $header;
}

sub HeaderSBATCH_Expanse {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --nodes __nodes__\n";
  return $header;
}

sub HeaderSBATCH_Ocean {
  my ($ppn) = @_;
  my $header = <<END;
#SBATCH -J __SpEC__                   # Job Name
#SBATCH -o SpEC.stdout                # Output file name
#SBATCH -e SpEC.stderr                # Error file name
#SBATCH -n __cores__                  # Number of cores
#SBATCH --ntasks-per-node $ppn        # number of MPI ranks per node
#SBATCH -p __Queue__                  # Queue name
#SBATCH -t __hours__:__minutes__:00   # Run time
END
  return $header;
}

sub HeaderSBATCH_Stampede2 {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --nodes __nodes__\n";
  return $header;
}

sub HeaderSBATCH_Frontera {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --nodes __nodes__\n";
  return $header;
}

sub HeaderSBATCH_Graham {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --mem-per-cpu=4096M\n";
  return $header;
}

sub HeaderSBATCH_Wheeler {
  my ($ppn, $preempt) = @_;
  my $header = HeaderSBATCH($ppn);
  if ($preempt == -1) {
    $header .= "#SBATCH -p sxs-low\n";
  } elsif ($preempt == 1) {
    $header .= "#SBATCH -p sxs-hi\n";
  } elsif ($preempt != 0) {
    Utils::Die("preempt=$preempt, but must be one of -1,0,1");
  }
  return $header;
}

sub HeaderSGE {
  my ($ppn) = @_;
  my $header = <<END;
#\$ -cwd                              # Start job in  submission directory
#\$ -l h_rt=__hours__:__minutes__:00  # Wall time (hh:mm:ss)
#\$ -pe ${ppn}way __cores__           # Number of cores.
#\$ -N __SpEC__                       # Job name
#\$ -o SpEC.stdout                    # Output file
#\$ -e SpEC.stderr                    # Output file
#\$ -q __Queue__                      # Which queue to run in.
#\$ -A __Account__                    # Which account to charge
#\$ -S /bin/bash                      # Run with bash
END
  return $header;
}

sub HeaderSGE_Datura {
  my ($ppn) = @_;
  my $header = HeaderSGE($ppn);
  $header =~ s/ \d+way / daturamon /;
  # qsub does not like comments on the option lines so we remove them.
  $header =~ s/^(#\$[^#]*[^#[:space:]])\s*#.*$/$1/mg;
  # ask for a reservation to avoid being constantly blocked by smaller jobs
  $header .= "#\$ -R yes\n";
  return $header;
}

sub HeaderPBS {
  my ($ppn) = @_;
  my $header = <<END;
#PBS -l nodes=__nodes__:ppn=$ppn
#PBS -l walltime=__hours__:__minutes__:00
#PBS -A __Account__
#PBS -q __Queue__
#PBS -N __SpEC__
#PBS -o SpEC.stdout
#PBS -e SpEC.stderr
#PBS -d .
#PBS -W umask=022
#PBS -S /bin/bash
END
  return $header;
}

sub HeaderPBS_Pleiades {
    my ($ppn) = @_;
    my $header = <<END;
#PBS -l select=__nodes__:ncpus=16:model=san
#PBS -l walltime=__hours__:__minutes__:00
#PBS -W group_list=__Account__
#PBS -q __Queue__
#PBS -N __SpEC__
#PBS -o SpEC.stdout
#PBS -e SpEC.stderr
#PBS -W umask=022
#PBS -S /bin/bash
END
  return $header;
}

sub HeaderPBS_ScinetGravity {
  my ($ppn) = @_;
  my $header = HeaderPBS($ppn);
  $header =~ s/(-l nodes=__nodes__:ppn=$ppn)/$1:gpus=2/;
  return $header;
}

sub HeaderPBS_BlueWaters {
  my ($ppn) = @_;
  my $header = HeaderPBS($ppn);
  $header =~ s/(-l nodes=__nodes__:ppn=[0-9]+)/$1:xe/;
  $header .= "source /etc/profile   # to define 'module'\n";
  return $header;
}

sub HeaderPBS_Trillian {
  my ($ppn) = @_;
  my $header = HeaderPBS($ppn);
  $header =~ s/(PBS -d)/##$1/;
  $header .= "source /etc/profile   # to define 'module'\n";
  $header .= "cd \"\$PBS_O_WORKDIR\" # move to submission directory\n";
  # The following two options provide a meaningful trace if the code crashes
  $header .= "export ATP_ENABLED=1\n";
  $header .= "export FOR_IGNORE_EXCEPTIONS=true\n";

  return $header;
}

sub HeaderPBS_Maple {
  my ($ppn) = @_;
  my $header = HeaderPBS($ppn);
  $header =~ s/(PBS -d)/##$1/;
  $header .= "source /etc/profile   # to define 'module'\n";
  $header .= "cd \"\$PBS_O_WORKDIR\" # move to submission directory\n";
  return $header;
}

#--------------------------------------------------------------------

# Get the number of cores on a local machine (i.e. 1 node machines)
sub CoresPerNodeLocal {
  my $ncpus = Utils::CoreCount();
  return ($ncpus);
}

sub RestrictByQueue_Hydra {
  my (%data) = @_;
  return DetermineQueue::RestrictForHydra(%data);
}

#--------------------------------------------------------------------

sub BackupStampede {
  my ($WRK, $SCR, $Email, $JobID) = @_;

  my $SCRATCH = Utils::GetEnvVar("SCRATCH");
  my $ARCHIVE = Utils::GetEnvVar("ARCHIVE");
  my $ARCHIVER = Utils::GetEnvVar("ARCHIVER");
  my $PATH = Utils::GetEnvVar("PATH");
  my $MACHINE = new Machines();
  my $ACCOUNT = $MACHINE->GetDefaultAccount();
  my $OPT_ACCOUNT = $ACCOUNT ? "-A $ACCOUNT" : "";

  # Backup dir is $SCR with the archive prefix
  (my $TarDir = $SCR) =~ s|$SCRATCH|$ARCHIVE|;
  my $TarCmd = "TarDirectory -a '$Email' -v2 -sz $SCR $ARCHIVER:$TarDir";
  Utils::OverwriteFile("$WRK/BACKUP.$JobID.sh",
                       "#!/bin/bash\nexport PATH=$PATH; $TarCmd");
  my $QsubCmd = "sbatch -D $WRK -J BACKUP.$JobID -n 1 -N 1 -p serial"
    . " $OPT_ACCOUNT -t 12:00:00 --begin=now+300 $WRK/BACKUP.$JobID.sh";
  system("ssh login1.stampede.tacc.utexas.edu '$QsubCmd'");
}

#--------------------------------------------------------------------

sub BackupStampede2 {
  my ($WRK, $SCR, $Email, $JobID) = @_;

  my $SCRATCH = Utils::GetEnvVar("SCRATCH");
  my $ARCHIVE = Utils::GetEnvVar("ARCHIVE");
  my $ARCHIVER = Utils::GetEnvVar("ARCHIVER");
  my $PATH = Utils::GetEnvVar("PATH");
  my $MACHINE = new Machines();
  my $ACCOUNT = $MACHINE->GetDefaultAccount();
  my $OPT_ACCOUNT = $ACCOUNT ? "-A $ACCOUNT" : "";

  # Backup dir is $SCR with the archive prefix
  (my $TarDir = $SCR) =~ s|$SCRATCH|$ARCHIVE|;
  my $TarCmd = "TarDirectory -a '$Email' -v2 -sz $SCR $ARCHIVER:$TarDir";
  Utils::OverwriteFile("$WRK/BACKUP.$JobID.sh",
                       "#!/bin/bash\nexport PATH=$PATH; $TarCmd");
  my $QsubCmd = "sbatch -D $WRK -J BACKUP.$JobID -n 1 -N 1 -p skx-normal"
    . " $OPT_ACCOUNT -t 12:00:00 --begin=now+300 $WRK/BACKUP.$JobID.sh";
  system("ssh login1.stampede2.tacc.utexas.edu '$QsubCmd'");
}

#--------------------------------------------------------------------

sub SimpleScratch {
  my ($Work,$JOBID) = @_;
  return "$Work/Run";
}

# Choose the scratch directory to be the same as the submission directory
# except for the beginning of the path, and with '$JOBID' on the end.
sub LinkedScratch_Hydra {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  if ($Work =~ m|^/hydra/u/|) {
    ($Scratch = $Work) =~ s|^/hydra/u/|/hydra/ptmp/|;
  } elsif ($Work =~ m|^/hydra/ptmp/|) {
    $Scratch = $Work;
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

# Choose the scratch directory to be the same as the submission directory
# except for the beginning of the path, and with '$JOBID' on the end.
sub LinkedScratch_Jureca {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  if ($Work =~ m|^/home./|) {
    ($Scratch = $Work) =~ s|^/home./|/work/|;
  } elsif ($Work =~ m|^/work/|) {
    $Scratch = $Work;
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}
sub HeaderSBATCH_Rusty {
  my ($ppn) = @_;
  my $header = HeaderSBATCH($ppn);
  $header .= "#SBATCH --constraint=ib\n";
  return $header;
}
sub LinkedScratch_Minerva {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  if ($Work =~ m|^/home/|) {
    ($Scratch = $Work) =~ s|^/home/|/scratch/|;
  } elsif ($Work =~ m|^/work/|) {
    ($Scratch = $Work) =~ s|^/work/|/scratch/|;
  } elsif ($Work =~ m|^/scratch/|) {
    $Scratch = $Work;
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

# Choose the scratch directory to be the same as the submission directory
# except for the beginning of the path, and with '$JOBID' on the end.
sub LinkedScratchTACC {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  if ($Work =~ m|^/work/|) {
    ($Scratch = $Work) =~ s|^/work/|/scratch/|;
  } elsif ($Work =~ m|^/home1/|) {
    ($Scratch = $Work) =~ s|^/home1/|/scratch/|;
  } elsif ($Work =~ m|^/scratch/|) {
    $Scratch = $Work;
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

sub LinkedScratchGPC {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  # this captures both scratch and project-space.  Note that the
  # evolution goes into the same space where the input files are
  # prepared.  (scratch -> scratch; project -> project).
  if ($Work =~ m|^/sgfs1/|) {
    $Scratch = $Work;
# This was an earlier variant from Mark Scheel, who prepared input
# files on project, but sent the Run/ directories onto scratch.
# Dangerous, as it opens the possibility that a user assumes the data
# is safe, because the input files were prepared on project-space.
# } elsif ($Work =~ m|^/project/pfeiffer/|) {
#    ($Scratch = $Work) =~ s|^/project/pfeiffer/|/oldscratch/|;
  } elsif ($Work =~ m|^/home|) {
    die "You are submitting a job from /home, which is not writable from "
       ."batch jobs on this system. Submit from scratch instead.\n";
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

sub LinkedScratchNiagara {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  # this captures scratch space.  Note that the
  # evolution goes into the same space where the input files are
  # prepared.  (scratch -> scratch).
  if ($Work =~ m|^/gpfs/fs0/|) {
    $Scratch = $Work;
  } elsif ($Work =~ m|^/home|) {
    die "You are submitting a job from /home, which is not writable from "
       ."batch jobs on this system. Submit from scratch instead.\n";
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

# Choose the scratch directory to be the same as the submission directory
# except for the beginning of the path, and with '$JOBID' on the end.
sub LinkedScratchBridges {
  my ($Work,$JOBID) = @_;

  my $Scratch;
  if ($Work =~ m|^/pylon2/|) {
    ($Scratch = $Work) =~ s|^/pylon2/|/pylon5/|;
  } elsif ($Work =~ m|^/pylon5/|) {
    $Scratch = $Work;
  } else {
    die "Cannot figure out scratch directory from $Work";
  }
  return "$Scratch/$JOBID";
}

sub MasterNodeName_Slurm {
  my $host = `hostname`;
  chomp $host;
  return $host;
}

sub NodeNameList_Slurm {
  my @nodenamelist;
  foreach my $host (`srun hostname | sort | uniq`) {
    chomp $host;
    push(@nodenamelist,$host);
  }
  return @nodenamelist;
}

#--------------------------------------------------------------------

sub Run {
  my ($Nprocs) = @_;
  warn "Only running on 1 proc" if ($Nprocs != 1);
  return "";
}

sub Ibrun {
  my ($Nprocs) = @_;
  return "ibrun -np $Nprocs";
}

sub Mpirun {
  my ($Nprocs) = @_;
  return "mpirun -np $Nprocs";
}
sub Mpirun_threads {
  my ($Nnodes,$NthreadsPerNode) = @_;
  my $pernode = "--map-by ppr:1:node:pe=$NthreadsPerNode";
  # --map-by ppr:1:node is the correct option since OpenMPI 1.7.5.
  # For older MPI these may be correct instead:
  my $help_opts = "--help >/dev/null 2>&1";
  $pernode = '--pernode' unless system("mpirun $pernode $help_opts") == 0;
  $pernode = '-ppn 1' unless system("mpirun $pernode $help_opts") == 0;
  return "mpirun -np $Nnodes $pernode";
}

# OpenMPI v3+ changed the default to disallow oversubscribing of processes on
# the machine. This makes it hard to use MPI in tests, because a user may have
# a machine with a low number of procs, or may want to run tests in parallel.
# This returns to the old OpenMPI behavior by overriding the new default.
sub MpirunLocal {
  my ($Nprocs) = @_;
  # Check if MPI is OpenMPI v3+
  my $VersionFirstLine = ( split /\n/, `mpirun --version`)[0];
  if ($VersionFirstLine =~ /mpirun \(Open MPI\)/) {
    my $MajorVersion = ($VersionFirstLine =~ /(\d+)/)[0];
    if ($MajorVersion >= 3) {
      return "mpirun --oversubscribe -np $Nprocs";
    }
  }
  # For OpenMPI v1 or v2, or for MPICH:
  return "mpirun -np $Nprocs";
}

sub Srun {
  my ($Nprocs) = @_;
  return "srun --ntasks $Nprocs";
}

sub Mpirun_Maple {
  my ($Nprocs) = @_;
  return "mpirun -genv MV2_ENABLE_AFFINITY 0 -n $Nprocs";
}

sub Mpirun_Pleiades{
    my ($Nprocs) = @_;
    return "mpiexec -n $Nprocs";
}

sub Mpirun_Minerva {
  my ($Nprocs) = @_;
  return "mpiexec.hydra -genv MALLOC_MMAP_MAX_=0 -genv MALLOC_TRIM_THRESHOLD_=536870912 -n $Nprocs";
}

sub Aprun {
  my ($Nprocs) = @_;
  # The -q flag suppresses output from aprun that we don't want (or need).
  # The "-j 1" flag activates one "processor" (i.e. integer core) per Bulldozer
  # floating-point unit. This means each MPI process has full use of a FP unit.
  return "aprun -q -j 1 -n $Nprocs";
}
sub Aprun_threads {
  my ($Nnodes,$NthreadsPerNode) = @_;
  return "aprun -q -j 1 -n $Nnodes -N 1 -d $NthreadsPerNode";
}

sub Aprun_Trillian {
  my ($Nprocs) = @_;
  return "aprun -q -n $Nprocs";
}

sub Poe {
  my ($Nprocs) = @_;
  return "MP_PROCS=$Nprocs /usr/bin/poe -procs";
}

#--------------------------------------------------------------------

sub NprocsLoadLeveler {
  my $NodeFile = Utils::GetEnvVar("LOADL_HOSTFILE");
  my $Nprocs = `cat $NodeFile | wc -l`;
  chomp $Nprocs;
  return $Nprocs;
}

sub NprocsSBATCH {
  return Utils::GetEnvVar("SLURM_NTASKS");
}

sub NprocsPBS {
  my $NodeFile = Utils::GetEnvVar("PBS_NODEFILE");
  my $Nprocs = `cat $NodeFile | wc -l`;
  chomp $Nprocs;
  return $Nprocs;
}

sub NprocsSGE {
  return Utils::GetEnvVar("NSLOTS");
}

#--------------------------------------------------------------------

sub JobidLoadLeveler {
  # from http://www-01.ibm.com/support/knowledgecenter/SSFJTW_5.1.0/com.ibm.cluster.loadl.v5r1.load100.doc/am2ug_envvars.htm?lang=en
  return Utils::GetEnvVar("LOADL_JOB_NAME");
}

sub JobidSBATCH {
  return Utils::GetEnvVar("SLURM_JOB_ID");
}

sub JobidSGE {
  return Utils::GetEnvVar("JOB_ID");
}

sub JobidPBS {
  my $id = Utils::GetEnvVar("PBS_JOBID");
  $id =~ s|\..*||;
  return $id;
}

# seconds since the last epoch
sub JobidDate {
  my $id = `date +%s`;
  chomp $id;
  return $id;
}

#--------------------------------------------------------------------

sub JobnameLoadLeveler {
  # from http://www-01.ibm.com/support/knowledgecenter/SSFJTW_5.1.0/com.ibm.cluster.loadl.v5r1.load100.doc/am2ug_envvars.htm?lang=en
  # this is actually the same as the job id
  return Utils::GetEnvVar("LOADL_JOB_NAME");
}

sub JobnameSBATCH {
  return Utils::GetEnvVar("SLURM_JOB_NAME");
}

sub JobnameSGE {
  return Utils::GetEnvVar("JOB_NAME");
}

sub JobnamePBS {
  return Utils::GetEnvVar("PBS_JOBNAME");
}

#--------------------------------------------------------------------
1; # This is how to end sources
