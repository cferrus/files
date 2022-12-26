#!/usr/bin/env perl
package Utils;  # Start the utils namespace

require 5;
use warnings FATAL => 'all';
use strict;
use File::Basename;
use File::Spec qw(rel2abs);
use Carp;
use File::Copy;
use Cwd;

#-------------------------------------------------
# Die with a backtrace
sub Die {
  my($error)=@_;
  print STDERR "\n*********** ERROR ***********\n";
  $Carp::MaxArgLen=0;      # These settings make all characters of all
  $Carp::MaxArgNums=undef; # function arguments show up in a backtrace.
  Carp::cluck();
  die "$error";
}

#-------------------------------------------------
# Return the hostname of the system
sub Hostname {
  my $hostname = lc(`hostname`);
  chomp $hostname;

  my ($name, $aliases, $addrtype,
          $length, @addrs) = gethostbyname $hostname;
  if($name) {
    # out of all our names, find one that has a dot, if possible.
    my @names = ($name, $hostname, split / /,$aliases);
    my @dottednames = grep {not m/^localhost(\.|$)/} ( grep {m/\./} @names );
    $hostname = $dottednames[0] if scalar @dottednames;
  }

  return $hostname;
}

#-------------------------------------------------
# Return ENV hash from a subshell executing a single-line command
# ENV inherits from the current ENV, but is overwritten/added to by the command
# e.g. %ENV = Utils::GetEnvFromSubshell(". this_machine.env");
sub GetEnvFromSubshell {
  my ($cmd, $opt_v) = @_;
  unless(defined $opt_v) { $opt_v=0; }   # Default verbosity

  # quote escaping gymnastics ensue...
  my $PrintEnv = 'while(($key,$value)=each(%ENV)){ print "key=$key\nvalue=$value\0"; }';

  my $FullCmd = "$cmd && perl -e '\\''$PrintEnv'\\''";
  # Must execute in bash (many commands assume it, e.g. module)
  my $env = SystemOutput("bash -c '$FullCmd'", $opt_v);

  my %NEWENV;
  my @env_array = split('\0', $env);
  foreach my $var (@env_array) {
    $var =~ m/key=(.*)\nvalue=(.*)$/s;
    my $key = $1;
    my $value = $2;
    Die("Regexp failed for $var\n") if (!defined($key) or !defined($value));
    $NEWENV{$key} = $value;
  }
  return %NEWENV;
}

#-------------------------------------------------
# Safely get a variable from the environment
sub GetEnvVar {
  my ($var) = @_;
  die "Environment variable $var undefined" unless (exists $ENV{$var});
  my $value = $ENV{$var}; chomp $value;
  return $value;
}

#-------------------------------------------------
# Runs a system command. 
sub System {
  my($command,$opt_v,$die_on_failure)=@_;
  unless(defined $opt_v) { $opt_v=0; }   # Default verbosity
  unless(defined $die_on_failure) { $die_on_failure=1; } # Default
  print STDERR "$command\n" if $opt_v>0;

  if(system($command) != 0) {
    my $errno = $!;
    my $FailInfo = "Failed system command '$command' ";
    if ($? == -1) {
      $FailInfo .= "...failed to execute: $errno\n";
    }
    elsif ($? & 127) {
      $FailInfo .= sprintf("...died with signal %d, %s coredump\n",
                           ($? & 127), ($? & 128) ? 'with' : 'without');
    }
    else {
      $FailInfo .= sprintf("...exited with value %d\n", $? >> 8);
    }
    if($die_on_failure) {
      Die($FailInfo);
    } else {
      print($FailInfo);
      return 0;
    }
  }
  return 1;
}

#----------------------------------------
# Tries a system command multiple times.
sub TrySystem {
  my($cmd,$num_iters,$sleep_value,$opt_v)=@_;
  foreach my $iter (0..$num_iters) {
    if(System($cmd,$opt_v,0)) {
      return 1;
    }
    print("Command $cmd failed iter $iter, sleeping and trying again.");
    sleep($sleep_value);
  }
  Die("Command $cmd failed $num_iters times. Aborting.\n");
}

#----------------------------------------
# Runs an ssh command on $host.
sub SystemSsh {
  my($cmd,$host,$opt_v,$die_on_failure)=@_;
  my $command = "ssh -o StrictHostKeyChecking=no $host \"bash -l -c '$cmd'\"";
  return System($command,$opt_v,$die_on_failure);
}

#----------------------------------------
# Tries an ssh command on $host multiple times.
sub TrySystemSsh {
  my($cmd,$host,$num_iters,$sleep_value,$opt_v)=@_;
  foreach my $iter (0..$num_iters) {
    if(SystemSsh($cmd,$host,$opt_v,0)) {
      return 1;
    }
    print("Command $cmd failed iter $iter, sleeping and trying again.");
    sleep($sleep_value);
  }
  Die("Command $cmd failed $num_iters times. Aborting.\n");
}

#-------------------------------------------------
# Runs a system command, capturing output. Dies on failure
sub SystemOutput {
  my($cmd,$opt_v) = @_;
  unless(defined $opt_v) { $opt_v=0; }   # Default verbosity
  print STDERR "$cmd\n" if $opt_v>0;

  open(SYS,"$cmd |") || Die("Cannot open command '$cmd': $!");
  if (wantarray) {
    # Put command output into separate elements of an array.
    my @text = <SYS>;
    close(SYS) || Die("Cannot close command '$cmd': $!. Command returned text '@text'");
    return @text;
  } else {
    # Put command output into one string.
    my $text = do {local $/; <SYS>}; # Slurp in the entire pipe.
    close(SYS) || Die("Cannot close command '$cmd': $!. Command returned text '$text'");
    return $text;
  }
}

#-------------------------------------------------
# Count the number of cores on this machine (if multi-node, will return
# the number of cores per node)
sub CoreCount {
  my $ncpus;
  if (-e "/proc/cpuinfo") {
    # Linux
    $ncpus = SystemOutput("cat /proc/cpuinfo|grep processor|wc -l");
  } elsif (-e "/usr/sbin/sysctl") {
    # Mac OSX
    my @ncpus = map{m/(\d+)/} SystemOutput("sysctl hw.ncpu");
    $ncpus = $ncpus[0];
  } else {
    Die("Don't know how to count cores on this machine!");
  }
  chomp($ncpus);
  return $ncpus;
}

#-------------------------------------------------
# This is the code that makes libs absolute in perl scripts.
# Might only be necessary to allow compiled script mobility.
sub MakeLibsPathAbsolute {
 for (@INC) {
   if (! ref && -d && !File::Spec->file_name_is_absolute($_)) {
     # Should be ok not to use tilde expansion here.
     $_ = File::Spec->rel2abs($_);
   }
 }
}

#-------------------------------------------------
# Given a set of references to strings that contain pathnames that are
# possibly relative, replace them with absolute pathnames.
# Does tilde expansion to usernames if the string begins with a tilde.
# Does not resolve symbolic links (due to use of rel2abs)
sub MakePathsAbsolute {
  my(@pathrefs)=@_;

  foreach my $pathref (@pathrefs) {
    next if (!defined($$pathref) || $$pathref eq '');

    if($$pathref =~ /^~/) {
      # glob does tilde expansion
      ( $$pathref ) = glob( $$pathref );
    } else {
      # Otherwise run rel2abs (which does not do tilde expansion)
      $$pathref = File::Spec->rel2abs($$pathref);
    }
    
  }
}

#-------------------------------------------------
# Link to the canonical path of $src (following symlinks), relative to $dest.
sub LinkToOriginal {
  my($src,$dest)=@_;
  my $abs_src = Cwd::abs_path($src);
  my $destdir = dirname($dest);
  $src = File::Spec->abs2rel($abs_src, $destdir);
  MySymlink($src, $dest);
}

#-------------------------------------------------
# Create a symlink or die trying
sub MySymlink {
  my($src,$dest)=@_;
  symlink($src, $dest) || Die("Cannot link '$src' to '$dest': $!");
}

#-------------------------------------------------
# Copies, preserving permissions.
sub MyCopy {
  my($src,$dest)=@_;

  # Because we change permissions, we prohibit the (file,dir) syntax
  Die("Destination must be a file, not a directory") if (-d $dest);

  copy($src,$dest) or Die("Cannot copy '$src' to '$dest', $!");
  # It is really annoying that File::copy doesn't preserve permissions!
  # So need to do this by hand...making sure that the files are at least user
  # writable
  my $mode = ((stat($src))[2] & 07777) | 00200;
  chmod($mode, $dest) || Die("Cannot chmod '$dest': $!");

  return 1;
}

#-------------------------------------------------
# Reads a file and returns contents (either in single string or array)
sub ReadFile {
  my($file)=@_;

  open(FILE,"<$file") or Die("Cannot open '$file': $!");
  if (wantarray) {
    # Reads each line of a file into separate elements of an array.
    my @text = <FILE>;
    close(FILE) or Die("Cannot close '$file': $!");
    return @text;
  } else {
    # Reads entire file into one string.
    my $text = do {local $/; <FILE>}; # Slurp in the entire file.
    close(FILE) or Die("Cannot close '$file': $!");
    return $text;
  }
}

#-------------------------------------------------
# Writes text to a file, overwriting current contents.
sub OverwriteFile {
  my($file,$text)=@_;

  open(FILE,">$file") || Die("Cannot open '$file' for writing: $!");
  print FILE $text;
  close(FILE) || Die("Cannot close '$file'");
}

#-------------------------------------------------
# Compares two files and returns if they are equal
sub CompareFiles {
  my ($f1, $f2) = @_;

  my $file_passed;
  my $DiffFile = "$f1.diff";  # default diff file name

  if ($f1 =~ /\.h5$/) {
    # Use 'h5diff' for .h5 files
    system("hash h5diff &>/dev/null")==0 || Die("Cannot find 'h5diff' in PATH.");
    # Ticket 448 workaround: check if file comparison fails, and if it does,
    # then write the h5diff report separately
    my $output = `h5diff $f1 $f2 2>&1`; my $result = $?;
    $file_passed = ($result==0 && $output eq "");
    system("h5diff -v -c -n 17 -r $f1 $f2 > $DiffFile 2>&1");
  } elsif ($f1 =~ /\.png$/ && system("hash compare 2>&-")==0) {
    # Use 'compare' from ImageMagick for .png files, if it is in PATH
    $DiffFile = dirname($f1)."/".basename($f1,".png").".diff.png";
    my $DiffCmd = '[ $(compare -metric AE -compose src %s %s %s 2>&1) -eq 0 ]';
    my $RealDiffCmd = sprintf($DiffCmd, $f1, $f2, $DiffFile);
    $file_passed = (system($RealDiffCmd) == 0);
  } else {
    # Use whitespace-insensitive diff as a fallthrough
    $file_passed = (system("diff -b -u $f1 $f2 > $DiffFile 2>&1") == 0);
  }

  if ($file_passed) {
    unlink($DiffFile) || Die("Cannot remove $DiffFile");
  }
  return $file_passed;
}

#-------------------------------------------------
# Writes text to a file, appending to current contents.
sub AppendToFile {
  my($file,$text)=@_;

  open(FILE,">>$file") || Die("Cannot open '$file' for appending: $!");
  print FILE $text;
  close(FILE) || Die("Cannot close '$file'");
}

#-------------------------------------------------
# Send e-mail. Dies on failure. Arguments are as follows:
#   subject     = subject of e-mail
#   body        = body of e-mail
#   to          = comma-delimited string of e-mail addresses
#   host        = ssh to this machine to mail (default undef)
#   mailer      = run this exec as the mailer (default "mail")
#   sleep       = sleep time in seconds after e-mail is sent (default 30)
#   v           = verbosity (default 0)
sub SendMail {
  my (%args) = @_;
  VerifyHashArgs(\%args,"subject","body","to","host","mailer","sleep","v");
  $args{v}=0 unless (exists $args{v});
  $args{mailer}="mail" unless (exists $args{mailer});
  $args{sleep}=30 unless (exists $args{sleep});
  System("which $args{mailer} >/dev/null 2>&1");   # Verify mailer is in $PATH

  # Check that email addresses are of the form 'string@string'
  my @ValidEmails;
  foreach my $address (split(/,/,$args{to})) {
    if ($address =~ /^\S+\@\S+$/) {
      push(@ValidEmails, $address);
    } else {
      warn "Email address '$address' is invalid!";
    }
  }
  Die("No valid email addresses specified!\n") unless (@ValidEmails);

  # Send the mail (always escape single quotes before single quoting!)
  (my $subject = $args{subject}) =~ s/'/'\\''/g;
  my $cmd = "$args{mailer} -s '$subject' @ValidEmails";
  if ($args{host}) {
    $cmd =~ s/'/'\\''/g;
    $cmd = "ssh $args{host} '$cmd'";
  }
  open(MAIL, "| $cmd") || Die("Cannot open pipe to mail: $!\n");
  print(MAIL $args{body});
  close(MAIL) || Die("Cannot send mail\n");

  print STDERR "Sent email to: @ValidEmails\n" if ($args{v} > 0);
  print STDERR "echo \$Body | $cmd\n" if ($args{v} > 1);

  # The mail may not be sent if the environment is wiped before the
  # mail daemon has time to send (since it occurs asynchronously).
  print STDERR "Sleeping for $args{sleep} seconds to give the mailer " .
               "time to send.\n" if ($args{v} > 1);
  sleep $args{sleep};
}

#-------------------------------------------------
# Copy source file to destination, making substitutions
# indicated. Optionally require all substitutions.
# This used to live in DoMultipleRuns.
sub CopyAndSubInFile {
    my($src,$dest,$vars,$require_replacement)=@_;

    # This variable is for a sanity check below
    # We will record whether each key was found and changed at least once.
    my @MissingKeys;

    # Read file and do modifications, all at once
    # Warn if a key is replaced by an undefined variable
    my $contents = ReadFile($src);
    while (my ($key, $value) = each %$vars) {
        unless(defined($value)) {
            die "Replacement for '$key' in '$src' is an undefined variable!\n";
        }
        unless($contents =~ s/$key/$value/g) { # Replace $key by $vars->{$key}
            push(@MissingKeys, $key); # Log lack of replacement
        }
    }

    # Now do sanity check
    # Make sure each key was replaced at least once in the file
    if ($require_replacement and @MissingKeys) {
      die "Keys '@MissingKeys' not found in $src\n";
    }

    OverwriteFile($dest, $contents);
    chmod((stat($src))[2] & 07777, $dest) or die "Cannot chmod '$dest'";
}

#-------------------------------------------------
# Returns the complement of array (elements in array, but NOT in exclude)
# Note that array is a ref because you can't have two '@' args to a sub
# WARNING: an element will be excluded if it's in array twice!
sub Complement {
  my ($array, @exclude) = @_;

  my %count;
  foreach my $e (@$array, @exclude) { $count{$e}++ }

  my @comp;
  foreach my $e (keys %count) {
    if ($count{$e} == 1) {
      # Element must be in @array to be in the result
      push(@comp, $e) if (grep {$_ eq $e} @$array);
    }
  }
  return @comp;
}

#-------------------------------------------------
# Check that the arguments passed to a function via a hash do not
# contain superfluous keys, which could imply a typo in the input.
sub VerifyHashArgs {
  my ($args, @AllowedKeys) = @_;
  Die("args must be a hash") if (not ref($args) eq "HASH");
  my @BadKeys;
  foreach my $key (keys %{$args}) {
    push(@BadKeys,$key) if (not grep {/^$key$/} @AllowedKeys);
  }
  if (@BadKeys) {
    my $sub = (caller(1))[3];
    $sub = $sub ? $sub : (caller(0))[3];
    Die("Unexpected args '@BadKeys' to $sub. Allowed args: @AllowedKeys.\n");
  }
}

#-------------------------------------------------
# Gets Perl code from a file and evaluates it.
# Note that any variables defined in the file become package
# variables in the Utils:: namespace.
sub EvalPerlCodeInFile {
  my($file)=@_;   # Gets the local variable $file from passed-in argument
  local $_;       # Makes $_ a local instead of a global var in this scope.

  # Test syntax of the file (e.g. DoMultipleRuns.input)
  Die("Cannot find $file\n") unless (-f $file);
  my $testsyntax=`perl -wc $file 2>&1`;
  Die("Perl does not understand $file:\n$testsyntax") if($?>>8);

  my $code="no strict;"; # Don't require strictness is user input code
  $code .= ReadFile($file);

  # Two-stage evaluation avoids absurdly verbose error messages
  my $sub = eval "sub { $code }"; # eval compiles and executes $code as perl.

  Die($@) if $@; # Catch error messages from $code
  $sub->();
}

#-------------------------------------------------
# Finds files matching filepattern recursively starting in dir.
# Returns an array of filepaths with "dir" prepended.
# filepattern can not contain "/" characters.
sub Find {
  my($filepattern,$dir)=@_;

  die "Filepattern should not contain '/' characters" if($filepattern =~m|/|);

  opendir(DIR,$dir) || die "Cannot opendir $dir\n";
  my(@entry)=readdir(DIR);
  closedir(DIR);

  my @result;
  foreach my $f (@entry) {
    next if($f eq '.' || $f eq '..');
    my $fullpath = "$dir/$f";
    if(-d $fullpath) { # does NOT skip links
      push(@result,Find($filepattern,$fullpath));
    } elsif(-f $fullpath && $f =~ m|^${filepattern}$|) {
      push(@result,$fullpath);
    }
  }
  return @result;
}

# Returns 1 if $file is a non-empty directory, 0 otherwise.
sub IsNonEmptyDirectory {
  my ($file)=@_;
  return 0 if not -e $file; # Does not exist.
  return 0 if not -d $file; # Is not a directory.
  opendir(my $dh, $file) or return 0; # Cannot open directory (permissions?)
  # Here we check if three successive 'readdir' calls succeed, indicating
  # that there is at least one nontrivial entry in the directory.
  readdir $dh; # Returns '.', which always exists, even in empty directories.
  readdir $dh; # Returns '..', which always exists, even in empty directories.
  return 1 if (readdir $dh); # Reads a file/directory if one exists.
  return 0; # No files in the directory.
}

# Returns 1 if $file is a readable directory, 0 otherwise.
sub IsReadableDirectory {
  my ($file)=@_;
  return 0 if not -e $file; # Does not exist.
  return 0 if not -d $file; # Is not a directory.
  opendir(my $dh, $file) or return 0; # Cannot open directory (permissions?)
  # Here we check if two successive 'readdir' calls succeed, indicating
  # that '.' and '..' exist in the directory. This is a sanity check.
  readdir $dh; # Returns '.', which always exists, even in empty directories.
  return 1 if (readdir $dh); # Reads '..' which always exists, even if empty.
  return 0; # Cannot read '.' or '..', so something is wrong.
}

#-------------------------------------------------
# Copied from List::MoreUtils (a cpan module)
sub uniq (@) {
  my %seen = ();
  grep { not $seen{$_}++ } @_;
}

#-------------------------------------------------
1; # This is how to end sources
