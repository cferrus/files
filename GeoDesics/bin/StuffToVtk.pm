#!/usr/bin/env perl
use warnings;
package StuffToVtk;  # Start the StuffToVtk namespace
use Exporter;        # load Exporter module
@ISA=qw(Exporter);   # Inherit from Exporter
@EXPORT=qw(CopyToRemoteDir FindMatchDirRecursive);

require 5;
use strict;
use Cwd;
use File::Basename;

# Load Utils module
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;

#-------------------------------------------------
# For remote copying.
# Assumes currentdir is of the form /some/path/dddd/more/path
# and puts the /dddd/more/path/ part after $remotebasedir.
# If it doesnt find this form, then it looks for
# /some/path/username/more/path, and puts the /more/path part
# after $remotebasedir
sub GetRemoteDir {
  my($currentdir,$remotedest)=@_;

  my ($REMOTEHOST,$remotebasedir) = GetRemoteHostAndBaseDir();
  if(defined $remotedest) {
    return $remotebasedir . "/" . $remotedest;
  }

  my $remotedir;
  my $user = `whoami`; chomp $user;
  if($currentdir =~ m|.*(/\d{4}/.*)$|) {
    $remotedir = $remotebasedir . $1;
  } elsif($currentdir =~ m|.*/$user(/.*)$|) {
    $remotedir = $remotebasedir . $1;
  } else {
    Utils::Die("Cannot construct remote dir for source dir " .
        "'$currentdir'\n" .
        "I assume currentdir is of the form /some/path/dddd/more/path\n".
        "or /some/path/username/more/path\n".
        "Perhaps the current dir does not match the expected form.\n");
  }
  return $remotedir;
}

#-------------------------------------------------
sub GetRemoteHostAndBaseDir {
  my $file = $ENV{'HOME'} . "/.specvtkremoteinfo";
  my $host = undef;
  my $dir  = undef;
  if(-f $file) {
    my $text = Utils::ReadFile($file);
    $text    =~ s/\n+//g;   # Remove newlines
    $text    =~ s/\s+//g;   # Remove whitespace
    if($text =~ m|Host=([^;]+);|) { $host=$1;}
    if($text =~ m|Dir=([^;]+);|)  { $dir =$1;}
    Utils::Die("Cannot find Host in $file") unless(defined $host);
    Utils::Die("Cannot find Dir  in $file") unless(defined $dir);
  } else {
    print STDOUT 
      "Please answer a few questions. The answers will be written\n",
      "to a file '$file' for subsequent uses of this script.\n";
    print STDOUT
      "Please enter the remote hostname to copy Vtk files: >";
    $host=<STDIN>; chomp $host;
    print STDOUT
      "Please enter the remote base directory on $host: >";
    $dir=<STDIN>; chomp $dir;
    Utils::OverwriteFile($file, "Host=$host;\nDir=$dir;\n");
  }
  return ($host,$dir);
}

#-------------------------------------------------
# Copies $setname directory and .visit,.pvd files to $remotedir on a
# remote host. Tar is invoked from basename($setname).
# Also, $setname must be relative to $BaseDir
sub CopyToRemoteDir {
  my($BaseDir,$SrcDir,$RemoteDest,$setname)=@_;
  local $_;

  # Get remote dir
  my $RemoteDir = GetRemoteDir($BaseDir,$RemoteDest);

  # Get Vtk files to copy
  my $FullSetname = "$SrcDir/$setname";
  my @files = ($FullSetname, glob "$FullSetname.*");
  Utils::Die("Empty file list for tar achive.") unless (@files);

  # Get the right path prefixes for files and directories
  my $TarDir = dirname($FullSetname);
  @files = grep {s|$TarDir/||} @files;

  my ($Host,$bla) = GetRemoteHostAndBaseDir();
  print "Copying to $Host:$RemoteDir\n";
  Utils::System("tar cf - -C $TarDir @files | ssh $Host " .
                "'mkdir -p $RemoteDir; cd $RemoteDir; rm -rf @files; tar xf -'");
}

#---------------------------------------------------------
sub FindMatchDirRecursive {
  my($label,$regexp,$dir,$opt_v)=@_;
  print STDERR "Looking for $label in $dir\n" if($opt_v);

  # Try current directory first
  return $dir if IsThisMatchDir($label,$regexp,$dir);
  print STDERR "Cannot find $label in top dir $dir\n" if($opt_v);

  # Not in current directory, try all subdirs.
  opendir(DIR,$dir) or die "Cannot opendir $dir";
  my(@entry)=readdir(DIR);
  closedir(DIR);

  foreach my $e (@entry) {
    next if ($e eq '..' || $e eq '.');
    next unless (-d "$dir/$e");
    print STDERR "Will try $dir/$e\n" if($opt_v);
    my $newdir = FindMatchDirRecursive($label,$regexp,"$dir/$e",$opt_v);
    return $newdir if(defined $newdir);
  }

  return undef;
}

#---------------------------------------------------------
sub IsThisMatchDir {
  my($label,$regexp,$dir)=@_;

  opendir(DIR,$dir) || die "Cannot opendir $dir\n";
  my(@entries) = readdir(DIR);
  closedir(DIR);

  foreach my $entry (@entries) {
    next if ($entry eq '.' || $entry eq '..');
    return 1 if($entry =~ m/$regexp/);
    if($entry =~ m/h5$/) {
      my $h5entries = Utils::SystemOutput("h5ls $dir/$entry");
      return 1 if($h5entries =~ m/$label/);
    }
  }
  return undef;
}

#-------------------------------------------------
1; # This is how to end sources
