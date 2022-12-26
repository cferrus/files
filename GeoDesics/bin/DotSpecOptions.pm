#!/usr/bin/env perl
package DotSpecOptions; # Start the DotSpecOptions namespace
use Exporter;           # load Exporter module
@ISA=qw(Exporter);      # Inherit from Exporter
@EXPORT=qw(DotSpecHelp GetHash GetValue);

require 5;
use strict;
use warnings FATAL => 'all';
use Cwd;
use File::Basename;
# Load Utils module
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;

#-------------------------------------------
# Help text to supplement any script utilizing DotSpecOptions
sub DotSpecHelp {
  return <<END;
  Syntax for specifying options in the ~/.SpEC file is:
      Key1=Value1; Key2=Value2; ...
  Please refrain from using quotes in the option-values, but all
  whitespace is removed, so newlines and spaces are fine. Lines
  beginning with '#' are ignored.

OPTION-KEYS:
  TensorYlmDataBaseDir      # Path to the TensorYlmDataBase directory
                            # Default=''
  EmailAddresses            # Comma-delimited list of e-mail addresses
                            # to send run termination info.
                            # Default=''
  TerminationInfoFile       # File path to append run termination info.
                            # Default=''
  MaxMakeThreads            # Maximum number of threads to use when calling
                            # 'make parallel*' targets. Ignored if >Nprocs.
                            # Default=4
  UseMakeColors             # Enable some colored output in our make system.
                            # Default=0
  UseCcache                 # Enable use of 'ccache'.
                            # Default=0
  DefaultAccount            # Override the DefaultAccount field in Machines.pm.
                            # Default=''
END
}

#-------------------------------------------
# Returns a hash with the option-key and option-value pairs
# @keys = list of option-keys to parse (empty => use all keys)
sub GetHash {
  my(@keys) = @_;

  # populate the database and set the array of keys to get
  my %database = (
    "TensorYlmDataBaseDir" => [\&GetTensorYlmDataBaseDir, ''],
    "EmailAddresses"       => [\&GetValueFromKey, ''],
    "TerminationInfoFile"  => [\&GetValueFromKey, ''],
    "MaxMakeThreads"       => [\&GetValueFromKey, 4],
    "UseCcache"            => [\&GetValueFromKey, 0],
    "UseMakeColors"        => [\&GetUseMakeColors, 0],
    "DefaultAccount"       => [\&GetValueFromKey, ''],
  );
  my @AllKeys = keys %database;
  @keys = @AllKeys if (not @keys);

  # Populate a dictionary with key-value pairs from DotSpecFile (if it exists).
  # There can be multiple semi-colon delimited pairs on each line.
  my $DotSpecFile = "$ENV{'HOME'}/.SpEC";
  my %file_options = ();
  if (-f $DotSpecFile) {
    my @text = Utils::ReadFile($DotSpecFile);
    foreach my $line (@text) {
      $line =~ s/\s+//g;   # Remove whitespace
      next if ($line =~ m/^#/);   # Skip lines starting with '#'
      my @options = split(/;/, $line);
      foreach my $option (@options) {
        if ($option =~ m/([^=]+)=(.*)/) {
          $file_options{$1} = $2;
        } else {
          Utils::Die("Wrong format for '$option'; expected 'key=value'.");
        }
      }
    }
  }

  # Make sure all specified options are valid
  Utils::VerifyHashArgs(\%file_options, @AllKeys);

  # Determine the value associated with each key
  my %result = ();
  foreach my $key (@keys) {
    Utils::Die("Key '$key' not valid.\n") if (not grep {/^$key$/} @AllKeys);
    my $Value = $database{$key}[0]->($key, %file_options);
    my $DefaultValue = $database{$key}[1];
    $Value = $DefaultValue if (not defined $Value);
    $result{$key} = $Value;
  }
  return %result;
}

# Return the option-value for a single key
sub GetValue {
  my($key) = @_;
  my %result = GetHash($key);
  return $result{$key};
}

############################################
# Option-key Get subfunctions
############################################

sub GetValueFromKey {
  my($key, %file_options) = @_;
  my $value = undef;
  if (grep(/^$key$/, keys(%file_options))) {
    $value = $file_options{$key};
  }
  return $value;
}

#-------------------------------------------

sub GetUseMakeColors {
  my($key, %file_options) = @_;
  my $Value = GetValueFromKey($key, %file_options);
  return $Value if (not defined $Value);
  # case-insensitive, and match various bool syntaxes
  if(lc($Value) =~ m/(1|yes|true)/) {
    return 1;
  } elsif(lc($Value) =~ m/(0|no|false)/) {
    return 0;
  } else {
    die "Do not understand option '$Value' for key '$key'.\n";
  }
}

#-------------------------------------------

sub GetTensorYlmDataBaseDir {
  my($key, %file_options) = @_;
  my $TensorYlmDataBaseDir = GetValueFromKey($key, %file_options);
  return $TensorYlmDataBaseDir;
}

#-------------------------------------------------
1; # This is how to end sources
