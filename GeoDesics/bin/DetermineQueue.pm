#!/usr/bin/env perl
package DetermineQueue;  # Start the DetermineQueue namespace

require 5;
use warnings FATAL => 'all';
use strict;
use List::Util 'max';
use Cwd;
use File::Basename;
# Load Utils module
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;

sub RestrictForHydra {
  my (%data) = @_;

  # http://www.rzg.mpg.de/services/computing/hydra/batch-system says that small
  # jobs with less than 16 nodes wil most likely run on the Sandy Bridge nodes
  # For now, respect this and ask for 16 cores per node if the job size is small.
  my %result = ();
  my $threshold = 16*16;
  foreach my $cpn (keys %data) {
    if ($cpn == 16) {
      $result{$cpn} = [ grep { $_ <  $threshold } @{$data{$cpn}} ];
    } elsif ($cpn == 20) {
      $result{$cpn} = [ grep { $_ >= $threshold } @{$data{$cpn}} ];
    }
  }
  return %result;
}

#--------------------------------------------------------------------
1; # This is how to end sources
