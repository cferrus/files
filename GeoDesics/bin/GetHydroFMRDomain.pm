#!/usr/bin/env perl
package GetHydroFMRDomain;

require 5;
use warnings FATAL => 'all';
use strict;
use Cwd qw(cwd abs_path);
use File::Basename;
use List::Util qw(min max);
use POSIX ();

# Build HyDomain.input with Fixed Mesh Refinement.
sub GetHydroFMRDomain {
  # Options:
  # Nlev : number of refinement levels.
  # Extents : number of cells per dimension at each level.
  # GridHalfLength : finest level of refinement covers [-L,L].
  # Segments : Number of subdomains in each dimension.
  # SdTag : Subdomains will be called Interval$SdTag-Lev...
  my($Nlev,$Extents,$GridHalfLength,$Segments,$SdTag) = @_;

  $Nlev--;
  my $newdomainstring = "HistoryFile=<<NONE>>;\n";
  $newdomainstring .= "\n";
  $newdomainstring .= "SubdomainStructure =\n";
  $newdomainstring .= "PerimBlocks(BaseName = Interval$SdTag-Lev$Nlev-;\n";
  foreach my $d ("x","y","z"){
      $newdomainstring .= "  $d-Axis = (\n";
      $newdomainstring .= "            Extents = $Extents;\n";
      $newdomainstring .= "            Bounds  = -$GridHalfLength,$GridHalfLength;\n";
      $newdomainstring .= "            SplitIntervals = $Segments;\n";
      $newdomainstring .= "            Maps     = Lin;\n";
      $newdomainstring .= "            IndexMap = Uniform;\n";
      $newdomainstring .= "            Topology = I1;\n";
      $newdomainstring .= "            GhostZones = 3;\n";
      $newdomainstring .= "            GhostZonesOnBoundary = true;\n";
      $newdomainstring .= "            Centering = Cell;\n";
      $newdomainstring .= "            );\n";
  }
  $newdomainstring .= "),\n";
  $Nlev--;
  for (;$Nlev>=0;$Nlev--) {
      $newdomainstring .= "PerimBlocks(BaseName = Interval$SdTag-Lev$Nlev-;\n";
      $GridHalfLength*=2;
      foreach my $d ("x","y","z"){
          $newdomainstring .= "  $d-Axis = (\n";
          $newdomainstring .= "            Extents = $Extents;\n";
          $newdomainstring .= "            Bounds  = -$GridHalfLength,$GridHalfLength;\n";
          $newdomainstring .= "            SplitIntervals = $Segments;\n";
          $newdomainstring .= "            Maps     = Lin;\n";
          $newdomainstring .= "            IndexMap = Uniform;\n";
          $newdomainstring .= "            Topology = I1;\n";
          $newdomainstring .= "            GhostZones = 3;\n";
          $newdomainstring .= "            GhostZonesOnBoundary = true;\n";
          $newdomainstring .= "            Centering = Cell;\n";
          $newdomainstring .= "            );\n";
      }
      $newdomainstring .= "  SkipTheseSubdomains=";
      foreach my $x (2..5){
        foreach my $y (2..5){
          foreach my $z (2..5){
            $newdomainstring .= "Interval$SdTag-Lev$Nlev-$x.$y.$z,"
          }
        }
      }
      $newdomainstring .= ";\n),\n";
  }
  $newdomainstring .= ";\n";
  $newdomainstring .= "FileBaseName=HyDomain;\n";
  return $newdomainstring;
}
