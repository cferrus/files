#!/usr/bin/env perl
package ComputeRequiredVout;  # Start the namespace
use Exporter;                 # load Exporter module
@ISA=qw(Exporter);            # Inherit from Exporter
@EXPORT=qw(ComputeRequiredVout);

require 5;
use warnings FATAL => 'all';
use strict;
use Cwd;
use File::Basename;
# Load SpEC perl modules
use lib dirname(Cwd::abs_path(__FILE__));
use Utils;
use SpEC;
use Machines;

my $ME    = __FILE__;
# We want to use $RealBin for finding executables to call using system
# calls.  But dirname(__FILE__) points to either a bin directory or a
# Support/Perl directory. If the former, then we use it.  If the
# latter, then there is always a Support/bin directory next to
# Support/Perl, and we want to use Support/bin because Support/bin
# contains things (like python scripts and SpEC executables) that are
# not in Support/Perl.  So "/../bin" does the trick for both cases.
my $RealBin = dirname(__FILE__) . "/../bin";

my $opt_v = 1;  # default verbosity
my $opt_d = 1;  # debugging information

my $MPI1 = new Machines()->GetMpiCmd(1);  #for single process execs
my $APPLYOBS1 = "$MPI1 $RealBin/ApplyObservers";

sub ComputeRequiredVout {
  my($CenterStrings,$time,$h5fileprefix,$MinAllowedCharSpeed,
     $olddomaindir,$oldspatialcoordmapdir,$datadir,
     $useinitialdatavariables,$domainhasnohistoryfile,
     @extrafilestocopy)=@_;

  # Label holding
  my $StartLabel = (scalar(@$CenterStrings)>1 ? "A" : "");
 
  # Temporary files
  my $tempdir  = "ComputeRequiredVout"; 
  # The following is just in case tempdir exists.
  if(-d $tempdir) {
    $tempdir .= "a";
  }
  while(-d $tempdir) {
    ++$tempdir;
  }
  my $tempfile = "ApplyObservers.$$.input"; 

  # Save original dir, cd into temp dir.
  my $origdir = `/bin/pwd`; chomp $origdir;
  mkdir($tempdir) || die "Cannot mkdir $tempdir";
  chdir $tempdir  || die "Cannot chdir $tempdir";

  # List of files that we will delete at the end.
  my @filestocleanup;
  push(@filestocleanup,$tempfile) unless ($opt_d);

  foreach my $file(@extrafilestocopy) {
    my $file2=$file; $file2 =~ s|.*/||;
    Utils::MyCopy($file,$file2);
    push(@filestocleanup,$file2);
  }

  # Construct ApplyObservers string.  We care about Vout on the gr domain.
  my $DomainInput = 'GrDomain.input';
  my $inertial    = SpEC::IsHydro($olddomaindir)
      ? "SpectralInertial" : "Inertial";
  
  {
    my $deltat        = $time*1.e-12+1.e-12;
    my $ComputeItems;
    my $vars        ;
    if($useinitialdatavariables) {
      $vars         = "g(Dim=3;Sym=11;Input=Nid_g),
                       Shift(Dim=3;Sym=1;Input=Nid_Shift),
                       Lapse(Dim=3;Sym=;Input=Nid_N)";
      $ComputeItems = "Subdomain(Items=
                       InvertMatrix(Input=g;Output=Invg);
                       ),";
    } else {
      $vars         = "psi(Dim=4;Sym=11;Input=psi)";
      $ComputeItems = "Subdomain(Items=
                       Add3Plus1ItemsFromGhPsi(psi=psi;OutputPrefix=;);
                       ),";
    }
    my $ApplyObsInput = 
        "DataBoxItems = 
    ReadFromFile(File=SpatialCoordMap.input),
    Domain
    (Items =
     ################################################################
     # 3+1 quantities in inertial frame
     ################################################################
     ReadTensorsFromDiskWithMap
     (Dir                        = $datadir;
      DomainDir                  = $olddomaindir;
      SpatialCoordMapDir         = $oldspatialcoordmapdir;
      Time                       = $time;
      MapPrefixGridToInertial    = GridTo${inertial};
      MapPrefixSrcGridToInertial = GridTo${inertial};
      DeltaT                     = $deltat;
      H5FilePrefix               = $h5fileprefix;
      DomainFile                 = $DomainInput;";
    if($domainhasnohistoryfile) {
      $ApplyObsInput .= " 
      DomainHasNoHistoryFile     = true;";
    }
    $ApplyObsInput .= " 
      Input                      = $vars;
      Interpolator=ParallelAdaptive(TopologicalInterpolator
                                    =CardinalInterpolator;);
      );
    ),
    $ComputeItems";
    if($opt_d) {
      $ApplyObsInput .= "
    Subdomain(Items=";
      my $Label=$StartLabel;
      foreach my $centerstring(@$CenterStrings) {
        $ApplyObsInput .= "
              MinIngoingRadialPhysCharSpeed
              (Output        = CharSpeedMin$Label;
               Center        = $centerstring;
               Invg          = Invg;
               Lapse         = Lapse;
               Shift         = Shift;
               FrameVelocity = GridTo${inertial}::FrameVelocity;
               InvJac        = GridTo${inertial}::InvJacobian;
               ),
              CharSpeedLam00Multiplier
              (Output              = CharSpeedLam00Multiplier$Label;
               Center              = $centerstring;
               Invg                = Invg;
               GridToInertial      = GridTo${inertial};
               DistortedToInertial = DistortedTo${inertial};
               ),";
        ++$Label;
      }
      $ApplyObsInput .= "
    ),";
    }
    $ApplyObsInput .= "
    Boundary(Items=
              GlobalSliceIntegrator(Integrator=Spectral),
              ExtractFromParent(Input=Invg),
              ExtractFromParent(Input=g),
              ExtractFromParent(Input=Lapse),
              ExtractFromParent(Input=Shift),
              EvaluateScalarFormula(Output=OneScalar;Formula=1.;),";
    my $Label=$StartLabel;
    foreach my $centerstring(@$CenterStrings) {
    $ApplyObsInput .= "
              MinIngoingRadialPhysCharSpeed
              (Output        = CharSpeedMin$Label;
               Center        = $centerstring;
               Invg          = Invg;
               Lapse         = Lapse;
               Shift         = Shift;
               FrameVelocity = GridTo${inertial}::FrameVelocity;
               InvJac        = GridTo${inertial}::InvJacobian;
               ),
              CharSpeedLam00Multiplier
              (Output              = CharSpeedLam00Multiplier$Label;
               Center              = $centerstring;
               Invg                = Invg;
               GridToInertial      = GridTo${inertial};
               DistortedToInertial = DistortedTo${inertial};
               ),";
    ++$Label;
    }
    $ApplyObsInput .= "
             );
     Observers =";
    $Label=$StartLabel;
    my $ending = "";
    foreach my $centerstring(@$CenterStrings) {
      $ApplyObsInput .= "
        $ending
        SliceIntegral(Input           = CharSpeedMin$Label;
                      SphericalSlices = 
                             TouchingFromOutside(Center  =$centerstring;
                                                 Globs   =*;
                                                 Eps     =1.e-10;
                                                 MinIndex=0;
                                                 MaxIndex=0;
                                                 );
                      BaseFileName    = CharSpeedMin$Label;
                      MinMaxInstead   = Min;
                      CombineSlices   = yes;
                      VolumeMetric    = None;
                     ),";
      if($opt_d) {
        $ApplyObsInput .= "
        DumpTensors(Input               = CharSpeedMin$Label;
                    OnlyTheseSubdomains = Sphere${Label}0;
                   ),
        ConvertToVtk(Input              = CharSpeedMin$Label;
                     Basename           = GridCharSpeedMin$Label;
                     Subdomains         = Sphere${Label}0;
                     Connectivity       = Auto;
                     FillPoles          = TriangularReduction;
                    ),
        ConvertToVtk(Input              = CharSpeedMin$Label;
                     Basename           = InertialCharSpeedMin$Label;
                     Subdomains         = Sphere${Label}0;
                     Connectivity       = Auto;
                     Coords             = GridTo${inertial}::MappedCoords;
                     FillPoles          = TriangularReduction;
                    ),
        MinMaxOverAngles(Input          = CharSpeedMin$Label;
                         OutputFileName = CharSpeedMin${Label}VsRadius;
                         SphereName     = Sphere${Label}0;
                         MinMax         = Min;
                    ),";
      }
      $ApplyObsInput .= "
        SliceIntegral(Input           = CharSpeedLam00Multiplier$Label,
                                        OneScalar;
                      SphericalSlices = 
                             TouchingFromOutside(Center  =$centerstring;
                                                 Globs   =*;
                                                 Eps     =1.e-10;
                                                 MinIndex=0;
                                                 MaxIndex=0;
                                                );
                      BaseFileName    = CharSpeedLam00Multiplier$Label;
                      MinMaxInstead   = ;
                      VolumeMetric    = g;
                      CombineSlices   = yes;";
      $ending = "),";
      push(@filestocleanup,"CharSpeedLam00Multiplier$Label.dat",
           "CharSpeedMin$Label.dat");
      ++$Label;
    }
    $ApplyObsInput .= "
             );";
    Utils::OverwriteFile($tempfile, $ApplyObsInput . "\n");
  }

  my $applyobserveropts = "-domaininput=$DomainInput -NoDomainHistory -E -B -UseTimes $time ";

  # Run ApplyObservers
  Utils::System("$APPLYOBS1 $applyobserveropts $tempfile", $opt_v);

  # Get output of ApplyObservers
  my @requiredvouts;
  my $Label=$StartLabel;
  foreach my $centerstring(@$CenterStrings) {
    my $RequiredVOut=undef;
    {
      my $MinSpeed = undef;
      {
        my $text =  Utils::ReadFile("CharSpeedMin$Label.dat");
        $text    =~ s/#.*$//m; # Remove comments
        if($text =~ m/^\s*(\S+)\s+(\S+)$/m) {
          $MinSpeed=$2;
        }
        unless(defined $MinSpeed) {die "Cannot read MinSpeed\n";}
      }
      my $FacInt  = undef;
      my $FacNorm = undef;
      {
        my $text =  Utils::ReadFile("CharSpeedLam00Multiplier$Label.dat");
        $text    =~ s/#.*$//m; # Remove comments
        if($text =~ m/^\s*(\S+)\s+(\S+)\s+(\S+)$/m) {
          $FacInt =$2;
          $FacNorm=$3;
        }
        unless(defined $FacInt)  {die "Cannot read FacInt\n";}
        unless(defined $FacNorm) {die "Cannot read FacNorm\n";}
        unless($FacNorm>0) {die "Bad Facnorm $FacNorm\n";}
      }
      my $Fac = $FacInt/$FacNorm;  # Note that $Fac is < 0
      print "Found MinSpeed = $MinSpeed, Fac = $Fac\n" if ($opt_v);
      # v_new = v_old + Fac ExtraLamdot
      #       = v_old - Fac Extrav/Y00
      # So Extrav = (v_new-v_old)*Y00/(-Fac)
      # Now let's demand that v_new > some limit
      my $Pi        = 4*atan2(1,1);
      my $Y00       = sqrt(0.25/$Pi);
      $RequiredVOut = ($MinAllowedCharSpeed-$MinSpeed)*$Y00/(-$Fac);
      print "Found RequiredVOut = $RequiredVOut\n" if ($opt_v);
    }
    push(@requiredvouts,$RequiredVOut);
    ++$Label;
  }

  # Clean up
  unlink(glob "*.txt");
  unlink(@filestocleanup,glob("D*.input")) || die;
  chdir $origdir || die "Cannot chdir $origdir";
  unless($opt_d) {
    # If we have debug on, then we don't want to wipe out the
    # debug directory!
    rmdir $tempdir || die "Cannot rmdir $tempdir";
  }

  return @requiredvouts;
}

#-------------------------------------------------
1; # This is how to end sources

