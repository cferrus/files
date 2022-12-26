#!/usr/bin/env python
from __future__ import print_function, division

# Input files written by initial data routines

#=============================================================================
def CreateTranslationCpFile(outerV):
    '''Construct an *NthDeriv* style checkpoint file to initialize the
    translation control system to account for the motion of the black holes'''
    temp = """
    Time=  0.0;Version=0;Tsaved=0.0;Nc=3;DerivOrder=2;Tx= 0.0;dTx= %e;""" \
    """d2Tx= 0.0;Ty=0.0;dTy= %e;d2Ty= 0.0;Tz= 0.0;dTz=%e;d2Tz=  0.0;
"""%(outerV[0],outerV[1],outerV[2]);
    with open('ID_Init_FuncTrans.txt', 'w') as f:
        f.write(temp)

#=============================================================================

def ComputeExtr(IDType, useNegExpBC):

    if useNegExpBC:         # doesn't need extrapolation
        Extrapolation_opts = ""
    else:
        Extrapolation_opts = """
        Extrapolator = SphereTaylor(
          Order = 8;
          Interpolator = Spectral;
        );
        """

    temp = """
    DataBoxItems =
      ReadFromFile(File=GlobalItems.input),
      ReadFromFile(File=OrbitalParamItems.input;),
      Domain(
      Items = ReadTensorsFromDiskWithMap(
        DomainHasNoHistoryFile = true;
        Input  =
          ConformalFactor(Dim=3; Sym=;),
          LapseTimesConformalFactor(Dim=3; Sym=;),
          Shift(Dim=3; Sym=1;),
        ;
        Time   = 0;
        DeltaT = 1.e-8;
        Dir    = ../../EllData;
        Interpolator = Simple(
           TopologicalInterpolator = Spectral;
        );
        %s
      );
    ),
    Subdomain(Items=FirstDeriv(Input=Shift;Output=dShift;);),
    """%(Extrapolation_opts)

    # Add conformal data (depends on ID type)
    if IDType != "CFMS":
        temp += """
    Subdomain(Items=SuperposedKerrConformalData(
                      IDType = %s;
                      ConformalMassA = __AM__;
                      ConformalSpinA = __ASpinx__,__ASpiny__,__ASpinz__;
                      BoostVelocityA = __Avx__, __Avy__, __Avz__;
                      CenterA = __Acx__,__Acy__,__Acz__;
                      ConformalMassB = __BM__;
                      ConformalSpinB = __BSpinx__,__BSpiny__,__BSpinz__;
                      BoostVelocityB = __Bvx__, __Bvy__, __Bvz__;
                      CenterB = __Bcx__,__Bcy__,__Bcz__;
                      GaussianFalloffWidthA = __AW__;
                      GaussianFalloffWidthB = __BW__;
                      ),
      );
        """%IDType
    else:
        temp += """
    Subdomain(Items=
      EvaluateMatrixFormula(
        Output=ConformalMetric; Dim=3; Symm=11;
        M[0,0]=1; M[1,1]=1; M[2,2]=1
      ),
      EvaluateMatrixFormula(Output=dtConformalMetric; Dim=3; Symm=11;),
      EvaluateScalarFormula(Output=TrExtrinsicCurvature; Formula=0),
      AddMetricItems(Metric=ConformalMetric; Prefix=Conformal;)
    );
        """

    temp += """
    Observers =
    DumpTensors(
      Input = ConformalFactor,LapseTimesConformalFactor,Shift;
    ),
    ObservePhysicalInitialData(
      Items=
      Subdomain(
      Items=EvaluateVectorFormula(Output=OuterV;
                                  V[0]=__OuterVx__;
                                  V[1]=__OuterVy__;
                                  V[2]=__OuterVz__;),
      EvaluateFormula(Output=NewShift;
                      A = Shift;
                      B = OuterV;
                      Formula=A-B;
                      ),
       );
      Observers=
      ObserveInSubdir(Subdir=ForEvID;Observers=
      #In next line, dump only g, K, lapse, shift...these are for EvID
      DumpTensors(Input=g, K, Lapse, NewShift;
                  FileNames=Nid_g, Nid_K, Nid_N, Nid_Shift);),
    );
        """

    with open('ComputeExtr.input', 'w') as f:
        f.write(temp)

#=============================================================================
from math import log10,sqrt
def computeL1Point(m1,m2,D):
    '''Distance of the L1 Lagrangian point from m1, in Newtonian gravity'''
    return D*(0.5-0.227*log10(m2/m1))

def Domain(Filename,ExcFrac,IDType,NrA=None,NrB=None,LA=None,LB=None,D=None,
    q=None,preamble='',domain_opts=None):

    # Set up the dimensions of the domain
    m1 = q/(1+q)
    m2 = 1-m1
    L1_distA=computeL1Point(m1,m2,D)
    L1_distB=D-L1_distA
    if domain_opts:
        rOutA = domain_opts['rOutA']
        rOutB = domain_opts['rOutB']
    else:
        rOutA = 3./5*L1_distA
        rOutB = 3./5*L1_distB
        #rOutA = D/3.
        #rOutB = D/3.

    BbhIDDomain="""
       SubdomainStructure=
BbhIdDomain(
  IDType = %s;
  CenterA = __Acx__,__Acy__,__Acz__;
  CenterB = __Bcx__,__Bcy__,__Bcz__;
  rA = %s*__Arexc__;
  rB = %s*__Brexc__;
"""%(IDType, ExcFrac, ExcFrac)

    # RA/RB should be Boyer-Linqduist (BL) radii as the
    # HarmonicKerrHorizonConforming expects BL radius.
    # The Boyer-Linqduist radius is related to the HarmonicKerr radius
    # by rBL = rH + M, where M is the mass of the component BH.
    if IDType == "SHK":
        BbhIDDomain +="""
  RA = %e+__AM__;
  RB = %e+__BM__;
  MassA = __AM__;
  MassB = __BM__;
"""%(rOutA, rOutB)
    else:
        BbhIDDomain +="""
  RA = %e;
  RB = %e;
"""%(rOutA, rOutB)

    if NrA is not None:
        BbhIDDomain += "  NrA = %d;\n" % NrA
    if LA is not None:
        BbhIDDomain += "  LA = %d;\n" % LA

    if NrB is not None:
        BbhIDDomain += "  NrB = %d;\n" % NrB
    if LB is not None:
        BbhIDDomain += "  LB = %d;\n" % LB

    # If CF, simply omit all the lines with the information
    # needed for KerrHorizonConforming and Boost maps.
    # ConstructBbhIdDomain will then not apply the maps, as wanted
    if IDType == "SKS":
      BbhIDDomain += """
  # For KerrHorizonConforming - this is "a" of Kerr metric
  SpinA = __ASpinx__,__ASpiny__,__ASpinz__;
  SpinB = __BSpinx__,__BSpiny__,__BSpinz__;

  # For boost map
  BoostA = __ABoostGammax__,__ABoostGammay__, __ABoostGammaz__;
  BoostB = __BBoostGammax__,__BBoostGammay__, __BBoostGammaz__;
"""
    elif IDType == "SHK":
      BbhIDDomain += """
  CoordMapA = HarmonicKerrHorizonConforming(
                  Spin = __ASpinx__,__ASpiny__,__ASpinz__;  # Dimensions of M
                  Mass = __AM__;
            )
            >> RescaleAxes(Dim=3;
                           scaleVec=__ABoostGammax__,
                                    __ABoostGammay__,
                                    __ABoostGammaz__
                          );

  CoordMapB = HarmonicKerrHorizonConforming(
                  Spin = __BSpinx__,__BSpiny__,__BSpinz__;  # Dimensions of M
                  Mass = __BM__;
            )
            >> RescaleAxes(Dim=3;
                           scaleVec=__BBoostGammax__,
                                    __BBoostGammay__,
                                    __BBoostGammaz__);
"""

    BbhIDDomain+="""
  # SphereC  (always centered at origin)
  rC=2*__Sep__;  RC=1e9*__Sep__;

  # geometry coefficients (modify only after talking to Harald)
  # The alpha is expected in degrees. 
  # The value of alpha corresponds to 45 rad, the original setting in SpEC.
  alpha=58.3100780887044394559;
  fcyl=0.95;
  fblock=1.05;
  fC=1.1;
 );

"""
    # [Extr]Domain.input (from SKS)

    # Preamble has any ResizeTheseSubdomains changes from AMR
    temp=preamble
    temp += """
HistoryFile=<<NONE>>;
BoundaryInfo = (WarningLevel=1);\n"""
    temp += BbhIDDomain
    temp += "OwnerPolicy=LoadBalanced(AssignFewestPointsToProc0=false;);\n"

    with open(Filename, 'w') as f:
        f.write(temp)

#=============================================================================
def petsc(IDType,absolute=True):
    if not absolute:
        temp="""
    ##### Options of nonlinear solver
    -snes_rtol __snesrtol__
    -snes_monitor
    -snes_max_it __snes_max_it__
    #-snes_max_it 0
    """
    else:
        temp="""
    ##### Options of nonlinear solver
    -snes_atol __snesatol__
    -snes_stol __snesatol__
    -snes_monitor
    -snes_max_it __snes_max_it__
    #-snes_max_it 0
    """
    if IDType == "CFMS":
        temp+="-snes_max_linear_solve_fail 2"


    temp+="""
    ##### Options for linear solver
    -ksp_type fgmres
    -ksp_gmres_modifiedgramschmidt
    -ksp_gmres_restart 300
    #-ksp_xmonitor
    -ksp_monitor"""
    if not absolute:
        temp+="""
    -ksp_rtol __ksprtol__"""
    else:
        temp+="""
    -ksp_atol __kspatol__"""
    temp+="""
    -ksp_max_it 300
    #-ksp_max_it 6

    ################
    #### PRECONDITIONING
    ################

    ###### No preconditioning
    #-spells_pc_type none

    ###### Invert A_FD exactly
    #-spells_pc_type lu

    ###### Perform ILU(k) on A_FD
    #-spells_pc_type ilu
    #-spells_pc_ilu_levels 2

    ####### Performs a LU decomp on each (processor)-block
    #-spells_pc_type asm
    #-spells_sub_pc_type lu
    #-spells_pc_asm_overlap 0

    ###### Approximate inversion of A_FD
    # (Note: We set the preconditioner to be bjacobi.  This is the default
    # in parallel anyway.  Using it in serial won't hurt, but will allow
    # us to use the same ilu_levels in serial and parallel.)
    -spells_pc_type ksp
    -spells_ksp_ksp_max_it 100
    -spells_ksp_ksp_gmres_restart 200
    -spells_ksp_pc_type bjacobi
    -spells_ksp_sub_pc_factor_levels 4
    # NOTE: replace "factor" by "ilu" for petsc 2.3.0 (from 2.3.1)
    """

    with open('petsc.input', 'w') as f:
        f.write(temp)

#=============================================================================
def Elliptic(IDType, useNegExpBC):
    if IDType != "CFMS":
      Precon_LFF_BC = "Identity(Vars=LapseTimesConformalFactor),"
      InnerLapseTimesConformalFactorBC = \
        "Dirichlet(Var=LapseTimesConformalFactor;TakeFromTensor=%sLapse;)," \
        %IDType
      OuterShiftBC  = "Dirichlet(Var=Shift; TakeFromTensor=OuterV;),"
    else:
      Precon_LFF_BC = "NeumannFlat(Var=LapseTimesConformalFactor; Value=1),"
      InnerLapseTimesConformalFactorBC = \
        "FlatVonNeumann(Var=LapseTimesConformalFactor;),"
      OuterShiftBC  = "Dirichlet(Var=Shift; TakeFromTensor=OuterV;),"

    if useNegExpBC:         # Need some single BH quantities for BC
      def InnerConformalFactorBC(Hole):     # Hole should be 'A' or 'B'
        return """ApparentHorizon(
            SBH_Invg=Bh{0}::Invg;
            SBH_dg=Bh{0}::dg;
            SBH_K=Bh{0}::K;
        ),""".format(Hole)
      def InnerShiftBC(Hole):
        return """InnerShift(Center=__{0}cx__,__{0}cy__,__{0}cz__;
            OmegaInnerShiftBC={0}OmegaInnerShiftBC;
            SBH_Invg=Bh{0}::Invg;
            SBH_Lapse=Bh{0}::Lapse;
            SBH_Shift=Bh{0}::Shift;
            SBH_Velocity=__{0}vx__,__{0}vy__,__{0}vz__;
        ),""".format(Hole)

    else:   # Apparent Horizon BCs
      def InnerConformalFactorBC(Hole):
        return """ApparentHorizon(),"""

      def InnerShiftBC(Hole):
        return """InnerShift(
            Center=__{0}cx__,__{0}cy__,__{0}cz__;
            OmegaInnerShiftBC={0}OmegaInnerShiftBC;
        ),""".format(Hole)


    temp = """
# -*- specinput -*-
Observers=ReadFromFile(File=EllipticObservers.input);
VarsStructure = ConformalFactor, LapseTimesConformalFactor, Shift(Symm=1;);
VolumeTerms   = ExtendedConformalThinSandwich();
FirstColumnInOutputFiles = 0;

Preconditioner = AllAtOnce(
  VolumeTerms = MinusFlatLaplacian(
    Vars = ConformalFactor,LapseTimesConformalFactor,Shift;
  );
  BCs =
    SliceLFF*(
      Identity(Vars=Shift),
      """ + Precon_LFF_BC + """
      NeumannFlat(Var=ConformalFactor; Value=1),
    ),
    SliceUFF.SphereC0(
      Identity(Vars=ConformalFactor,LapseTimesConformalFactor,Shift;),
    ),
  ;
);

ResidualOptions=(DomainInterpolator=
                 ParallelAdaptive(TopologicalInterpolator=CardinalInterpolator;
                                  SendingSubdomainCutoff=1.1;
                                  ReceivingProcCutoff=0.9;););
PetscResolution = Sphere*(Scale=1,1,1; Offset=0,-2,-4),
                  Cyl*(Scale=1,1,1; Offset=0,-4,0),
                  Perim*(Scale=1,1,1; Offset=0,0,0);
CoreOptions = (SnesMonitorFile = snes.dat; KspMonitorFile = ksp.dat;);
ObserveEvenIfNotConverged = __Obs__;

BCs=SliceLFF.SphereA0(
        %s
        %s
        %s
    ),
    SliceLFF.SphereB0(
        %s
        %s
        %s
    ),
    """%(InnerConformalFactorBC('A'), InnerShiftBC('A'), \
        InnerLapseTimesConformalFactorBC, InnerConformalFactorBC('B'), \
        InnerShiftBC('B'), InnerLapseTimesConformalFactorBC)

    temp += """
    SliceUFF.SphereC0(
        Dirichlet(Var=ConformalFactor; Value=1;),
        Dirichlet(Var=LapseTimesConformalFactor; Value=1;),
        """+OuterShiftBC+"""
    );
    """

    with open('Elliptic.input', 'w') as f:
        f.write(temp)

#=============================================================================
def OrbitalParamItems():
    temp="""
    DataBoxItems = Domain(Items=
      ConstVectorOfDouble(Output=OmegaOrbit; Components=0,0,__OmegaOrbit__),
      ConstVectorOfDouble(Output=AOmegaInnerShiftBC;
                          Components=__AOmegarx__, __AOmegary__, __AOmegarz__),
      ConstVectorOfDouble(Output=BOmegaInnerShiftBC;
                          Components=__BOmegarx__, __BOmegary__, __BOmegarz__),
      ConstDouble(Output=dtExpansionFactor; Value=__adot0__),
    );
    """
    with open('OrbitalParamItems.input', 'w') as f:
        f.write(temp)

#=============================================================================
def GlobalItems():
    temp="""
    DataBoxItems =
    Boundary(Items=GlobalSliceIntegrator(Integrator=Spectral);),
    Subdomain(
      Items =
        GlobalIntegrator(Integrator=Spectral),
        GlobalDifferentiator(
          GlobalDifferentiator = Standard(TopologicalDifferentiator=Spectral);
        ),
      ;
    );
    """
    with open('GlobalItems.input', 'w') as f:
        f.write(temp)

#=============================================================================
def SpatialCoordMap():
    temp="""
    DataBoxItems = Domain(
        Items = SpatialCoordMapItems(
            InputCoords = <<Grid>>;
            Prefix = GridToInertial;
            Map = Scaling(SpatialDim=3; a=1; dta=0; dt2a=0);
        );
     );
    """
    with open("SpatialCoordMap.input", 'w') as f:
        f.write(temp)

#=============================================================================

def EllipticObservers(useNegExpBC, debug):

    # if useNegExpBC == True, the apparent horizon is not the inner boundary,
    # so we need to do an apparent horizon find
    def ApparentHorizonFinderDomainItems(Hole):    # Hole should be 'A' or 'B'
        return """
        #================================================================
        #  Ah{0} -- Only needed when using NegativeExpansionBC
        #================================================================
        ParallelStrahlkorperFinder
          (Output    = Ah{0};
           MapPrefix =;
           CheckpointFile = Cp-Strahlkorper{0}.dat;
           Center = __{0}cx__,__{0}cy__,__{0}cz__;
           InitialGuess  = Sphere(Radius=3*__{0}rexc__);
           TopologicalInterpolator=
                Radial(FallbackInterpolator=CardinalInterpolator);
           TerminateOnFailure=yes;
           Adaptive=yes;
           #LmaxFromSphereResolution = yes;
           MaxItsFailureRefinementIncrement = 1;
           MaxItsFailureMaxRefinementFactor = 3.0;
           MaxItsFailureMaxLmesh            = 80;
           OptsToStrahlkorperFinder=
           (MinL = __MinLAh{0}__;
            InitL = __InitLAh{0}__;
            NthetaFormula = floor(max(1.5*L,L+2));
            AdaptiveLSelectors=
            Residual(MinValue=__AhMinRes__;MaxValue=__AhMaxRes__;
                     UseNormalizedResidual=yes;),
            SurfaceAreaElement(SpatialMetric=g;
                               TakeSquareRoot = yes;
                               MinTruncationError = __AhMinTrunc__;
                               MaxTruncationError = __AhMaxTrunc__;),
            SurfaceShape(MinTruncationError=__AhMinTrunc__;
                         MaxTruncationError = __AhMaxTrunc__;);

            StrahlkorperResidual
            =ConstantExpansion(Expansion=0.0;
                               InvMetric=Invg;
                               ExCurv   =K;
                               Gamma    =Christoffel2ndKind);
            Prefix = Ah{0};
            StrahlkorperFinderAlgorithm=
            FastFlow(InvMetric=Invg;
                     Verbosity=0;
                     MaxInterpRetries=3;
            );
            SpatialMetric=g;
           )
        ),""".format(Hole)



    temp="""
Observers =
  DumpTensors(Input = ConformalFactor, LapseTimesConformalFactor, Shift;),
  ObserveMvDoubles(Input=PADM_conformal; Filename=PADM_c.dat;),
  ObserveMvDoubles(Input=JADM_conformal; Filename=JADM_c.dat;),

  # DEBUG
  NormOfTensor(
    Input = ResidualConformalFactor,
            ResidualLapseTimesConformalFactor,
            ResidualShift;
    MetricForTensors=ComponentsSeparately;
    Filename=Residual.dat;
    Op=L2;
  ),
  NormOfTensor(
    Input = DiffConformalFactor,
            DiffLapseTimesConformalFactor,
            DiffShift;
    MetricForTensors=ComponentsSeparately;
    Filename=Diff.dat;
    Op=L2;
  ),
  NormOfTensor(
    Input=Diff;
    MetricForTensors=None;
    Filename=TotalDiff.dat;
    Op=L2;
  ),
  NormOfTensor(
    Input=Diff;
    EachSubdomainSeparately=yes;
    Op=L2;
    Filename=DiffSubdomain.dat;
  ),
  PowerDiagnostics(
    PowerMonitors=PowerConformalFactor,PowerShift,
                  PowerLapseTimesConformalFactor;
    GridDiagnostics=;
    Subdomains=*;
    PowerMonitorOutputFormat=grace;
  ),


  PowerDiagnostics(
    GridDiagnostics=GridDiagPowerShift;
    PowerMonitors=PowerShift;
    Subdomains=*;
    PowerMonitorOutputFormat=dat;
  ),
  ChangeGridResolutionInitialData(
    GridDiagnostics=GridDiagPowerShift;
    OutputFile=curResolutionChanges;
    TruncationImprovement=__curTruncationImprovement__;
    TruncationTarget=__curTruncationTarget__;
    TruncationStop=__curTruncationStop__;
    ConformalBoundaries=(
        (Cyl0, Cyl1, Cyl2, Cyl3, Cyl4),
        (Cyl0, Cyl1, Cyl2, Cyl3, Cyl4));
    ConformalIndex=0,1;
    MinResolutionChange=0;
    ForceRadialResolutionInSpheres=false;
  ),
  ChangeGridResolutionInitialData(
    GridDiagnostics=GridDiagPowerShift;
    OutputFile=nextResolutionChanges;
    TruncationImprovement=__nextTruncationImprovement__;
    TruncationTarget=__nextTruncationTarget__;
    TruncationStop=__nextTruncationStop__;
    ConformalBoundaries=(
        (Cyl0, Cyl1, Cyl2, Cyl3, Cyl4),
        (Cyl0, Cyl1, Cyl2, Cyl3, Cyl4));
    ConformalIndex=0,1;
    MinResolutionChange=0;
    ForceRadialResolutionInSpheres=false;
  ),"""


    SurfacesList = "AhA, AhB"
    ShiftPerpBCItems = ""
    extraTensors = ""
    if useNegExpBC:   # Need to find apparent horizons
        if debug:     # Also get inner boundary surfaces for some diagnostics
            extraSurfacesToCopy = "SurfacesToCopy = InnerBdryA, InnerBdryB;"
            SurfacesList = "AhA, AhB, InnerBdryA, InnerBdryB"
            ShiftPerpBCItems = """
          # The following are only required for ShiftPerpBCOnStrahlkorper
          # observer

          # Assuming only z-rotation
          EvaluateVectorFormula(Output=OmegaxrHi; Dim=3;
                V[0] = 0*x2 - __OmegaOrbit__*x1;
                V[1] = __OmegaOrbit__*x0 - 0*x2;
                V[2] = 0*x1 - 0*x0;
          ),
          EvaluateVectorFormula(Output=dtExpansionFactorTimesrHi; Dim=3;
                V[0] = __adot0__*x0;
                V[1] = __adot0__*x1;
                V[2] = __adot0__*x2;
          ),
          # Co-rotating frame shift
          EvaluateFormula(Output=CoRotShift;
                O=OmegaxrHi; A=dtExpansionFactorTimesrHi;
                S=Shift;
                Formula= S+O+A;
          ),
            """
            extraTensors = "ExtraTensors=Lapse, CoRotShift;"
            extraDomainItems = ""
        else:
            extraSurfacesToCopy = ""
            extraDomainItems= "ConstDouble(Output=Time; Value=0),"

        temp += """
  ObservePhysicalInitialData(
    %s
    Items =
      Domain(
        Items =
          %s
          %s
          %s
          ;
      ),  # end Domain"""%(extraSurfacesToCopy, \
                           extraDomainItems, \
                           ApparentHorizonFinderDomainItems('A'), \
                           ApparentHorizonFinderDomainItems('B'))

    else:            # The apparent horizons are at the inner boundary itself
        temp += """
  ObservePhysicalInitialData(

    SurfacesToCopy = AhA,AhB;

    Items =
    """


    temp += """

      Subdomain(
        Items =
          Trace(
            Input=K;
            Output=TrK;
            Indices=0,1;
            PositionOfIndices=l,l;
            InvMetric=Invg;
          ),
          EvaluateScalarFormula(
            Output=LapseTimesConformalFactor;
            A = Lapse;
            B = ConformalFactor;
            Formula=A*B;
          ),
          UnaryOp(Op=Sqrt; Input=Detg; Output=SqrtDetg;),
          WeylMagnetic(
            Output = WeylB;
            CdK = CdK;
            g = g;
            SqrtDetg = SqrtDetg;
          ),

          %s

        ;
      ),    # end Subdomain
    ; #end Items

    EventAndItemAdders =
      AddStrahlkorperPhysicalQuantities(
        ObservationTrigger = Never;
        BaseName = Surfaces;
        Surfaces = %s;
        AddRedshift=no;
        Metric = g;
        InvMetric = Invg;
        ExtrinsicCurvature = K;
        Christoffel2ndKind = Christoffel2ndKind;
        Ricci = Ricci;
        RicciScalar = RicciScalar;
        WeylMagnetic = WeylB;
        MapPrefixFromGridToAHFrame = <<Identity>>;
        MapFromAHToMeasurementFrame = <<Identity>>;
        %s
      )
    ;
    """%(ShiftPerpBCItems, SurfacesList, extraTensors)


    temp += """

    Observers =

      SurfaceInfo(Input=AhA; OutputCoefficients=yes;),
      SurfaceInfo(Input=AhB; OutputCoefficients=yes;),

      DumpTensors(
        Input=g,RicciScalar,K,Lapse,Ham,Mom,dtg,dtK;
      ),
      NormOfTensor(
        Input=Ham,Mom;
        MetricForTensors=ComponentsSeparately;
        Filename=Constraint.dat;
        Op=L2;
      ),
      NormOfTensor(
        Input=Ham,Mom;
        MetricForTensors=ComponentsSeparately;
        Filename=Constraints_All.dat;
        Op=L2;
        EachSubdomainSeparately = yes;
      ),
      AdmIntegrals(
        Sphere=SphereC0;
        MonopoleTerms=ConformalFactor,LapseTimesConformalFactor;
        Invg=Invg; dg=dg;
        K=K; TrK=TrK;
        ConformalFactor=ConformalFactor;
        g=g;
        CenterOfMassFile=CoM.dat;
        Verbose=yes;
      ),
      StrahlkorperCoordinateSpin(
        Surface=AhA;
        Metric=g;
        ExtrinsicCurvature=K;
        Interpolator=Spectral;
      ),
      StrahlkorperCoordinateSpin(
        Surface=AhB;
        Metric=g;
        ExtrinsicCurvature=K;
        Interpolator=Spectral;
      ),
      """

    # The following give diagnostics for ApparentHorizonBC and InnerShiftBC
    # Epsilon is defined in InnerShift.cpp
    if useNegExpBC and debug:
        temp += """
      ExpansionOnStrahlkorper(
        StrahlkorperDataBoxBaseName = Surfaces;
        FileName = ExpansionOnSurfaces.dat;
        ExpansionType=Outgoing;
        Invg=Invg;
        K=K;
        Christoffel2ndKind=Christoffel2ndKind;
      ),
      ShiftPerpBCOnStrahlkorper(
        StrahlkorperDataBoxBaseName = Surfaces;
        FileName = ShiftPerpBCOnSurfaces.dat;
        Invg=Invg;
        Lapse=Lapse;
        CoRotShift=CoRotShift;
      ),
        """

    temp += """
    ;
  ), #end ObservePhysicalInitialData
;
"""

    with open('EllipticObservers.input', 'w') as f:
        f.write(temp)

#=============================================================================
def EllipticItems(useNegExpBC, debug):

    def SurfaceInfoFromSlice(tag):
        return """
      SurfaceInfoFromSlice(
        Center=Center_A;
        Subdomain=SphereA0;
        SlicePos=0;
        Output=%sA;
      ),
      SurfaceInfoFromSlice(
        Center=Center_B;
        Subdomain=SphereB0;
        SlicePos=0;
        Output=%sB;
      ),
        """%(tag, tag)

    # If using a NegativeExpansionBC, the apparent horizon is within the
    # domain and we need to do an apparent horizon find. So, we don't need
    # SurfaceInfo from the inner boundary. However, if debug==True, we will
    # also get the SurfaceInfo from the inner boundary for some diagnostics.
    if useNegExpBC:
        if debug:
            SurfaceInfoItems = SurfaceInfoFromSlice('InnerBdry')
        else:
            SurfaceInfoItems = ""

    # If using an apparent horizon BC, we don't do an apparent horizon find
    # but instead get the apparent horizon from the inner boundary.
    else:
        SurfaceInfoItems = SurfaceInfoFromSlice('Ah')


    temp="""
# Elliptic items common to both CF and SuperposedKerr
Items =
  ReadFromFile(File=ExtraEllipticItems.input),
  ReadFromFile(File=GlobalItems.input),
  ReadFromFile(File=OrbitalParamItems.input;),

  Domain(
    Items =
      # For BCIMSS
      ConstDouble(Output=Time;Value=0),
      ConstVectorOfDouble(
        Output=Center_A;
        Components=__Acx__,__Acy__,__Acz__;
      ),
      ConstVectorOfDouble(
        Output=Center_B;
        Components=__Bcx__,__Bcy__,__Bcz__;
      ),

      %s
      """%(SurfaceInfoItems)

    temp +="""

      AdmLinearMomentumFromConformalData(
        Output=PADM_conformal;
        AdmConformalIntegrand=AdmConformalIntegrand;
        ConformalFactor=ConformalFactor;
        ConformalMetric=ConformalMetric;
        InverseConformalMetric=InvConformalMetric;
        ConformalChristoffel2ndKind=ConformalChristoffel2ndKind;
        OuterSphere=SphereC0;
        MaxRad=__MaxRad__;
      ),
      AdmAngularMomentumFromConformalData(
        Output=JADM_conformal;
        AdmConformalIntegrand=AdmConformalIntegrand;
        ConformalFactor=ConformalFactor;
        ConformalMetric=ConformalMetric;
        InverseConformalMetric=InvConformalMetric;
        ConformalChristoffel2ndKind=ConformalChristoffel2ndKind;
        OuterSphere=SphereC0;
        IncludeVolumeTerm=True;
        MaxRad=__MaxRad__;
      ),
      AdmLinearMomentumVolumeTerm(
        Output=PADM_VT;
        AdmConformalIntegrand=AdmConformalIntegrand;
        ConformalFactor=ConformalFactor;
        ConformalMetric=ConformalMetric;
        InverseConformalMetric=InvConformalMetric;
        ConformalChristoffel2ndKind=ConformalChristoffel2ndKind;
        OuterSphere=SphereC0;
      ),
    ;
  ), #end Domain

  Subdomain(
    Items =
      NumberOfFilteredModesFromOptions(
        NumberOfFilteredModesForI1=0;
        NumberOfFilteredModesForS1=2;
        NumberOfFilteredModesForB2=0;
        NumberOfFilteredModesForB3=0;
        NumberOfFilteredModesForS2=2;
        NumberOfFilteredModesForB2Radial=0;
        NumberOfFilteredModesForB3Radial=0;
      ),
      EvaluateVectorFormula(
        Output=OuterV;
        V[0]=__OuterVx__;
        V[1]=__OuterVy__;
        V[2]=__OuterVz__;
      ),
      ComputeAdmConformalIntegrand(
        Output=AdmConformalIntegrand;
        Shift=Shift;ConformalFactor=ConformalFactor;
        LapseTimesConformalFactor=LapseTimesConformalFactor;
        InverseConformalMetric=InvConformalMetric;
        dtConformalMetric=dtConformalMetric;
        ConformalChristoffel2ndKind=ConformalChristoffel2ndKind;
        TraceExtrinsicCurvature=TrExtrinsicCurvature;
      ),
      #### For Observers ####
      GridDiagnostics(PowerMonitorItems=PowerShift;),
      PowerMonitor(
        Inputs=ConformalFactor,Shift,LapseTimesConformalFactor;
        UseTensorYlmForS2=no;
        AddRadialMonitorForB2B3=yes;
      ),
    ;
  ), #end Subdomain

  Boundary(
    Items =
      ExtractFromParent(
        Input = ConformalMetric,
                DetConformalMetric,InvConformalMetric,
                ConformalChristoffel2ndKind,
                ConformalRicciScalar, TrExtrinsicCurvature,
                dtTrExtrinsicCurvature,
                OuterV,
                dtConformalMetric
        ;
      ),
      ExtractDerivFromParent(
        Input = dConformalMetric,
                dConformalChristoffel2ndKind,
                dTrExtrinsicCurvature,
                ddtTrExtrinsicCurvature,
                ddtConformalMetric
        ;
      ),
    ;
  ), #end Boundary

; #end Items
"""

    with open('EllipticItems.input', 'w') as f:
        f.write(temp)

#=============================================================================

def CF_ExtraEllipticItems():
    temp="""
DataBoxItems =
  Subdomain(
    Items =

      #Conformal metric
      EvaluateMatrixFormula(
        Output=ConformalMetric;
        Dim=3; Symm=11;
        M[0,0]=1; M[1,1]=1; M[2,2]=1;
      ),
      AddMetricItems(Metric=ConformalMetric; Prefix=Conformal;),

      # set \partial_t \\tilde g_{ij}=0
      EvaluateFormula(
        B=ConformalMetric;
        Output=dtConformalMetric;
        Formula=B*0;
      ),
      FirstDeriv(Input=dtConformalMetric; Output=ddtConformalMetric;),

      #TrK
      EvaluateScalarFormula(Output=TrExtrinsicCurvature;Formula=0.;),
      EvaluateScalarFormula(Output=dtTrExtrinsicCurvature;Formula=0.;),
      FirstDeriv(Input=TrExtrinsicCurvature; Output=dTrExtrinsicCurvature;),

      ################################################################
      # INIITIA GUESS
      ################################################################
      EllipticInitialGuess(
        ImportDir=__LASTDIR__;
        SubdomainItems=
          #### Black hole A ####
          AnalyticEinsteinSolution(
            Solution = IsotropicMaxSliceSchwarzschild(
              C=__C__;
              Mass=__rA__/__RadiusOfBH__;
              Center=__Acx__, __Acy__,__Acz__;
            );
            Output=SolutionA;
          ),
          # compute relevant stuff from first analytic solution
          AnalyticEinstein::g(Input=SolutionA; Output=SolutionA_g),
          Determinant(Input=SolutionA_g; Output=SolutionA_Detg),
          EvaluateFormula(
            D=SolutionA_Detg;
            Output=SolutionA_Psi;
            Formula=D^(1./12);
          ),
          AnalyticEinstein::N(Input=SolutionA; Output=SolutionA_Lapse),
          AnalyticEinstein::Shift(Input=SolutionA; Output=SolutionA_Shift),

          #### Black hole B ####
          AnalyticEinsteinSolution(
            Solution = IsotropicMaxSliceSchwarzschild(
              C=__C__;
              Mass=__rB__/__RadiusOfBH__;
              Center=__Bcx__,__Bcy__,__Bcz__;
            );
            Output=SolutionB;
          ),
          AnalyticEinstein::g(Input=SolutionB; Output=SolutionB_g),
          Determinant(Input=SolutionB_g; Output=SolutionB_Detg),
          EvaluateFormula(
            D=SolutionB_Detg;
            Output=SolutionB_Psi;
            Formula=D^(1./12);
          ),
          AnalyticEinstein::N(Input=SolutionB; Output=SolutionB_Lapse),
          AnalyticEinstein::Shift(Input=SolutionB; Output=SolutionB_Shift),

          # Compose the intiial guess from SolutionA and SolutionB
          EvaluateFormula(
            Output=IG-ConformalFactor;
            A=SolutionA_Psi;B=SolutionB_Psi;
            Formula=A+B-1;
          ),
          EvaluateFormula(
            Output=IG-Shift;
            A=SolutionA_Shift; B=SolutionB_Shift;
            Formula=A+B;
          ),
          EvaluateFormula(
            Output=IG-LapseTimesConformalFactor;
            A=SolutionA_Lapse;    B=SolutionB_Lapse;
            X=SolutionA_Psi;      Y=SolutionB_Psi;
            Formula=A*X+B*Y-1;
          ),
        ;
      ), #end EllipticInitialGuess
    ;
  ), #end Subdomain
;

"""
    with open('ExtraEllipticItems.input', 'w') as f:
        f.write(temp)


#=============================================================================
def SBHSolItems(useNegExpBC):

    # Get the required single BH quantities for NegativeExpansionBC
    if useNegExpBC:
        SbhBdryItems = """
                BhA::Invg,
                BhA::K,
                BhA::Lapse,
                BhA::Shift,
                BhB::Invg,
                BhB::K,
                BhB::Lapse,
                BhB::Shift,
        """
        SbhBdryDerivItems = """
      ExtractDerivFromParent(
        Input = BhA::dg,
                BhB::dg,
      ),
      """
    else:
        SbhBdryItems = ""
        SbhBdryDerivItems = ""
    return SbhBdryItems, SbhBdryDerivItems


#=============================================================================

def SuperposedKerr_ExtraEllipticItems(IDType, useNegExpBC):


    SbhBdryItems, SbhBdryDerivItems = SBHSolItems(useNegExpBC)

    temp="""
DataBoxItems =
  Subdomain(
    Items =

      SuperposedKerrConformalData(
        IDType=%s;
        ConformalMassA = __AM__;
        ConformalSpinA = __ASpinx__,__ASpiny__,__ASpinz__;
        BoostVelocityA = __Avx__, __Avy__, __Avz__;
        CenterA = __Acx__,__Acy__,__Acz__;
        ConformalMassB = __BM__;
        ConformalSpinB = __BSpinx__,__BSpiny__,__BSpinz__;
        BoostVelocityB = __Bvx__, __Bvy__, __Bvz__;
        CenterB = __Bcx__,__Bcy__,__Bcz__;
        GaussianFalloffWidthA = __AW__;
        GaussianFalloffWidthB = __BW__;
      ),"""%IDType

    temp +="""

      # ********************* Conformal Metric Related ******************
      #Actual metric = physical metric using whatever conformal factor
      #spells finds.
      EvaluateScalarFormula(A=ConformalFactor;Output=Psi4;Formula=A*A*A*A;),
      BinaryOp(Output=gActual; A=ConformalMetric; B=Psi4; Op=A*B;),
      BinaryOp(Output=InvgActual; A=InvConformalMetric; B=Psi4; Op=A/B;),

      #Dump residual to file
      TensorNorm(Output=ShiftNorm; Input=Shift; Metric=gActual;),
      TensorNorm(Output=ResShiftNorm; Input=ResidualShift; Metric=gActual;),

      #For physical observers
      FirstDeriv(Input=gActual;Output=dgActual),

      EllipticInitialGuess(
        ImportDir=__LASTDIR__;
        SubdomainItems=
          #The initial guess
          EvaluateScalarFormula(Output=IG-ConformalFactor; Formula=1),
          EvaluateFormula(Output=IG-Shift; A={0}Shift; Formula=A),
          EvaluateFormula(
            Output=IG-LapseTimesConformalFactor;
            A={0}Lapse;
            Formula=A;
          ),
        ;
      ),
    ;
  ), #end Subdomain

  Boundary(
    Items =
      ExtractFromParent(
        Input = gActual,
                IG-ConformalFactor,
                {0}Lapse,
                {0}Shift,
                IG-LapseTimesConformalFactor,IG-Shift,
                InvgActual,
                {1}
        ;
      ),

      {2}
    ;
  ), #end Boundary
;
""".format(IDType, SbhBdryItems, SbhBdryDerivItems)

    with open('ExtraEllipticItems.input', 'w') as f:
        f.write(temp)

#=============================================================================
# Expects all kwargs to be strings, floats, or ints
def Metadata(**kwargs):
    from datetime import datetime
    from socket import gethostname
    import sys
    # Convert all values to strings
    for k,v in kwargs.items():
        if isinstance(v,float) or isinstance(v,int):
            kwargs[k] = "%20.16f" % v
        # Now check for strings.  Python2 and Python3 treat strings
        # differently, and some of the inputs come from functions that
        # in python2 may be unicode or may be str objects.
        elif sys.version_info[0] < 3:
            # python 2: allow 'str' and 'unicode'.
            if type(v) != type(u'Hi') and type(v) != type('Hi'):
                msg = "kwarg ({},{}) should not be a {}".format(k,v,type(v))
                raise TypeError(msg)
        else:
            # python 3: allow only 'str'.
            if type(v) != type('Hi'):
                msg = "kwarg ({},{}) should not be a {}".format(k,v,type(v))
                raise TypeError(msg)

    kwargs["Date"] = datetime.utcnow().strftime('%Y-%m-%d %X UTC')
    kwargs["Loc"]  = gethostname()+":"+kwargs["ID_Dir"]

    # Construct contents of ID_Params.perl
    ID_Params = """
# Approx. radii of apparent horizons in initial data
$ID_rA = {ID_rA};
$ID_rB = {ID_rB};

# Excision radii in extrapolated initial data
$ID_rExcA = {ID_rExcA};
$ID_rExcB = {ID_rExcB};

# Centers of holes
@ID_cA = ({cAx},{cAy},{cAz});
@ID_cB = ({cBx},{cBy},{cBz});
# Coordinate distance between centers of holes
$ID_d = {ID_d};

# Orbital Parameters
$ID_Omega0 = {ID_Omega0};
$ID_adot0  = {ID_adot0};

# ADM-energy of initial data
$ID_Eadm = {ID_Eadm};

# Apparent horizon masses of AH's in initial data
$ID_MAhA = {ID_MAhA};
$ID_MAhB = {ID_MAhB};

# Christoudoulou masses of BHs in ID:
$ID_MA = {ID_MA};
$ID_MB = {ID_MB};

# Dimensionless spins in initial data
@ID_chiA = ({chiAx},{chiAy},{chiAz});
@ID_chiB = ({chiBx},{chiBy},{chiBz});
$ID_chiAMagnitude = {ID_chiAMagnitude};
$ID_chiBMagnitude = {ID_chiBMagnitude};

# ADM linear momentum
@ID_Padm = ({Px},{Py},{Pz});

# ADM angular momentum
@ID_Jadm = ({Jx},{Jy},{Jz});

$ID_Type = {ID_Type};

$ID_Origin="{ID_Type} constructed on {Date}, using SpEC {Revision}\\n".
           "ID location {Loc}";
""".format(**kwargs)

    with open('ID_Params.perl', 'w') as f:
        f.write(ID_Params)
