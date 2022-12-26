#!/usr/bin/env python
from __future__ import print_function, division

import sys, os
from math import acos, pi
import argparse as ap
from datetime import datetime
from time import sleep
from shutil import copy, copytree
from glob import glob
import re  # For subn
from numpy import array, dot, sqrt, genfromtxt, cross
sys.path.insert(0, os.path.realpath(__file__ + '/../../../Support/Python'))
from Utils import System, SystemOutput, error, warning, norm, call_perl, ReadH5
import BBH_ID_WriteInputFiles
import numpy as np
from numpy.linalg import inv
import h5py

# This command forces print to flush on newline
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)

# this number should be chosen consistently with the inner lapse-BC.
# It is the radius of an isolated non-spinning, non-boosted black hole
# of unit-mass, in those coordinates.
RADIUS_OF_BH = 0.85949977036002983

#=============================================================================


def computeL1Point(m1, m2, D):
    '''Distance of the L1 Lagrangian point from m1, in Newtonian gravity'''
    return D * (0.5 - 0.227 * np.log10(m2 / m1))


def NegativeExpansionBCSettings(debug):

    if useNegExpBC:

        ## If useNegExpBC==True, we never extrapolate, but still use ExtrFrac
        ## to place the inner boundary inside the apparent horizon. This means
        ## that Domain.input and ExtrDomain.input are the same for this case.
        InnerRadiusFrac = "__ExtrFrac__"

        # These are used by ApparentHorizonFinder, since the inner boundary
        # is no longer an apparent horizon, we need to do a search for it.
        ApparentHorizonFinderSettings = """
# Apparent horizon finder options
# Get L from Domain.input extents for SphereA0 and SphereB0
$MinLAhA = max(`more Domain.input | grep SphereA0 | cut -d '=' -f 2 |""" \
    """ cut -d ',' -f 2`, 10);
$MinLAhB = max(`more Domain.input | grep SphereB0 | cut -d '=' -f 2 |""" \
    """ cut -d ',' -f 2`, 10);
$InitLAhA = $MinLAhA;
$InitLAhB = $MinLAhB;
$AhMinRes=10**($curTruncationTarget);
$AhMaxRes=10**($curTruncationTarget+1);
$AhMinTrunc=$AhMinRes;
$AhMaxTrunc=$AhMaxRes;
        """
        AdditionalEllipObsSettings = """
                         "__Acx__" => $Acx,
                         "__Acy__" => $Acy,
                         "__Acz__" => $Acz,
                         "__Bcx__" => $Bcx,
                         "__Bcy__" => $Bcy,
                         "__Bcz__" => $Bcz,
                         "__Arexc__" => $Arexc,
                         "__MinLAhA__" => $MinLAhA,
                         "__InitLAhA__" => $InitLAhA,
                         "__Brexc__" => $Brexc,
                         "__MinLAhB__" => $MinLAhB,
                         "__InitLAhB__" => $InitLAhB,
                         "__AhMinRes__" => $AhMinRes,
                         "__AhMaxRes__" => $AhMaxRes,
                         "__AhMinTrunc__" => $AhMinTrunc,
                         "__AhMaxTrunc__" => $AhMaxTrunc,
        """
        if debug:
            AdditionalEllipObsSettings += """
                         "__OmegaOrbit__" => $OmegaOrbit,
                         "__adot0__" => $adot0,
            """

        AdditionalDomainSettings = """
                                        "__ExtrFrac__" => $ExtrFrac,
        """

        AdditionalEllipticSettings = """
                                           "__Avx__" => $Avx,
                                           "__Avy__" => $Avy,
                                           "__Avz__" => $Avz,
                                           "__Bvx__" => $Bvx,
                                           "__Bvy__" => $Bvy,
                                           "__Bvz__" => $Bvz,
        """

    else:  # These are for an Apparent horizon BC
        InnerRadiusFrac = "1.0"
        AdditionalEllipObsSettings = ""
        ApparentHorizonFinderSettings = ""
        AdditionalDomainSettings = ""
        AdditionalEllipticSettings = ""

    return InnerRadiusFrac, AdditionalEllipObsSettings, \
        AdditionalDomainSettings, AdditionalEllipticSettings, \
        ApparentHorizonFinderSettings


def PrepareSuperposedKerrInputFiles(initial_guess_dir,
                                    Lev,
                                    D,
                                    Omega0,
                                    adot,
                                    chiA,
                                    chiB,
                                    outerV,
                                    rA,
                                    rB,
                                    OmegaA,
                                    OmegaB,
                                    cxA,
                                    cyA,
                                    czA,
                                    debug,
                                    q=None,
                                    preamble='',
                                    elliptic_opts=None,
                                    amr_opts=None,
                                    domain_opts=None,
                                    SuperposedKerr_opts=None,
                                    rad_factor=None):

    ####PASS CONFORMAL MASS AND CONFORMAL SPIN INSTEAD OF rA,rB
    ####But should we change that later?? (Use R instead of M as free param)
    ####(Easy linear relation)

    InnerRadiusFrac, AdditionalEllipObsSettings, AdditionalDomainSettings, \
        AdditionalEllipticSettings, ApparentHorizonFinderSettings \
        = NegativeExpansionBCSettings(debug)

    if IDType == "SHK":
        # Domain.input for SHK also needs component masses
        AdditionalDomainSettings += """
                                        "__AM__" => $AM,
                                        "__BM__" => $BM,
        """
        AdditionalExtrDomainSettings = """
                                        "__AM__" => $AM,
                                        "__BM__" => $BM,
        """
    else:
        AdditionalExtrDomainSettings = ""

    # InnerRadiusFrac = '__ExtrFrac__' when useNegExpBC == True, thus
    # Domain.input and ExtrDomain.input are the same and we don't extrapolate.
    # InnerRadiusFrac = '1.0' if useNegExpBC == False
    BBH_ID_WriteInputFiles.Domain("Domain.input",
                                  InnerRadiusFrac,
                                  IDType,
                                  D=D,
                                  q=q,
                                  preamble=preamble,
                                  domain_opts=domain_opts)
    BBH_ID_WriteInputFiles.Domain("ExtrDomain.input",
                                  "__ExtrFrac__",
                                  IDType,
                                  D=D,
                                  q=q)

    # Update the domain with the correct resolution from previous iteration
    BBH_ID_WriteInputFiles.GlobalItems()
    BBH_ID_WriteInputFiles.OrbitalParamItems()
    BBH_ID_WriteInputFiles.SpatialCoordMap()
    BBH_ID_WriteInputFiles.Elliptic(IDType, useNegExpBC)
    BBH_ID_WriteInputFiles.EllipticItems(useNegExpBC, debug)
    BBH_ID_WriteInputFiles.SuperposedKerr_ExtraEllipticItems(
        IDType, useNegExpBC)
    BBH_ID_WriteInputFiles.EllipticObservers(useNegExpBC, debug)

    # Elliptic options
    # Be default adaptive *absolute* tolerances are used
    if "ksp_rtol" in elliptic_opts:
        ksp_rtol = float(elliptic_opts['ksp_rtol'])
        snes_rtol = float(elliptic_opts['snes_rtol'])
        absolute = False
    else:
        ksp_atol = 1. / 5 * 10**(-float(elliptic_opts['ksp_atol_base']) - Lev)
        snes_atol = 10**(-float(elliptic_opts['snes_atol_base']) - Lev)
        absolute = True

    BBH_ID_WriteInputFiles.petsc(IDType, absolute=absolute)
    # Adaptive mesh refinement options
    base_lev = amr_opts['base_lev']
    trunc_impr = 0

    # Set up AMR observer for the *current* Lev
    cur_AMR_tol = 10**(-(base_lev + Lev))
    cur_trunc_targ = np.log10(cur_AMR_tol)
    cur_trunc_stop = np.log10(cur_AMR_tol)

    # Set up AMR observer for the *next* Lev
    next_AMR_tol = 10**(-(base_lev + Lev + 1))
    next_trunc_targ = np.log10(next_AMR_tol)
    next_trunc_stop = np.log10(next_AMR_tol)

    chiAtargetSq = chiAtarget[0]**2 + chiAtarget[1]**2 + chiAtarget[2]**2
    chiBtargetSq = chiBtarget[0]**2 + chiBtarget[1]**2 + chiBtarget[2]**2
    if chiAtargetSq > 1:
        error("You have set SpinA with magnitude greater than 1 in your\n"
              "Params.input file. This does not yield a black hole! Please\n"
              "fix your SpinA value, which is currently\n"
              "SpinA =" + str(chiAtarget))
    if chiBtargetSq > 1:
        error("You have set SpinB with magnitude greater than 1 in your\n"
              "Params.input file. This does not yield a black hole! Please\n"
              "fix your SpinB value, which is currently\n"
              "SpinB =" + str(chiBtarget))

    if IDType == "SKS":
        AM = rA / (1 + sqrt(1 - chiAtargetSq))
        BM = rB / (1 + sqrt(1 - chiBtargetSq))
    elif IDType == "SHK":
        # The mass has to computed assuming rA is a Harmonic-Kerr
        # horizon radius. This is related to Boyer-Lindquist
        # coordinates as: rH = rBL - M
        AM = rA / (sqrt(1 - chiAtargetSq))
        BM = rB / (sqrt(1 - chiBtargetSq))

    ASpinx = chiA[0] * AM
    ASpiny = chiA[1] * AM
    ASpinz = chiA[2] * AM
    Avx = cxA * adot
    Avy = cxA * Omega0
    Avz = 0

    if SuperposedKerr_opts:
        AW = SuperposedKerr_opts['wA']
        BW = SuperposedKerr_opts['wB']
    else:
        L1_distA = computeL1Point(q / (1.0 + q), 1.0 / (1.0 + q), D)
        L1_distB = D - L1_distA
        AW = 3. / 5 * L1_distA
        BW = 3. / 5 * L1_distB

    ABoostGammax = 1 / sqrt(1 - Avx * Avx)
    ABoostGammay = 1 / sqrt(1 - Avy * Avy)
    ABoostGammaz = 1 / sqrt(1 - Avz * Avz)
    BSpinx = chiB[0] * BM
    BSpiny = chiB[1] * BM
    BSpinz = chiB[2] * BM
    Bcx = cxA - D
    Bcy = cyA
    Bcz  = czA
    Bvx = Bcx * adot
    Bvy = Bcx * Omega0
    Bvz = 0

    BBoostGammax = 1 / sqrt(1 - Bvx * Bvx)
    BBoostGammay = 1 / sqrt(1 - Bvy * Bvy)
    BBoostGammaz = 1 / sqrt(1 - Bvz * Bvz)

    if IDType == "SKS":
        Arexc = rA
        Brexc = rB
    elif IDType == "SHK":
        # These have to be in BL radius, not Harmonic-Kerr because
        # the HarmonicKerrHorizonConforming map expects BL radius.
        Arexc = rA + AM
        Brexc = rB + BM

    MaxRad = max(AW, BW) * rad_factor

    if AW < 2 * rA:
        error("For your configuration, the width of the Gaussian weight \
 function around BH A is less than twice its radius.  WA=%g, rA=%g.\n" %
              (AW, rA) +
              "Please adjust the width, choose other parameters, or remove"
              " this test and proceed at your own risk.\n")
    if BW < 2 * rB:
        error("SuperposedKerr WARNING\n"
              "For your configuration, the width of the Gaussian weight \n"
              "function around BH B is less than twice its radius.  \n"
              "WB=%g, rB=%g.\n" % (BW, rB) +
              "Please adjust the width, choose other parameters, or remove \n"
              "this test and proceed at your own risk.\n")

    # DoMultipleRuns.input
    temp = """
# -*- perl -*-
use List::Util qw[min max];

$Level = %i;

$Sep = %20.16f;

$OmegaOrbit = %20.16f;  # Assumed only in z-direction


###Change to adot0
###Python variable as well
$adot0=%20.16f;   """ % (Lev, D, Omega0, adot)

    temp += """

$AOmegarx = %20.16f;
$AOmegary = %20.16f;
$AOmegarz = %20.16f;
$BOmegarx = %20.16f;
$BOmegary = %20.16f;
$BOmegarz = %20.16f;   """ % (OmegaA[0], OmegaA[1], OmegaA[2] + Omega0,
                              OmegaB[0], OmegaB[1], OmegaB[2] + Omega0)

    temp += """

$AM = %20.16f;
$BM = %20.16f;  """ % (AM, BM)

    temp += """

#ASpin is angular momentum per unit mass!

$ASpinx = %20.16f;
$ASpiny = %20.16f;
$ASpinz = %20.16f;
$Avx = %20.16f;
$Avy = %20.16f;
$Avz = %20.16f;
$Acx = %20.16f;
$Acy = %20.16f;
$Acz = %20.16f;
$AW = %20.16f;

##Gammax, Gammaz can be hardcoded to 0

$ABoostGammax = %20.16f;
$ABoostGammay = %20.16f;
$ABoostGammaz = %20.16f; """ % (ASpinx, ASpiny, ASpinz, Avx, Avy, Avz, cxA,
                                cyA, czA, AW, ABoostGammax, ABoostGammay,
                                ABoostGammaz)

    temp += """

#BSpin is angular momentum per unit mass!

$BSpinx = %20.16f;
$BSpiny = %20.16f;
$BSpinz = %20.16f;
$Bvx = %20.16f;
$Bvy = %20.16f;
$Bvz = %20.16f;
$Bcx = %20.16f;
$Bcy = %20.16f;
$Bcz = %20.16f;
$BW = %20.16f;
$BBoostGammax = %20.16f;
$BBoostGammay = %20.16f;
$BBoostGammaz = %20.16f; """ % (BSpinx, BSpiny, BSpinz, Bvx, Bvy, Bvz, Bcx,
                                Bcy, Bcz, BW, BBoostGammax, BBoostGammay,
                                BBoostGammaz)

    temp += """

$Arexc = %20.16f;
$Brexc = %20.16f; """ % (Arexc, Brexc)

    temp += """


$OuterVx = %e ;
$OuterVy = %e;
$OuterVz = %e;

$ExtrFrac = %20.16f;

""" % (outerV[0], outerV[1], outerV[2], opts.ExtrFrac)

    temp += """
# ADAPTIVE MESH REFINEMENT
# We have 2 AMR drivers, which differ in the desired precision by
# a factor of 10.
# Note that TruncationStop = log10(10^{-3-Lev})
$TruncationImprovement=%.2f;
$curTruncationTarget=%.2f;
$curTruncationStop=%.2f;
$nextTruncationTarget=%.2f;
$nextTruncationStop=%.2f;
""" % (trunc_impr, cur_trunc_targ, cur_trunc_stop, next_trunc_targ,
       next_trunc_stop)

    # if this is the second or later iteration,
    # use the output of the previous iteration as initial guess
    if initial_guess_dir == "":
        temp += "$last_dir='';"
    else:
        temp += "$last_dir='%s';" % (initial_guess_dir)

    temp += """
$Prefix="Lev";

"""

    temp += """
$dir = $Prefix.$Level;

$THISDIR = `/bin/pwd`; chomp $THISDIR;

# Petsc options
$snes_max_it = """
    temp += SnesIts0
    if not absolute:
        temp += """;
$ksprtol = "%e";""" % ksp_rtol
        temp += """
$snesrtol = "%e";""" % snes_rtol
    else:
        temp += """;
$kspatol = "%e";""" % ksp_atol
        temp += """
$snesatol = "%e";""" % snes_atol

    temp += """
# Cut-off radius for PADM calculation
$MaxRad = %e;

%s

""" % (MaxRad, ApparentHorizonFinderSettings)

    temp += """
RunInDirectory($dir,
                    {"Domain.input" => {
                                        "__ABoostGammax__" => $ABoostGammax,
                                        "__ABoostGammay__" => $ABoostGammay,
                                        "__ABoostGammaz__" => $ABoostGammaz,
                                        "__ASpinx__" => $ASpinx,
                                        "__ASpiny__" => $ASpiny,
                                        "__ASpinz__" => $ASpinz,
                                        "__Acx__" => $Acx,
                                        "__Acy__" => $Acy,
                                        "__Acz__" => $Acz,
                                        "__BBoostGammax__" => $BBoostGammax,
                                        "__BBoostGammay__" => $BBoostGammay,
                                        "__BBoostGammaz__" => $BBoostGammaz,
                                        "__BSpinx__" => $BSpinx,
                                        "__BSpiny__" => $BSpiny,
                                        "__BSpinz__" => $BSpinz,
                                        "__Bcx__" => $Bcx,
                                        "__Bcy__" => $Bcy,
                                        "__Bcz__" => $Bcz,


                                        "__Arexc__" => $Arexc,
                                        "__Brexc__" => $Brexc,
                                        "__Sep__" => $Sep,

                                        %s
                     },""" % (AdditionalDomainSettings)

    temp += """
                     "ExtrDomain.input" => {

                                        "__ABoostGammax__" => $ABoostGammax,
                                        "__ABoostGammay__" => $ABoostGammay,
                                        "__ABoostGammaz__" => $ABoostGammaz,
                                        "__ASpinx__" => $ASpinx,
                                        "__ASpiny__" => $ASpiny,
                                        "__ASpinz__" => $ASpinz,
                                        "__Acx__" => $Acx,
                                        "__Acy__" => $Acy,
                                        "__Acz__" => $Acz,
                                        "__BBoostGammax__" => $BBoostGammax,
                                        "__BBoostGammay__" => $BBoostGammay,
                                        "__BBoostGammaz__" => $BBoostGammaz,
                                        "__BSpinx__" => $BSpinx,
                                        "__BSpiny__" => $BSpiny,
                                        "__BSpinz__" => $BSpinz,
                                        "__Bcx__" => $Bcx,
                                        "__Bcy__" => $Bcy,
                                        "__Bcz__" => $Bcz,


                                        "__Arexc__" => $Arexc,
                                        "__Brexc__" => $Brexc,
                                        "__ExtrFrac__" => $ExtrFrac,
                                        "__Sep__" => $Sep,
                                        %s
                                        },
    """ % AdditionalExtrDomainSettings

    temp+="""
                     "EllipticItems.input" => {
                                               "__Acx__" => $Acx,
                                               "__Acy__" => $Acy,
                                               "__Acz__" => $Acz,
                                               "__Bcx__" => $Bcx,
                                               "__Bcy__" => $Bcy,
                                               "__Bcz__" => $Bcz,
                                               "__OuterVx__" => $OuterVx,
                                               "__OuterVy__" => $OuterVy,
                                               "__OuterVz__" => $OuterVz,
                                               "__MaxRad__" => $MaxRad,
                                               },
                     "ExtraEllipticItems.input" => {
                         "__LASTDIR__" => $last_dir,
                         "__AM__" => $AM,
                         "__BM__" => $BM,
                         "__Acx__" => $Acx,
                         "__Acy__" => $Acy,
                         "__Acz__" => $Acz,
                         "__Bcx__" => $Bcx,
                         "__Bcy__" => $Bcy,
                         "__Bcz__" => $Bcz,
                         "__Avx__" => $Avx,
                         "__Avy__" => $Avy,
                         "__Avz__" => $Avz,
                         "__Bvx__" => $Bvx,
                         "__Bvy__" => $Bvy,
                         "__Bvz__" => $Bvz,
                         "__ASpinx__" => $ASpinx,
                         "__ASpiny__" => $ASpiny,
                         "__ASpinz__" => $ASpinz,
                         "__BSpinx__" => $BSpinx,
                         "__BSpiny__" => $BSpiny,
                         "__BSpinz__" => $BSpinz,
                         "__AW__" => $AW,
                         "__BW__" => $BW,
                     },


                     "Elliptic.input" => { "__Acx__" => $Acx,
                                           "__Acy__" => $Acy,
                                           "__Acz__" => $Acz,
                                           "__Bcx__" => $Bcx,
                                           "__Bcy__" => $Bcy,
                                           "__Bcz__" => $Bcz,
                                           "__Obs__" => '%s',
                                           %s
                                           },
                                           """%(ObserveIfNotConverged, \
                                                AdditionalEllipticSettings)

    temp += """
                 "EllipticObservers.input" => {
                     "__curTruncationImprovement__" => $TruncationImprovement,
                     "__curTruncationTarget__" => $curTruncationTarget,
                     "__curTruncationStop__" => $curTruncationStop,
                     "__nextTruncationImprovement__" => $TruncationImprovement,
                     "__nextTruncationTarget__" => $nextTruncationTarget,
                     "__nextTruncationStop__" => $nextTruncationStop,
                     %s
                  },""" % AdditionalEllipObsSettings
    if not absolute:
        temp += """

                      "petsc.input" => {"__snes_max_it__" => $snes_max_it,
                                        "__ksprtol__" => $ksprtol,
                                        "__snesrtol__" => $snesrtol,
                                       },"""
    else:
        temp += """
                       "petsc.input" => {"__snes_max_it__" => $snes_max_it,
                                        "__kspatol__" => $kspatol,
                                        "__snesatol__" => $snesatol,
                                       },"""
    temp += """
            "OrbitalParamItems.input" => {
              "__OmegaOrbit__" => $OmegaOrbit,
              "__AOmegarx__" => $AOmegarx,
              "__AOmegary__" => $AOmegary,
              "__AOmegarz__" => $AOmegarz,
              "__BOmegarx__" => $BOmegarx,
              "__BOmegary__" => $BOmegary,
              "__BOmegarz__" => $BOmegarz,
              "__adot0__" => $adot0,
            },

      }
    );
    $last_dir="../".$dir;
    """

    f = open('DoMultipleRuns.input', 'w')
    f.write(temp)
    f.close()


#=============================================================================
def PrepareCFInputFiles(initial_guess_dir,
                        Lev,
                        D,
                        Omega0,
                        adot,
                        outerV,
                        debug,
                        rA,
                        rB,
                        OmegaA,
                        OmegaB,
                        cxA,
                        cyA,
                        czA,
                        q=None,
                        preamble='',
                        elliptic_opts=None,
                        amr_opts=None,
                        domain_opts=None):

    ############################################################
    # Elliptic solver options
    ############################################################
    # Be default adaptive *absolute* tolerances are used
    if "ksp_rtol" in elliptic_opts:
        ksp_rtol = float(elliptic_opts['ksp_rtol'])
        snes_rtol = float(elliptic_opts['snes_rtol'])
        absolute = False
    else:
        ksp_atol = 1. / 5 * 10**(-float(elliptic_opts['ksp_atol_base']) - Lev)
        snes_atol = 10**(-float(elliptic_opts['snes_atol_base']) - Lev)
        absolute = True

    BBH_ID_WriteInputFiles.petsc(IDType, absolute=absolute)

    ############################################################

    #############################################################
    # Adaptive mesh refinement options
    #############################################################

    # Truncation improvement is always set to 0 for BBH runs!

    if amr_opts is not None:
        base_lev = amr_opts['base_lev']
    else:
        base_lev = 3

    trunc_impr = 0.0
    # Set up AMR observer for the *current* Lev
    cur_AMR_tol = 10**(-(base_lev + Lev))
    cur_trunc_targ = np.log10(cur_AMR_tol)
    cur_trunc_stop = np.log10(cur_AMR_tol)

    # Set up AMR observer for the *next* Lev
    next_AMR_tol = 10**(-(base_lev + Lev + 1))
    next_trunc_targ = np.log10(next_AMR_tol)
    next_trunc_stop = np.log10(next_AMR_tol)
    ############################################################

    BBH_ID_WriteInputFiles.SpatialCoordMap()
    BBH_ID_WriteInputFiles.OrbitalParamItems()
    BBH_ID_WriteInputFiles.petsc(IDType)
    BBH_ID_WriteInputFiles.Domain("Domain.input","0001.000000000",IDType, \
                                 D=D,q=q,preamble=preamble)
    BBH_ID_WriteInputFiles.Elliptic(IDType, useNegExpBC)
    BBH_ID_WriteInputFiles.EllipticItems(useNegExpBC, debug)
    BBH_ID_WriteInputFiles.CF_ExtraEllipticItems()
    BBH_ID_WriteInputFiles.GlobalItems()
    BBH_ID_WriteInputFiles.EllipticObservers(useNegExpBC, debug)

    Arexc = rA
    Brexc = rB

    # DoMultipleRuns.input
    temp = """# -*- perl -*-
$Level = %d;
# Radii of excision spheres
$RadiusOfBH = '%20.16f';
$rA=%20.16f;
$rB=%20.16f;               """ % (Lev, RADIUS_OF_BH, rA, rB)

    temp += """
# Petsc options
$snes_max_it = """
    temp += SnesIts0

    if not absolute:
        temp += """;
$ksprtol = "%e";""" % ksp_rtol
        temp += """
$snesrtol = "%e";""" % snes_rtol
    else:
        temp += """;
$kspatol = "%e";""" % ksp_atol
        temp += """
$snesatol = "%e";""" % snes_atol
    temp += """
# Coordinate separation between centers of excision spheres
# Centers of holes
$D  = %20.16f;
$Acx= %20.16f;
$Acy= %20.16f;
$Acz= %20.16f;
$Bcz= $Acz;
$Bcx= $Acx-$D;
$Bcy= $Acy;             """ % (D, cxA, cyA, czA)

    temp += """
# Orbital Parameters
$Omega0 = %20.16f;
$adot   = %20.16f;        """ % (Omega0, adot)

    temp += """
# Angular frequencies of horizons
$OmegaAx = %20.16f;
$OmegaAy = %20.16f;
$OmegaAz = %20.16f;
$OmegaBx = %20.16f;
$OmegaBy = %20.16f;
$OmegaBz = %20.16f;       """ % (OmegaA[0], OmegaA[1], OmegaA[2], OmegaB[0],
                                 OmegaB[1], OmegaB[2])
    temp += """
$Arexc = %20.16f;
$Brexc = %20.16f; """ % (Arexc, Brexc)
    temp += """
$OuterVx = %e ;
$OuterVy = %e;
$OuterVz = %e;
""" % (outerV[0], outerV[1], outerV[2])

    temp += """
# ADAPTIVE MESH REFINEMENT
# We have 2 AMR drivers, which differ in the desired precision by
# a factor of 10.
# Note that TruncationStop = log10(10^{-3-Lev})
$TruncationImprovement=%.2f;
$curTruncationTarget=%.2f;
$curTruncationStop=%.2f;
$nextTruncationTarget=%.2f;
$nextTruncationStop=%.2f;
""" % (trunc_impr, cur_trunc_targ, cur_trunc_stop, next_trunc_targ,
       next_trunc_stop)

    temp += """
# For initial guess of elliptic solver
$C = 1.7370341836426595287;

# MaxRad for CF should use the default value in the Adm classes
$MaxRad = 1.e5;

"""

    # if this is the second or later iteration,
    # use the output of the previous iteration as initial guess
    if initial_guess_dir == "":
        temp += "$last_dir='';"
    else:
        temp += "$last_dir='%s';" % (initial_guess_dir)

    temp += """
$Prefix="Lev";
"""

    temp += """
$dir = $Prefix.$Level;

RunInDirectory($dir,
                    {"Domain.input" => {
                                        "__Acx__" => $Acx,
                                        "__Acy__" => $Acy,
                                        "__Acz__" => $Acz,
                                        "__Bcx__" => $Bcx,
                                        "__Bcy__" => $Bcy,
                                        "__Bcz__" => $Bcz,
                                        "__Arexc__" => $Arexc,
                                        "__Brexc__" => $Brexc,
                                        "__Sep__" => $D
                     },
                     "EllipticItems.input" => {
                         "__Acx__" => $Acx,
                         "__Acy__" => $Acy,
                         "__Acz__" => $Acz,
                         "__Bcx__" => $Bcx,
                         "__Bcy__" => $Bcy,
                         "__Bcz__" => $Bcz,
                         "__OuterVx__" => $OuterVx,
                         "__OuterVy__" => $OuterVy,
                         "__OuterVz__" => $OuterVz,
                         "__MaxRad__" => $MaxRad,
                     },
                     "ExtraEllipticItems.input" => {
                         "__Acx__" => $Acx,
                         "__Acy__" => $Acy,
                         "__Acz__" => $Acz,
                         "__Bcx__" => $Bcx,
                         "__Bcy__" => $Bcy,
                         "__Bcz__" => $Bcz,
                         "__RadiusOfBH__" => $RadiusOfBH,
                         "__LASTDIR__" => $last_dir,
                         "__C__" => $C,
                         "__rA__" => $rA,
                         "__rB__" => $rB,
                     },
                     "EllipticObservers.input" => {
                     "__curTruncationImprovement__" => $TruncationImprovement,
                     "__curTruncationTarget__" => $curTruncationTarget,
                     "__curTruncationStop__" => $curTruncationStop,
                     "__nextTruncationImprovement__" => $TruncationImprovement,
                     "__nextTruncationTarget__" => $nextTruncationTarget,
                     "__nextTruncationStop__" => $nextTruncationStop,
                      },"""
    if not absolute:
        temp += """

                      "petsc.input" => {"__snes_max_it__" => $snes_max_it,
                                        "__ksprtol__" => $ksprtol,
                                        "__snesrtol__" => $snesrtol,
                                       },"""
    else:
        temp += """
                       "petsc.input" => {"__snes_max_it__" => $snes_max_it,
                                        "__kspatol__" => $kspatol,
                                        "__snesatol__" => $snesatol,
                                       },"""
    temp += """
                     "OrbitalParamItems.input" => {
                         "__OmegaOrbit__" => $Omega0,
                         "__AOmegarx__" => $OmegaAx,
                         "__AOmegary__" => $OmegaAy,
                         "__AOmegarz__" => $OmegaAz,
                         "__BOmegarx__" => $OmegaBx,
                         "__BOmegary__" => $OmegaBy,
                         "__BOmegarz__" => $OmegaBz,
                         "__adot0__" => $adot,
                     },
                     "Elliptic.input" => {
                         "__Obs__" => '""" + ObserveIfNotConverged + """',
                         "__Acx__" => $Acx,
                         "__Acy__" => $Acy,
                         "__Acz__" => $Acz,
                         "__Bcx__" => $Bcx,
                         "__Bcy__" => $Bcy,
                         "__Bcz__" => $Bcz,
                     },

                    }
        );
        $last_dir="../".$dir;

"""

    f = open('DoMultipleRuns.input', 'w')
    f.write(temp)
    f.close()


#=============================================================================
def ReadLastLineAsFloats(filename):
    if not os.path.exists(filename):
        error("File %s not found, terminating.\n"
              "Perhaps the EllipticSolver failed to converge and did not "
              "generate data?" % filename)

    f = open(filename)
    p = [float(x) for x in f.readlines()[-1].split()]
    f.close()
    return p


#=============================================================================
def ParseSpinMass_Ah(filename, res, surface):
    # parse spins and masses from file 'filename', and add them
    # into the dictionary 'res'
    h5dir = "{}.dir".format(surface)
    with ReadH5(filename) as f:
        res[surface + "-chi"] = f["{}/chiInertial.dat".format(h5dir)][0][1:]
        res[surface + "-M"] = f["{}/ChristodoulouMass.dat".format(h5dir)][0][1]
        res[surface + "-Mirr"] = f["{}/ArealMass.dat".format(h5dir)][0][1]


#=============================================================================
def ParseAdmIntegrals(filename, res):
    f = ReadLastLineAsFloats(filename)
    # NOTE: if no iterations, initial guess => Monopole(ConformalFactor)=0
    '''
    if not opts.no_iterations and abs(f[1]-f[8])/abs(f[1]+f[8]) > 0.2:
        error("E_adm and Monopole(ConformalFactor) are different\n"
              "For non-flat backgrounds, use E_ADM.\n"
              "E_ADM=%g, MonopoleConformalFactor=%g"%(f[1],f[8]))
    '''
    res["Eadm"] = f[8]
    res["Mk"] = 0.5 * (f[8] - f[9])  # because f[9] = Eadm-2Mk
    res["Padm"] = array([f[2], f[3], f[4]])
    res["Jadm"] = array([f[5], f[6], f[7]])


def SurfToCmd(Nprocs):
    return MpiCmd("SurfaceToSpatialCoordMapFiles", Nprocs)


def ApplyObsCmd(Nprocs):
    return MpiCmd("ApplyObservers", Nprocs)


def MpiCmd(Exec, Nprocs):
    """ApplyObservers command with correct bin dir and MPI wrapper"""
    ExecPath = "%s/%s" % (BinDir, Exec)
    if not os.path.exists(ExecPath):
        error("Could not find executable " + ExecPath)
    MPI = SystemOutput("%s/GetMpiCmd -p %d" % (BinDir, Nprocs))
    return "%s %s" % (MPI, ExecPath)


#=============================================================================
# Function used in WriteOutputIntoEvID
# Generates checkpoint files for initialization of Shape/Size maps
# Returns the horizon radii in the mapped frame
# INPUT: LAh = max L used in the evolution control systems
def FindAhCoefs(LAh, LBh):
    Ls = [LAh, LBh]
    print("# Entering FindAhCoefs")

    # Get the horizon coefficients: in evolution frame, horizons will be
    # spheres
    MappedHorizonRadii = {}
    i = 0
    for X in ["A", "B"]:
        OrigAhCoefs = "../private/EllData/Ah%s_L%02dCoefs.dat" % (X, Ls[i])
        copy(OrigAhCoefs, "ID_Ah%sCoefs.dat" % (X))

        # NOTE:
        # The value of rminfac here must agree with $rExcA/$rInitAhA
        # in DoMultipleRuns.input in the evolution.  Otherwise shape
        # control will not have Q=0 at the start of the evolution.
        rminfac = sqrt(opts.ExtrFrac)

        # Here we use the -p option, we do not specify rAH, and we
        # set the deriv order to zero.  This means that all of the resulting
        # coefficients are independent of -rmaxfac.  So we choose -rmaxfac
        # to 10 arbitrarily.  The coefficients *do* depend on rminfac,
        # and rminfac should be the same in the subsequent evolution as it
        # is here, or else the initial Q for shape control will not be zero.

        opt = " -l 6 -s 0 -1 -p -T 0 -rmaxfac 10 -rminfac %f " % rminfac
        cmd = SurfToCmd(1) + opt + OrigAhCoefs

        SystemOutput(cmd)

        infile = open("Init-GridParams.txt", 'r')
        text = infile.readlines()
        NumLines = len(text)
        if (NumLines != 1):
            error("Init-GridParams.txt should be only one line")
        match = re.search("rAH\s*=\s*([^;]*);", text[0])
        MappedHorizonRadii[X] = match.group(1)
        print("MappedHorizonRadii = ", MappedHorizonRadii[X], "\n")

        os.rename("Init-FuncLambdaFactor0.txt",
                  "ID_Init_FuncLambdaFactor%s0.txt" % X)
        os.rename("Init-FuncLambdaFactor.txt",
                  "ID_Init_FuncLambdaFactor%s.txt" % X)
        os.rename("Init-SmoothAhRadius.txt",
                  "ID_Init_SmoothAhRadius%s.txt" % X)
        os.remove("Init-GridParams.txt")
        os.remove("Init-ShapeTimeScales.txt")
        i += 1
    print("# Leaving FindAhCoefs")
    return MappedHorizonRadii


################################################################
def WriteOutputIntoEvID(IDType):
    ################################################################$
    print("# Entering WriteOutputIntoEvID")

    copy("../private/Extr/Domain.input", "./GrDomain.input")
    copy("../private/EllData/SpatialCoordMap.input", ".")
    for filename in glob("../private/Extr/Extrapolation/ForEvID/*.h5"):
        copy(filename, ".")

    # With AMR, the L of the 2 horizons might be different,
    # so search for the highest L for each horizon
    def GetLmax(H):
        Ah_files = glob("../private/EllData/Ah{}_L*.dat".format(H))
        Ah_Lmax = sorted(Ah_files)[-1]
        Lmax = int(re.findall('\d+', Ah_Lmax)[0])
        return Lmax

    LmaxA = GetLmax('A')
    LmaxB = GetLmax('B')

    # copy Ah_Coefs
    MappedHorizonRadii = FindAhCoefs(LmaxA, LmaxB)

    if IDType != "CFMS":
        # The radii in ThisItParams are wrong for SuperposedKerr since they
        # assume that the radius in the distorted frame and in the grid frame
        # are the same. The mapped radii are the correct values
        rA = MappedHorizonRadii["A"]
        rB = MappedHorizonRadii["B"]
    else:
        rA = ThisItParams["rA"]
        rB = ThisItParams["rB"]

    # Generate metadata from the final iteration (using this iterate's params)
    cxA = ThisItParams["cxA"]
    cyA = ThisItParams["cyA"]
    czA = ThisItParams["czA"]
    WriteMetadata(D, Omega0, adot, Eadm, MirrA, MirrB, MA, MB, chiA, chiB,
                  Padm, Jadm, rA, rB, cxA, cyA, czA)

    if IDType != "CFMS":
        print("#Computing vout")
        System("%s/ComputeIDVout --mapdir=%s" %(BinDir,EvID_Dir) +\
               " --domaindir=%s --datadir=%s" %(EvID_Dir,EvID_Dir) +\
               " --inputfilesdir=%s" %EvInputPath +\
               " --cA='%20.16f,%20.16f,%20.16f'" %(cxA,cyA,czA) +\
               " --cB='%20.16f,%20.16f,%20.16f'" %(cxA-D,cyA,czA))

    print("#Returning from WriteOutputIntoEvID")
    return


#=============================================================================
def WriteMetadata(D, Omega0, adot, Eadm, MirrA, MirrB, MA, MB, chiA, chiB,
                  Padm, Jadm, rA, rB, cxA, cyA, czA):

    rExcA = opts.ExtrFrac * float(rA)
    rExcB = opts.ExtrFrac * float(rB)

    cAx, cAy, cAz = [cxA, cyA, czA]
    cBx, cBy, cBz = [cxA - D, cyA, czA]
    chiAx, chiAy, chiAz = chiA
    chiBx, chiBy, chiBz = chiB
    Px, Py, Pz = Padm
    Jx, Jy, Jz = Jadm

    # Construct ID_Params.perl
    BBH_ID_WriteInputFiles.Metadata(
        cAx=cAx,
        cAy=cAy,
        cAz=cAz,
        cBx=cBx,
        cBy=cBy,
        cBz=cBz,
        Px=Px,
        Py=Py,
        Pz=Pz,
        Jx=Jx,
        Jy=Jy,
        Jz=Jz,
        chiAx=chiAx,
        chiAy=chiAy,
        chiAz=chiAz,
        chiBx=chiBx,
        chiBy=chiBy,
        chiBz=chiBz,
        ID_chiAMagnitude=norm(chiA),
        ID_chiBMagnitude=norm(chiB),
        ID_rA=rA,
        ID_rB=rB,
        ID_rExcA=rExcA,
        ID_rExcB=rExcB,
        ID_d=D,
        ID_Omega0=Omega0,
        ID_adot0=adot,
        ID_Eadm=Eadm,
        ID_MAhA=MirrA,
        ID_MAhB=MirrB,
        ID_MA=MA,
        ID_MB=MB,
        ID_Type='BBH_%s' % IDType,
        ID_Dir=ID_Dir,
        Revision=Revision,
    )


def measure_improvement(residual_hist):
    '''Given a record of residuals, compute the root-finding improvement'''
    sz = len(residual_hist)
    if len(residual_hist) < 4:
        return np.zeros(len(residual_hist[0]))

    improvement = np.zeros(len(residual_hist[0]))
    # Loop over every quantity, e.g. mA
    for i in range(len(residual_hist[0])):
        improvement[i] = sqrt( \
            sqrt(residual_hist[-4][i]*residual_hist[-3][i]) \
            /sqrt(residual_hist[-2][i]*residual_hist[-1][i]))

    return improvement


def parseAMRResults(target_dir,
                    preamble,
                    residual,
                    k,
                    impr,
                    base_lev=3,
                    eps=100.0,
                    max_lev=6):
    '''Check the output of the AMR driver and return the *absolute* path of
    the resolution changes to be used for the next iteration as well as
    the truncation error estimate'''

    pat = re.compile("(-\d*\.?\d*)")
    # First, check if the current changes exist
    curdir = os.getcwd()
    os.chdir(target_dir)
    curChanges = "curResolutionChanges.dat"
    curTermination = "curResolutionChangesTermination.dat"
    nextChanges = "nextResolutionChanges.dat"
    nextTermination = "nextResolutionChangesTermination.dat"

    if os.path.isfile(curChanges):
        # This means that the current resolution goal was *not* achieved
        # We do not need to change k
        fp = open(curChanges, 'r')
        preamble = "".join(fp.readlines())
        fp.close()

        with open(curTermination, "r") as fp:
            for line in fp:
                if "Estimated truncation error" in line:
                    trunc = re.findall(pat, line)

                    trunc = float(trunc[0])
                    break
        fp.close()

    else:
        # We have achieved the desired resolution
        # Is the root finding error small enough?
        if k < max_lev and (residual < eps * 10**(-(base_lev + k))
                            or impr < sqrt(1.5)):
            try:
                fp = open(nextChanges, 'r')
                preamble = "".join(fp.readlines())
                fp.close()
            except IOError:

                pass
            k = k + 1

        with open(nextTermination, "r") as fp:
            for line in fp:
                if "Estimated truncation error" in line:
                    trunc = re.findall(pat, line)
                    trunc = float(trunc[0])
                    break
        fp.close()

    # Always update params!
    update_params = True
    os.chdir(curdir)
    return trunc, k, preamble, update_params


#=============================================================================
def ExecuteSpells(lev, BinDir, cmd):
    '''Run Spells and look at the ouput of the 2 AMR observers
    to decide the next resolution'''

    # Copy and replace all the parameters in all the input files
    os.system("%s/DoMultipleRuns -n -B -e '' -1" % BinDir)
    os.chdir("Lev%d" % lev)
    # Run spells

    System("%s %s/Spells > spells.out" % (cmd, BinDir))
    target_dir = os.getcwd()
    return target_dir


def updateDomain(fnewDomain, oldDir):
    '''Update the domain for extrapolation with the changes in the last
    iteration'''
    fp = open(oldDir + "/Domain.input")
    txt = fp.read()

    changes = re.findall(r'(ResizeTheseSubdomains =.*;).*<<NONE>>',
                         txt,
                         flags=re.DOTALL)
    if changes:
        changes = changes[0]
    fp.close()

    # If there are changes, write them:

    if changes:
        with open(fnewDomain, 'r') as fw:
            newDomain = fw.read()

            # Reaplce the old ReizeTheseSubdomains if possible
        temp = re.sub('ResizeTheseSubdomains =.*HistoryFile=<<NONE>>',
                      changes + "\nHistoryFile=<<NONE>>",
                      newDomain,
                      flags=re.DOTALL)

        with open(fnewDomain, 'w') as fw:
            if temp == newDomain:
                fw.write(changes + newDomain)
            else:
                fw.write(temp)


def check_status(target_dir):
    '''Check if there are any Error*.txt files in target_dir. Then try to
    diagnoze some common problems:
    i) if ksp linear solve did not converge, increase resolution
    ii) if snes solve did not converge, ask for more non-linear iterations'''
    cur_dir = os.getcwd()
    os.chdir(target_dir)
    error_list = glob("Error*.txt")

    os.chdir(cur_dir)
    if error_list:

        return "Failure"
    else:
        return "Success"


def iterateDomain(repeat):
    '''Attempt to recover a failing run by increasing starting resolution'''
    # For spheres
    Nr = 10 + 4 * repeat
    L = 11 + 4 * repeat
    # For cylinders
    Nrho = 9 + 4 * repeat
    Nphi = 16 + 4 * repeat
    Nz = 9 + 4 * repeat
    # For perim blocks
    Ns = 9 + 4 * repeat
    skeleton = """
ResizeTheseSubdomains =
SphereA0(Extents=%d,%d,%d),
SphereB0(Extents=%d,%d,%d),
SphereC0(Extents=%d,%d,%d),
Cyl0(Extents=%d,%d,%d),
Cyl1(Extents=%d,%d,%d),
Cyl2(Extents=%d,%d,%d),
Cyl3(Extents=%d,%d,%d),
Cyl4(Extents=%d,%d,%d),
Perim0(Extents=%d,%d,%d),
Perim2(Extents=%d,%d,%d),
Perim4(Extents=%d,%d,%d),
;
""" % (Nr, L, 2 * L, Nr, L, 2 * L, Nr, L, 2 * L, Nrho, Nphi, Nz, Nrho, Nphi,
       Nz, Nrho, Nphi, Nz, Nrho, Nphi, Nz, Nrho, Nphi, Nz, Ns, Ns, Ns, Ns, Ns,
       Ns, Ns, Ns, Ns)
    return skeleton


#=============================================================================
def DoInitialDataSolve(Params=None,
                       outerV=None,
                       nprocs=None,
                       initial_guess_dir="",
                       dirname="bla",
                       D=None,
                       Omega0=None,
                       adot=None,
                       chiA=None,
                       chiB=None,
                       specID="__specID__",
                       debug=False,
                       IDType="CFMS",
                       lev=0,
                       q=1,
                       preamble='',
                       elliptic_opts=None,
                       amr_opts=None,
                       domain_opts=None,
                       SuperposedKerr_opts=None,
                       rad_factor=None,
                       ObserveIfNotConverged=False):
    """
    This function writes a set of input files for the elliptic
    solver, and then performs a DoMultipleRuns-call to run the
    elliptic solver.  It returns a number of useful things in
    a dictionary:
      Jadm, Eadm, Padm, Mk  -- ADM quantites
      Ah{,B}-{Mirr, M, chi}  -- data about AhA and AhB


    The input arguments are:
      dirname - directory in which the elliptic solve is run (if it exists,
                it will be removed)
      MaxLevels - run this many levels of refinement
    """

    ### SANITY CHECKS
    # All parameters should be specified
    for k, v in locals().items():
        if v is None:
            error("please specify " + k)

    original_directory = os.getcwd()

    if os.path.exists(dirname):
        error("Directory " + dirname + " already exists")

    os.mkdir(dirname)
    os.chdir(dirname)

    rA = Params['rA']
    rB = Params['rB']
    OmegaA = Params['OmegaA']
    OmegaB = Params['OmegaB']
    cxA = Params['cxA']
    cyA = Params['cyA']
    czA = Params['czA']

    if IDType != "CFMS":
        PrepareSuperposedKerrInputFiles(
            initial_guess_dir,
            lev,
            D,
            Omega0,
            adot,
            chiA,
            chiB,
            outerV,
            rA,
            rB,
            OmegaA,
            OmegaB,
            cxA,
            cyA,
            czA,
            debug,
            q=q,
            preamble=preamble,
            elliptic_opts=elliptic_opts,
            amr_opts=amr_opts,
            domain_opts=domain_opts,
            SuperposedKerr_opts=SuperposedKerr_opts,
            rad_factor=rad_factor)
    else:
        PrepareCFInputFiles(initial_guess_dir,
                            lev,
                            D,
                            Omega0,
                            adot,
                            outerV,
                            debug=debug,
                            q=q,
                            preamble=preamble,
                            elliptic_opts=elliptic_opts,
                            amr_opts=amr_opts,
                            domain_opts=domain_opts,
                            **Params)

    sleep(2)  # wait to allow FS to catch up

    # Run Spells
    cmd = SystemOutput(BinDir + "/GetMpiCmd -p %i" % (nprocs))

    target_dir = ExecuteSpells(lev, BinDir, cmd)
    os.chdir(target_dir)


    ###### EXTRACT DATA FROM OUTPUT OF ELLIPTIC SOLVER
    sleep(2)  # wait to allow FS to catch up
    status = check_status(target_dir)
    # parse output into dictionary 'res'
    res = {}
    if status == "Success":
        ParseAdmIntegrals("AdmIntegrals.dat", res)
        ParseSpinMass_Ah("Horizons.h5", res, "AhA")
        ParseSpinMass_Ah("Horizons.h5", res, "AhB")
        res['Padm'] = genfromtxt("PADM_c.dat")[1:4]
        res['CoM'] = genfromtxt("CoM.dat")[1:4]
        res['Jadm'] = genfromtxt("JADM_c.dat")[1:4]
    os.chdir(original_directory)
    return res, target_dir, status


################################################################
################################################################


# Initialization for iteration 0
def Initialize():
    # Make 'private' directory and chdir there
    System("rm -rf private")
    os.mkdir("private")
    os.chdir("private")

    # Make bin directory inside 'private' for SuperposedKerr runs
    System("%s/MakeBinDirectory -e %s -b %s" %
           (SourceBinDir, SourceSpells, SourceBinDir))

    # Write Params.out file
    paramsout = open('Parameters.dat', 'w')
    paramsout.write("# Beginning ID-rootfinding for q=%f, D=%f, Omega0=%f, "
                    "adot=%f \n" % (q, D, Omega0, adot))
    paramsout.write("""# [1] = its
# [2] = refinement level
# [3] = rA
# [4] = rB
# [5] = OmegaAx
# [6] = OmegaAy
# [7] = OmegaAz
# [8] = OmegaBx
# [9] = OmegaBy
# [10] = OmegaBz
# [11] = cxA
# [12] = cyA
# [13] = czA
# [14] = NumEs
""")
    paramsout.close()

    # Write initial physout file
    physout = open('PhysValues.dat', 'w')
    physout.write("# Beginning ID-rootfinding for q=%f, D=%f, Omega0=%f, "
                  "adot=%f \n" % (q, D, Omega0, adot))
    physout.write("""# [1] = its
# [2] = MA
# [3] = MB
# [4] = chiAx
# [5] = chiAy
# [6] = chiAz
# [7] = chiBx
# [8] = chiBy
# [9] = chiBz
# [10] = Px
# [11] = Py
# [12] = Pz
# [13] = Eadm
# [14] = Mk
# [15] = outerVx
# [16] = outerVy
# [17] = outerVz
""")
    physout.close()

    # Write initial errout file
    errout = open('Errors.dat', 'w')
    errout.write("""# Residuals of root-finder
# [1] = its
# [2] = Number of AMR its.
# [3] = err-massA
# [4] = err-massB
# [5] = err-chiA
# [6] = err-chiB
# [7] = err-Padm
# [8] = est. AMR trunc. error
# [9] = residual = rms(err-M{A,B},err-chi{A,B}, P)
"""
    )
    errout.close()

    # Set initial rA and rB
    if opts.rA is None:
        if IDType != "CFMS":
            rA = 0.86 * rPlusAFull
        else:
            rA = q / (1 + q) * RADIUS_OF_BH
    else:
        rA = opts.rA
    if opts.rB is None:
        if IDType != "CFMS":
            rB = 0.86 * rPlusBFull
        else:
            rB = 1 / (1 + q) * RADIUS_OF_BH
    else:
        rB = opts.rB

    # Set up initial OmegaA and OmegaB
    # This works for both small and large spins
    if IDType == "SKS":
        # See Eqs. (5) - (7) of 1506.01689
        rPlusA = (q / (1 + q)) * (1 + sqrt(1 - dot(chiAtarget, chiAtarget)))
        rPlusB = (1 / (1 + q)) * (1 + sqrt(1 - dot(chiBtarget, chiBtarget)))
        OmegaA = chiAtarget / (-2 * rPlusA)
        OmegaB = chiBtarget / (-2 * rPlusB)
    elif IDType == "SHK":
        # See Eqs. (5) - (7) of 1506.01689, but here rPlusA/rPlusB are in the
        # Cook-Scheel Harmonic coordinates
        rPlusA = (q / (1 + q)) * sqrt(1 - dot(chiAtarget, chiAtarget))
        rPlusB = (1 / (1 + q)) * sqrt(1 - dot(chiBtarget, chiBtarget))
        OmegaA = chiAtarget/(-2*(q/(1+q))*(1+ sqrt(1 \
                   - dot(chiAtarget,chiAtarget))))
        OmegaB = chiBtarget/(-2*(1/(1+q))*(1+ sqrt(1 \
            - dot(chiBtarget,chiBtarget))))
    else:
        OmegaA = -0.25 / (q / (1 + q)) * chiAtarget + array([0, 0, Omega0])
        OmegaB = -0.25 / (1 / (1 + q)) * chiBtarget + array([0, 0, Omega0])

    # But...if the spin is *very* large, you must be careful
    # not to guess too high an initial OmegaA or OmegaB
    # So, if spin is bigger than bigSpin,
    # then multiply the initial OmegaA or OmegaB by
    # rescaleFactor
    chiAtargetMag = sqrt(chiAtargetSq)
    chiBtargetMag = sqrt(chiBtargetSq)
    bigSpin = 0.98
    rescaleFactor = 0.95
    if chiAtargetMag > bigSpin:
        OmegaA *= rescaleFactor
    if chiBtargetMag > bigSpin:
        OmegaB *= rescaleFactor

    # Initial cxA and cyA
    cyA = 0
    cxA = D / (1 + q)
    czA = 0
    # Initialize history of rA and rB
    rAhistory = []
    rBhistory = []
    OmegaAhistory = []
    OmegaBhistory = []
    MAhistory = []
    MBhistory = []
    chiAhistory = []
    chiBhistory = []
    errMAhistory = []
    errMBhistory = []
    errchiAhistory = []
    errchiBhistory = []
    # Initialize CurrLevels
    CurrLevels = NextLevels

    # Initialize NumES
    NumES = 0  # number of elliptic solves at current resolution

    Params = PackParams(rA=rA,
                        rB=rB,
                        OmegaA=OmegaA,
                        OmegaB=OmegaB,
                        cxA=cxA,
                        cyA=cyA,
                        czA=czA)

    # Update Parameters.dat before first iteration with initial values
    WriteParamsDat(0, CurrLevels, NumES, **Params)

    # write initial J_old.h5 file
    Jaco = h5py.File('J_old.h5', 'w')
    Jaco.close()

    # write preamble.dat file
    prewrite = open('preamble.dat', 'w')
    prewrite.close()

    # write Residuals.dat file
    resiwrite = open('Residuals.dat', 'w')
    resiwrite.close()

    return Params,rAhistory,rBhistory,CurrLevels,NumES,OmegaAhistory, \
        OmegaBhistory,MAhistory,MBhistory,chiAhistory,chiBhistory, \
        errMAhistory,errMBhistory,errchiAhistory,errchiBhistory


#=============================================================================
def Restart(StartingIt):
    # StartingIt means the last iteration that has been done. The code
    # will restart from StartingIt+1 Chdir to 'private' directory.
    os.chdir("private")
    # Read the Parameters.dat file to get the params
    rAhistory = []
    rBhistory = []
    cxAhistory = []
    cyAhistory = []
    czAhistory = []
    OmegaAhistory = []
    OmegaBhistory = []
    itshistory = []
    Levelshistory = []
    NumEShistory = []

    MAhistory = []
    MBhistory = []
    chiAhistory = []
    chiBhistory = []
    errMAhistory = []
    errMBhistory = []
    errchiAhistory = []
    errchiBhistory = []
    outerVhistory = []

    paramsin = open('Parameters.dat', 'r')
    count = 0
    for line in paramsin:
        # Count from 0 to StartingIt. Because Parameters.dat
        # record data from its=0
        if line[0] != "#" and count <= StartingIt:
            cols = line.split()
            # Read its and Levels
            itshistory.append(int(cols[0]))
            LastPreviousIter = int(cols[0])
            Levelshistory.append(int(cols[1]))
            # Read initial rA and rB
            rAhistory.append(float(cols[2]))
            rBhistory.append(float(cols[3]))
            # Read OmegaA and OmegaB
            OmegaAtmp = [float(cols[4]), float(cols[5]), float(cols[6])]
            OmegaBtmp = [float(cols[7]), float(cols[8]), float(cols[9])]
            OmegaAhistory.append(OmegaAtmp)
            OmegaBhistory.append(OmegaBtmp)
            # Read cxA and cyA
            cxAhistory.append(float(cols[10]))
            cyAhistory.append(float(cols[11]))
            czAhistory.append(float(cols[12]))
            # Read NumES
            NumEShistory.append(int(cols[13]))
            count += 1
    paramsin.close()
    count = 0
    physin = open('PhysValues.dat', 'r')
    for line in physin:
        # Count from 0 to StartingIt-1. Because PhysValues.dat record
        # data from its=1
        if line[0] != "#" and count <= StartingIt - 1:
            cols = line.split()
            # Read MA errMA and MB errMB
            MAhistory.append(float(cols[1]))
            errMAhistory.append(float(cols[1]) - q / (1. + q))
            MBhistory.append(float(cols[2]))
            errMBhistory.append(float(cols[2]) - 1. / (1. + q))
            # Read chiA and chiB
            chiAtmp = [float(cols[3]), float(cols[4]), float(cols[5])]
            chiBtmp = [float(cols[6]), float(cols[7]), float(cols[8])]
            errchiAtmp=[float(cols[3])-chiAtarget[0],float(cols[4])-chiAtarget[1],\
             float(cols[5])-chiAtarget[2]]
            errchiBtmp=[float(cols[6])-chiBtarget[0],float(cols[7])-chiBtarget[1],\
            float(cols[8])-chiBtarget[2]]
            chiAhistory.append(chiAtmp)
            chiBhistory.append(chiBtmp)
            errchiAhistory.append(errchiAtmp)
            errchiBhistory.append(errchiBtmp)
            # Read outerV
            outerVtmp = [float(cols[14]), float(cols[15]), float(cols[16])]
            outerVhistory.append(outerVtmp)
            count += 1
    physin.close()
    print(StartingIt, len(rAhistory))
    # Compute initial rA and rB
    rA = rAhistory[StartingIt]
    rB = rBhistory[StartingIt]

    # Compute initial OmegaA and OmegaB
    OmegaA = OmegaAhistory[StartingIt]
    OmegaB = OmegaBhistory[StartingIt]

    # Compute initial cxA and cyA
    cxA = cxAhistory[StartingIt]
    cyA = cyAhistory[StartingIt]
    czA = czAhistory[StartingIt]
    # Initialize CurrLevels
    CurrLevels = Levelshistory[StartingIt]

    # Initialize NumES
    NumES = NumEShistory[StartingIt]

    # Initialize OuterV
    # Because PhysValues.dat record data from its=1, the index of
    # outerVini should be StartingIt-1 instead of StartingIt
    outerVini = outerVhistory[StartingIt - 1]

    # Initialize residuals
    if StartingIt == 1:
        residualslist = []
        residuals = np.loadtxt('Residuals.dat').tolist()
        residualslist.append(residuals)
    else:
        residuals = np.loadtxt('Residuals.dat')[0:StartingIt]
        residualslist = residuals.tolist()

    Params = PackParams(rA=rA,
                        rB=rB,
                        OmegaA=OmegaA,
                        OmegaB=OmegaB,
                        cxA=cxA,
                        cyA=cyA,
                        czA=czA)
    # rAhistory, rBhistory, OmegaAhistory, OmegaBhistory are updated in
    # the current iteration (in UpdatingFormulae). So they should take
    # data until the last iteration
    return Params,rAhistory[:-1],rBhistory[:-1],CurrLevels,NumES,\
 LastPreviousIter,OmegaAhistory[:-1],OmegaBhistory[:-1],\
 MAhistory,MBhistory,chiAhistory,chiBhistory,errMAhistory,\
 errMBhistory,errchiAhistory,errchiBhistory,outerVini,residualslist


def compute_outerV(outerV, Padm, CoM, Omega0, MA, MB, errMA, errMB, D, cxA,
                   cyA, czA):
    '''Update outerV such that Padm->0'''
    errM = errMA + errMB
    M = MA + MB
    Omega_orb = np.array([0.0, 0.0, Omega0])
    Sep = array([D, 0.0, 0.0])
    cA = array([cxA, cyA, czA])
    dcxA, dcyA, dczA = compute_COM_update(CoM, cxA, cyA, czA, MA, MB, errMA,
                                          errMB, D)
    delta_cA = array([dcxA, dcyA, dczA])
    set1 = array([
        cxA, cyA, czA, dcxA, dcyA, dczA, Padm[0], Padm[1], Padm[2], errMA,
        errMB
    ])
    fp = open("OuterV_Diag.dat", "a")
    np.savetxt(fp, set1)
    fp.close()
    return outerV - Padm/M + errM/M*(outerV+ cross(Omega_orb,(cA-CoM))) \
        - cross(Omega_orb,CoM) - cross(Omega_orb,delta_cA) \
        - errMB/M*cross(Omega_orb,Sep)


def fitFuncR(p, x):
    M = x[:, 0]
    chi = x[:, 1:]

    chisq = np.array([dot(chi[i], chi[i]) for i in range(len(chi))])
    return p[0] * M * (1 + sqrt(1 - chisq))


def fitFuncOmega(p, x):
    r = x[:, 0]
    chi = x[:, 1:]
    chisq = np.array([norm(o) for o in chi])
    return p[0] * chisq / (2 * r)


def errFunc(p, x, y, tp):
    if tp == "R":
        return fitFuncR(p, x) - y
    else:
        return fitFuncOmega(p, x) - y


def compute_COM_update(CoM, cxA, cyA, czA, MA, MB, errMA, errMB, D):
    dcxA = -CoM[0] + (-MA * errMB + MB * errMA) / (MA + MB)**2 * D
    dcyA = -CoM[1]
    dczA = -CoM[2]

    return dcxA, dcyA, dczA


def Jacobian_SKS(rA, OmegaA, rB, OmegaB):
    '''Given the radius and the angular velocity of the horizon, compute the
    Jacobian, assuming an isolated Kerr-Schild hole'''
    J = np.zeros((8, 8))
    tempA = 1 + sqrt(1 - 4 * rA**2 * dot(OmegaA, OmegaA))
    dMAdrA = 1/tempA + (4*rA**2*dot(OmegaA,OmegaA)) \
            /(tempA**2*sqrt(1-4*rA**2*dot(OmegaA,OmegaA)))
    dMAdrB = 0
    dMAdOmegaA_1 = (4*rA**3*OmegaA[0]) \
        /(tempA**2*sqrt(1-4*rA**2*dot(OmegaA,OmegaA)))
    dMAdOmegaA_2 = (4*rA**3*OmegaA[1]) \
        /(tempA**2*sqrt(1-4*rA**2*dot(OmegaA,OmegaA)))
    dMAdOmegaA_3 = (4*rA**3*OmegaA[2]) \
        /(tempA**2*sqrt(1-4*rA**2*dot(OmegaA,OmegaA)))
    dMAdOmegaB_1 = 0
    dMAdOmegaB_2 = 0
    dMAdOmegaB_3 = 0

    dChiA1drA = -2 * OmegaA[0]
    dChiA2drA = -2 * OmegaA[1]
    dChiA3drA = -2 * OmegaA[2]
    dChiA1drB = 0
    dChiA2drB = 0
    dChiA3drB = 0
    dChiA1dOmegaA_1 = -2 * rA
    dChiA2dOmegaA_2 = -2 * rA
    dChiA3dOmegaA_3 = -2 * rA
    dChiA1dOmegaB_1 = 0
    dChiA2dOmegaB_2 = 0
    dChiA3dOmegaB_3 = 0

    tempB = 1 + sqrt(1 - 4 * rB**2 * dot(OmegaB, OmegaB))
    dMBdrB = 1/tempB + (4*rB**2*dot(OmegaB,OmegaB)) \
        /(tempB**2*sqrt(1-4*rB**2*dot(OmegaB,OmegaB)))

    dMBdOmegaB_1 = (4*rB**3*OmegaB[0]) \
        /(tempB**2*sqrt(1-4*rB**2*dot(OmegaB,OmegaB)))
    dMBdOmegaB_2 = (4*rB**3*OmegaB[1]) \
        /(tempB**2*sqrt(1-4*rB**2*dot(OmegaB,OmegaB)))
    dMBdOmegaB_3 = (4*rB**3*OmegaB[2]) \
        /(tempB**2*sqrt(1-4*rB**2*dot(OmegaB,OmegaB)))

    dChiB1drB = -2 * OmegaB[0]
    dChiB2drB = -2 * OmegaB[1]
    dChiB3drB = -2 * OmegaB[2]

    dChiB1dOmegaB_1 = -2 * rB
    dChiB2dOmegaB_2 = -2 * rB
    dChiB3dOmegaB_3 = -2 * rB

    J[0, 0] = dMAdrA
    #J[0,1] = dMAdrB
    J[0, 2] = dMAdOmegaA_1
    J[0, 3] = dMAdOmegaA_2
    J[0, 4] = dMAdOmegaA_3
    #J[0,5] = dMAdOmegaB_1
    #J[0,6] = dMAdOmegaB_2
    #J[0,7] = dMAdOmegaB_3

    #J[1,0] = dMBdrA
    J[1, 1] = dMBdrB
    #J[1,2] = dMBdOmegaA_1
    #J[1,3] = dMBdOmegaA_2
    #J[1,4] = dMBdOmegaA_3
    J[1, 5] = dMBdOmegaB_1
    J[1, 6] = dMBdOmegaB_2
    J[1, 7] = dMBdOmegaB_3

    J[2, 0] = dChiA1drA
    J[2, 2] = dChiA1dOmegaA_1

    J[3, 0] = dChiA2drA
    J[3, 3] = dChiA2dOmegaA_2

    J[4, 0] = dChiA3drA
    J[4, 4] = dChiA3dOmegaA_3

    J[5, 1] = dChiB1drB
    J[5, 5] = dChiB1dOmegaB_1

    J[6, 1] = dChiB2drB
    J[6, 6] = dChiB2dOmegaB_2

    J[7, 1] = dChiB3drB
    J[7, 7] = dChiB3dOmegaB_3
    return J


def Cubic(a, b, c):
    """ Solves the cubic equation x^3 + a x^2 + b x + c = 0.
    Follows Chapter 5.6 of Numerical recipes (3rd Edition).
    Assumes a,b,c are real.
    Returns only the real roots.
    """

    Q = (a**2. - 3. * b) / 9.
    R = (2. * a**3. - 9. * a * b + 27. * c) / 54.

    ## Will have three real roots
    if R**2 < Q**3:
        theta = np.arccos(R / np.sqrt(Q**3))
        x1 = -2 * np.sqrt(Q) * np.cos(theta / 3.) - a / 3.
        x2 = -2 * np.sqrt(Q) * np.cos(theta / 3. + 2 * np.pi / 3.) - a / 3.
        x3 = -2 * np.sqrt(Q) * np.cos(theta / 3. - 2 * np.pi / 3.) - a / 3.
        real_roots = [x1, x2, x3]
    ## Only one real root
    else:
        A = -np.sign(R) * (np.abs(R) + np.sqrt(R**2 - Q**3))**(1. / 3)
        if A == 0:
            B = 0
        else:
            B = Q / A
        x1 = A + B - a / 3.
        real_roots = [x1]

    return np.array(real_roots)


def HarmonicKerr_MassSpinFromROmega(r, Omega, M_guess):
    """ Computes Mass and dimensionless spin of a Harmonic Kerr BH
    given its horizon radius and horizon frequency analogous to
    Eqs.(7-8) of arxiv:1506.01689. In the case of Harmonic Kerr
    we need to solve a cubic equation to obtain the Mass, this
    can have multiple roots (2 valid typically), but one of them
    tends to have a high mass and spin, so we choose the one
    that is closest to some initial guess. This guess should be
    the desired mass of the BH.
    """

    OmegaSq = np.sum(np.array(Omega)**2)
    if OmegaSq == 0:
        # this is trivial, but will lead to 1/0 in the general case
        M = r
        chi = 0 * Omega
    else:
        M_roots = Cubic(r, -1. / 4. / OmegaSq, r / 4. / OmegaSq)
        M_roots = M_roots[M_roots > 0]  # pick only positive roots
        # pick root closest to initial guess
        idx = np.argmin(np.abs(M_roots - M_guess))
        M = M_roots[idx]
        chi = -2 * (r + M) * np.array(Omega)

    return M, chi


def Jacobian_SBH_Harmonic(r, Omega, M_guess):
    '''Given the radius and the angular velocity of the horizon, compute the
    Jacobian, assuming an isolated Harmonic Kerr hole.
    Here Jac^a_b = dX^a/dx^b where X = {M, chi} and x = {r, Omega}

    NOTE: Within this function, r corresponds to the outer horizon radius in
    Cook-Scheel Harmonic coordinates. The Kerr-Schild outer horizon radius
    will be explicitly identified as r_{KS}.

    We start with the Kerr-Schild relations between {M, chi} and
    {r^{KS}_{+}, Omega}, see Eqs. (7) and (8) of 1506.01689
        M = r_{KS}/(1 + \sqrt{1 - 4 r^2_{KS} Omega^2 })
        chi = - 2 r_{KS} Omega

    Now, Harmonic r is related to Kerr-Schild r_{KS} by r = r_{KS} - M.
    So, we get
        M = r/(\sqrt{1 - 4 (r+M)^2 Omega^2 })       --- (1)
        chi = - 2 (r + M) Omega

    Unlike the Kerr-Schild case, it's not straightforward to write {M,chi}
    as explicit expressions of {r, Omega}, making the Jacobian computation
    messy. Instead, it is easy to write {r, Omega} as explicit expressions of
    {M, chi}. Rewriting Eq (1) above:
        M = r/(\sqrt{1-chi^2})
        chi = -2 (r + M) Omega
    Finally,
        r = M \sqrt{1-chi^2}
        Omega = -chi/2M/(1+\sqrt{1-chi^2})

    So to compute the Jacobian, Jac^a_b = dX^a/dx^b where X = {M, chi} and
    x = {r, Omega}, we first compute the inverse Jacobian, which is
    straightforward. Then we invert the inverse Jacobian.
    '''

    # The inverse Jacobian is much easier to compute, so we compute that first
    # and invert it
    InvJac = np.zeros((4, 4))
    M, chi = HarmonicKerr_MassSpinFromROmega(r, Omega, M_guess)

    ## Square Root of One Minus Chi Squared
    SROMCS = np.sqrt(1 - np.sum(np.array(chi)**2))

    # dr/dM
    InvJac[0, 0] = SROMCS
    # dr/dchix
    InvJac[0, 1] = -M * chi[0] / SROMCS
    # dr/dchiy
    InvJac[0, 2] = -M * chi[1] / SROMCS
    # dr/dchiz
    InvJac[0, 3] = -M * chi[2] / SROMCS

    # dOmegax/dM
    InvJac[1, 0] = chi[0] / (2. * M**2 * (1 + SROMCS))
    # dOmegax/dchix
    InvJac[1,1] = (-1 - SROMCS + chi[1]**2 + chi[2]**2) \
        /(2.*M*SROMCS*(1+SROMCS)**2)
    # dOmegax/dchiy
    InvJac[1, 2] = -(chi[0] * chi[1]) / (2. * M * SROMCS * (1 + SROMCS)**2)
    # dOmegax/dchiz
    InvJac[1, 3] = -(chi[0] * chi[2]) / (2. * M * SROMCS * (1 + SROMCS)**2)

    # dOmegay/dM
    InvJac[2, 0] = chi[1] / (2. * M**2 * (1 + SROMCS))
    # dOmegay/dchix
    InvJac[2, 1] = -(chi[0] * chi[1]) / (2. * M * SROMCS * (1 + SROMCS)**2)
    # dOmegay/dchiy
    InvJac[2,2] = (-1 - SROMCS + chi[0]**2 +chi[2]**2) \
        /(2.*M*SROMCS*(1+SROMCS)**2)
    # dOmegay/dchiz
    InvJac[2, 3] = -(chi[1] * chi[2]) / (2. * M * SROMCS * (1 + SROMCS)**2)

    # dOmegaz/dM
    InvJac[3, 0] = chi[2] / (2. * M**2 * (1 + SROMCS))
    # dOmegaz/dchix
    InvJac[3, 1] = -(chi[0] * chi[2]) / (2. * M * SROMCS * (1 + SROMCS)**2)
    # dOmegaz/dchiy
    InvJac[3, 2] = -(chi[1] * chi[2]) / (2. * M * SROMCS * (1 + SROMCS)**2)
    # dOmegaz/dchiz
    InvJac[3,3] = (-1 - SROMCS + chi[0]**2 + chi[1]**2) \
        /(2.*M*SROMCS*(1+SROMCS)**2)

    # Invert to get the Jacobian
    return inv(InvJac)


def Jacobian_SHK(rA, OmegaA, rB, OmegaB, mA_target, mB_target):
    '''Given the radii and the horizon angular velocities of the two BHs,
    compute the Jacobian, assuming isolated Harmonic Kerr holes. Here, the
    Jacobian is for the transformation from {r,Omega} to {M,chi}.'''

    JacA = Jacobian_SBH_Harmonic(rA, OmegaA, mA_target)
    JacB = Jacobian_SBH_Harmonic(rB, OmegaB, mB_target)
    J = np.zeros((8, 8))

    # dMAdrA
    J[0, 0] = JacA[0, 0]
    # dMAdOmegaA_i
    J[0, 2] = JacA[0, 1]
    J[0, 3] = JacA[0, 2]
    J[0, 4] = JacA[0, 3]

    # dChiA1drA
    J[2, 0] = JacA[1, 0]
    # dChiA1dOmegaA_i
    J[2, 2] = JacA[1, 1]
    J[2, 3] = JacA[1, 2]
    J[2, 4] = JacA[1, 3]

    # dChiA2drA
    J[3, 0] = JacA[2, 0]
    # dChiA2dOmegaA_i
    J[3, 2] = JacA[2, 1]
    J[3, 3] = JacA[2, 2]
    J[3, 4] = JacA[2, 3]

    # dChiA2drA
    J[4, 0] = JacA[3, 0]
    # dChiA2dOmegaA_i
    J[4, 2] = JacA[3, 1]
    J[4, 3] = JacA[3, 2]
    J[4, 4] = JacA[3, 3]

    # dMBdrB
    J[1, 1] = JacB[0, 0]
    # dMBdOmegaB_i
    J[1, 5] = JacB[0, 1]
    J[1, 6] = JacB[0, 2]
    J[1, 7] = JacB[0, 3]

    # dChiB1drB
    J[5, 1] = JacB[1, 0]
    # dChiB1dOmegaB_i
    J[5, 5] = JacB[1, 1]
    J[5, 6] = JacB[1, 2]
    J[5, 7] = JacB[1, 3]

    # dChiB2drB
    J[6, 1] = JacB[2, 0]
    # dChiB2dOmegaB_i
    J[6, 5] = JacB[2, 1]
    J[6, 6] = JacB[2, 2]
    J[6, 7] = JacB[2, 3]

    # dChiB2drB
    J[7, 1] = JacB[3, 0]
    # dChiB2dOmegaB_i
    J[7, 5] = JacB[3, 1]
    J[7, 6] = JacB[3, 2]
    J[7, 7] = JacB[3, 3]

    return J


#=============================================================================
def UpdatingFormulae(rA,
                     rB,
                     OmegaA,
                     OmegaB,
                     cxA,
                     cyA,
                     czA,
                     CoM=None,
                     OmegaAhistory=None,
                     OmegaBhistory=None,
                     MAhistory=None,
                     MBhistory=None,
                     chiAhistory=None,
                     chiBhistory=None,
                     errMAhistory=None,
                     errMBhistory=None,
                     errchiAhistory=None,
                     errchiBhistory=None,
                     mA_target=None,
                     mB_target=None,
                     chiA_target=None,
                     chiB_target=None,
                     high_spin_update=True,
                     J_old=None,
                     currLevels=0):

    rAhistory.append(rA)
    rBhistory.append(rB)
    OmegaAhistory.append(OmegaA)
    OmegaBhistory.append(OmegaB)
    MAhistory.append(MA)
    MBhistory.append(MB)
    chiAhistory.append(chiA)
    chiBhistory.append(chiB)

    #errMAhistory.append(errMA)
    #errMBhistory.append(errMB)
    #errchiAhistory.append(errchiA)
    #errchiBhistory.append(errchiB)

    if IDType != "CFMS":
        if not high_spin_update:

            # By default, high_spin_update is always True. This needs
            # to get updated if that changes.
            if IDType == "SHK":
                error("Not implemented for SHK")

            sqrtA = sqrt(1 - dot(chiA, chiA))

            rA = rA - errMA * (1 + sqrtA) + MA * dot(chiA, errchiA) / sqrtA
            OmegaA = OmegaA  - errMA*chiA/(2*MA*MA*(1+sqrtA)) \
                + errchiA / (2*MA*(1+sqrtA)) \
                + chiA*dot(chiA, errchiA) / (2*MA*(1+sqrtA)**2*sqrtA)

            sqrtB = sqrt(1 - dot(chiB, chiB))
            rB = rB - errMB * (1 + sqrtB) + MB * dot(chiB, errchiB) / sqrtB
            OmegaB = OmegaB - errMB*chiB/(2*MB*MB*(1+sqrtB)) \
                + errchiB / (2*MB*(1+sqrtB)) \
                + chiB*dot(chiB, errchiB) / (2*MB*(1+sqrtB)**2*sqrtB)
        else:

            F = np.array([
                errMA, errMB, errchiA[0], errchiA[1], errchiA[2], errchiB[0],
                errchiB[1], errchiB[2]
            ])

            if len(rAhistory) < 2:
                if IDType == "SKS":
                    # If we are just starting, use analytic Kerr-Schild
                    # Jacobian
                    J = Jacobian_SKS(rA, OmegaA, rB, OmegaB)
                elif IDType == "SHK":
                    # If we are just starting, use analytic Harmonic-Kerr
                    # Jacobian
                    J = Jacobian_SHK(rA, OmegaA, rB, OmegaB, mA_target,
                                     mB_target)

                J_inverse = inv(J)
            else:
                # Update via the Sherman-Morrison formula
                J_inverse = 1.0 * J_old
                derrMA = errMA - errMAhistory[-2]
                derrMB = errMB - errMBhistory[-2]
                derrchiA = errchiA - errchiAhistory[-2]
                derrchiB = errchiB - errchiBhistory[-2]
                deltaF = array([
                    derrMA, derrMB, derrchiA[0], derrchiA[1], derrchiA[2],
                    derrchiB[0], derrchiB[1], derrchiB[2]
                ])
                delta_x = array([
                    rA - rAhistory[-2], rB - rBhistory[-2],
                    OmegaA[0] - OmegaAhistory[-2][0],
                    OmegaA[1] - OmegaAhistory[-2][1],
                    OmegaA[2] - OmegaAhistory[-2][2],
                    OmegaB[0] - OmegaBhistory[-2][0],
                    OmegaB[1] - OmegaBhistory[-2][1],
                    OmegaB[2] - OmegaBhistory[-2][2]
                ])

                u = dot(J_inverse, F)
                c = dot(delta_x, (delta_x + u))
                J_inverse = J_inverse - 1. / c * dot(np.outer(u, delta_x),
                                                     J_inverse)

            delta_param = -1. * dot(J_inverse, F)
            drA = delta_param[0]
            drB = delta_param[1]
            dOmegaA = delta_param[2:5]
            dOmegaB = delta_param[5:]

            rA = rA + drA
            rB = rB + drB
            OmegaA = OmegaA + dOmegaA
            OmegaB = OmegaB + dOmegaB

            J_old = J_inverse

    else:
        ### CONFORMALLY FLAT
        if opts.old_radius_update:
            rA = rA - errMA * RADIUS_OF_BH
            rB = rB - errMB * RADIUS_OF_BH
        else:
            rA = rA - errMA * rA / MA
            rB = rB - errMB * rB / MB
            # horizon frequency update (vectorial equations)
        OmegaA = OmegaA + errchiA / (4 * MA) - errMA / (4 * MA * MA) * chiA
        OmegaB = OmegaB + errchiB / (4 * MB) - errMB / (4 * MB * MB) * chiB

    if Omega0 > 1.e-10:

        dcxA, dcyA, dczA = compute_COM_update(CoM, cxA, cyA, czA, MA, MB,
                                              errMA, errMB, D)
        cxAt = cxA + dcxA
        cyAt = cyA + dcyA
        czAt = czA + dczA

        # If we predict the shift to be wrong, DO NOT move COM
        if cxAt <= 0 or cxAt - D >= 0:
            pass
        else:
            cxA = 1.0 * cxAt
            cyA = 1.0 * cyAt
            czA = 1.0 * czAt
    else:
        cxA = cxA + (cxA * errMA + (cxA - D) * errMB) / (MA + MB)
        cyA = cyA + cyA * (errMA + errMB) / (MA + MB)

    if opts.shanks_transform_for_mass and NumES == 4:
        temp = rAhistory
        rA = (temp[-3] * temp[-1] - temp[-2]**2) / (temp[-3] - 2. * temp[-2] +
                                                    temp[-1])
        temp = rBhistory
        rB = (temp[-3] * temp[-1] - temp[-2]**2) / (temp[-3] - 2. * temp[-2] +
                                                    temp[-1])
        print("Applying Shanks transformation to excision radii.")
        print("Please verify enhanced convergence by plotting " 
              "private/Errors.dat.\n")

    if rA < 0 or rB < 0 or cxA < 0 or cxA - D > 0:
        error("""
Iter%02i resulted in update that caused wrong signs in some parameters:
   rA=%g  (should be positive)
   rB=%g  (should be positive)
   cxA=%g (should be positive)
   cxB=cxA-D=%g (should be negative)
Check the previous Iter?? for convergence, in particular check that the
AH-masses are reasonable.  Perhaps it helps to increase the resolution
demanded of AMR by using the option --base_lev and setting to >3.
            """ % (its, rA, rB, cxA, cxA - D))

    # Params for the next iteration
    Params = PackParams(rA=rA,
                        rB=rB,
                        OmegaA=OmegaA,
                        OmegaB=OmegaB,
                        cxA=cxA,
                        cyA=cyA,
                        czA=czA)

    # Update Parameters.dat for new parameters, before next iteration.
    WriteParamsDat(its, currLevels, NumES, **Params)

    return Params, J_old


#=============================================================================
# Used for packing up (rA,rB,OmegaA,OmegaB,cxA,cyA) into a dictionary
def PackParams(**kwargs):
    return kwargs


def WriteParamsDat(its, CurrLevels, NumES, rA, rB, OmegaA, OmegaB, cxA, cyA,
                   czA):
    paramsout = open('Parameters.dat', 'a')
    paramsout.write(" %i %i %.15g %.15g %.15g %.15g %.15g %.15g"\
                    " %.15g %.15g %.15g %.15g %.15g %i\n"%
                    (its, CurrLevels, rA, rB,
                     OmegaA[0], OmegaA[1], OmegaA[2],
                     OmegaB[0], OmegaB[1], OmegaB[2],
                     cxA, cyA,czA, NumES))
    paramsout.close()


# Subroutine to safely compute acos, keeping acosarg within the valid range
def SafeaCos(acosarg):
    if acosarg > 1:
        return 0.0
    elif acosarg < -1:
        return 180.0
    else:
        return acos(acosarg) * 180 / pi


################################################################
################################################################
################################################################
#### MAIN ROUTINE
################################################################
################################################################
################################################################

if __name__ == "__main__":
    StartTime = datetime.now()

    # Program options

    usage = """%(prog)s [opts]
Runs several Conformal Thin-Sandwich initial data solves to perform
root-finding to generate a final initial data set which has a total
mass (sum of Christoudoulou masses) of unity, a certain mass-ratio,
coordinate separation 'D', orbital frequency 'Omega0' and radial
expansion factor 'adot'.

This script assumes that PrepareID was used to set up the initial
data directory (i.e. that 'Ev' and 'bin' directories exist in the
appropriate place; see PrepareID for details).
"""

    p = ap.ArgumentParser(usage=usage)

    # Physical parameter options
    pp = p.add_argument_group("Physical parameter options")
    pp.add_argument("--q",
                    type=float,
                    required=True,
                    help="Desired mass-ratio, MA/MB")
    pp.add_argument("--chiA",
                    required=True,
                    help="Desired S/M^2 on BH A.  Specified as 'Sx,Sy,Sz'")
    pp.add_argument("--chiB",
                    required=True,
                    help="Desired S/M^2 on BH B.  Specified as 'Sx,Sy,Sz'")
    pp.add_argument("--D",
                    type=float,
                    required=True,
                    help="Desired coordinate separation D")
    pp.add_argument("--Omega0",
                    type=float,
                    required=True,
                    help="Desired orbital frequency")
    pp.add_argument("--adot", type=float, required=True, help="Desired adot")
    pp.add_argument(
        "--rA",
        type=float,
        help="Initial radius of sphere A. Default = q/(1+q)*RADIUS_OF_BH")
    pp.add_argument(
        "--rB",
        type=float,
        help="Initial radius of sphere B. Default = 1/(1+q)*RADIUS_OF_BH")
    pp.add_argument(
        "--IDType",
        default="CFMS",
        type=str,
        help="Initial Data type: 'SKS', 'SHK' or 'CFMS'. Default 'CFMS'")
    pp.add_argument(
        "--useNegExpBC",
        action='store_true',
        default=False,
        help="If used, place the boundary inside the horizon and avoid "
        "extrapolation.")
    pp.add_argument(
        "--ExtrFrac",
        type=float,
        default=0.89,
        help="Relative distance to extrapolate inward to produce an "
        "excision boundary for evolution (default: %(default)s). "
        "If useNegExpBC is True, will set the inner boundary "
        "with same ExtrFrac relative to apparent horizon radius, "
        "but will not need to extrapolate.")

    # Elliptic solver options
    pe = p.add_argument_group("Elliptic solver options")
    pe.add_argument("--ksp_rtol",
                    type=float,
                    help="The tolerance for the"
                    "petsc linear solver")
    pe.add_argument("--snes_rtol",
                    type=float,
                    help="The tolerance for the"
                    " petsc nonlinear solver")

    pe.add_argument(
        "--ksp_atol_base",
        default=4,
        type=float,
        help="Determines ksp_atol via ksp_atol = 1/5 * 10^{-(kspt_atol_base+k)}")
    pe.add_argument(
        "--snes_atol_base",
        default=4,
        type=float,
        help="Determines snes_atol via snes_atol = 10^{-(snes_atol_base+k)}")

    # AMR options
    pa = p.add_argument_group("AMR options")

    pa.add_argument(
        "--base_lev",
        default=3,
        type=int,
        help="Determines the initial expectation for AMR error via "
        "trunc_err = 10^{-(base_lev+k)}")
    pa.add_argument("--improvement_indicator", default="max", type=str,
        help="Determines the overall residual improvement indicator used for AMR. "
        "For each of the relevant physical parameters the improvement of the residual "
        "in the last iteration is computed. The overall improvement indicator is constructed "
        "from the values in this list. "
        "Possible options: 'max' uses the largest value in the list of improvements. "
        "'rms' uses the root mean square of all values in the list of improvments.")

    # Parameter updating options
    pu = p.add_argument_group("Parameter updating options")
    pu.add_argument(
        "--shanks_transform_for_mass",
        default=False,
        action="store_true",
        help="If specified, after four iterations on "
        "same Lev, apply a Shanks transform to the excision radii. This "
        "should enhance convergence.")
    pu.add_argument(
        "--old_radius_update",
        default=False,
        action="store_true",
        help="If specified, use the ratio of mass/radius from the last "
        "iteration when determining the update of the radius.  If not "
        "specified use the single BH value 0.85.")

    pu.add_argument(
        "--high_spin_update",
        default=True,
        action="store_true",
        help="If True, use updating formulas that may accelerate convergence "
        "for high spins/mass ratios. WARNING: not tested!")

    pu.add_argument(
        "--rad_factor",
        type=float,
        default=1e3,
        help="The radius of the cut-off for PADM calculation, *in units of "
        "the larger of SuperposedKerr Gaussian widths*")

    # Logistical options
    pl = p.add_argument_group("Logistical options")
    pl.add_argument(
        "--nprocs",
        type=int,
        required=True,
        help="Number of processors to run elliptic solver (1, 6, or 11 "
        "recommended).")
    pl.add_argument("--levels",
                    type=int,
                    required=True,
                    help="Finest refinement level")
    pl.add_argument("--its",
                    type=int,
                    default=30,
                    help="Maximum number of iterations (default: %(default)s)")
    pl.add_argument("--target_residual",type=float, default=1e-7,
                    help="Target residual of the parameter-adjustment. "
                    "Will try to reach |M1-M1_target|+|M2-M2_target|+... "
                    "smaller than this value")

    pl.add_argument(
        "--evolve",
        default=False,
        action="store_true",
        help="If specified, start evolution using this initial data.")
    pl.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="If specified, keep all output of the EllipticSolver. Otherwise "
        "erase most files for all but the final Lev of the final Iter.")
    pl.add_argument(
        "--no_iterations",
        default=False,
        action="store_true",
        help="If specified, don't do any elliptic iterations. For testing "
        "only!")

    # Restarting and convergence run options
    pr = p.add_argument_group("Restarting options")
    pr.add_argument(
        "--initial_data_from_this_iteration",
        type=int,
        default=0,
        help="To continue an aborted BBH_ID.py that has successfully "
        "completed iterations up to some number N, set this option to N. "
        "Similarly if you set --postprocess_only, set this option to the N of"
        " the iteration you want to use for postprocess. (default: "
        "%(default)s)")
    pr.add_argument(
        "--convergence_test_at_this_iter",
        type=int,
        default=None,
        help="Taking the parameters from iteration N, create a new dir "
        "ConvergenceTest and run the elliptic solve from Lev0 to Lev5. By "
        "default it is not run. To use the *last* iteration of the elliptic"
        "solver, set N=-1 ")

    pr.add_argument(
        "--postprocess_only",
        action="store_true",
        dest="postprocess_only",
        help="Do not do any elliptic solves. "
        "Instead, use the last successful elliptic solve (specified "
        "via --initial_data_from_this_iteration) and compute the "
        "postprocessed data from it.")

    # Domain adjustment
    pd = p.add_argument_group("Domain options")
    pd.add_argument(
        "--rOutA",
        type=float,
        default=None,
        help="The outer radius of the inner spherical shell around black hole "
        "A")
    pd.add_argument(
        "--rOutB",
        type=float,
        default=None,
        help="The outer radius of the inner spherical shell around black hole "
        "B")

    pd.add_argument(
        "--domainDefFromFile",
        type=str,
        default=None,
        help="A file that contains 'ResizeTheseSubdomains' from an AMR output "
        "to use for the first iteration instead of the defaults")

    # Superposition nuisance parameters
    pS = p.add_argument_group("Superposition options")
    pS.add_argument("--wA",
                    type=float,
                    default=None,
                    help="The width of the Gaussian around black hole A")
    pS.add_argument("--wB",
                    type=float,
                    default=None,
                    help="The width of the Gaussian around black hole B")

    # Parse all the arguments and do sanity checks
    opts = p.parse_args()

    q = opts.q
    D = opts.D

    chiAtarget = array([float(x) for x in opts.chiA.split(",")])
    if len(chiAtarget) != 3:
        error("chiA=%s does not have three components" % opts.chiA)

    chiBtarget = array([float(x) for x in opts.chiB.split(",")])
    if len(chiBtarget) != 3:
        error("chiB=%s does not have three components" % opts.chiB)

    Omega0 = opts.Omega0
    adot = opts.adot

    MaxIts = opts.its
    MaxLevel = opts.levels
    target_residual = opts.target_residual
    NextLevels = 0
    nprocs = opts.nprocs
    debug = opts.debug
    IDType = opts.IDType
    useNegExpBC = opts.useNegExpBC

    if IDType == "CFMS" and useNegExpBC:
        error("Negative expansion BC not implemented for CFMS.")

    if IDType == "SKS" or IDType == "SHK":
        print("# constructing %s initial data for MA/MB=%g" % (IDType, q))
    elif IDType == "CFMS":
        print("# constructing conformally flat initial data for MA/MB=%g" % q)
    else:
        error("Error: Invalid IDType, use 'SKS', 'SHK' or 'CFMS'")

    # Check that if we are continuing an abandoned run, the convergence
    # test will be performed after the *last* iteration
    ConvergenceIt = opts.convergence_test_at_this_iter
    StartingIt = opts.initial_data_from_this_iteration

    if ConvergenceIt is not None and StartingIt > 0:
        error("If you are doing initial data solve from scratch or continuing \
          a solution from an abandoned one, the convergence test *must be \
         run after the last iteration")

    if opts.no_iterations:
        SnesIts0 = "0"
        SnesIts1 = "0"
    else:
        SnesIts0 = "50"
        SnesIts1 = "10"

    # Get several paths.
    ID_Dir = os.getcwd()
    BinDir = ID_Dir + "/private/bin"
    Spells = ID_Dir + "/private/SpEC"

    EvDir = ID_Dir + "/../Ev/"
    SourceBinDir = ID_Dir + "/../bin"
    SourceSpells = ID_Dir + "/../bin/Spells"
    EvInputPath = EvDir

    # Set bool for optional evolution
    EvolveAfterID = opts.evolve

    #Elliptic solver options

    # If we have specified ksp and snes *rtol* then we use relative
    # tolerances, otherwise use absolute!
    elliptic_opts = {}

    if opts.ksp_rtol is not None and opts.snes_rtol is not None:
        elliptic_opts['ksp_rtol'] = opts.ksp_rtol
        elliptic_opts['snes_rtol'] = opts.snes_rtol
    else:
        elliptic_opts['ksp_atol_base'] = opts.ksp_atol_base
        elliptic_opts['snes_atol_base'] = opts.snes_atol_base

    amr_opts = {}
    amr_opts['base_lev'] = opts.base_lev

    domain_opts = {}
    if opts.rOutA and opts.rOutB:
        domain_opts['rOutA'] = float(opts.rOutA)
        domain_opts['rOutB'] = float(opts.rOutB)

    SuperposedKerr_opts = {}
    if opts.wA and opts.wB:
        SuperposedKerr_opts['wA'] = float(opts.wA)
        SuperposedKerr_opts['wB'] = float(opts.wB)

    # Constants
    chiAtargetSq = chiAtarget[0]**2 + chiAtarget[1]**2 + chiAtarget[2]**2
    chiBtargetSq = chiBtarget[0]**2 + chiBtarget[1]**2 + chiBtarget[2]**2
    if IDType == "SHK":
        # Harmonic Kerr radius: r = rBL - M
        rPlusAFull = (q / (1 + q)) * (sqrt(1 - chiAtargetSq))
        rPlusBFull = (1 / (1 + q)) * (sqrt(1 - chiBtargetSq))
    else:
        rPlusAFull = (q / (1 + q)) * (1 + sqrt(1 - chiAtargetSq))
        rPlusBFull = (1 / (1 + q)) * (1 + sqrt(1 - chiBtargetSq))

    # Initialize
    # (python2 note: the line below, when originally written for python2,
    #  was 'if StartingIt > 0 and ConvergenceIt <=0:'.
    #  If ConvergenceIt is None, then ConvergenceIt <=0 is true in python2.
    #  In python3, comparing None vs any other type is an error.  So here
    #  we keep the same behavior in python3 and python2.)
    if StartingIt > 0 and (ConvergenceIt is None or ConvergenceIt <= 0):

        if opts.postprocess_only:
            # In this case, the parameters here are used only for creating
            # the metadata files.  Because these metadata files contain the
            # *INPUT* parameters to each elliptic solve, the correct values
            # here are the ones at the *output* of the  *second to last*
            # elliptic solve. Thus the 'StartingIt-1'.
            RestartIt = StartingIt - 1
        else:
            RestartIt = StartingIt
        initial_guess_dir = os.path.join(os.getcwd(), "private",
                                         "Iter%02i" % StartingIt)
        NextItParams,rAhistory,rBhistory,\
            CurrLevels,NumES,LastPreviousIter,OmegaAhistory, \
            OmegaBhistory,MAhistory,MBhistory,chiAhistory,chiBhistory, \
            errMAhistory,errMBhistory,errchiAhistory,errchiBhistory,outerVini, \
     residualslist = Restart(RestartIt)
        levs = glob("%s/Lev*" % initial_guess_dir)
        highest_lev = levs[-1]
        initial_guess_dir = highest_lev

        AMR_tol = 1e-10

    elif ConvergenceIt is not None and ConvergenceIt > 0:

        bdir = os.getcwd()
        RestartIt = ConvergenceIt
        NextItParams,rAhistory,rBhistory,\
            CurrLevels,NumES,LastPreviousIter = Restart(RestartIt)
        initial_guess_dir = bdir + "/private/Iter%02i" % ConvergenceIt
        levs = glob("%s/Lev*" % initial_guess_dir)
        highest_lev = levs[-1]
        initial_guess_dir = highest_lev
    else:

        NextItParams,rAhistory,rBhistory,CurrLevels,NumES,OmegaAhistory, \
            OmegaBhistory,MAhistory,MBhistory,chiAhistory,chiBhistory, \
            errMAhistory,errMBhistory,errchiAhistory,errchiBhistory \
            = Initialize()

        LastPreviousIter = 0
        initial_guess_dir = ""

    # Check that Spells exec exists
    if not os.path.exists(Spells):
        error("Could not find executable " + Spells)

    # Get 'Revision'
    AOText = SystemOutput("%s -h" % ApplyObsCmd(1))
    FoundMatch = re.search('Code Revision ([^\s]+)\s+', AOText)
    if FoundMatch:
        Revision = FoundMatch.group(1)
    else:
        Revision = "UNKNOWN"

    print("# BH parameters: MA=%g, chiA=(%g, %g, %g); "
          "MB=%g, chiB=(%g, %g, %g)" \
          %(q/(1+q), chiAtarget[0], chiAtarget[1], 
            chiAtarget[2], 1/(1+q), chiBtarget[0], 
            chiBtarget[1], chiBtarget[2]))

    print("# Distance D=%g, Omega0=%g, adot=%g" % (D, Omega0, adot))
    print("# (see Parameters.dat, PhysValues.dat and Errors.dat for"
          " parameters used in")
    print("#  each iteration, as well as physical parameters and errors.)")

    # Dirname and ThisItParams will be overwritten in the iteration.
    # They are set here in case there are no iterations.
    dirname = "Iter%02i" % StartingIt
    ThisItParams = NextItParams

    # First real iteration is 1. Iteration 0 exists to write initial values
    # to file.
    its = StartingIt + 1
    tolerances_met = 0

    if opts.postprocess_only:
        # We need ID_Params and Metadata, without doing another elliptic solve.
        # Read data from last iteration
        rres = {}
        FinalDirr = dirname + "/Lev" + '%(MaxLevel)g' % vars()
        ParseAdmIntegrals(FinalDirr + "/AdmIntegrals.dat", rres)
        rres["Padm"] = genfromtxt(FinalDirr + "/PADM_c.dat")[1:4]
        ParseSpinMass_Ah(FinalDirr + "/Horizons.h5", rres, "AhA")
        ParseSpinMass_Ah(FinalDirr + "/Horizons.h5", rres, "AhB")
        Eadm = rres["Eadm"]
        Padm = rres["Padm"]
        Jadm = rres["Jadm"]
        chiA = rres["AhA-chi"]
        chiB = rres["AhB-chi"]
        MA = rres["AhA-M"]
        MB = rres["AhB-M"]
        MirrA = rres["AhA-Mirr"]
        MirrB = rres["AhB-Mirr"]

    # The intial lev is by definition k=0
    # The initial desired resolution is 10^{-(3+k)}
    k = 0
    # However, we may be restarting via
    # the --initial_data_from_this_iteration option.
    # If this is the case, we want to start from the same lev
    # as the previous directory.
    if StartingIt > 0:
        k = CurrLevels
    if not opts.domainDefFromFile:
        preamble = ""
    else:
        fp = open(opts.domainDefFromFile)
        preamble = "".join(fp.readlines())
        fp.close()

    if StartingIt > 0:
        with open("preamble.dat", "r") as prefi:
            prelines = prefi.readlines()
            index = [preind for preind, preval in enumerate(prelines) \
            if preval==str(StartingIt)+'\n'][0]
            preamble = "".join(prelines[index + 1:index + 14])

    improvements = []

    if StartingIt > 0:
        residuals = residualslist
    elif StartingIt == 0:
        residuals = []
    # The factor to decide when to change resolution.
    # This is \epsilon_{R} in arXiv:1506.01689
    eps = 100

    # Velocity of the outer boundary is initially 0
    outerV = np.zeros(3)

    if StartingIt > 0:
        outerV = outerVini

    # Inverse of the Jacobian
    J_old = None
    if StartingIt > 0 and IDType != "CFMS":
        with h5py.File('J_old.h5', 'r') as jfile:
            J_old = np.array(jfile['its' + str(StartingIt)])

    # Number of times we have tried to succeed on the first
    # iterations
    repeats = 0
    MaxRepeats = 4

    # Begin iterations

    while  (its < MaxIts and not  tolerances_met  \
               and not opts.postprocess_only and  \
                (ConvergenceIt is None or ConvergenceIt <0) ):

        ThisItParams = NextItParams  # use updated parameters

        print("") # some whitespace for clarity
        print("Iter%02i (up to Lev%i) starting at %s"\
              %(its, k, datetime.now().strftime("%H:%M:%S")))
        sys.stdout.flush()

        ObserveIfNotConverged = "yes" if (opts.no_iterations) else "no"

        specID = "q%4.1f,Iter%02i" % (q, its)
        dirname = "Iter%02i" % its
        # The following writes the input files, and then runs spells on them

        res, target_dir, status = DoInitialDataSolve(
            nprocs=nprocs,
            outerV=outerV,
            initial_guess_dir=initial_guess_dir,
            dirname=dirname,
            lev=k,
            D=D,
            Omega0=Omega0,
            Params=ThisItParams,
            chiA=chiAtarget,
            chiB=chiBtarget,
            adot=adot,
            specID=specID,
            debug=debug,
            IDType=IDType,
            q=q,
            preamble=preamble,
            elliptic_opts=elliptic_opts,
            amr_opts=amr_opts,
            domain_opts=domain_opts,
            SuperposedKerr_opts=SuperposedKerr_opts,
            ObserveIfNotConverged=ObserveIfNotConverged,
            rad_factor=opts.rad_factor)

        if status != "Success":
            if its != 1:
                error("ID solve failed at iteration %d. Check the output"
                      " directory %s" % (its, target_dir))
            elif its == 1 and repeats < MaxRepeats:
                its = 1
                preamble = iterateDomain(repeats + 1)
                os.system("mv Iter01 Iter01_Failed_%d" % repeats)
                repeats += 1

                continue
            else:
                error("ID solve failed at first iteration after %d  attempts. "
                      "Check the output directory %s" % (repeats, target_dir))

        # The output of the elliptic solver
        MA = res["AhA-M"]
        MB = res["AhB-M"]
        MirrA = res["AhA-Mirr"]
        MirrB = res["AhB-Mirr"]
        Padm = res["Padm"]
        Eadm = res["Eadm"]
        Mk = res["Mk"]
        Jadm = res["Jadm"]
        chiA = res["AhA-chi"]
        chiB = res["AhB-chi"]
        CoM = res["CoM"]

        ################################################################
        # Errors in quantities
        ################################################################
        errMA = MA - q / (1. + q)
        errMB = MB - 1. / (1. + q)

        errchiA = chiA - chiAtarget
        errchiB = chiB - chiBtarget

        errMAhistory.append(errMA)
        errMBhistory.append(errMB)
        errchiAhistory.append(errchiA)
        errchiBhistory.append(errchiB)
        cxA = ThisItParams["cxA"]
        cyA = ThisItParams["cyA"]
        czA = ThisItParams["czA"]

        # Compute the new velocity of the outer boundary
        outerV = compute_outerV(outerV, Padm, CoM, Omega0, MA, MB, errMA,
                                errMB, D, cxA, cyA, czA)

        physout = open('PhysValues.dat', 'a')
        physout.write("%i %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g "
                      "%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n" %
                      (its, MA, MB, chiA[0], chiA[1], chiA[2], chiB[0],
                       chiB[1], chiB[2], Padm[0], Padm[1], Padm[2], Eadm, Mk,
                       outerV[0], outerV[1], outerV[2]))
        physout.close()
        # Compute total residual - this is the RMS error in quantities
        # of interest.
        # Note that because CoM is convenience only it is *not* included
        residual = sqrt(1/5.*((errMA**2+errMB**2)/(MA+MB)**2 \
            +dot(errchiA,errchiA)+ dot(errchiB,errchiB) \
            +dot(Padm,Padm)/(MA+MB)**2))
        residuals.append(
            [abs(errMA),
             abs(errMB),
             norm(errchiA),
             norm(errchiB),
             norm(Padm)])
        improvement = measure_improvement(residuals)
        improvements.append(improvement)

        # write residuals
        resi = open('Residuals.dat', 'a')
        resi.write(
            "%.15g %.15g %.15g %.15g %.15g\n" %
            (abs(errMA), abs(errMB), norm(errchiA), norm(errchiB), norm(Padm)))
        resi.close()

        improvement_indicator = 0
        if opts.improvement_indicator == "rms":
            improvement_indicator=sqrt(dot(improvement,improvement)/len(improvement))
        elif opts.improvement_indicator == "max":
            improvement_indicator=improvement.max()
        else:
            error("Unknown improvement indicator: %s"%opts.improvement_indicator)

        # Look at AMR output and decide if the domain should be updated
        AMR_error, k_next, preamble, update_params = parseAMRResults(
            target_dir,
            preamble,
            residual,
            k,
            improvement_indicator,
            eps=eps,
            max_lev=MaxLevel)
        if (improvement.max()<1.1 and AMR_error<=-1*MaxLevel-3) \
            or its>MaxIts or (residual<=target_residual
                              and AMR_error<=-1*MaxLevel-3):
            tolerances_met = 1

        # diagnostic output into file and to screen
        errout = open('Errors.dat', 'a')
        errout.write("%i  %i %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n"%
                     (its, k,
                      abs(errMA),  abs(errMB),
                      sqrt(dot(errchiA,errchiA)),
                      sqrt(dot(errchiB,errchiB)),
                      sqrt(dot(Padm,Padm)),
                      10**(AMR_error),
                      residual
                      ))
        errout.close()

        print("Spells completed (%4i seconds)" % ((datetime.now() - StartTime).seconds))
        print("            MA        MB   |             chiA           |" 
              "             chiB           |            P_ADM")
        print("Values: %8.6f  %8.6f | %8.5f %8.5f %8.5f | %8.5f %8.5f %8.5f" 
              " | %8.5f %8.5f %8.5f "%(MA, MB, chiA[0], chiA[1], chiA[2], 
                                       chiB[0], chiB[1], chiB[2], 
                                       Padm[0], Padm[1], Padm[2]))
        print("Errors: %8.1e  %8.1e | %8.1e %8.1e %8.1e | %8.1e %8.1e %8.1e"
              " | %8.1e %8.1e %8.1e"%(errMA, errMB, errchiA[0], errchiA[1],
                                      errchiA[2], errchiB[0], errchiB[1],
                                      errchiB[2], Padm[0],
                                      Padm[1],Padm[2]))
        # decide whether we are done
        if residual<=target_residual and AMR_error<=-1*MaxLevel-3:
            print("residual {:.3e} < target_residual={:.3e}, terminating"
                  .format(residual,target_residual))
            tolerances_met=1
        elif improvement_indicator<1.1 and AMR_error<=-1*MaxLevel-3:
            print("%s improvement small, and AMR is resolved.  Terminating" 
                   % opts.improvement_indicator)
            tolerances_met=1
        elif its>MaxIts:
            print("reached MaxIts={}.  Terminating".format(MaxIts))
            tolerances_met=1
        else:
            print("residual {:.3e} (target_residual={:.1e}), "
                  "improvement_indicator={:.3f} promising, "
                  "AMR_error {:.3f} -- continue iterating"
                  .format(residual, target_residual,
                          improvement_indicator, AMR_error))
        if k==k_next and tolerances_met==0:
            print("keep Lev{}, b/c residual {:.3e}>10^(-1-k)={:.0e} and impr {:.3f}>1.224".format(
                k, residual, 10**(-1-k), improvement_indicator))
        if k_next>k and tolerances_met==0:
            print("Lev{} -> Lev{}, b/c residual {:.3e}<10^(-1-k)={:.0e} or impr {:.3f}<1.224".format(
                k, k_next, residual, 10**(-1-k), improvement_indicator))

        ################################################################
        ## updating formulae
        ################################################################

        if update_params:
            NextItParams, J_old = UpdatingFormulae(
                CoM=CoM,
                OmegaAhistory=OmegaAhistory,
                OmegaBhistory=OmegaBhistory,
                MAhistory=MAhistory,
                MBhistory=MBhistory,
                chiAhistory=chiAhistory,
                chiBhistory=chiBhistory,
                errMAhistory=errMAhistory,
                errMBhistory=errMBhistory,
                errchiAhistory=errchiAhistory,
                errchiBhistory=errchiBhistory,
                mA_target=q / (1. + q),
                mB_target=1. / (1 + q),
                chiA_target=chiAtarget,
                chiB_target=chiBtarget,
                high_spin_update=opts.high_spin_update,
                J_old=J_old,
                currLevels=k,
                **ThisItParams)

        else:
            NextItParams = ThisItParams

# write J_old
        if IDType != "CFMS":
            Jaco = h5py.File('J_old.h5', 'a')
            if 'its' + str(its) in Jaco.keys():
                # Overwrite iterations that are already in the file but
                # after the restart location
                del Jaco['its' + str(its)]
            Jaco.create_dataset('its' + str(its), data=J_old)
            Jaco.close()

# write preamble
        prewrite = open('preamble.dat', 'a')
        prewrite.write("%i\n" % its)
        prewrite.write(preamble)
        prewrite.write("\n")
        prewrite.close()

        initial_guess_dir = target_dir
        its += 1
        k=k_next

# Perform the convergence test if we are supposed to
    AMR_tol = 1e-10
    if ConvergenceIt is not None:
        # Reset everything for convergence test
        k = 0
        preamble = ''

        # Set the residual to 0 so that the next domain
        # is used as long as AMR tolerance is low enough

        if ConvergenceIt == 0:
            initial_guess_dir = ''

        residual = 0.0
        if ConvergenceIt > 0:
            ThisItParams = NextItParams
        print("# Convergence test  starting from Iteration %d"%(ConvergenceIt))
        sys.stdout.flush()

        specID = "q%4.1f,Iter%02i" % (q, its)
        os.mkdir("ConvergenceTest")
        fp = open("AMR_error.dat", "w")
        target_dir2 = initial_guess_dir
        i = 0
        while True:

            dirname = "ConvergenceTest/Iter%02i" % i
            res2, target_dir2, status = DoInitialDataSolve(
                nprocs=nprocs,
                outerV=outerV,
                initial_guess_dir=target_dir2,
                dirname=dirname,
                lev=k,
                D=D,
                Omega0=Omega0,
                Params=ThisItParams,
                chiA=chiAtarget,
                chiB=chiBtarget,
                adot=adot,
                specID=specID,
                debug=debug,
                IDType=IDType,
                q=q,
                preamble=preamble,
                elliptic_opts=elliptic_opts,
                amr_opts=amr_opts,
                domain_opts=domain_opts,
                SuperposedKerr_opts=SuperposedKerr_opts,
                ObserveIfNotConverged=ObserveIfNotConverged,
                rad_factor=opts.rad_factor)
            MA = res2["AhA-M"]
            MB = res2["AhB-M"]
            MirrA = res2["AhA-Mirr"]
            MirrB = res2["AhB-Mirr"]
            Padm = res2["Padm"]
            Eadm = res2["Eadm"]
            Mk = res2["Mk"]
            Jadm = res2["Jadm"]
            chiA = res2["AhA-chi"]
            chiB = res2["AhB-chi"]
            improvement = 3.0
            fp.write("%d %d\t" % (i, k))
            AMR_error, k, preamble, update_params = parseAMRResults(
                target_dir2,
                preamble,
                residual,
                k,
                improvement,
                eps=eps,
                max_lev=MaxLevel)
            fp.write("%.7f\n" % AMR_error)
            if AMR_error <= -1 * MaxLevel - 3:
                initial_guess_dir2 = target_dir2
                break
            i += 1

        print("Done convergence test")
        fp.close()

    # Done iterations
    # Write improvements and residuals to disk
    improvements = array(improvements)
    np.savetxt("Improvements.dat", improvements)

    FinalDir = initial_guess_dir

    print("The final Directory is %s" % FinalDir)
    print("# Finished rootfinding at %s (%i seconds = %4.2f hours)"\
        %(datetime.now().strftime("%H:%M:%S"),
              (datetime.now()-StartTime).seconds,
              (datetime.now()-StartTime).seconds/3600.))
    print("# best Elliptic data is in", FinalDir)

    ################################
    # more sanity checks
    # compute EvID only if those checks worked
    ################################
    if its > MaxIts:
        error("Limit its=%i reached!!" % MaxIts)

################################################################
    if not (os.path.isdir("EllData")):
        if os.path.exists("EllData"):
            error("Trying to copy to EllData, which is not a directory")
        else:
            print("# Copying %s into EllData/" % FinalDir)
            copytree(FinalDir, "EllData")

################################################################

    if not (os.path.isdir("Extr")):
        if os.path.exists("Extr"):
            error("Trying to copy to Extr, which is not a directory")

        print("# Extrapolating into horizons at %s (%i seconds)"\
              %(datetime.now().strftime("%H:%M:%S"),
              (datetime.now()-StartTime).seconds))

        os.mkdir("Extr")
        os.chdir("Extr")

        copy("../EllData/GlobalItems.input", ".")
        copy("../EllData/OrbitalParamItems.input", ".")

        if IDType != "CFMS":
            copy("../EllData/ExtrDomain.input", "./Domain.input")
        else:
            copy("../EllData/Domain.input", "./Domain.input")
            System("sed -e \"s/001.000000000/%s/\" ./Domain.input "
                   "> DomainNew.input" % opts.ExtrFrac)
            System("mv DomainNew.input Domain.input")

        updateDomain("Domain.input", "../EllData/")

        fp = open("../EllData/DoMultipleRuns.input", "r")
        temp = ""
        for line in fp:
            if "ADAPTIVE" not in line:
                temp += line
            else:
                break
        fp.close()

        if IDType != "CFMS":
            temp += """
RunInDirectory("Extrapolation",
         {"ComputeExtr.input" => {
                                  "__AM__" => $AM,
                                  "__BM__" => $BM,
                                  "__ASpinx__" => $ASpinx,
                                  "__ASpiny__" => $ASpiny,
                                  "__ASpinz__" => $ASpinz,
                                  "__Avx__" => $Avx,
                                  "__Avy__" => $Avy,
                                  "__Avz__" => $Avz,
                                  "__Acx__" => $Acx,
                                  "__Acy__" => $Acy,
                                  "__Acz__" => $Acz,

                                  "__BSpinx__" => $BSpinx,
                                  "__BSpiny__" => $BSpiny,
                                  "__BSpinz__" => $BSpinz,
                                  "__Bvx__" => $Bvx,
                                  "__Bvy__" => $Bvy,
                                  "__Bvz__" => $Bvz,
                                  "__Bcx__" => $Bcx,
                                  "__Bcy__" => $Bcy,
                                  "__Bcz__" => $Bcz,
                                  "__OuterVx__" => $OuterVx,
                                  "__OuterVy__" => $OuterVy,
                                  "__OuterVz__" => $OuterVz,
                                   "__AW__" => $AW,
                                   "__BW__" => $BW,
                                   },
         }
  );
            """
        else:
            temp += """
RunInDirectory("Extrapolation",
               {"ComputeExtr.input" => {
                                  "__OuterVx__" => $OuterVx,
                                  "__OuterVy__" => $OuterVy,
                                  "__OuterVz__" => $OuterVz,

             },
        }
 );
            """

        fw = open("DoMultipleRuns.input", "w")
        fw.write(temp)
        fw.close()

        BBH_ID_WriteInputFiles.ComputeExtr(IDType, useNegExpBC)
        System(BinDir + "/DoMultipleRuns -B -e '' -c '%s -v -domaininput "
               "Domain.input  -NoDomainHistory -UseTimes 0 ComputeExtr.input "
               "> ComputeExtr.out 2>&1'" % (ApplyObsCmd(nprocs)))

################################################################
    print("# Computing EvID at %s (%i seconds)"\
        %(datetime.now().strftime("%H:%M:%S"),
          (datetime.now()-StartTime).seconds))

    EvID_Dir = os.path.join(ID_Dir, "EvID")
    System("rm -rf %s" % EvID_Dir)
    os.mkdir(EvID_Dir)
    os.chdir(EvID_Dir)

    BBH_ID_WriteInputFiles.CreateTranslationCpFile(outerV)
    WriteOutputIntoEvID(IDType)

    ################################################################

    if not debug:
        os.chdir(os.path.join(ID_Dir, "private"))
        print("# Cleanup: removing Iter??, and tarring EllData and Extr")
        System("rm -rf Iter??", verbose=True)
        if os.path.isdir("EllData"):
            System("tar -czf EllData.tgz EllData; rm -rf EllData",
                   verbose=True)
        if os.path.isdir("Extr"):
            System("tar -czf Extr.tgz Extr; rm -rf Extr", verbose=True)

################################################################
    print("-------------------------------------------------------------------")
    print("BBH_ID.py finished at %s (%i seconds)"\
          %(datetime.now().strftime("%H:%M:%S"),
          (datetime.now()-StartTime).seconds))

    print("Parameters of final initial data set:")
    print("MA=%12.10f, MB=%12.10f, M1+M2=%f12.10, M1/M2=%10.8f" \
          % (MA, MB, MA + MB, MA / MB))
    print("chiA=", chiA, "  chiB=", chiB)
    print("E_ADM=%f12.10" % Eadm)
    print("P_ADM=", Padm)
    print("J_ADM=", Jadm)

    LJangle = SafeaCos(Jadm[2] / sqrt(dot(Jadm, Jadm)))
    chiAJangle = SafeaCos(
        dot(chiA, Jadm) / sqrt(dot(Jadm, Jadm) * dot(chiA, chiA)))
    chiBJangle = SafeaCos(
        dot(chiB, Jadm) / sqrt(dot(Jadm, Jadm) * dot(chiB, chiB)))
    chiAchiBangle = SafeaCos(
        dot(chiA, chiB) / sqrt(dot(chiA, chiA) * dot(chiB, chiB)))
    print("angle(L,J)=%4.2f, angle(chiA, J)=%4.2f, angle(chiB,J)=%4.2f " \
          "angle(chiA,chiB)=%4.2f"%(LJangle,chiAJangle,
                                    chiBJangle,chiAchiBangle))

    ################################################################
    # Optionally evolve the initial data
    ################################################################
    def call_machines(func):
        return call_perl("Machines", "new Machines()->%s" %(func), \
            want_array=False)

    if EvolveAfterID:
        print("Attempting to submit evolution.")
        ResubCmd = call_machines(
            "GetResubCmd(('WorkDir'=>'%s','SubCmd'=>'./StartJob.sh'))" % EvDir)
        System(ResubCmd, verbose=True)

    print("BBH_ID.py complete.")
