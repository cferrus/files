#!/usr/bin/env python
from __future__ import division
import matplotlib
matplotlib.use('Agg')
import sys, os
sys.path.insert(0,os.path.realpath(__file__+'/../../Python'))
from BbhDiagnosticsImpl import Compute_OrbitalFrequency
from Utils import ReadH5, ReadDat, norm
from Utils import System, error, warning, SystemOutput
from ParseIdParams import ParseIdParams
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, size, sin, cos, array, pi, cross, dot, mean
from numpy.random import rand
import argparse
import re
from scipy import optimize
from scipy.interpolate import InterpolatedUnivariateSpline


"""
Performs eccentricity fitting, and outputs the estimated eccentricity
and updated initial data parameters.
"""

#==============================================================================

# Choose the earliest tmin **after** the junk radiation.
# This function computes the absolute value of the second derivative of
# omega, ODblDot; it then takes a running average of this abs. value
# over avgp points; separately, it averages ODblDot from time 500 until the end
# and calls this yfin; with this, it calculates drop which is the value of
# ODblDot @t=0 - yfin; a threshold of ODblDot, ythr, is the sum of yfin and
# threshp*drop; then tmin is the time at which ODblDot first goes under this
# threshold, ythr, + tshiftp(set @200); however, tmin is capped at 500;
# in other words, tmin=min(tmin, 500). (Author: Robert McGehee)
def FindTmin(t_temp2, dOmegadt, max_tmin):
  # Compute the 2nd deriv. of Omega, take its abs. value
  d2Omegadt = (dOmegadt[2:]-dOmegadt[0:-2])/(t_temp2[2:]-t_temp2[0:-2])
  t=t_temp2[1:-1]
  ODblDot=abs(d2Omegadt)

  # Compute the running average of ODblDot with a given avg param
  avgp=10         #How many points in the running average (has to be an even #)
  y= [0.00]*size(t)
  for i in range(0,avgp//2):
      for j in range(0,avgp+1):
          if j!=i:
              y[i]=y[i]+ODblDot[j]
      y[i]=y[i]/avgp
  for i in range(avgp//2,size(t)-avgp//2):
      for j in range(i-avgp//2,i+avgp//2+1):
          if j!=i:
              y[i]=y[i]+ODblDot[j]
      y[i]=y[i]/avgp
  for i in range(size(t)-avgp//2,size(t)):
      for j in range(size(t)-avgp-1,size(t)):
          if j!=i:
              y[i]=y[i]+ODblDot[j]
      y[i]=y[i]/avgp

  # Compute the average,yfin, of ODblDot from t=500 until the end
  sum=0.0000000
  for i in range(500,size(t)):
      sum=sum+ODblDot[i]
  yfin=sum/(size(t)-500)
  drop=ODblDot[0]-yfin #the drop in magnitude determines a threshold
  threshp=0.001        #parameter which controls how low the threshold ythr is
  ythr=yfin+threshp*drop

  # Loop through y to find first time where y < ythr
  k=0
  while (k<size(t)) and (y[k]>=ythr):
      k=k+1
  if k==size(t):
      k=k-1
  tshiftp=200          #parameter controls how much "safety" is added to tmin
  tmin=t[k]+tshiftp

  # Don't allow tmin to be arbitrarily large
  return min(tmin, max_tmin)

#==============================================================================

def old_fit(x,y,F, inparams):
    # x, y -- arrays of equal length, having x and y values of data
    # F -- function: F(p,x) taking parameters p and set of x-values x
    # inparam -- starting values of parameters
    errfunc=lambda p,x,y: F(p,x)-y
    p,success=optimize.leastsq(errfunc, inparams[:], args=(x,y))

    # compute rms error of fit
    e2=(errfunc(p,x,y))**2
    rms=sqrt(sum(e2)/size(e2))
    return p,rms,success


# fit the data (t,y) to the model F[p,t], by least-squares fitting the params p.
def fit(t, y, F, p0, bounds, jac, name):
    if args.unbounded_fits:
      return old_fit(t, y, F, p0)

    # t,y -- arrays of equal length, having t and y values of data
    # F   -- function: F(p,t) taking parameters p and set of t-values
    # p0  -- starting values of parameters
    errfunc = lambda p,t,y: sqrt(mean((F(p,t)-y)**2))
    jacfunc = lambda p,t,y: array([ 1/errfunc(p,t,y) * mean((F(p,t)-y)*dFdp(p,t)) for dFdp in jac ])
    res = optimize.minimize(
      errfunc,
      method = 'TNC',  #options: SLSQP, L-BFGS-B, TNC
      x0 = p0,
      bounds = bounds,
      jac = jacfunc,
      tol = 1e-15,
      options = {
        'disp': False,
        'maxiter': 3000,
        'eta': 0.5,
      },
      args = (t,y),
    )

    if not res.success:
      if name == "F2cos2_SS":
        error("minimize failed in {}: {}".format(name, res.message))
      else:
        warning("minimize failed in {}: {}".format(name, res.message))

    return res.x, res.fun, res.success

def plot_fitted_function(p, F, x, y, idxFit, idxPlot, idxZoom, name, style):
    # add the plot to all four panels
    # p, F: params and fitting function
    # x, y - complete data
    # idxFit -- indices used in fit, for indicating fit-interval
    # idxPlot -- indices to be plotted in left windows
    # idxZoom -- indices to be plotted in right windows
    # name -- string placed into legend
    # style - plot-style

    xBdry=array([x[idxFit][0], x[idxFit][-1]])
    yBdry=array([y[idxFit][0], y[idxFit][-1]])

    ## Top-left plot: fit
    plt.subplot(2,2,1)
    v=plt.axis()
    plt.plot(x[idxPlot],1e6*F(p,x[idxPlot]),style, label=name)

    # add point at begin and end of fit interval:
    plt.plot(xBdry,1e6*F(p,xBdry),'o')
    plt.axis(v)

    # bottom-left plot: residual
    plt.subplot(2,2,3)
    q1=plt.axis()
    data=1e6*(F(p,x)-y)
    plt.plot(x[idxPlot],data[idxPlot],style,label=name)
    plt.plot(xBdry,1e6*(F(p,xBdry)-yBdry),'o')
    q2=plt.axis()
    miny=min(q1[2], min(data[idxFit]))
    maxy=max(q1[3], max(data[idxFit]))
    plt.axis([q2[0], q2[1], miny, maxy ] )
    plt.title("1e6 residual")


    ## Top-right plot: zoom of fit
    plt.subplot(2,2,2)
    plt.plot(x[idxZoom],1e6*F(p,x[idxZoom]),style, label=name)
    plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.3)
           ,labelspacing=0.25,handletextpad=0.0,fancybox=True
           )

    # add point at begin of fit interval:
    plt.plot([xBdry[0]],[1e6*F(p,xBdry)[0]],'o')

    # bottom-right plot: zoom of residual
    plt.subplot(2,2,4)
    plt.plot(x[idxZoom],1e6*(F(p,x[idxZoom])-y[idxZoom]),style,label=name)
    plt.plot([xBdry[0]],[(1e6*(F(p,xBdry)-yBdry))[0]],'o')
    plt.title("1e6 residual")

    return

################################################################

#TODO: this should be updated now that we bound omega
def CheckPeriastronAdvance(omega,Omega0,name,ecc,ecc_str):
  if omega/Omega0<0.5 or omega/Omega0>1.2:
    ecc = 9.99999
    ecc_str = str(ecc)
    err="OmegaDotEccRemoval: {name}-fit resulted in large (>1.2) or negative (<0.5) periastron advance:\nomega/Omega0={omegafrac}. This is likely wrong (or a bad omega fit from too small ecc oscillation amplitudes),\nso eccentricity was set to {ecc}\n".format(name=name,omegafrac=omega/Omega0,ecc=ecc)
    warning(err)
    summary.write(err)
  return ecc,ecc_str

################################################################

def ComputeUpdate(Omega0, adot0, D0,
                  Tc, B, omega,
                  phi, phi_tmin,
                  name, tmin, tmax, rms,
                  Improved_Omega0_update, no_check, Source):
    delta_adot0=B/(2.*Omega0)*cos(phi)
    delta_Omega0=-B*omega/(4.*Omega0*Omega0)*sin(phi)
    if(Improved_Omega0_update):
        # extra factor Omega0/omega in delta_Omega0
        delta_Omega0=-B/(4.*Omega0)*sin(phi)

    delta_D0=-B*D0*omega*sin(phi)/(2*Omega0*(Omega0**2+2/D0**3))
    ecc=B/(2.*Omega0*omega)
    ecc_str="{:7.7f}".format(ecc)

    if rms/B>0.4:
      summary.write("Large residual of ecc-fit, report bound on ecc\n")
      # for a sine-wave, the rms is 1/2 its amplitude.  Therefore, assume
      # that the amplitude of a bad-fit is 2*rms.  Double that, for safety
      # and because we have to disregard the term omega/Omega0
      ecc=4.*rms/(2.*Omega0*Omega0)
      ecc_str="<{:.1e}".format(ecc)
    elif not no_check:
      ecc,ecc_str = CheckPeriastronAdvance(omega,Omega0,name,ecc,ecc_str)
    pi=np.pi
    summary.write("%s:  %+11.8f    %+11.8f     %+9.6f   %s    %9.6f       %9.6f\n"\
        %(name,delta_Omega0,delta_adot0,delta_D0,ecc_str,
          (phi-pi/2.)%(2.*pi), (phi_tmin-pi/2.)%(2.*pi)  # mean anomaly, in [0, 2pi]
      ))

    f = open("Params_%s.dat"%name, 'w')
    f.write("# EccRemoval.py utilizing orbital frequency Omega, fit %s\n"%name)
    f.write("# Source file=%s\n" % Source)
    f.write("# Omega0=%10.8g, adot0=%10.8g, D0=%10.8g\n"%(Omega0,adot0,D0))
    f.write("# Fitting interval [tmin,tmax]=[%g,%g]\n"%(tmin,tmax))
    f.write("# oscillatory part of fit: (B,omega,phi)=(%g,%g,%g)\n"%(B,omega,phi))
    f.write("# ImprovedOmega0Update=%s\n"%Improved_Omega0_update)
    f.write("# [1] = Omega0\n")
    f.write("# [2] = 1e4 adot0\n")
    f.write("# [3] = D0\n")
    f.write("# [4] = ecc\n")
    f.write("%10.12f\t%10.12f\t%10.12f\t%10.12f\n"%(Omega0, 1e4*adot0, D0, ecc))
    f.write("%10.12f\t%10.12f\t%10.12f\t%10.12f\n"%(Omega0+delta_Omega0, 1e4*(adot0+delta_adot0), D0+delta_D0, 0.))
    f.close()

    f = open("Fit_%s.dat"%name, 'w')
    f.write("# EccRemoval.py utilizing orbital frequency Omega, fit %s\n"%name)
    f.write("# Source file=%s\n" % Source)
    f.write("# Omega0=%10.8g, adot0=%10.8g\n"%(Omega0,adot0))
    f.write("# Fitting interval [tmin,tmax]=[%g,%g]\n"%(tmin,tmax))
    f.write("# oscillatory part of fit: (B,omega,phi)=(%g,%g,%g)\n"%(B,omega,phi))
    f.write("# [1] = Tstart\n")
    f.write("# [2] = Tend\n")
    f.write("# [3] = Tc\n")
    f.write("# [4] = B\n")
    f.write("# [5] = omega\n")
    f.write("# [6] = sin(phi)\n")
    f.write("# [7] = rms residual\n")
    f.write("# [8] = rms residual/B\n")
    f.write("# [9] = omega/Omega0\n")
    f.write("# [10] = ecc\n")

    f.write("%g %g %g %g %g %g %g %g %g %g\n"%(tmin, tmax, Tc, B, omega, sin(phi), rms, rms/B, omega/Omega0, ecc))
    f.close()

################################################################

def ComputeOmega500Update(Omega0, adot0, D0,
                          Tc, B, omega, phi, name, tmin, tmax, rms,
                          deltaOmega0, Source):
    """Compute updating formulae that change Omega by deltaOmega, and adjust
    D accordingly.  This is useful to achieve Omega(t=500)=Omega500,
    in addition to e=0."""
    delta_adot0=B/(2.*Omega0)*cos(phi)
    #delta_Omega0=-B*omega/(4.*Omega0*Omega0)*sin(phi)
    #if(Improved_Omega0_update):
    #    # extra factor Omega0/omega in delta_Omega0
    #    delta_Omega0=-B/(4.*Omega0)*sin(phi)

    delta_D0=-1./6*B*D0/Omega0/omega*sin(phi) - 2./3. *deltaOmega0/Omega0*D0
    ecc=B/(2.*Omega0*omega)
    summary.write("%s:  %7.5f   %+11.8f    %+11.8f     %+9.6f\n"\
        %(name,ecc,delta_adot0,deltaOmega0,delta_D0))

    f = open("Params_%s.dat"%name, 'w')
    f.write("# EccRemoval.py utilizing orbital frequency Omega, fit %s\n"%name)
    f.write("# Source file=%s\n" % Source)
    f.write("# Omega0=%10.8g, adot0=%10.8g, D0=%10.8g\n"%(Omega0,adot0,D0))
    f.write("# Fitting interval [tmin,tmax]=[%g,%g]\n"%(tmin,tmax))
    f.write("# oscillatory part of fit: (B,omega,phi)=(%g,%g,%g)\n"%(B,omega,phi))
    f.write("# NOTE:  This adjusts Omega0 to desired deltaOmega0=%f.\n"%deltaOmega0)
    f.write("# NOTE:  CHANGE ALL THREE OF adot0, Omega0, adot0!\n")
    f.write("# [1] = Omega0\n")
    f.write("# [2] = 1e4 adot0\n")
    f.write("# [3] = D0\n")
    f.write("# [4] = ecc\n")
    f.write("%10.8f\t%10.8f\t%10.8f\t%10.8f\n"%(Omega0, 1e4*adot0, D0, ecc))
    f.write("%10.8f\t%10.8f\t%10.8f\t%10.8f\n"%(Omega0+deltaOmega0, 1e4*(adot0+delta_adot0), D0+delta_D0, 0.))
    f.close()


#==== logic for idperl option
def ParseIDParams(path):
    File = os.path.basename(os.path.realpath(path))
    Dir  = os.path.dirname(os.path.realpath(path))
    D = ParseIdParams(Dir, file=File, array=True)

    Omega0 = float(D['ID_Omega0'][0])
    adot0  = float(D['ID_adot0'][0])
    D0     = float(D['ID_d'][0])

    return Omega0,adot0,D0

def GetTrajectories(Dir, Type, tmax_fit):
  if Type == "bbh":
    Horizons = ReadH5(os.path.join(Dir,"Horizons.h5"))
    TrajA = Horizons["AhA.dir/CoordCenterInertial.dat"]
    TrajB = Horizons["AhB.dir/CoordCenterInertial.dat"]
  elif Type == "bhns":
    Horizons = ReadH5(os.path.join(Dir,"Horizons.h5"))
    TrajA = Horizons["AhA.dir/CoordCenterInertial.dat"]
    TrajB = ReadDat(os.path.join(Dir,"InertialCenterOfMass.dat"))
  elif Type == "bhnsH5":
    Horizons = ReadH5(os.path.join(Dir,"Horizons.h5"))
    Matter = ReadH5(os.path.join(Dir,"Matter.h5"))
    TrajA = Horizons["AhA.dir/CoordCenterInertial.dat"]
    TrajB = Matter["InertialCenterOfMass.dat"]
  elif Type == "nsns":
    TrajA = ReadDat(os.path.join(Dir,"CoM-NSA-InertialFrame.dat"))
    TrajB = ReadDat(os.path.join(Dir,"CoM-NSB-InertialFrame.dat"))
  elif Type == "nsnsH5":
    Matter = ReadH5(os.path.join(Dir,"Matter.h5"))
    TrajA = Matter["CoM-NSA-InertialFrame.dat"]
    TrajB = Matter["CoM-NSB-InertialFrame.dat"]

  # Truncate the trajectories to a region around the fit interval
  tmax = 1.5 * tmax_fit
  tA = TrajA[:,0]
  tB = TrajB[:,0]
  return TrajA[tA<tmax,:], TrajB[tB<tmax,:]



def ComputeOmegaAndDerivsFromFile(TrajA, TrajB):
  '''Compute Omega, dOmega, ellHat and nHat'''

  # While method="Smooth" appears to be more accurate, "Fit"
  # does a slightly better job of reducing eccentricity.
  t_raw, OmegaVec = Compute_OrbitalFrequency(TrajA, TrajB, 10, method='Fit')
  Omega_raw = norm(OmegaVec, axis=1)

  # Compute derivative
  dOmega    = (Omega_raw[2:]-Omega_raw[0:-2])/(t_raw[2:]-t_raw[0:-2])
  Omega     = Omega_raw[1:-1]
  OmegaVec  = OmegaVec[1:-1]
  t = t_raw[1:-1]

  return t,Omega,dOmega,OmegaVec

def interpolateVector(t_old,v,t_new):
  assert len(v[0])==3, \
    "This function expects a list of 3d vectors, got {}d".format(len(v[0]))
  res = np.zeros((len(t_new),3))

  for i in range(3):
    intp =InterpolatedUnivariateSpline(t_old,v[:,i],k=5)
    res[:,i]=intp(t_new)

  return res

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def GetVarsFromSpinData(Dir, XA, XB, OmegaVec, t, tmin):
  '''Compute quantities that depends on the mass and spin'''
  Horizons = ReadH5(os.path.join(Dir, "Horizons.h5"))
  sA = Horizons["AhA.dir/DimensionfulInertialSpin.dat"]
  sB = Horizons["AhB.dir/DimensionfulInertialSpin.dat"]
  MA = Horizons["AhA.dir/ChristodoulouMass.dat"]
  MB = Horizons["AhB.dir/ChristodoulouMass.dat"]

  m_idx = find_nearest(MA[:,0],tmin)
  mA = MA[m_idx,1]
  mB = MB[m_idx,1]

  t_idx = find_nearest(t, tmin)
  Omega_tmin = norm(OmegaVec[t_idx,:])

  alpha_intrp, S_0_perp_n = ComputeSpinTerms(XA,XB,sA,sB,mA,mB,OmegaVec,t)

  nu = mA * mB / (mA + mB) ** 2
  m = mA + mB

  x_i = (m * Omega_tmin) ** (2. / 3)
  T_merge = tmin + 5 / 256.*m ** 2 / nu * x_i ** (-4)
  Amp = (5. / 256.*m ** 2 / nu) ** (3. / 8) * m


  return alpha_intrp, S_0_perp_n, T_merge, Amp


def ComputeSpinTerms(TrajA, TrajB, sA, sB, mA, mB, OmegaVec, t):
    t_s = sA[:,0]

    def PointwiseDot(A,B):
      '''Dot product for each element in two lists of vectors'''
      return np.array([np.dot(a,b) for a,b in zip(A,B)])
    def PointwiseNormalize(A):
      '''Normalization for each element in a list of vectors'''
      return A/norm(A, axis=1)[:,np.newaxis]

    ellHat    = PointwiseNormalize(OmegaVec)
    nHat_raw  = PointwiseNormalize(TrajA[:,1:4] - TrajB[:,1:4])
    nHat      = PointwiseNormalize(interpolateVector(TrajA[:, 0], nHat_raw, t))
    lambdaHat = PointwiseNormalize(cross(ellHat, nHat, axis=1))

    S_0 = (1 + mB / mA) * sA[:,1:4] + (1 + mA / mB) * sB[:,1:4]
    S_0 = interpolateVector(t_s, S_0, t)
    S_0_perp = S_0 - PointwiseDot(S_0, ellHat)[:,np.newaxis] * ellHat
    S_0_perp_n = norm(S_0_perp, axis=1)

    alpha = np.arctan2(PointwiseDot(S_0_perp, lambdaHat), \
                       PointwiseDot(S_0_perp, nHat))
    alpha_intrp = InterpolatedUnivariateSpline(t, alpha, k=5)

    # Sanity check
    got = sin(2 * alpha) * S_0_perp_n ** 2 * 0.5
    expected = PointwiseDot(S_0_perp, nHat) * PointwiseDot(S_0_perp, lambdaHat)
  
    # We *fit* OmegaVec, so it fails to be perpendicular to nHat by about 1e-6
    np.testing.assert_allclose(got,expected,rtol=1e-5,atol=1e-8)

    return alpha_intrp, S_0_perp_n

class FitBounds:
  """Set bounds for some of the variables being fit"""
  def __init__(self, tmax, theOmega0):
    # Time to coalescence must be greater than the fit interval
    self.Tc    = [tmax, None]
    # Frequency should be positive and similar to the initial frequency
    self.omega = [0.6*theOmega0, 1.4*theOmega0]
    # Amplitude of the cos term is chosen to be positive (negative is phi+=pi)
    self.B     = [1e-16, None]
    # phi is a phase, trivially bounded
    self.phi   = [0, 2*pi]
    # For various variables that have no limits
    self.none  = [None, None]


def perform_fits(t,dOmegadt,idxFit,tFit,dOmegadtFit,tmin,tmax,theOmega0,theD0,theadot0,opts,Source):
  '''Given the data, perform all the fits and output the data'''

  # Set bounds for some of the variables
  lim = FitBounds(tmax, theOmega0)

  # fit a0(T-t)^(-11./8)
  Tmerger = t[-1]+1000

  F1 = lambda p,t: p[1]*(p[0]-t)**(-11/8)
  p0 = [Tmerger, 1e-5]
  pBounds = [lim.Tc, lim.none]
  jac = [ lambda p,t: (-11/8)*p[1]*(p[0]-t)**(-19/8),
          lambda p,t:              (p[0]-t)**(-11/8) ]
  pF1, rmsF1, F1_status = fit(tFit, dOmegadtFit, F1, p0, pBounds, jac, "F1")

  #================

  F1cos1 = lambda p,t: p[1]*(p[0]-t)**(-11/8) + p[2]*cos(p[3]*t+p[4])
  jac = [ lambda p,t: (-11/8)*p[1]*(p[0]-t)**(-19/8),
          lambda p,t:              (p[0]-t)**(-11/8),
          lambda p,t:       cos(p[3]*t + p[4]),
          lambda p,t: -p[2]*sin(p[3]*t + p[4])*t,
          lambda p,t: -p[2]*sin(p[3]*t + p[4]) ]

  rmsF1cos1 = 2*rmsF1
  pF1cos1 = pF1
  # try a few initial guesses for the phase to ensure good convergence
  for phi in range(0,6):
    p0 = [pF1[0], pF1[1], rmsF1, 0.8*theOmega0, phi]
    pBounds = [lim.Tc, lim.none, lim.B, lim.omega, lim.phi]
    ptemp, rmstemp, F1cos1_status = fit(tFit, dOmegadtFit, F1cos1,
                                        p0, pBounds, jac, "F1cos1")
    if(rmstemp<rmsF1cos1):
        rmsF1cos1=rmstemp
        pF1cos1=ptemp

  #================

  F1cos2 = lambda p,t: p[1]*(p[0]-t)**(-11/8) + p[2]*cos(p[3]*t+p[4]+p[5]*t*t)
  jac = [ lambda p,t: (-11/8)*p[1]*(p[0]-t)**(-19/8),
          lambda p,t:              (p[0]-t)**(-11/8),
          lambda p,t:       cos(p[3]*t + p[4] + p[5]*t*t),
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t)*t,
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t),
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t)*t*t ]
  p0 = [pF1cos1[0], pF1cos1[1], pF1cos1[2], pF1cos1[3], pF1cos1[4], 0]
  pBounds = [lim.Tc, lim.none, lim.B, lim.omega, lim.phi, lim.none]
  pF1cos2, rmsF1cos2, F1cos2_status = fit(tFit, dOmegadtFit, F1cos2,
                                          p0, pBounds, jac, "F1cos2")

  #================

  F2cos1 = lambda p,t: p[1]*(p[0]-t)**(-11/8) + p[2]*(p[0]-t)**(-13/8) + p[3]*cos(p[4]*t+p[5])
  jac = [ lambda p,t: (-11/8)*p[1]*(p[0]-t)**(-19/8) + (-13/8)*p[2]*(p[0]-t)**(-21/8),
          lambda p,t: (p[0]-t)**(-11/8),
          lambda p,t: (p[0]-t)**(-13/8),
          lambda p,t:       cos(p[4]*t + p[5]),
          lambda p,t: -p[3]*sin(p[4]*t + p[5])*t,
          lambda p,t: -p[3]*sin(p[4]*t + p[5]) ]
  p0 = [pF1cos1[0], pF1cos1[1], 0., pF1cos1[2], pF1cos1[3], pF1cos1[4]]
  pBounds = [lim.Tc, lim.none, lim.none, lim.B, lim.omega, lim.phi]
  pF2cos1, rmsF2cos1, F2cos1_status = fit(tFit, dOmegadtFit, F2cos1,
                                          p0, pBounds, jac, "F2cos1")

  #================

  # F2cos2 = a0(Tc-t)^(-11/8)+a1(Tc-t)^(-13/8)+Bcos(omega t+phi+b t^2)
  F2cos2 = lambda p,t: p[1]*(p[0]-t)**(-11/8) + p[2]*(p[0]-t)**(-13/8) + p[3]*cos(p[4]*t+p[5]+p[6]*t*t)
  jac = [ lambda p,t: (-11/8)*p[1]*(p[0]-t)**(-19/8) + (-13/8)*p[2]*(p[0]-t)**(-21/8),
          lambda p,t: (p[0]-t)**(-11/8),
          lambda p,t: (p[0]-t)**(-13/8),
          lambda p,t:           cos(p[4]*t + p[5] + p[6]*t*t),
          lambda p,t: -p[3]*sin(p[4]*t + p[5] + p[6]*t*t)*t,
          lambda p,t: -p[3]*sin(p[4]*t + p[5] + p[6]*t*t),
          lambda p,t: -p[3]*sin(p[4]*t + p[5] + p[6]*t*t)*t*t ]
  p0 = [pF1cos2[0], pF1cos2[1], 0., pF1cos2[2], pF1cos2[3], pF1cos2[4], pF1cos2[5]]
  pBounds = [lim.Tc, lim.none, lim.none, lim.B, lim.omega, lim.phi, lim.none]
  pF2cos2, rmsF2cos2, F2cos2_status = fit(tFit, dOmegadtFit, F2cos2,
                                          p0, pBounds, jac, "F2cos2")

  fp = open("FitSuccess.txt","w")
  fp.write("F1  %d\n" % F1_status)
  fp.write("F1cos1  %d\n" % F1cos1_status)
  fp.write("F1cos2  %d\n" % F1cos2_status)
  fp.write("F2cos1  %d\n" % F2cos1_status)
  fp.write("F2cos2  %d\n" % F2cos2_status)
  fp.close()

  #==== output all four residuals
  summary.write("tmin=%f  (determined by "%tmin)
  if(opts.tmin==None):
    summary.write(" fit)\n")
  else:
    summary.write(" option)\n")

  summary.write("""
RESIDUALS
    F1cos1 rms=%g   \tF1cos2 rms=%g
    F2cos1 rms=%g   \tF2cos2 rms=%g

DIAGNOSTICS
                Tc         B     sin(phi)    rms/B   omega/Omega0
""" %(rmsF1cos1,rmsF1cos2,rmsF2cos1,rmsF2cos2,))

  # add first two fields for remaining fits
  tmp="%-8s %10.1f   %6.2e    %6.4f     %5.3f      %5.3f\n"
  summary.write(tmp%("F1cos1", pF1cos1[0], pF1cos1[2], sin(pF1cos1[4]),
                               rmsF1cos1/pF1cos1[2], pF1cos1[3]/theOmega0))
  summary.write(tmp%("F1cos2", pF1cos2[0], pF1cos2[2], sin(pF1cos2[4]),
                               rmsF1cos2/pF1cos2[2], pF1cos2[3]/theOmega0))
  summary.write(tmp%("F2cos1", pF2cos1[0], pF2cos1[3], sin(pF2cos1[5]),
                               rmsF2cos1/pF2cos1[3], pF2cos1[4]/theOmega0))
  summary.write(tmp%("F2cos2", pF2cos2[0], pF2cos2[3], sin(pF2cos2[5]),
                               rmsF2cos2/pF2cos2[3], pF2cos2[4]/theOmega0))

  #==== compute updates and generate EccRemoval_FIT.dat files
  summary.write("""
ECCENTRICITY AND UPDATES
         delta_Omega0   delta_adot0     delta_D0     ecc      mean_anomaly(0)  mean_anomaly(tmin)
""")

  ComputeUpdate(theOmega0, theadot0, theD0,
                pF1cos1[0], pF1cos1[2], pF1cos1[3],
                pF1cos1[4], pF1cos1[3]*tmin+pF1cos1[4],  # arg of cos(..) at t=0 and at tmin
                "F1cos1",tmin,tmax,rmsF1cos1,
                opts.improved_Omega0_update, opts.no_check, Source
                )

  ComputeUpdate(theOmega0, theadot0, theD0,
                pF1cos2[0], pF1cos2[2], pF1cos2[3],
                pF1cos2[4], pF1cos2[4]+pF1cos2[3]*tmin+pF1cos2[5]*tmin*tmin,
                "F1cos2",tmin,tmax,rmsF1cos2,
                opts.improved_Omega0_update, opts.no_check, Source
                )

  ComputeUpdate(theOmega0, theadot0, theD0,
                pF2cos1[0], pF2cos1[3], pF2cos1[4],
                pF2cos1[5], pF2cos1[5]+pF2cos1[4]*tmin,
                "F2cos1",tmin,tmax,rmsF2cos1,
                opts.improved_Omega0_update, opts.no_check, Source
                )

  ComputeUpdate(theOmega0, theadot0, theD0,
                pF2cos2[0], pF2cos2[3], pF2cos2[4],
                pF2cos2[5], pF2cos2[5]+pF2cos2[4]*tmin+pF2cos2[6]*tmin*tmin,
                "F2cos2",tmin,tmax,rmsF2cos2,
                opts.improved_Omega0_update, opts.no_check, Source
                )

  if(opts.Omega500!=None):
    Omega500=opts.Omega500

    # perform fit to Omega itself, so we can compue Omega(t=500)
    # fit a0(T-t)^(-11./8)
    Tmerger=t[-1]+1000

    Fomg = lambda p,t: p[1]*(p[0]-t)**(-3./8.)
    jac = [ lambda p,t: (-3/8)*p[1]*(p[0]-t)**(-11/8),
            lambda p,t:             (p[0]-t)**(-3/8) ]
    p0 = [Tmerger, Omega[0]*(Tmerger-t[0])**(3/8)]
    pBounds = [lim.Tc, lim.none]
    pFomg,rmsFomg=fit(tFit, OmegaFit, Fomg, p0, pBounds, jac, "Omega500")

    Omega500fit=Fomg(pFomg,500)

    summary.write("""
OMEGA t=500 ADJUSTMENT
  fit Omega(t) = %f*(%f-t)^(-3/8)
  from fit: Omega(t=500)=%f
"""%(pFomg[1], pFomg[0], Omega500fit)
    )

    ComputeOmega500Update(theOmega0, theadot0, theD0,
                          pF2cos2[0], pF2cos2[3], pF2cos2[4], pF2cos2[5],
                          "Omega500_F2cos2",tmin,tmax,rmsF2cos2,
                          Omega500-Omega500fit, Source
                          )



  if not opts.no_plot:
    #=== set plotting intervals
    idxPlot= (t> tmin - 0.2*(tmax-tmin)) & (t< tmax + 0.35*(tmax-tmin) )
    idxZoom=(t< tmin + 0.2*(tmax-tmin))
    tPlot=t[idxPlot]
    dOmegadtPlot=dOmegadt[idxPlot]


    plt.figtext(0.5,0.95,"%s [%g, %g]"%(os.path.basename(Source), tmin,tmax),color='b',size='large',ha='center')
    #=== Top left plot dOmega/dt ===
    plt.subplot(2,2,1)
    plt.plot(tPlot[tPlot>=tmin], 1e6*dOmegadtPlot[tPlot>=tmin], 'k', label="dOmega/dt", linewidth=2)
    QQ1=plt.axis()
    plt.cla()
    plt.plot(tPlot, 1e6*dOmegadtPlot, 'k', label="dOmega/dt", linewidth=2)
    QQ2=plt.axis()
    plt.axis([QQ2[0], QQ2[1], QQ1[2], QQ1[3] ])
    plt.title("1e6 dOmega/dt")


    #==== bottom left ====
    # set x-axes to top-left scale, and y-axes to something small
    # as initial conditions for adding line by line above
    plt.subplot(2,2,3)
    plt.axis([QQ1[0], QQ1[1], -1e-10, 1e-10])

    #==== Top right plot -- zoom of dOmega/dt ====
    plt.subplot(2,2,2)
    plt.plot(t[idxZoom],1e6*dOmegadt[idxZoom], 'k',linewidth=2) #  label="dOmega/dt",
    plt.title("1e6 dOmega/dt           .")  # extra space to avoid legend

    # Plot the individual fits
    plot_fitted_function(pF1cos1,F1cos1,t,dOmegadt,idxFit,idxPlot,idxZoom,
                         'F1cos1', ':')
    plot_fitted_function(pF1cos2,F1cos2,t,dOmegadt,idxFit,idxPlot,idxZoom,
                         'F1cos2', ':')
    plot_fitted_function(pF2cos1,F2cos1,t,dOmegadt,idxFit,idxPlot,idxZoom,
                         'F2cos1', '--')
    plot_fitted_function(pF2cos2,F2cos2,t,dOmegadt,idxFit,idxPlot,idxZoom,
                         'F2cos2', '--')

    # zoom out of the y-axis in the lower left panel by 15%
    plt.subplot(2,2,3)
    q1=plt.axis()
    Deltay=q1[3]-q1[2]
    plt.axis([q1[0], q1[1], q1[2]-0.15*Deltay, q1[3]+0.15*Deltay ] )

    # adjust margins of figure
    plt.subplots_adjust(left=0.09,right=0.95,bottom=0.07,hspace=0.25)

    plt.savefig("FigureEccRemoval.pdf")


def perform_fits_SS(t,dOmegadt,idxFit,tFit,OmegaFit,dOmegadtFit,T_merge, Amp,alpha_intrp,S_0_perp_n,
                    tmin,tmax,theOmega0,theD0,theadot0,opts,Source):
  '''Same as perform_fits but uses functions including spin-spin   terms'''

  # Set bounds for some of the variables
  lim = FitBounds(tmax, theOmega0)

  # Fit *Omega*
  OmegaFunc = lambda p,t: p[1]*(p[0]-t)**(-3/8)   # power law form for Omega
  jac = [ lambda p,t: (-3/8)*p[1]*(p[0]-t)**(-11/8),
          lambda p,t:             (p[0]-t)**(-3/8) ]
  p0 = [T_merge, Amp]
  pBounds = [lim.Tc, lim.none]
  pOmega, rmsOmega, Omega_status = fit(tFit, OmegaFit, OmegaFunc,
                                       p0, pBounds, jac, "OmegaFunc")
  FitTc = pOmega[0]  # Fit time to coalescence

  #================

  F1_SS = lambda p,t: p[0]*(FitTc-t)**(-11/8)
  jac = [ lambda p,t: (FitTc-t)**(-11/8) ]
  p0 = [1e-5]
  pBounds = [lim.none]
  pF1_SS, rmsF1_SS, F1_SS_status = fit(tFit, dOmegadtFit, F1_SS,
                                       p0, pBounds, jac, "F1_SS")

  #================

  # This function is *only* used to provide an initial guess for F1cos1_SS below
  # We use it because generally the amplitude and time to merger one obtains with
  # the "usual" F1cos1 is highly degenerate and can give strange values.
  F1cos1_helper = lambda p,t: p[0]*(FitTc-t)**(-11/8) + p[1]*cos(p[2]*t + p[3])
  jac = [ lambda p,t: (FitTc-t)**(-11/8),
          lambda p,t:       cos(p[2]*t + p[3]),
          lambda p,t: -p[1]*sin(p[2]*t + p[3])*t,
          lambda p,t: -p[1]*sin(p[2]*t + p[3]) ]

  rmsF1cos1_helper = 2*rmsF1_SS
  # try a few initial guesses for the phase to ensure good convergence
  for phi in range(0, 6):
    p0 = [pF1_SS[0], rmsF1_SS, 0.8*theOmega0, phi]
    pBounds = [lim.none, lim.B, lim.omega, lim.phi]
    ptemp, rmstemp,status = fit(tFit, dOmegadtFit, F1cos1_helper,
                                p0, pBounds, jac, "F1cos1_helper")

    if(rmstemp < rmsF1cos1_helper):
      rmsF1cos1_helper = rmstemp
      pF1cos1_helper = ptemp

  #================

  rmsF1cos1_SS = 3*rmsF1_SS

  F1cos1_SS = lambda p,t : p[0]*(FitTc-t)**(-11/8) + p[1]*cos(p[2]*t + p[3]) - p[4]*sin(2*alpha_intrp(t) + p[5])
  jac = [ lambda p,t: (FitTc-t)**(-11/8),
          lambda p,t:       cos(p[2]*t + p[3]),
          lambda p,t: -p[1]*sin(p[2]*t + p[3])*t,
          lambda p,t: -p[1]*sin(p[2]*t + p[3]),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[5]),
          lambda p,t: -p[4]*cos(2*alpha_intrp(t) + p[5]) ]
  for phi in range(0, 6):
    p0 = [pF1cos1_helper[0], pF1cos1_helper[1],
          pF1cos1_helper[2], pF1cos1_helper[3],
          0.5*S_0_perp_n[0]**2 * (theOmega0/theD0)**2, phi]
    pBounds = [lim.none, lim.B, lim.omega, lim.phi, lim.none, lim.phi]
    ptemp, rmstemp, F1cos1_SS_status = fit(tFit, dOmegadtFit, F1cos1_SS,
                                           p0, pBounds, jac, "F1cos1_SS")

    if(rmstemp <= rmsF1cos1_SS):
      rmsF1cos1_SS = rmstemp
      pF1cos1_SS = ptemp

  #================

  F1cos2_SS = lambda p,t: p[0]*(FitTc-t)**(-11/8) + p[1]*cos(p[2]*t + p[3] + p[4]*t*t) - p[5]*sin(2*alpha_intrp(t) + p[6])
  jac = [ lambda p,t: (FitTc-t)**(-11/8),
          lambda p,t:       cos(p[2]*t + p[3] + p[4]*t*t),
          lambda p,t: -p[1]*sin(p[2]*t + p[3] + p[4]*t*t)*t,
          lambda p,t: -p[1]*sin(p[2]*t + p[3] + p[4]*t*t),
          lambda p,t: -p[1]*sin(p[2]*t + p[3] + p[4]*t*t)*t*t,
          lambda p,t:      -sin(2*alpha_intrp(t) + p[6]),
          lambda p,t: -p[5]*cos(2*alpha_intrp(t) + p[6]) ]
  p0 = [pF1cos1_SS[0], pF1cos1_SS[1], pF1cos1_SS[2], pF1cos1_SS[3],
        0, pF1cos1_SS[4],pF1cos1_SS[5]]
  pBounds = [lim.none, lim.B, lim.omega, lim.phi, lim.none, lim.none, lim.phi]
  pF1cos2_SS, rmsF1cos2_SS, F1cos2_SS_status = fit(tFit, dOmegadtFit, F1cos2_SS,
                                                   p0, pBounds, jac, "F1cos2_SS")

  #================

  F2cos1_SS = lambda p,t: p[0]*(FitTc-t)**(-11/8) + p[1]*(FitTc-t)**(-13/8) + p[2]*cos(p[3]*t + p[4]) - p[5]*sin(2*alpha_intrp(t) + p[6])

  jac = [ lambda p,t: (FitTc-t)**(-11/8),
          lambda p,t: (FitTc-t)**(-13/8),
          lambda p,t:       cos(p[3]*t + p[4]),
          lambda p,t: -p[2]*sin(p[3]*t + p[4])*t,
          lambda p,t: -p[2]*sin(p[3]*t + p[4]),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[6]),
          lambda p,t: -p[5]*cos(2*alpha_intrp(t) + p[6]) ]
  p0 = [pF1cos1_SS[0]*1.01, pF1cos1_SS[0]/100, pF1cos1_SS[1]*1.01, pF1cos1_SS[2], pF1cos1_SS[3], pF1cos1_SS[4], pF1cos1_SS[5]]
  pBounds = [lim.none, lim.none, lim.B, lim.omega, lim.phi, lim.none, lim.phi]
  pF2cos1_SS, rmsF2cos1_SS, F2cos1_SS_status = fit(tFit, dOmegadtFit, F2cos1_SS,
                                                   p0, pBounds, jac, "F2cos1_SS")

  #================

  F2cos2_SS = lambda p,t: p[0]*(FitTc-t)**(-11/8) + p[1]*(FitTc-t)**(-13/8) + p[2]*cos(p[3]*t + p[4] + p[5]*t*t) - p[6]*sin(2*alpha_intrp(t) + p[7])

  jac = [ lambda p,t: (FitTc-t)**(-11/8),
          lambda p,t: (FitTc-t)**(-13/8),
          lambda p,t:       cos(p[3]*t + p[4] + p[5]*t*t),
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t)*t,
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t),
          lambda p,t: -p[2]*sin(p[3]*t + p[4] + p[5]*t*t)*t*t,
          lambda p,t:      -sin(2*alpha_intrp(t) + p[7]),
          lambda p,t: -p[6]*cos(2*alpha_intrp(t) + p[7]) ]
  p0 = [pF1cos2_SS[0], pF1cos2_SS[0]/100, pF1cos2_SS[1], pF1cos2_SS[2], pF1cos2_SS[3],
        pF1cos2_SS[4], pF1cos2_SS[5], pF1cos2_SS[6]]
  pBounds = [lim.none, lim.none, lim.B, lim.omega, lim.phi, lim.none, lim.none, lim.phi]

  pF2cos2_SS, rmsF2cos2_SS, F2cos2_SS_status = fit(tFit, dOmegadtFit, F2cos2_SS,
                                                   p0, pBounds, jac, "F2cos2_SS")


  fp = open("FitSuccess_SS.txt","w")
  fp.write("FOmega  %d\n" % Omega_status)
  fp.write("F1  %d\n" % F1_SS_status)
  fp.write("F1cos1_SS  %d\n" % F1cos1_SS_status)
  fp.write("F1cos2_SS  %d\n" % F1cos2_SS_status)
  fp.write("F2cos1_SS  %d\n" % F2cos1_SS_status)
  fp.write("F2cos2_SS  %d\n" % F2cos2_SS_status)
  fp.close()

  #==== output all four residuals
  summary.write("tmin=%f  (determined by " % tmin)
  if(opts.tmin == None):
    summary.write(" fit)\n")
  else:
    summary.write(" option)\n")

  summary.write("""
RESIDUALS
    F1cos1 rms=%g   \tF1cos2 rms=%g
    F2cos1 rms=%g   \tF2cos2 rms=%g

DIAGNOSTICS
                Tc         B     sin(phi)    rms/B   omega/Omega0
""" % (rmsF1cos1_SS, rmsF1cos2_SS, rmsF2cos1_SS, rmsF2cos2_SS))

  # add first two fields for remaining fits
  tmp = "%-8s %10.1f   %6.2e    %6.4f     %5.3f      %5.3f\n"

  summary.write(tmp % ("F1cos1_SS", FitTc, pF1cos1_SS[1], sin(pF1cos1_SS[3]),
                       rmsF1cos1_SS / pF1cos1_SS[1], pF1cos1_SS[2] / theOmega0))
  summary.write(tmp % ("F1cos2_SS", FitTc, pF1cos2_SS[1], sin(pF1cos2_SS[3]),
                       rmsF1cos2_SS / pF1cos2_SS[1], pF1cos2_SS[2] / theOmega0))
  summary.write(tmp%("F2cos1_SS", FitTc, pF2cos1_SS[2], sin(pF2cos1_SS[4]),
                     rmsF2cos1_SS/pF2cos1_SS[2], pF2cos1_SS[3]/theOmega0))
  summary.write(tmp%("F2cos2_SS", FitTc, pF2cos2_SS[2], sin(pF2cos2_SS[4]),
                     rmsF2cos2_SS/pF2cos2_SS[2], pF2cos2_SS[3]/theOmega0))

  #==== compute updates and generate EccRemoval_FIT.dat files
  summary.write("""
ECCENTRICITY AND UPDATES
             delta_Omega0   delta_adot0     delta_D0     ecc     mean_anomaly(0)  mean_anomaly(tmin)
""")

  ComputeUpdate(theOmega0, theadot0, theD0,
              FitTc, pF1cos1_SS[1], pF1cos1_SS[2],
              pF1cos1_SS[3], pF1cos1_SS[2]*tmin+pF1cos1_SS[3],  # arg of cos(..) at t=0 and at tmin
              "F1cos1_SS", tmin, tmax, rmsF1cos1_SS,
              opts.improved_Omega0_update, opts.no_check, Source
              )

  ComputeUpdate(theOmega0, theadot0, theD0,
                FitTc, pF1cos2_SS[1], pF1cos2_SS[2],
                pF1cos2_SS[3], pF1cos2_SS[3]+pF1cos2_SS[2]*tmin+pF1cos2_SS[4]*tmin*tmin,
                "F1cos2_SS", tmin, tmax, rmsF1cos2_SS,
                opts.improved_Omega0_update, opts.no_check, Source
              )

  ComputeUpdate(theOmega0, theadot0, theD0,
              FitTc, pF2cos1_SS[2], pF2cos1_SS[3],
                pF2cos1_SS[4], pF2cos1_SS[4]+pF2cos1_SS[3]*tmin,
              "F2cos1_SS", tmin, tmax, rmsF2cos1_SS,
              opts.improved_Omega0_update, opts.no_check, Source
              )


  ComputeUpdate(theOmega0, theadot0, theD0,
              FitTc, pF2cos2_SS[2], pF2cos2_SS[3],
                pF2cos2_SS[4], pF2cos2_SS[4]+pF2cos2_SS[3]*tmin+pF2cos2_SS[5]*tmin*tmin,
              "F2cos2_SS", tmin, tmax, rmsF2cos2_SS,
              opts.improved_Omega0_update, opts.no_check, Source
              )


  if not opts.no_plot:
     #=== set plotting intervals
    idxPlot= (t> tmin - 0.2*(tmax-tmin)) & (t< tmax + 0.35*(tmax-tmin) )
    idxZoom=(t< tmin + 0.2*(tmax-tmin))
    tPlot=t[idxPlot]
    dOmegadtPlot=dOmegadt[idxPlot]
    plt.clf()
    plt.figtext(0.5, 0.95, "%s [%g, %g]" % (os.path.basename(Source), tmin, tmax), color='b', size='large', ha='center')
    #=== Top left plot dOmega/dt === 
    plt.subplot(2, 2, 1)
    plt.plot(tPlot[tPlot >= tmin], 1e6 * dOmegadtPlot[tPlot >= tmin], 'k', label="dOmega/dt", linewidth=2)
    QQ1 = plt.axis()
    plt.cla()
    plt.plot(tPlot, 1e6 * dOmegadtPlot, 'k', label="dOmega/dt", linewidth=2)
    QQ2 = plt.axis()
    plt.axis([QQ2[0], QQ2[1], QQ1[2], QQ1[3] ])
    plt.title("1e6 dOmega/dt")

    #==== bottom left ====
    # set x-axes to top-left scale, and y-axes to something small
    # as initial conditions for adding line by line above
    plt.subplot(2, 2, 3)
    plt.axis([QQ1[0], QQ1[1], -1e-10, 1e-10])

    #==== Top right plot -- zoom of dOmega/dt ====
    plt.subplot(2, 2, 2)
    plt.plot(t[idxZoom], 1e6 * dOmegadt[idxZoom], 'k', linewidth=2)  #  label="dOmega/dt", 
    plt.title("1e6 dOmega/dt           .")  # extra space to avoid legend

    # Plot individual fits
    plot_fitted_function(pF1cos1_SS, F1cos1_SS, t, dOmegadt, idxFit, idxPlot, idxZoom,
                         'F1cos1_SS', ':')
    plot_fitted_function(pF1cos2_SS, F1cos2_SS, t, dOmegadt, idxFit, idxPlot, idxZoom,
                         'F1cos2_SS', ':')
    plot_fitted_function(pF2cos1_SS, F2cos1_SS, t, dOmegadt, idxFit, idxPlot, idxZoom,
                         'F2cos1_SS', '--')
    plot_fitted_function(pF2cos2_SS, F2cos2_SS, t, dOmegadt, idxFit, idxPlot, idxZoom,
                         'F2cos2_SS', '--')

    # zoom out of the y-axis in the lower left panel by 15%
    plt.subplot(2, 2, 3)
    q1 = plt.axis()
    Deltay = q1[3] - q1[2]
    plt.axis([q1[0], q1[1], q1[2] - 0.15 * Deltay, q1[3] + 0.15 * Deltay ])

    # adjust margins of figure
    plt.subplots_adjust(left=0.09, right=0.95, bottom=0.07, hspace=0.25)

    plt.savefig("FigureEccRemoval_SS.pdf")



#################################################################

def main(opts):
  
  Source = opts.d
  theOmega0,theadot0,theD0 = ParseIDParams(opts.idperl)

  # Determine the latest possible fit time so we can truncate data
  max_tmin_fit = 500
  max_tmax_fit = max_tmin_fit + 5*pi/theOmega0

  XA,XB = GetTrajectories(Source, opts.t, max_tmax_fit)
  t,Omega,dOmegadt,OmegaVec = ComputeOmegaAndDerivsFromFile(XA,XB)

  #==== Determine bounds of fit
  if opts.tmin is None:
    tmin = FindTmin(t, dOmegadt, max_tmin_fit)
  else:
    tmin = opts.tmin

  if(opts.tmax==None):
    tmax=min(t[-1], tmin+5*pi/theOmega0)
  else:
    tmax=opts.tmax

  idxFit=  (t>=tmin) & (t<=tmax)
  tFit=t[idxFit]
  OmegaFit=Omega[idxFit]
  dOmegadtFit=dOmegadt[idxFit]
  
  # First, do the usual fits
  perform_fits(t,dOmegadt,idxFit,tFit,dOmegadtFit,tmin,tmax,theOmega0,
               theD0,theadot0,opts,Source)
  
  # Now, do fits that include spin-spin interactions,
  # see Buonnano et al., 2010 (arXiv 1012.1549v2).

  # Everything below is only for BBH (since NS have negligible spins)
  if opts.t=="bbh":
    alpha_intrp, S_0_perp_n, T_merge, Amp = GetVarsFromSpinData(opts.d, XA, XB,
                                                                OmegaVec, t, tmin)
    
    perform_fits_SS(t,dOmegadt,idxFit,tFit,OmegaFit,dOmegadtFit,T_merge,
                    Amp,alpha_intrp,S_0_perp_n,tmin,tmax,theOmega0,
                    theD0,theadot0,opts,Source)

#=============================================================================
if __name__ == "__main__":

  p = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter
  )

  p1 = p.add_argument_group("required arguments")
  p1.add_argument("--idperl", type=str, required=True, metavar="FILE",
    help="File with initial data variables ID_{Omega0,adot0,D0},\ne.g. ID_Params.perl")

  p1.add_argument("-d", type=str, required=True, metavar="DIR",
    help="Directory containing the evolution trajectory file(s)")
  p1.add_argument("-t", type=str, required=True, choices=['bbh','bhns','bhnsH5','nsns','nsnsH5'],
    help="""Type of binary determines which files are expected in DIR.
bbh:
  Horizons.h5/Ah{A,B}.dir/CoordCenterInertial.dat
  Horizons.h5/Ah{A,B}.dir/DimensionfulInertialSpin.dat
  Horizons.h5/Ah{A,B}.dir/ChristodoulouMass.dat
bhns:
  Horizons.h5/AhA.dir/CoordCenterInertial.dat
  InertialCenterOfMass.dat
bhnsH5:
  Horizons.h5/AhA.dir/CoordCenterInertial.dat
  Matter.h5/MatterCenterOfMass.dat
nsns:
  CoM-NS{A,B}-InertialFrame.dat
nsnsH5:
  Matter.h5/CoM-NS{A,B}-InertialFrame.dat
""")

  p.add_argument("--tmin",type=float, metavar="FLOAT",
    help="""Fit points with t>tmin.
Default determined by FindTmin estimation scheme.""")
  p.add_argument("--tmax",type=float, metavar="FLOAT",
    help="""Fit points with t<tmax.
Defaults to min(t.back(), 5pi/Omega0+tmin).""")

  p.add_argument("--improved_Omega0_update", action="store_true",
    help="""If specified, multiply Omega0-update by factor Omega0/omega_r.
Empirically, this is recommended to improve convergence.
See ticket #305 for details.""")
  p.add_argument("--Omega500", type=float, metavar="FLOAT",
    help="""If specified, compute updating formulae, such that: (1) e=0
(2) Omega(t=500)=Omega500. This will give changes to Omega0,
adot0 and D0, and all three must be applied. The updating
formulae are only computed for F2cos2.""")
  p.add_argument("--no_check", action="store_true",
    help="""If specified, do not check for large or negative periastron
advances. This will prevent the eccentricity from being
set to 9.9999""")
  p.add_argument("--no_plot", action="store_true",
    help="If specified, do not generate plots.")
  p.add_argument("--unbounded_fits", action="store_true",
    help="Use old unbounded fits.")

  #==== parse the options
  args = p.parse_args()
  summary = open("summary.txt","w")
  main(args)
  summary.close()
  System("cat summary.txt")
