#!/usr/bin/env python
from __future__ import division,print_function
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
from scipy.interpolate import UnivariateSpline
from numpy.linalg import inv

"""
Performs eccentricity fitting, and outputs the estimated eccentricity and updated initial data parameters.
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
  # (previously this was the average from 500 TIME SAMPLES until the end,
  # and it failed if you had less than 500 time samples.
  restricted_oddot = [oddot for (tee, oddot) in zip(t,ODblDot) if tee > 500.0]
  yfin=sum(restricted_oddot)/len(restricted_oddot)
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
  #return 200

#==============================================================================
# fit the data (t,y) to the model F[p,t], by least-squares fitting the params p.
def fit(t, y, F, p0, bounds, jac, name):
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
        'maxiter': 10000,
        'eta': 0.5,
      },
      args = (t,y),
    )

    if not res.success:
        warning("minimize failed in {}: {}".format(name, res.message))

    return res.x, res.fun, res.success


# fit the data (t,y) to the model F[p,t], by least-squares fitting the params p.
def fitKepler(t, y, F, p0, bounds, jac, name):
    # t,y -- arrays of equal length, having t and y values of data
    # F   -- function: F(p,t) taking parameters p and set of t-values
    # p0  -- starting values of parameters
    errfunc = lambda p : sqrt(sum((F(p,t[i])-y[i])**2 for i in range(0,t.size))/t.size)
    jacfunc = lambda p : array([ 1/errfunc(p) * (sum((F(p,t[i])-y[i])*dFdp(p,t[i]) for i in range(0,t.size))/t.size) for dFdp in jac ])

    res = optimize.minimize(
      errfunc,
      method = 'TNC',  #options: SLSQP, L-BFGS-B, TNC
      x0 = p0,
      bounds = bounds,
      jac = jacfunc,
      tol = 1e-5,
      options = {
        'disp': False,
        'maxiter': 10000,
        'eta': 0.5,
      },
      args = (),
    )

    if not res.success:
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
    plt.legend(loc='upper right', bbox_to_anchor=(1.15,1.3),labelspacing=0.25,handletextpad=0.0,fancybox=True)

    # add point at begin of fit interval:
    plt.plot([xBdry[0]],[1e6*F(p,xBdry)[0]],'o')

    # bottom-right plot: zoom of residual
    plt.subplot(2,2,4)
    plt.plot(x[idxZoom],1e6*(F(p,x[idxZoom])-y[idxZoom]),style,label=name)
    plt.plot([xBdry[0]],[(1e6*(F(p,xBdry)-yBdry))[0]],'o')
    plt.title("1e6 residual")

    return

################################################################
def ComputeUpdate(InParTrue,Tref, OrElT,eta, afit, efit, dfit, name, tmin, tmax, rms, Source):
    Amp = lambda a,e: -2*e*sqrt(1-e*e)*a**(-3)+e*(4*(4-eta)+e*e*(5*eta-22))/sqrt(1-e*e)*a**(-4)
    Om = lambda a,e: 1/a/sqrt(a)+(eta-9)*0.5/a/a/sqrt(a)
    Ct = lambda a,e: e-e*(8-3*eta)*0.5/a
    Cf = lambda a,e: e+e*eta*0.5/a

    uoft = lambda a, e, d, t : Om(a,e)*t + d

    rmodel = lambda a,e,d,t: a*(1-e*cos(uoft(a,e,d,t)))
    rdotmodel = lambda a,e,d,t:  Om(a,e)*a*e*sin(uoft(a,e,d,t))/(1-Ct(a,e)*cos(uoft(a,e,d,t)))
    Omegamodel = lambda a,e,d,t: (sqrt(1-e*e)*a**(-1.5)-0.5*((3-eta)+e*e*(2*eta-9))*a**(-2.5)/sqrt(1-e*e))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-1)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-1)

    etfit = Ct(afit,efit)
    ephifit = Cf(afit,efit)
    OmR = Om(afit,efit)/InParTrue[0]
    rmsR = rms/Amp(afit,efit)
    aT = OrElT[0]
    eT = OrElT[1]
    dT = OrElT[2]
    q=2
    a0=afit
    aold = np.power(-64./5./q*(1-1./q)*4*(1+73./24.*efit*efit+37./96.*efit*efit*efit*efit)*(-Tref)+a0*a0*a0*a0,1./4.) 

    a0=aT
    aTold = np.power(-64./5./q*(1-1./q)*4*(1+73./24.*eT*eT+37./96.*eT*eT*eT*eT)*(-Tref)+a0*a0*a0*a0,1./4.) 

    InParFit=[Omegamodel(afit,efit,dfit,0),rdotmodel(afit,efit,dfit,0)/rmodel(afit,efit,dfit,0),rmodel(afit,efit,dfit,0)]
    InParT=[Omegamodel(aT,eT,dT,0),rdotmodel(aT,eT,dT,0)/rmodel(aT,eT,dT,0),rmodel(aT,eT,dT,0)]

    #InParFit=[Omegamodel(aold,efit,dfit,0),rdotmodel(aold,efit,dfit,0)/rmodel(aold,efit,dfit,0),rmodel(aold,efit,dfit,0)]
    #InParT=[Omegamodel(aTold,eT,dT,0),rdotmodel(aTold,eT,dT,0)/rmodel(aTold,eT,dT,0),rmodel(aTold,eT,dT,0)]


    print(InParFit)
    dOmega = InParT[0]-InParFit[0]
    drdot = InParT[1]-InParFit[1]
    dr = InParT[2]-InParFit[2]

    if rmsR>0.4:
      summary.write("Large residual of ecc-fit\n")

    summary.write("%s:  %+11.8f    %+11.8f     %+9.6f\n"\
        %(name,dOmega,drdot,dr))

    f = open("Params_%s.dat"%name, 'w')
    f.write("# InitialDataAdjustment.py utilizing orbital frequency Omega, fit %s\n"%name)
    f.write("# Source file=%s\n" % Source)
    f.write("# Omega0=%10.8g, adot0=%10.8g, D0=%10.8g\n"%(InParTrue[0],InParTrue[1],InParTrue[2]))
    f.write("# Fitting interval [tmin,tmax]=[%g,%g]\n"%(tmin,tmax))
    f.write("# Parameters of the orbit: (a,e,delta)=(%g,%g,%g)\n"%(afit,efit,dfit))
    f.write("# Eccentricity parameters: (et,epfi)=(%g,%g)\n"%(etfit,ephifit))
    f.write("# [1] = Omega0\n")
    f.write("# [2] = 1e4 adot0\n")
    f.write("# [3] = D0\n")
    f.write("# [4] = e0\n")
    f.write("%10.12f\t%10.12f\t%10.12f\t%10.12f\t\n"%(InParTrue[0], 1e4*InParTrue[1], InParTrue[2], efit))
    f.write("%10.12f\t%10.12f\t%10.12f\t%10.12f\t\n"%(InParTrue[0]+dOmega, 1e4*(InParTrue[1]+drdot), InParTrue[2]+dr, eT))
    f.close()

    f = open("Fit_%s.dat"%name, 'w')
    f.write("# InitialDataAdjustment.py utilizing orbital frequency Omega, fit %s\n"%name)
    f.write("# Source file=%s\n" % Source)
    f.write("# Omega0=%10.8g, adot0=%10.8g, D0=%10.8g\n"%(InParTrue[0],InParTrue[1],InParTrue[2]))
    f.write("# Fitting interval [tmin,tmax]=[%g,%g]\n"%(tmin,tmax))
    f.write("# Parameters of the orbit: (a,e,delta)=(%g,%g,%g)\n"%(afit,efit,dfit))
    f.write("# Eccentricity parameters: (et,epfi)=(%g,%g)\n"%(etfit,ephifit))
    f.write("# [1] = Tstart\n")
    f.write("# [2] = Tend\n")
    f.write("# [4] = a (PN semi major axis)\n")
    f.write("# [5] = e (PN eccentricity)\n")
    f.write("# [6] = delta (PN eccentric anomaly constant) \n")
    f.write("# [7] = rms residual\n")
    f.write("# [8] = rms residual/Amp\n")
    f.write("# [9] = omega/Omega0\n")

    f.write("%g %g %g %g %g %g %g %g\n"%(tmin, tmax, afit, efit, dfit, rms, rmsR, OmR))
    f.close()

    return InParTrue[0]+dOmega, InParTrue[1]+drdot, InParTrue[2]+dr

def ComputeSpinUpdate(opts, XA, XB, max_tmax_fit, OmegaVec, t, tmin, Tref, SpinATref, SpinBTref, theOmega0):
  Horizons = ReadH5(os.path.join(opts.d, "Horizons.h5"))
  sA = Horizons["AhA.dir/DimensionfulInertialSpin.dat"]
  sB = Horizons["AhB.dir/DimensionfulInertialSpin.dat"]
  MA = Horizons["AhA.dir/ChristodoulouMass.dat"]
  MB = Horizons["AhB.dir/ChristodoulouMass.dat"]

  m_idx = find_nearest(MA[:,0],tmin)
  mA = MA[m_idx,1]
  mB = MB[m_idx,1]
  m = mA + mB
  tmax = 1.5 * max_tmax_fit
  time = XA[:,0]
  sAvec = sA[time<tmax,1:4]
  sBvec = sB[time<tmax,1:4]

  #SpinATref and SpinBTref are defined in a system where the BHs are along the x axis
  nvec = (XA[:,1:4]-XB[:,1:4])
  ellvec = OmegaVec

  nvecinterp = [InterpolatedUnivariateSpline(time,nvec[:,0],k=5),InterpolatedUnivariateSpline(time,nvec[:,1],k=5),InterpolatedUnivariateSpline(time,nvec[:,2],k=5)]
  ellvecinterp = [InterpolatedUnivariateSpline(t,ellvec[:,0],k=5),InterpolatedUnivariateSpline(t,ellvec[:,1],k=5),InterpolatedUnivariateSpline(t,ellvec[:,2],k=5)]


  nvecTref = np.array([nvecinterp[0](Tref),nvecinterp[1](Tref),nvecinterp[2](Tref)])
  nvecT0 = np.array([nvecinterp[0](0),nvecinterp[1](0),nvecinterp[2](0)])
  nhatTref= nvecTref/norm(nvecTref)
  ellvecTref = np.array([ellvecinterp[0](Tref),ellvecinterp[1](Tref),ellvecinterp[2](Tref)])
  ellhatTref = ellvecTref/norm(ellvecTref)
  lambdahatTref = np.cross(ellhatTref,nhatTref)


  #transform the target spins from the BH frame to the cartesian frame of the simulation
  T = coordtrans(nhatTref,lambdahatTref,ellhatTref)
  Tinv = coordtransback(nhatTref,lambdahatTref,ellhatTref)
  sAcart = list(np.squeeze(np.asarray(dot(T,SpinATref))))
  sBcart = list(np.squeeze(np.asarray(dot(T,SpinBTref))))


  sAvecinterp = [InterpolatedUnivariateSpline(time,sAvec[:,0],k=5),InterpolatedUnivariateSpline(time,sAvec[:,1],k=5),InterpolatedUnivariateSpline(time,sAvec[:,2],k=5)]
  sBvecinterp = [InterpolatedUnivariateSpline(time,sBvec[:,0],k=5),InterpolatedUnivariateSpline(time,sBvec[:,1],k=5),InterpolatedUnivariateSpline(time,sBvec[:,2],k=5)]
  # Jinitial=[sAvecinterp[0](0)+sBvecinterp[0](0),sAvecinterp[1](0)+sBvecinterp[1](0),sAvecinterp[2](0)+sBvecinterp[2](0)+mA*mB*(m*theOmega0)**(-1/3)]
  # Jinitialunit=Jinitial/norm(Jinitial)
  
  sAvecinitial = np.array([sAvecinterp[0](0),sAvecinterp[1](0),sAvecinterp[2](0)])
  sAvecinitialunit = sAvecinitial/norm(sAvecinitial)
  sAvecTref = np.array([sAvecinterp[0](Tref),sAvecinterp[1](Tref),sAvecinterp[2](Tref)])
  sAvecTrefunit = sAvecTref/norm(sAvecTref)

  #Rotation perpendicular to the S1(0), S1(Tref) plane
  cossAinTref = np.inner(sAvecinitialunit,sAvecTrefunit)
  crosssAinTref = np.cross(sAvecinitialunit,sAvecTrefunit)
  RA = np.identity(3) + skew(crosssAinTref) + 1/(1+cossAinTref)*skew(crosssAinTref)**2 
  RAinv = inv(RA)
  newsA = dot(RAinv,sAcart)
  # print sAvecinitial/mA/mA
  # print newsA

  # #Rotation perpendicular to the S1(0), S1(Tref) plane
  # cossAinTref = np.inner(sAvecTrefunit,sAcart/norm(sAcart))
  # crosssAinTref = np.cross(sAvecTrefunit,sAcart/norm(sAcart))
  # RA = np.identity(3) + skew(crosssAinTref) + 1/(1+cossAinTref)*skew(crosssAinTref)**2 
  # RAinv = inv(RA)
  # newsA = dot(RA,sAvecinitialunit)*norm(sAcart)

  # print newsA

  if (abs(norm(sAvecTref/(mA*mA))- norm(sAcart)) > 1e-3):
    summary.write("""The magnitude of the target and the initial spin A are not equal (target sA = %g, initial sA = %g). Rescaling the initial spin magnitude to match the initial spin magnitude. New initial sA = %g\n"""
      %(norm(sAcart),norm(sAvecTref/(mA*mA)), norm(newsA)))

  # #Rotation around J
  # sAinitialJcomp=(sAvecinitial-np.inner(sAvecinitial,Jinitialunit)*Jinitialunit)/norm(sAvecinitial-np.inner(sAvecinitial,Jinitialunit)*Jinitialunit)
  # sATrefJcomp=(sAvecTref-np.inner(sAvecTref,Jinitialunit)*Jinitialunit)/norm(sAvecTref-np.inner(sAvecTref,Jinitialunit)*Jinitialunit)
  # costheta = np.inner(sAinitialJcomp,sATrefJcomp)
  # RA = costheta*np.identity(3) + sqrt(1-costheta*costheta)*skew(Jinitialunit) + (1-costheta)*outer(Jinitialunit) 
  # RAinv = inv(RA)
  # newsA = dot(RAinv,sAcart)

  sAcart_str = str(sAcart)
  sAvecTref_str=str(list(sAvecTref/(mA*mA)))
  newsA_str=str(list(np.squeeze(np.asarray(newsA))))
  sAdiff_str=str(list(SpinATref-list(np.squeeze(np.asarray(dot(Tinv,sAvecTref))))/(mA*mA)))
  summary.write("""Spin A adjustment: updated spin: %s, difference: %s\n\n"""
    %(newsA_str,sAdiff_str))

  sBvecinitial = np.array([sBvecinterp[0](0),sBvecinterp[1](0),sBvecinterp[2](0)])
  sBvecinitialunit = sBvecinitial/norm(sBvecinitial)
  sBvecTref = np.array([sBvecinterp[0](Tref),sBvecinterp[1](Tref),sBvecinterp[2](Tref)])
  sBvecTrefunit = sBvecTref/norm(sBvecTref)

  #Rotation perpendicular to the S2(0), S2(Tref) plane
  cossBinTref = np.inner(sBvecinitialunit,sBvecTrefunit)
  crosssBinTref = np.cross(sBvecinitialunit,sBvecTrefunit)
  RB = np.identity(3) + skew(crosssBinTref) + 1/(1+cossBinTref)*skew(crosssBinTref)**2 
  RBinv = inv(RB)
  newsB = dot(RBinv,sBcart)
  # print sBvecinitial/mB/mB
  # print newsB

  # #Rotation perpendicular to the S1(0), S1(Tref) plane
  # cossBinTref = np.inner(sBvecTrefunit,sBcart/norm(sBcart))
  # crosssBinTref = np.cross(sBvecTrefunit,sBcart/norm(sBcart))
  # RB = np.identity(3) + skew(crosssBinTref) + 1/(1+cossBinTref)*skew(crosssBinTref)**2 
  # RBinv = inv(RB)
  # newsB = dot(RB,sBvecinitialunit)*norm(sBcart)

  # print newsB

  if (abs(norm(sBvecTref/(mB*mB))- norm(sBcart)) > 1e-3):
    summary.write("""The magnitude of the target and the initial spin B are not equal (target sB = %g, initial sB = %g). Rescaling the initial spin magnitude to match the initial spin magnitude. New initial sB = %g\n"""
      %(norm(sBcart),norm(sBvecTref/(mB*mB)), norm(newsB)))

  # #Rotation around J
  # sBinitialJcomp=(sBvecinitial-np.inner(sBvecinitial,Jinitialunit)*Jinitialunit)/norm(sBvecinitial-np.inner(sBvecinitial,Jinitialunit)*Jinitialunit)
  # sBTrefJcomp=(sBvecTref-np.inner(sBvecTref,Jinitialunit)*Jinitialunit)/norm(sBvecTref-np.inner(sBvecTref,Jinitialunit)*Jinitialunit)
  # costheta = np.inner(sBinitialJcomp,sBTrefJcomp)
  # RB = costheta*np.identity(3) + sqrt(1-costheta*costheta)*skew(Jinitialunit) + (1-costheta)*outer(Jinitialunit) 
  # RBinv = inv(RB)
  # newsB = dot(RBinv,sBcart)

  sBcart_str = str(sBcart)
  sBvecTref_str = str(list(sBvecTref/(mB*mB)))
  newsB_str = str(list(np.squeeze(np.asarray(newsB))))
  sBdiff_str = str(list(SpinBTref-list(np.squeeze(np.asarray(dot(Tinv,sBvecTref))))/(mB*mB)))
  summary.write("""Spin B adjustment: updated spin: %s, difference: %s\n"""
    %(newsB_str,sBdiff_str))


  f = open("SpinUpdates.dat", 'w')
  f.write("# EccRemoval.py utilizing data from source file=%s and target spin values from file=%s\n" % (opts.d,opts.p))
  f.write("# Rotate spins around the plane perpendicular to Si(0) and Si(Tref)")
  f.write("Spin A adjustment: target spin at %d : %s,\ncurrect spin at %d: %s,\nupdated spin: %s,\ndifference: %s\n\n"
    %(Tref, sAcart_str, Tref,sAvecTref_str,newsA_str,sAdiff_str))
  f.write("Spin B adjustment: target spin at %d : %s,\ncurrect spin at %d: %s,\nupdated spin: %s,\ndifference: %s\n\n"
    %(Tref, sBcart_str, Tref,sBvecTref_str,newsB_str,sBdiff_str))
  f.close() 

  return list(np.squeeze(np.asarray(newsA))),list(np.squeeze(np.asarray(newsB)))

################################################################

def ParseParams(path):
    File = os.path.basename(os.path.realpath(path))
    Dir  = os.path.dirname(os.path.realpath(path))
    D = ParseIdParams(Dir, file=File, array=True)

    Omega0 = float(D['Omega0'][0])
    adot0  = float(D['adot0'][0])
    D0     = float(D['D0'][0])

    return Omega0,adot0,D0

def ParseTargetParams(path):
    File = os.path.basename(os.path.realpath(path))
    Dir  = os.path.dirname(os.path.realpath(path))
    D = ParseIdParams(Dir, file=File, array=True)

    Tref = float(D['ReferenceTime'][0])
    eT  = float(D['Eccentricity'][0])
    aT  = float(D['SemiMajorAxis'][0])
    dT  = float(D['AnomalyAngle'][0])
    
    MassRatio = float(D['MassRatio'][0])
    SpinATref = [float(D['SpinA'][0]),float(D['SpinA'][1]),float(D['SpinA'][2])]
    SpinBTref = [float(D['SpinB'][0]),float(D['SpinB'][1]),float(D['SpinB'][2])]
    IDType = D['IDType'][0]

    return Tref,eT,aT,dT,MassRatio,SpinATref,SpinBTref,IDType

def PointwiseDot(A,B):
      '''Dot product for each element in two lists of vectors'''
      return np.array([np.dot(a,b) for a,b in zip(A,B)])

def PointwiseNormalize(A):
      '''Normalization for each element in a list of vectors'''
      return A/norm(A, axis=1)[:,np.newaxis]

def GetTrajectories(Dir, tmax_fit):
  Horizons = ReadH5(os.path.join(Dir,"Horizons.h5"))
  TrajA = Horizons["AhA.dir/CoordCenterInertial.dat"]
  TrajB = Horizons["AhB.dir/CoordCenterInertial.dat"]


  # Truncate the trajectories to a region around the fit interval
  tmax = 1.5 * tmax_fit
  tA = TrajA[:,0]
  tB = TrajB[:,0]
  return TrajA[tA<tmax,:], TrajB[tB<tmax,:]

def ComputeOmegaAndDerivsFromFile(TrajA, TrajB):
  '''Compute Omega, dOmega, ellHat and nHat'''

  # While method="Smooth" appears to be more accurate, "Fit"
  # does a slightly better job of reducing eccentricity.
  t_raw, OmegaVec = Compute_OrbitalFrequency(TrajA, TrajB, 10, method='Smooth')
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

def skew(x):
    return np.matrix([[0, -x[2], x[1]],
                     [x[2], 0, -x[0]],
                     [-x[1], x[0], 0]])

def outer(x):
    return np.matrix([[x[0]*x[0], x[0]*x[1], x[0]*x[2]],
                     [x[1]*x[0], x[1]*x[1], x[1]*x[2]],
                     [x[2]*x[0], x[2]*x[1], x[2]*x[2]]])

def coordtrans(nhat,lambdahat,ellhat):
    return np.matrix([[nhat[0], lambdahat[0], ellhat[0]],
                      [nhat[1], lambdahat[1], ellhat[1]],
                      [nhat[2], lambdahat[2], ellhat[2]]])

def coordtransback(nhat,lambdahat,ellhat):
    return np.matrix([[nhat[0], nhat[1], nhat[2]],
                      [lambdahat[0], lambdahat[1], lambdahat[2]],
                      [ellhat[0], ellhat[1], ellhat[2]]])


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
    self.phi   = [-pi, pi]
    # eccentricity is between 0 and 1
    self.ecc   = [0, 1]
    # For various variables that have no limits
    self.none  = [None, None]
    # Separation
    self.a  = [2, 100]



################################################################

def perform_fits_SS(t,Omega,dOmegadt,idxFit,tFit,OmegaFit,dOmegadtFit,T_merge,Amp,alpha_intrp,S_0_perp_n,tmin,tmax,Tref,eT,aT,dT,MassRatio,theOmega0,theD0,theadot0,opts,Source):
  
  # Set bounds for some of the variables
  lim = FitBounds(tmax, theOmega0)

  eta = MassRatio/(1+MassRatio)/(1+MassRatio)

  # Fit *Omega*
  OmegaFunc = lambda p,t: p[1]*(p[0]-t)**(-3/8)   # power law form for Omega
  jac = [ lambda p,t: (-3/8)*p[1]*(p[0]-t)**(-11/8),
          lambda p,t:             (p[0]-t)**(-3/8) ]
  p0 = [T_merge, Amp]
  pBounds = [lim.Tc, lim.none]
  pOmega, rmsOmega, Omega_status = fit(tFit, OmegaFit, OmegaFunc,p0, pBounds, jac, "OmegaFunc")
  FitTc = pOmega[0]  # Fit time to coalescence

  #================
  F1_SS = lambda p,t: p[0]*(FitTc-t)**(-11/8)
  jac = [ lambda p,t: (FitTc-t)**(-11/8) ]
  p0 = [1e-5]
  pBounds = [lim.none]
  pF1_SS, rmsF1_SS, F1_SS_status = fit(tFit, dOmegadtFit, F1_SS, p0, pBounds, jac, "F1_SS")
  #import datetime


  Amp = lambda a,e: -2*e*sqrt(1-e*e)*a**(-3)+e*(4*(4-eta)+e*e*(5*eta-22))/sqrt(1-e*e)*a**(-4)
  B = lambda a,e: e+e/a*(eta-2)
  Om = lambda a,e: 1/a/sqrt(a)+(eta-9)*0.5/a/a/sqrt(a)
  Ct = lambda a,e: e-e*(8-3*eta)*0.5/a
  Cf = lambda a,e: e+e*eta*0.5/a

  # Keplereq = lambda u, a, e, d, t : u-Ct(a,e)*sin(u)-Om(a,e)*t-d
  # Keplereqprime = lambda u, a, e, d, t: 1-Ct(a,e)*cos(u)
  # Keplereqprime2 = lambda u, a, e, d, t: Ct(a,e)*sin(u)
  # uoft = lambda a, e, d, t : optimize.newton(Keplereq, Om(a,e)*t+d, fprime=Keplereqprime, args=(a,e,d,t), maxiter=10,fprime2=Keplereqprime2) 

  uoft = lambda a, e, d, t : Om(a,e)*t + d

  Omegadotmodel = lambda a,e,d,t:  Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)

  dAmpde = lambda a,e: 2/sqrt(1-e*e)*(2*e*e-1)*a**(-3)+(4*(4-eta)+e*e*(5*eta-22)*(3-2*e*e))/sqrt(1-e*e)/(1-e*e)*a**(-4)
  dAmpda = lambda a,e: 6*e*sqrt(1-e*e)*a**(-4)-4*e*(4*(4-eta)+e*e*(5*eta-22))/sqrt(1-e*e)*a**(-5)
  dBde = lambda a,e: 1+1/a*(eta-2)
  dBda = lambda a,e: -e/a/a*(eta-2)
  dOmda = lambda a,e: -1.5/a/a/sqrt(a)-2.5*(eta-9)*0.5/a/a/a/sqrt(a)
  dCtde = lambda a,e: 1-(8-3*eta)*0.5/a
  dCtda = lambda a,e: e*(8-3*eta)*0.5/a/a
  dCfde = lambda a,e: 1+eta*0.5/a
  dCfda = lambda a,e: -e*eta*0.5/a/a

  dOmegadotmodeldAmp= lambda a,e,d,t:  (1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)
  dOmegadotmodeldB= lambda a,e,d,t: -Amp(a,e)*cos(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)
  # dOmegadotmodeldOm = lambda a,e,d,t: t*Amp(a,e)*B(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      + t*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 2*t*Cf(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 3*t*Ct(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-5)
  # dOmegadotmodeldd = lambda a,e,d,t: Amp(a,e)*B(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      + Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 2*Cf(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 3*Ct(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-5)
  # dOmegadotmodeldCt = lambda a,e,d,t:  3*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      + Amp(a,e)*B(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      + Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*cos(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 2*Cf(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)\
  #                                      - 3*Ct(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-5)
  
  dOmegadotmodeldOm = lambda a,e,d,t: t*Amp(a,e)*B(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       + t*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       - 2*t*Cf(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       - 3*t*Ct(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)
  dOmegadotmodeldd = lambda a,e,d,t: Amp(a,e)*B(a,e)*sin(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       + Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       - 2*Cf(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)\
                                       - 3*Ct(a,e)*sin(uoft(a,e,d,t))*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)
  dOmegadotmodeldCt = lambda a,e,d,t:  3*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-2)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-4)


  dOmegadotmodeldCf = lambda a,e,d,t:  2*Amp(a,e)*(1-B(a,e)*cos(uoft(a,e,d,t)))*cos(uoft(a,e,d,t))*sin(uoft(a,e,d,t))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-3)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-3)

  #print datetime.datetime.now()
  #================
  rmsF1cos1_SS = 3*rmsF1_SS
  F1cos1_SS = lambda p,t : p[3]*(FitTc-t)**(-11/8) +Omegadotmodel(p[0],p[1],p[2],t) - p[4]*sin(2*alpha_intrp(t) + p[5])

  jac = [ lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2],t)*dAmpda(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2],t)*dBda(p[0],p[1])\
  +dOmegadotmodeldOm(p[0],p[1],p[2],t)*dOmda(p[0],p[1])+dOmegadotmodeldCf(p[0],p[1],p[2],t)*dCfda(p[0],p[1])\
  +dOmegadotmodeldCt(p[0],p[1],p[2],t)*dCtda(p[0],p[1]),
            lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2],t)*dAmpde(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2],t)*dBde(p[0],p[1])\
            +dOmegadotmodeldCf(p[0],p[1],p[2],t)*dCfde(p[0],p[1])+dOmegadotmodeldCt(p[0],p[1],p[2],t)*dCtde(p[0],p[1]),
            lambda p,t:   dOmegadotmodeldd(p[0],p[1],p[2],t),
            lambda p,t:   (FitTc-t)**(-11/8),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[5]),
          lambda p,t: -p[4]*cos(2*alpha_intrp(t) + p[5]) ]

  for phi in range(-3, 4):
    p0 = [aT, eT, phi, pF1_SS[0], 0.5*S_0_perp_n[0]**2 * (theOmega0/theD0)**2, dT]
    pBounds = [lim.a, lim.ecc, lim.phi, lim.none, lim.none, lim.phi]
    #print "here"
    ptemp, rmstemp, F1cos1_SS_status = fit(tFit, dOmegadtFit, F1cos1_SS, p0, pBounds, jac, "F1cos1_SS")
    
    if(rmstemp <= rmsF1cos1_SS):
      rmsF1cos1_SS = rmstemp
      pF1cos1_SS = ptemp

  #print "here1"

  #================
  F1cos2_SS = lambda p,t: p[3]*(FitTc-t)**(-11/8) +Omegadotmodel(p[0],p[1],p[2]+p[6]*t*t,t) - p[4]*sin(2*alpha_intrp(t) + p[5])

  jac = [ lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2]+p[6]*t*t,t)*dAmpda(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2]+p[6]*t*t,t)*dBda(p[0],p[1])\
  +dOmegadotmodeldOm(p[0],p[1],p[2]+p[6]*t*t,t)*dOmda(p[0],p[1])+dOmegadotmodeldCf(p[0],p[1],p[2]+p[6]*t*t,t)*dCfda(p[0],p[1])\
  +dOmegadotmodeldCt(p[0],p[1],p[2]+p[6]*t*t,t)*dCtda(p[0],p[1]),
            lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2]+p[6]*t*t,t)*dAmpde(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2]+p[6]*t*t,t)*dBde(p[0],p[1])\
            +dOmegadotmodeldCf(p[0],p[1],p[2]+p[6]*t*t,t)*dCfde(p[0],p[1])+dOmegadotmodeldCt(p[0],p[1],p[2]+p[6]*t*t,t)*dCtde(p[0],p[1]),
            lambda p,t:   dOmegadotmodeldd(p[0],p[1],p[2]+p[6]*t*t,t),
            lambda p,t:   (FitTc-t)**(-11/8),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[5]),
          lambda p,t: -p[4]*cos(2*alpha_intrp(t) + p[5]),
           lambda p,t:   t*t*dOmegadotmodeldd(p[0],p[1],p[2]+p[6]*t*t,t)]

  p0 = [pF1cos1_SS[0], pF1cos1_SS[1], pF1cos1_SS[2], pF1cos1_SS[3], pF1cos1_SS[4], pF1cos1_SS[5],0]
  pBounds = [lim.a, lim.ecc, lim.phi, lim.none, lim.none, lim.phi, lim.none]
  pF1cos2_SS, rmsF1cos2_SS, F1cos2_SS_status = fit(tFit, dOmegadtFit, F1cos2_SS, p0, pBounds, jac, "F1cos2_SS")

  #print "here2"
  #================
  F2cos1_SS = lambda p,t : p[3]*(FitTc-t)**(-11/8) + p[4]*(FitTc-t)**(-13/8) +Omegadotmodel(p[0],p[1],p[2],t) - p[5]*sin(2*alpha_intrp(t) + p[6])
  jac = [ lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2],t)*dAmpda(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2],t)*dBda(p[0],p[1])\
  +dOmegadotmodeldOm(p[0],p[1],p[2],t)*dOmda(p[0],p[1])+dOmegadotmodeldCf(p[0],p[1],p[2],t)*dCfda(p[0],p[1])\
  +dOmegadotmodeldCt(p[0],p[1],p[2],t)*dCtda(p[0],p[1]),
            lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2],t)*dAmpde(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2],t)*dBde(p[0],p[1])\
            +dOmegadotmodeldCf(p[0],p[1],p[2],t)*dCfde(p[0],p[1])+dOmegadotmodeldCt(p[0],p[1],p[2],t)*dCtde(p[0],p[1]),
            lambda p,t:   dOmegadotmodeldd(p[0],p[1],p[2],t),
            lambda p,t:   (FitTc-t)**(-11/8),
            lambda p,t:   (FitTc-t)**(-13/8),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[6]),
          lambda p,t: -p[5]*cos(2*alpha_intrp(t) + p[6]) ]


  p0 = [pF1cos1_SS[0], pF1cos1_SS[1], pF1cos1_SS[2], pF1cos1_SS[3], pF1cos1_SS[3]/100, pF1cos1_SS[4], pF1cos1_SS[5]]     
  pBounds = [lim.a, lim.ecc, lim.phi, lim.none,lim.none, lim.none, lim.phi]
  pF2cos1_SS, rmsF2cos1_SS, F2cos1_SS_status = fit(tFit, dOmegadtFit, F2cos1_SS, p0, pBounds, jac, "F2cos1_SS")

  #print "here3"
  #================
  F2cos2_SS = lambda p,t: p[3]*(FitTc-t)**(-11/8) + p[4]*(FitTc-t)**(-13/8) +Omegadotmodel(p[0],p[1],p[2]+p[7]*t*t,t) - p[5]*sin(2*alpha_intrp(t) + p[6])
  jac = [ lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2]+p[7]*t*t,t)*dAmpda(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2]+p[7]*t*t,t)*dBda(p[0],p[1])\
  +dOmegadotmodeldOm(p[0],p[1],p[2]+p[7]*t*t,t)*dOmda(p[0],p[1])+dOmegadotmodeldCf(p[0],p[1],p[2]+p[7]*t*t,t)*dCfda(p[0],p[1])\
  +dOmegadotmodeldCt(p[0],p[1],p[2]+p[7]*t*t,t)*dCtda(p[0],p[1]),
            lambda p,t: dOmegadotmodeldAmp(p[0],p[1],p[2]+p[7]*t*t,t)*dAmpde(p[0],p[1])+dOmegadotmodeldB(p[0],p[1],p[2]+p[7]*t*t,t)*dBde(p[0],p[1])\
            +dOmegadotmodeldCf(p[0],p[1],p[2]+p[7]*t*t,t)*dCfde(p[0],p[1])+dOmegadotmodeldCt(p[0],p[1],p[2]+p[7]*t*t,t)*dCtde(p[0],p[1]),
            lambda p,t:   dOmegadotmodeldd(p[0],p[1],p[2]+p[7]*t*t,t),
            lambda p,t:   (FitTc-t)**(-11/8),
            lambda p,t:   (FitTc-t)**(-13/8),
          lambda p,t:      -sin(2*alpha_intrp(t) + p[6]),
          lambda p,t: -p[5]*cos(2*alpha_intrp(t) + p[6]),
           lambda p,t:   t*t*dOmegadotmodeldd(p[0],p[1],p[2]+p[7]*t*t,t)]
  
  p0 = [pF1cos2_SS[0], pF1cos2_SS[1], pF1cos2_SS[2], pF1cos2_SS[3], pF1cos2_SS[3]/100, pF1cos2_SS[4], pF1cos2_SS[5], pF1cos2_SS[6]] 
  pBounds = [lim.a, lim.ecc, lim.phi, lim.none,lim.none, lim.none, lim.phi, lim.none]
  pF2cos2_SS, rmsF2cos2_SS, F2cos2_SS_status = fit(tFit, dOmegadtFit, F2cos2_SS,p0, pBounds, jac, "F2cos2_SS")
  
  #print "here4"
  #print datetime.datetime.now()
  fp = open("FitSuccess_SS.txt","w")
  fp.write("FOmegalin  %d\n" % Omega_status)
  fp.write("F1  %d\n" % F1_SS_status)
  fp.write("F1cos1_SS  %d\n" % F1cos1_SS_status)
  fp.write("F1cos2_SS  %d\n" % F1cos2_SS_status)
  fp.write("F2cos1_SS  %d\n" % F2cos1_SS_status)
  fp.write("F2cos2_SS  %d\n" % F2cos2_SS_status)
  fp.close()

  OrElT=[aT,eT,dT]
  InParTrue = [theOmega0,theadot0,theD0]

  #==== output all four residuals
  summary.write("tmin=%f  (determined by " % tmin)
  if(opts.tmin == None):
    summary.write(" fit)\n")
  else:
    summary.write(" option)\n")

  summary.write("""\nRESIDUALS
    F1cos1_SS rms=%g   \tF1cos2_SS rms=%g
    F2cos1_SS rms=%g   \tF2cos2_SS rms=%g
   \nDIAGNOSTICS
                a         ecc     d    rms/Amp   omega/Omega0    et     ephi\n""" 
  % (rmsF1cos1_SS, rmsF1cos2_SS, rmsF2cos1_SS, rmsF2cos2_SS))

  # add first two fields for remaining fits
  tmp = "%-8s %10.1f   %6.2e    %6.4f     %5.3f      %5.3f     %5.3f      %5.3f\n"

  summary.write(tmp%("F1cos1_SS", pF1cos1_SS[0], pF1cos1_SS[1], pF1cos1_SS[2], rmsF1cos1_SS/Amp(pF1cos1_SS[0],pF1cos1_SS[1]), Om(pF1cos1_SS[0],pF1cos1_SS[1])/theOmega0, Ct(pF1cos1_SS[0],pF1cos1_SS[1]),Cf(pF1cos1_SS[0],pF1cos1_SS[1])))
  summary.write(tmp%("F1cos2_SS", pF1cos2_SS[0], pF1cos2_SS[1], pF1cos2_SS[2], rmsF1cos2_SS/Amp(pF1cos2_SS[0],pF1cos2_SS[1]), Om(pF1cos2_SS[0],pF1cos2_SS[1])/theOmega0, Ct(pF1cos2_SS[0],pF1cos2_SS[1]),Cf(pF1cos2_SS[0],pF1cos2_SS[1])))
  summary.write(tmp%("F2cos1_SS", pF2cos1_SS[0], pF2cos1_SS[1], pF2cos1_SS[2], rmsF2cos1_SS/Amp(pF2cos1_SS[0],pF2cos1_SS[1]), Om(pF2cos1_SS[0],pF2cos1_SS[1])/theOmega0, Ct(pF2cos1_SS[0],pF2cos1_SS[1]),Cf(pF2cos1_SS[0],pF2cos1_SS[1])))
  summary.write(tmp%("F2cos2_SS", pF2cos2_SS[0], pF2cos2_SS[1], pF2cos2_SS[2], rmsF2cos2_SS/Amp(pF2cos2_SS[0],pF2cos2_SS[1]), Om(pF2cos2_SS[0],pF2cos2_SS[1])/theOmega0, Ct(pF2cos2_SS[0],pF2cos2_SS[1]),Cf(pF2cos2_SS[0],pF2cos2_SS[1])))

  #==== compute updates and generate EccRemoval_FIT.dat files
  summary.write("""\nECCENTRICITY AND UPDATES
             delta_Omega0   delta_adot0     delta_D0\n""")

  Omega0new,adot0new,D0new = ComputeUpdate(InParTrue, Tref, OrElT, eta, pF1cos1_SS[0], pF1cos1_SS[1], pF1cos1_SS[2], "F1cos1_SS", tmin, tmax, rmsF1cos1_SS, Source)
  Omega0new,adot0new,D0new = ComputeUpdate(InParTrue, Tref, OrElT, eta, pF1cos2_SS[0], pF1cos2_SS[1], pF1cos2_SS[2], "F1cos2_SS", tmin, tmax, rmsF1cos2_SS, Source)
  Omega0new,adot0new,D0new = ComputeUpdate(InParTrue, Tref, OrElT, eta, pF2cos1_SS[0], pF2cos1_SS[1], pF2cos1_SS[2], "F2cos1_SS", tmin, tmax, rmsF2cos1_SS, Source)
  Omega0new,adot0new,D0new = ComputeUpdate(InParTrue, Tref, OrElT, eta, pF2cos2_SS[0], pF2cos2_SS[1], pF2cos2_SS[2], "F2cos2_SS", tmin, tmax, rmsF2cos2_SS, Source)


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
  plot_fitted_function(pF1cos1_SS, F1cos1_SS, t, dOmegadt, idxFit, idxPlot, idxZoom, 'F1cos1_SS', ':')
  plot_fitted_function(pF1cos2_SS, F1cos2_SS, t, dOmegadt, idxFit, idxPlot, idxZoom, 'F1cos2_SS', ':')
  plot_fitted_function(pF2cos1_SS, F2cos1_SS, t, dOmegadt, idxFit, idxPlot, idxZoom, 'F2cos1_SS', '--')
  plot_fitted_function(pF2cos2_SS, F2cos2_SS, t, dOmegadt, idxFit, idxPlot, idxZoom, 'F2cos2_SS', '--')

  # zoom out of the y-axis in the lower left panel by 15%
  plt.subplot(2, 2, 3)
  q1 = plt.axis()
  Deltay = q1[3] - q1[2]
  plt.axis([q1[0], q1[1], q1[2] - 0.15 * Deltay, q1[3] + 0.15 * Deltay ])

  # adjust margins of figure
  plt.subplots_adjust(left=0.09, right=0.95, bottom=0.07, hspace=0.25)

  plt.savefig("FigureOmegadotFit_SS.pdf")

  return Omega0new,adot0new,D0new 

#################################################################

def main(opts):
  Source = opts.d
  theOmega0,theadot0,theD0 = ParseParams(opts.p)
  Tref,eT,aT,dT,MassRatio,SpinATref,SpinBTref,IDType \
    = ParseTargetParams(opts.t)

  # Determine the latest possible fit time so we can truncate data
  max_tmin_fit = Tref
  max_tmax_fit = max_tmin_fit + 5*pi/theOmega0
  #max_tmax_fit = max_tmin_fit + 4*2*3.14159*thea0*sqrt(thea0)/sqrt(1+theMassRatio)

  XA,XB = GetTrajectories(Source, max_tmax_fit)
  t,Omega,dOmegadt,OmegaVec = ComputeOmegaAndDerivsFromFile(XA,XB)

  print("Length of t is ",len(t))
  
  plt.close()
  plt.plot(t,Omega,label="omega")
  plt.legend()
  plt.show()

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

  # Do fits that include spin-spin interactions, see Buonnano et al., 2010 (arXiv 1012.1549v2).
  alpha_intrp, S_0_perp_n, T_merge, Amp = GetVarsFromSpinData(opts.d, XA, XB,OmegaVec, t, tmin)

  Omega0new,adot0new,D0new = perform_fits_SS(t,Omega,dOmegadt,idxFit,tFit,OmegaFit,dOmegadtFit,T_merge,Amp,alpha_intrp,\
    S_0_perp_n,tmin,tmax,Tref,eT,aT,dT,MassRatio,theOmega0,theD0,theadot0,opts,Source)

  # Now, compute the new values for the spins
  newsA,newsB=ComputeSpinUpdate(opts,XA,XB,max_tmax_fit,OmegaVec,t,tmin,Tref,SpinATref,SpinBTref,theOmega0)


  f = open("Params.input", 'w')
  f.write("# Set the initial data parameters\n")
  f.write("\n")
  f.write("# Orbital parameters\n")
  f.write("$Omega0 = %11.8f;\n" % Omega0new)
  f.write("$adot0 = %g;\n" % adot0new)
  f.write("$D0 = %11.8f;\n" % D0new)
  f.write("\n")
  f.write("# Physical parameters (spins are dimensionless)\n")
  f.write("$MassRatio = %11.8f;\n" % MassRatio)
  f.write("@SpinA = (%11.8f, %11.8f, %11.8f);\n" % (newsA[0],newsA[1],newsA[2]))
  f.write("@SpinB = (%11.8f, %11.8f, %11.8f);\n" % (newsB[0],newsB[1],newsB[2]))
  f.write("\n")
  f.write("# Evolve after initial data completes?\n")
  f.write("$Evolve = 1;\n")
  f.write('# IDType: "SKS", "SHK" or "CFMS".\n')
  f.write('$IDType = "%s";\n' % IDType)
  f.close()

#=============================================================================
if __name__ == "__main__":

  p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawTextHelpFormatter)

  p1 = p.add_argument_group("required arguments")

  p1.add_argument("-p", type=str, required=True, metavar="FILE",help="File with initial data variables,\ne.g. Params.input")
  p1.add_argument("-t", type=str, required=True, metavar="FILE",help="File with target parameters,\ne.g. TargetParams.input")
  p1.add_argument("-d", type=str, required=True, metavar="DIR",help="Directory containing the evolution trajectory file(s)")
  p.add_argument("--tmin",type=float, metavar="FLOAT",help="""Fit points with t>tmin. Default determined by FindTmin estimation scheme.""")
  p.add_argument("--tmax",type=float, metavar="FLOAT",help="""Fit points with t<tmax. Defaults to min(t.back(), 5pi/Omega0+tmin).""")

  args = p.parse_args()  #parse the options
  summary = open("summary.txt","w")
  main(args)
  summary.close()
  System("cat summary.txt")




# def perform_fits_with_u(t,dOmegadt,idxFit,tFit,dOmegadtFit,tmin,tmax,theOmega0,theD0,theadot0,thee0,opts,Source):

#   dOmegaFitoft = InterpolatedUnivariateSpline(tFit,dOmegadtFit,k=5) # the data here is a function of t. we'll use this function to convert it into a function of u, the anomaly
#   ddOmegaFitoft = dOmegaFitoft.derivative()


#   F1cos1 = lambda p,u: p[1]*(p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-11/8) - p[2]*sin(u)*(1-p[5]*cos(u))**(-4) - dOmegaFitoft((u-p[4]-p[5]*sin(u))/p[3])
#   jac = [ lambda p,u: (-11/8)*p[1]*(p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-19/8),
#           lambda p,u:              (p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-11/8),
#           lambda p,u:       -sin(u)*(1-p[5]*cos(u))**(-4),
#           lambda p,u: (-11/8)*p[1]/p[3]/p[3]*(u-p[4]-p[5]*sin(u))*(p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-19/8)+(u-p[4]-p[5]*sin(u))/p[3]*ddOmegaFitoft((u-p[4]-p[5]*sin(u))/p[3])/p[3],
#           lambda p,u: (-11/8)*p[1]/p[3]*(p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-19/8)+ddOmegaFitoft((u-p[4]-p[5]*sin(u))/p[3])/p[3],
#           lambda p,u: -4*p[2]*cos(u)*sin(u)*(1-p[5]*cos(u))**(-5) - 11/8*p[1]/p[3]*sin(u)*(p[0]-(u-p[4]-p[5]*sin(u))/p[3])**(-19/8) + sin(u)*ddOmegaFitoft((u-p[4]-p[5]*sin(u))/p[3])/p[3]  ]
                
#     ptemp, rmstemp, F1cos1_status = fitwithu(F1cos1, p0, pBounds, jac, "F1cos1")


# # fit the data (t,y) to the model F[p,t], by least-squares fitting the params p.
# def fitwithu(F, p0, bounds, jac, name):
#     # F   -- function: F(p,u) taking parameters p and u, the eccentric anomaly
#     # p0  -- starting values of parameters
#     errfunc = lambda p,t,y: sqrt(mean((F(p,t)-y)**2))
#     jacfunc = lambda p,t,y: array([ 1/errfunc(p,t,y) * mean((F(p,t)-y)*dFdp(p,t)) for dFdp in jac ])

#     return res.x, res.fun, res.success




  
