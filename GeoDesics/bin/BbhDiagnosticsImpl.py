#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import os, re
import numpy as np
import collections
from Utils import ReadH5,error,norm

"""
SpEC Python module 'BbhDiagnosticsImpl'

Contains functions for consistency checking of standardized SpEC
waveform catalogs.
"""


def ListToString(lst):
    """
    if lst contains <=4 elements, just return a string of them.
    if lst contains >=5 elements, print the first 3, "...", and the last
    """

    comma=""
    out="["
    Nfront=3  # keep three elements in front
    for l in lst[0:Nfront]: # print at most three elements (*)
        out+=comma+l
        comma=", "
    if(len(lst)>Nfront+2): # elipsis
        out+=comma+"<<"+str(len(lst)-Nfront-1)+">>"
        out+=comma+lst[-1]
    else:
        for l in lst[Nfront:]: # print rest
            out+=comma+l
    out+="]"
    return out

################################################################

def FindGaps(tsteps):
    """
Given a 1-d ndarray (i.e. the time-column of a dat-set),
compute the following tuple  

  (tmin, tmax, Nsteps, MinDeltaT, Ngaps, t_1stgap, Length_1stgap, t_2nd_gap, Length_2ndgap, t_3rd_gap, Length_3rdgap)

and return it.  The definition of 'gap' is that the time-step increases by more
than a factor of 2 from last time-step.
"""

    tmin=min(tsteps)
    tmax=max(tsteps)
    Nsteps=len(tsteps)
    Ngaps=0
    t_gaps=[]

    DeltaT = np.diff(tsteps)
    MinDeltaT = min(DeltaT)

    RapidIncrease = DeltaT[1:] > 2.*DeltaT[:-1]

    gaps=[idx for (idx,gap) in enumerate(RapidIncrease) if gap]
    gap_list=[]
    if(gaps):
        #print "gaps found at idx=", gaps
        for idx in gaps:
            #print "gap {}: DeltaT={}, tsteps={}".format(idx,DeltaT[idx-2:idx+2],tsteps[idx:idx+4])
            gap_list += [tsteps[idx+1], tsteps[idx+2]-tsteps[idx+1]]
    gap_list+=[-1, -1, -1, -1, -1, -1]
    
    return [tmin, tmax, Nsteps, MinDeltaT, len(gaps),]+ gap_list[:7]
        

################################################################

def Compute_OrbitalFrequency(xA, xB, N, method="Fit", NSamples=None):
    """Given numpy arrays for the black hole locations xA, xB,
    perform fits to 2N+1 data-points around each point, and from
    the fit compute the instantaneous orbital frequency,
            Omega = r\times \dot{r} / r^2
    return t,Omega
    """

    def FitData(data, N, NSamples=None):
      """given a numpy array data with time as first column, perform fits 
        covering N points before and after each data-point for each column.
        return the fitted values and their first time-derivatives as a 
        numpy arrays, with first column time"""
      # collect output data, first column being time
      last_idx=len(data)-N-1

      if NSamples==None:
        step=1
      else:
        step=max(int(last_idx/NSamples),1)

      # The output times
      t_final = data[N:last_idx:step,0]

      x_tmp = []
      v_tmp = []
      for idx in range(N, last_idx, step):
        # Time at which we want the result
        T = data[idx,0]

        x = data[idx-N:idx+N+1,0]-T # Shift back to t=0
        y = data[idx-N:idx+N+1,1:]  # Fit all the columns at once!
        p0 = np.polyfit(x,y,2)
        x_tmp.append(p0[2]) # p0[2] is the constant part of fit
        v_tmp.append(p0[1]) # p0[1] is the linear part

      return np.column_stack((t_final,x_tmp)), np.column_stack((t_final,v_tmp))

    def FilterData(data,N):
      assert False, "The FilterData function is poorly tested"
      from scipy.ndimage import gaussian_filter1d
      mode = 'reflect'
      return gaussian_filter1d(data,N,mode=mode),\
             gaussian_filter1d(data,N,mode=mode,order=1)

    def SmoothData(data,N):
      from Utils import SmoothData
      t   = data[:,0]
      dat = data[:,1:]
      dt_interp = None
      return SmoothData(t,dat,N,DT=dt_interp), \
             SmoothData(t,dat,N,DT=dt_interp,Deriv=1)

    def SplineData(data):
      """Interpolating spline (no smoothing)"""
      from scipy.interpolate import splrep,splev
      t   = data[:,0]
      spline_x = splrep(t, data[:,1])
      spline_y = splrep(t, data[:,2])
      spline_z = splrep(t, data[:,3])
      dx = splev(t, spline_x, der=1)
      dy = splev(t, spline_y, der=1)
      dz = splev(t, spline_z, der=1)
      v = np.vstack((t,dx,dy,dz)).T
      return data, v

    if NSamples is not None and method!='Fit':
        error("NSamples only works with 'Fit'")

    if method=="Fit":
      data = np.column_stack((xA,xB[:,1:4]))
      xs_fit,vs = FitData(data,N,NSamples=NSamples)
      xA_fit = xs_fit[:,0:4]
      xB_fit = xs_fit[:,[0,4,5,6]]
      vA = vs[:,0:4]
      vB = vs[:,[0,4,5,6]]
    elif method=="Filter":
      xA_fit,vA=FilterData(xA,N)
      xB_fit,vB=FilterData(xB,N)
    elif method=="Smooth":
      xA_fit,vA=SmoothData(xA,N)
      xB_fit,vB=SmoothData(xB,N)
    elif method=="Spline":
      xA_fit,vA=SplineData(xA)
      xB_fit,vB=SplineData(xB)
    else:
      error("Don't know method '{}'".format(method))

    # Compute Orbital frequency (r x dr/dt)/r^2
    t  = xA_fit[:,0]
    dr = xA_fit[:,1:] - xB_fit[:,1:]
    dr2 = norm(dr,axis=1)**2   #slightly inefficient
    dv = vA[:,1:] - vB[:,1:]
    Omega = np.cross(dr,dv)/dr2[:,np.newaxis]
    return t,Omega



################################################################
def VerifyTimeStepsInGroup(h5group, exclude):
    """
VerifyTimeStepsInGroup
Given an h5py object representing an opened h5-group, process all 
"*.dat" groups in this h5group.  For each such group, consider the
first column [:,0]:
(a) ensure thay all these time-columns are identical, and keep one
    in time_column.  If time-columns are not identical, record an
    error in err_str. 

(b) examine all *.dat groups and for each one, compute a tuple
  (group-name, tmin, tmax, Nsteps, MinDeltaT, Ngaps, t_1stgap, Length_1stgap, t_2nd_gap, Length_2ndgap, t_3rd_gap, Length_3rdgap)


Returns (time_column, err_str, L)

'exclude' gives a list of dat-groups to *exclude* from this comparison
"""
    match=re.compile("(.*)\.dat$")

    #----------------
    # first pass: compute gaps
    #----------------
    err_str=""
    TstepInfo=[]
    last_tsteps=None
    for setname in h5group:
        M=match.match(setname)
        if(setname in exclude or not M): continue
        dataset=h5group[setname]
        tsteps=dataset[:,0] 

        if(last_tsteps is not None 
           and np.array_equal(last_tsteps,tsteps)
           ):
            # identical data to before, no need to call FillGaps
            info=TstepInfo[-1][:-1];
        else:
            # don't have analyzed this tsteps yet
            info=FindGaps(tsteps)

            # issue error only if tstep different from before
            if(last_tsteps is not None):
                err_str +="  {} differs from preceeding groups\n".format(setname)
            last_tsteps=tsteps
        info.append(h5group.name+"/"+setname)
        
        TstepInfo.append(info)
            
    #----------------
    # second pass: check lengths only, to generate simpler error message
    #----------------

    # dictionary from length -> [dsets, ...]
    LengthsD=collections.defaultdict(list) 
    for info in TstepInfo:
        LengthsD[info[2]].append(info[-1])

    if len(LengthsD)>1:
        # DATA INCONSISTENT
        err_str="DATA SIZE INCONSISTENT IN {:s}\n".format(h5group.name)
        for it in LengthsD:
            err_str+="   size {:5d} {}\n".format(it,ListToString(LengthsD[it]))
        last_tsteps=np.zeros(1)

    return (last_tsteps, err_str,TstepInfo)
        
        

################################################################
def VerifyTimeStepsInFile(h5file,regexp=".*",
                          exclude=[]):
    """
Given a h5py.file h5file, iterate through top-level directories of
this file.  Within each directory with name matching 'regexp' call
VarifyTimeStepsInGroup.

Return a dictionary:
  dir-name -> tsteps
where tsteps is the numpy array containing all TimeSteps in said directory
"""

    match=re.compile(regexp)
    D={}
    err_str=""
    TstepInfo=[]
    for grp in h5file:
        M=match.match(grp)
        if not M: continue
        #print(grp, M.group(1))
        tmp,err,tstepInfo=VerifyTimeStepsInGroup(h5file[grp], exclude)
        D[grp]=tmp
        err_str+=err
        for i in tstepInfo:
            TstepInfo.append(i+[h5file.filename])
    return D,err_str,TstepInfo


################################################################
def ParseRadii(D):
    """
    Takes a dictionary D = {"R0104.dir" ->  ndarray}
    Extracts the radius from the name 'R0104.dir' -> 104
    Returns a ndarray
    res[:,0] = radius
    res[:,1] = Nsteps at that radius
"""
    res=np.zeros([len(D),2])
    matchR=re.compile("R0*([0-9]+).dir")
    for idx in enumerate(D):
        matched_radius=matchR.match(idx[1])
        if(not matched_radius):
            print("ERROR encountered wrong format for radius")
        res[idx[0],0]=float(matched_radius.group(1))
        res[idx[0],1]=D[idx[1]].shape[0]
    return res

################################################################
def VerifyTimeStepsInDirectory(dir):
    """
Given a path to a directory that contains standardized BBH output, check
a) Presence of Horizons.h5, {rh,rPsi4}_FiniteRadii_CodeUnits.h5
b) consistency of timesteps within each {rh,rPsi}/R????.dir/ 
c) consistency of timesteps accross each {rh,rPsi}/R????.dir/
d) consistency of timesteps within each Horizons.h5/Ah?.dir/

Returns  tsteps, steps_h, steps_psi4, err
  tsteps     -- a list of all distinct tsteps inside rh, rPsi4, Horizons
  steps_h    -- 1/Rextr vs Nsteps at that extraction radius
  steps_psi4 -- 1/Rextr vs Nsteps at that extraction radius
  err -- string listing encountered errors.
  Tstep_info -- a list of tuples (tmin, tmax, Nsteps, MinDeltaT, 
                                  Ngaps, t_1stgap, Length_1stgap, 
                                  t_2nd_gap, Length_2ndgap, 
                                  t_3rd_gap, Length_3rdgap,
                                  .dat-name, filename)
"""

    tsteps=[]
    err=""
    steps_h=None
    steps_psi4=None
    TstepInfo=[]
    #--------------------------------
    # Wave-Extraction files
    #--------------------------------
    h5 = ReadH5(dir+"/rh_FiniteRadii_CodeUnits.h5")
    hD,h_err,h_tstep_info=VerifyTimeStepsInFile(
        h5, regexp=".*dir", exclude=["InitialAdmEnergy.dat"] ) 
    h5.close()
    steps_h=ParseRadii(hD)

    h5 = ReadH5(dir+"/rh_FiniteRadii_CodeUnits.h5")
    psi4D,psi4_err,psi4_tstep_info=VerifyTimeStepsInFile(
        h5, regexp=".*dir", exclude=["InitialAdmEnergy.dat"] ) 
    h5.close()
    steps_psi4=ParseRadii(psi4D)
    
    err=h_err+psi4_err
    # test that both Psi4 and h are at same time-steps
    if h_err=="" and psi4_err=="":
        if sorted(hD.keys())!=sorted(psi4D.keys()):
            err += "RADII MISMATCH BETWEEN rh AND rPsi4\n"
            hset=set(hD.keys())
            psi4set=set(psi4D.keys())
            if len(hset-psi4set) > 0:
                err += "   missing in rPsi4: "+str(hset-psi4set)+"\n"
            if len(psi4set-hset) > 0:
                err += "   missing in rh:    "+str(psi4set-hset)+"\n"
            print(err)


    #--------------------------------
    # Horizons.h5
    #--------------------------------

    # for Horizons.h5, currently employ workaround 
    # for differing lengths
    h5 = ReadH5(dir+"/Horizons.h5")
    Ah1,Ah1_err,Ah1_tstep_info=VerifyTimeStepsInFile(h5, 
                                      exclude=["CoordCenterInertial.dat",
                                               "DimensionfulInertialSpin.dat"])
    Ah2,Ah2_err,Ah2_tstep_info=VerifyTimeStepsInFile(
        h5, 
        exclude=["ArealMass.dat", "ChristodoulouMass.dat", 
                 "DimensionfulInertialSpinMag.dat",
                 "chiInertial.dat", "chiMagInertial.dat"])        

    for ts in list(hD.items())+list(psi4D.items())+list(Ah1.items())+list(Ah2.items()):
        for prev_t in tsteps:
            if np.array_equal(prev_t[0],ts[1]): 
                break
        else:
            tsteps.append( (ts[1], ts[0]) )


    err+=Ah1_err+Ah2_err
    Tstep_info=h_tstep_info+psi4_tstep_info+Ah1_tstep_info+Ah2_tstep_info
    return tsteps,steps_h,steps_psi4,err,Tstep_info

################################################################
def Compute_h_SkyDirection(Hgrp, theta_deg, phi_deg, l_max=None):
    from sYlm import sYlm

    th = theta_deg/180.*np.pi   # convert to radians
    ph = phi_deg/180.*np.pi     # convert to radians

    match = re.compile("Y_l(.*)_m(.*).dat")

    # dictionaries (l,m) -> hlm
    hlm = {}
    t = None
    for dset in Hgrp:
        M=match.match(dset)
        if not M: continue

        l = int(M.group(1))
        m = int(M.group(2))

        # compute complex coefficient
        hlm[(l,m)] = Hgrp[dset][:,1] + 1j * Hgrp[dset][:,2]
        if t is None:
            t = Hgrp[dset][:,0]

    summed_h = 0j*t  # initialize
    for l,m in hlm:
        if (l_max is not None) and (l > l_max): continue
        coef = sYlm(-2,l,m,th,ph)
        summed_h += hlm[(l,m)] * coef

    return t, summed_h
