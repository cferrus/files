#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import re,sys
import os
import subprocess
from subprocess import call
import os

def get_value_from_string(text,key,allow_none=False):
    match = re.search("^"+key+"\s*=\s*([^;]*);",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def get_EOS_string(text,key,allow_none=False):
    match = re.search("^"+key+"\s*=\s*\"(.+?)\";",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def get_value(text,key,allow_none=False):
    match = re.search("^"+key+"\s*=\s*([0-9a-f\.\-]*)",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def read_target_params_input(filename):
  with open(filename,'r') as f:
    text   = f.read()
    y      = get_value_from_string(text,"\$BHMass")
    q      = get_value_from_string(text,"\$ADMmassNS")
    EOS = get_EOS_string(text,"\$EOS")
    m      = float(y)/float(q) #Mass ratio
    Mass   = float(y)+float(q) #Mass of system
    spin_bhx = get_value_from_string(text,"\$BHSpinx")
    spin_bhy = get_value_from_string(text,"\$BHSpiny")
    spin_bhz = get_value_from_string(text,"\$BHSpinz")
    spin_nsx = get_value_from_string(text,"\$NSSpinx")
    spin_nsy = get_value_from_string(text,"\$NSSpiny")
    spin_nsz = get_value_from_string(text,"\$NSSpinz")
    aT     = get_value_from_string(text,"\$SemiMajorAxis")
    eT     = get_value_from_string(text,"\$Eccentricity")
    dT     = get_value_from_string(text,"\$AnomalyAngle")
    D0     = (2.0*float(aT))/float(Mass)
    rc = call("./bin/dmrgen_lessfluff -e '{}' -d {} -q m_adm {}>output.dat".format(EOS,aT,q), shell=True)
    path = os.getcwd()
    ac = call("python ./bin/ZeroEccParamsFromPN --q {} --chiA {},{},{} --chiB {},{},{} --D0 {}>outputZeroEccParams.dat".format(float(m),float(spin_bhx),float(spin_bhy),float(spin_bhz),float(spin_nsx),float(spin_nsy),float(spin_nsz),float(D0)),shell=True)
    jobname = get_value_from_string(text,"\$JOBNAME")
    path_parent = os.path.dirname(os.getcwd())
    IDDIR = os.getcwd() + "/EvID"
    EOSNAME = get_value_from_string(text,"\$EOSNAME")
    EOSDIR = get_value_from_string(text,"\$EOSDIR")
    useLeakage = get_value_from_string(text,"\$useLeakage")
    useM1 = get_value_from_string(text,"\$useM1")
    useMC = get_value_from_string(text,"\$useMC")
    NuLibFile = get_value_from_string(text,"\$NuLibFile")
    useYe = get_value_from_string(text,"\$useYe")
    HighPhaseAccuracyOpts = get_value_from_string(text,"\$HighPhaseAccuracyOpts")
    LEV = get_value_from_string(text,"\$LEV")
    return(m,D0,y,float(q),float(Mass),spin_bhx,spin_bhy,spin_bhz,spin_nsx,spin_nsy,spin_nsz,float(aT),float(eT),float(dT),EOS,jobname,IDDIR,EOSNAME,EOSDIR,useLeakage,useM1,useMC,NuLibFile,useYe,HighPhaseAccuracyOpts,LEV)
