#!/usr/bin/env python
from __future__ import print_function
from MakeSubmit import suggested_cpn, _call_machines
from Utils import call_perl
from BHNS_script_functions import get_value_from_string, get_EOS_string, read_target_params_input
import re,sys
import subprocess
from subprocess import call
import os

def main():
  (m,D0,y,q,Mass,spin_bhx,spin_bhy,spin_bhz,spin_nsx,spin_nsy,spin_nsz,a0,e,d,EOS, \
jobname,IDDIR,EOSNAME,EOSDIR,useLeakage,useM1,useMC,NuLibFile,useYe, \
HighPhaseAccuracyOpts,LEV) =  read_target_params_input("./BHNSTargetParams.input")

#moving into the Ev directory:

  os.chdir(os.getcwd() + "/Ev")

#Roughly estimating CPN:

  AllowedCPN = [int(y) for y in _call_machines('GetAllowedCPN()', want_array=True)]
  NPROCS = (int(float(AllowedCPN[0]))+(int(float(AllowedCPN[0])*int(LEV))))
  GrProcs = int(0.4*NPROCS)

#Replacing default parameters in Ev_Params
#with user-chosen values from TargetParams:

  replacements = {'$JOBNAME = "BhNs-M1L2"':'$JOBNAME = '+jobname,
                  '$IDDIR = undef':'$IDDIR = "' +IDDIR+'"',
                  '$HighPhaseAccuracyOpts = undef':'$HighPhaseAccuracyOpts = '+HighPhaseAccuracyOpts,
                  '$useLeakage = 0':'$useLeakage = ' + useLeakage,
                  '$useM1 = 0':'$useM1 = '+useM1,
                  '$useMC = 0':'$useMC = '+useMC,
                  '$useYe = undef':'$useYe = '+useYe,
                  '$EOSDIR = undef':'$EOSDIR = ' +EOSDIR,
                  '$EOSNAME = "HempDD2"':'$EOSNAME = '+EOSNAME,
                  'NuLibFile = undef':'NuLibFile = '+NuLibFile}

  f = open('Ev_Params.input','rt')
  filedata = f.read()
  for src, target in replacements.iteritems():
    filedata = filedata.replace(src,target)
  f = open('Ev_Params.input','w')
  f.write(filedata)
  f.close()

# Same for DoMultipleRuns.input

  replacements = {'$MinLev = undef':'$MinLev = '+LEV,
                  '$MaxLev = undef':'$MaxLev = '+LEV,
                  '$GrProcs = 60+10*$LEV':'$GrProcs = ' + str(GrProcs),
                  '$NPROCS = 64 + 32*$LEV':'$NPROCS = ' + str(NPROCS)}

  f = open('DoMultipleRuns.input','rt')
  filedata = f.read()
  for src, target in replacements.iteritems():
    filedata = filedata.replace(src,target)
  f = open('DoMultipleRuns.input','w')
  f.write(filedata)
  f.close()

if __name__ == "__main__":
  main()

