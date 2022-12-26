#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import re,sys
import os
from numpy import sqrt, size, sin, cos, pi
import subprocess
from subprocess import call
import numpy as np

def get_value_from_string(text,key,allow_none=False):
    match = re.search("^"+key+"\s*=\s*([^;]*);",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def get_value_string(text,key,allow_none=False):
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

def get_string(text,key,allow_none=False):
    match = re.search(key+"\s*=\s*([0-9\.0-9]*)",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def read_target_params_input(filename):
  with open(filename,'r') as f:
    text   = f.read()
    y      = get_value_from_string(text,"\$BHMass")
    q      = get_value_from_string(text,"\$ADMmassNS")
    EOS = get_value_string(text,"\$EOS")
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
  return(m,D0,y,float(q),float(Mass),spin_bhx,spin_bhy,spin_bhz,spin_nsx,spin_nsy,spin_nsz,float(aT),float(eT),float(dT),EOS)

def read_output_dat(filename):
  with open(filename,'r') as t:
    text1 = t.read()
    rhoc = get_value_from_string(text1,"\$rhocOne")
    massNS = get_value_from_string(text1,"\$mass1")
  return(rhoc,massNS)

def outputZeroEccParams(filename):
  with open(filename,'r') as s:
    text2 = s.read()
    Omega0 = get_value(text2,"Omega0")
    adot = get_value(text2,"adot0")
  return(Omega0,adot)

def main():
  (m,D0,y,q,Mass,spin_bhx,spin_bhy,spin_bhz,spin_nsx,spin_nsy,spin_nsz,a0,e,d,EOS) \
    = read_target_params_input("TargetParams.input")

  (rhoc,massNS) \
    = read_output_dat("output.dat")

  (Omega0,adot)  \
    = outputZeroEccParams("outputZeroEccParams.dat")

  dist = 2.0*float(a0)
  Omega = float(Omega0)/float(Mass)

  # Write params.input
  with open("Params.input","w") as f:
    f.write("""# Type of solve required                                                      
# Irrotational : use quasi-circular formalism and solve for Omega             
#                 Flow in NS assumed to be irrotational.                      
# LowEllipticity : use Omega and dta as given in initial guess                
#                Flow in NS assumed to be irrotational                        
#                WARNING: Eccentricity reduction is supposed to be            
#                          performed at CONSTANT SEPARATION. Make sure         
#                          that you set $d below to the separation between     
#                          the compact objects at the previous iteration       
#                          (this may differ from the separation provided       
#                          in Params.input, as the center of the NS is        
#                          moved by the ID).                                  
#                          Use $ID_d from the perl file generated at the      
#                          END of the ID solve for the previous iteration.    
# Corotational : same as Irrotational, but with corotational flow
# HeadOn       : Omega=0, dta=0                                                                              
# SingleStar   : Does not solve BH part\n""")
    f.write("$Type = \"LowEllipticity\";\n")
    f.write("# Set the initial data parameters\n")
    f.write("\n")
    f.write("# Orbital angular velocity\n")
    f.write("$Omega = {};\n".format(Omega))
    f.write("# Initial radial velocity\n")
    f.write("$dta = {};\n".format(adot))
    f.write("# Initial separation (in units of solar masses)\n")
    f.write("$d = {};\n".format(dist))
    f.write("\n")
    f.write("# BH info\n")
    f.write("$BHMass = {}; #(in units of solar masses)\n".format(y))
    f.write("# BH spin in cartesian coords\n")
    f.write("$BHSpinx = {}*$BHMass;\n".format(spin_bhx))
    f.write("$BHSpiny = {}*$BHMass;\n".format(spin_bhy))
    f.write("$BHSpinz = {}*$BHMass;\n".format(spin_bhz))
    f.write("\n")   
    f.write("# NS Info\n")
    f.write("$ADMmassNS = {}; #(in units of solar masses)\n".format(q))
    f.write('$massNS = {};\n'.format(massNS))
    f.write('$EOS = "{}";\n'.format(EOS))
    f.write('$rhoc = {};\n'.format(rhoc))
    f.write("# NS Spin (units of 1/Msun)\n")
    f.write("$NSSpinx = {};\n".format(spin_nsx))
    f.write("$NSSpiny = {};\n".format(spin_nsy))
    f.write("$NSSpinz = {};\n".format(spin_nsz))
    f.write("\n")
    f.write("# Surface of the star is defined by h = targeth\n")
    f.write("# Use 1 for analytic EoS, and h such that rho(h) is small for tabulated EoS\n")
    f.write("$TargetH=1.0;\n")
    f.write("\n")
    f.write("# Evolve after initial data completes?\n")
    f.write("$Evolve = 1;\n")
    f.write("# To modify variables $CF, $BgPsi, $TwoShell, see DoMultipleRuns.input \n")

if __name__ == "__main__":
  main()
