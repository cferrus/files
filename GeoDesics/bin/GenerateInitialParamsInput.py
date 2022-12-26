#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import re,sys
from numpy import sqrt, size, sin, cos, pi
import numpy as np
from scipy import optimize

def get_value_from_string(text,key,allow_none=False):
    match = re.search("^"+key+"\s*=\s*([^;]*);",text,re.MULTILINE)
    if match:
      return(match.groups()[0])
    elif not allow_none:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))
  
def get_vector_from_string(text,key):
    match = re.search("^"+key+
                      "\s*=\s*\(\s*(\S+)\s*,\s*(\S+)\s*,(\S+)\s*\)\s*;",
                      text,re.MULTILINE)
    if match:
      return(match.groups())
    else:
      sys.exit("Error: Cannot find {} in the text\n{}".format(key,text))

def read_target_params_input(filename):
  with open(filename,'r') as f:
    text   = f.read()
    q      = get_value_from_string(text,"\$MassRatio")
    id_type= get_value_from_string(text,"\$IDType")
    spin_a = get_vector_from_string(text,"\@SpinA")
    spin_b = get_vector_from_string(text,"\@SpinB")
    aT     = get_value_from_string(text,"\$SemiMajorAxis")
    eT     = get_value_from_string(text,"\$Eccentricity")
    dT     = get_value_from_string(text,"\$AnomalyAngle")
    Tref   = get_value_from_string(text,"\$ReferenceTime",allow_none=True)
  return(float(q),spin_a,spin_b,float(aT),float(eT),float(dT),Tref,id_type)

def main():
  (q,spin_a,spin_b,a0,e,d,Tref,id_type) \
    = read_target_params_input("TargetParams.input")

  eta = q/(1+q)/(1+q)

  Amp = lambda a,e: -2*e*sqrt(1-e*e)*a**(-3)+e*(4*(4-eta)+e*e*(5*eta-22))/sqrt(1-e*e)*a**(-4)
  B = lambda a,e: e+e/a*(eta-2)
  Om = lambda a,e: 1/a/sqrt(a)+(eta-9)*0.5/a/a/sqrt(a)
  Ct = lambda a,e: e-e*(8-3*eta)*0.5/a
  Cf = lambda a,e: e+e*eta*0.5/a

  #uoft = lambda a, e, d, t : Om(a,e)*t + d

  Keplereq = lambda u, a, e, d, t : u-Ct(a,e)*sin(u)-Om(a,e)*t-d
  Keplereqprime = lambda u, a, e, d, t: 1-Ct(a,e)*cos(u)
  Keplereqprime2 = lambda u, a, e, d, t: Ct(a,e)*sin(u)
  uoft = lambda a, e, d, t : optimize.newton(Keplereq, Om(a,e)*t+d, fprime=Keplereqprime, args=(a,e,d,t), maxiter=10,fprime2=Keplereqprime2) 

  r = lambda a,e,d,t: a*(1-e*cos(uoft(a,e,d,t)))
  rdot = lambda a,e,d,t:  Om(a,e)*a*e*sin(uoft(a,e,d,t))/(1-Ct(a,e)*cos(uoft(a,e,d,t)))
  Omega = lambda a,e,d,t: (sqrt(1-e*e)*a**(-1.5)-0.5*((3-eta)+e*e*(2*eta-9))*a**(-2.5)/sqrt(1-e*e))*(1-Cf(a,e)*cos(uoft(a,e,d,t)))**(-1)*(1-Ct(a,e)*cos(uoft(a,e,d,t)))**(-1)

  #Tref = 2*pi/Omega(a0,e,d,0) #Throw away a full cycle due to junk radiation
  #print Tref

  #Tref = 2*pi*np.power(a0,3./2.)
  #print Tref

  need_to_write_Tref = False
  if not Tref:
    Tref = 2*pi*np.power(a0,3./2.)*(1+0.5*(9-eta)/a0)
    need_to_write_Tref = True
  else:
    Tref = float(Tref)
  #print Tref
  #print 0.5*(9-eta)/a0

  aT = np.power(-64./5./q*(1-1./q)*4*(1+73./24.*e*e+37./96.*e*e*e*e)*(-Tref)+a0*a0*a0*a0,1./4.) #Peters' formula
  #eT=e*np.power(aT/a0,19./12.)
  #print eT, e

  #print r(a0,e,d,0), rdot(a0,e,d,0)/r(a0,e,d,0), Omega(a0,e,d,0), Tref
  #print r(aT,e,d,0), rdot(aT,e,d,0)/r(aT,e,d,0), Omega(aT,e,d,0), Tref
  #print r(aT,eT,d,0), rdot(aT,eT,d,0)/r(aT,eT,d,0), Omega(aT,eT,d,0), Tref
  Omega0=Omega(aT,e,d,0)
  adot0=rdot(aT,e,d,0)/r(aT,e,d,0)
  D0=r(aT,e,d,0)

  if need_to_write_Tref:
    with open("TargetParams.input") as f:
      text   = f.read()

    text = re.sub("(\$ReferenceTime\s*=\s*);","\g<1>{};".format(Tref),text)
    
    with open("TargetParams.input","w") as f:
      f.write(text)

  # Write params.input
  with open("Params.input","w") as f:
    f.write("# Set the initial data parameters\n")
    f.write("\n")
    f.write("# Orbital parameters\n")
    f.write("$Omega0 = {};\n".format(Omega0))
    f.write("$adot0  = {};\n".format(adot0))
    f.write("$D0     = {};\n".format(D0))
    f.write("\n")
    f.write("# Physical parameters (spins are dimensionless)\n")
    f.write("$MassRatio = {};\n".format(q))
    f.write("@SpinA = ({},{},{});\n".format(*spin_a))
    f.write("@SpinB = ({},{},{});\n".format(*spin_b))
    f.write("\n")
    f.write("# Evolve after initial data completes?\n")
    f.write("$Evolve = 1;\n")
    f.write('# IDType: "SKS", "SHK" or "CFMS".\n')
    f.write('$IDType = {};\n'.format(id_type))

if __name__ == "__main__":
  main()
