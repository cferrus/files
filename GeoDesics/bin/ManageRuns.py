from __future__ import print_function

import os
import re
from Utils import error, warning
from functools import cmp_to_key

def GetTerminationReason(dir):
  """Given Run directory dir, find the termination reason
for the SpEC simulation.  This function parses the fle
TerminationReason.txt, if not present, returns None
"""
  fname=os.path.join(dir,"TerminationReason.txt")
  if not os.path.isfile(fname):
    return None

  regex=re.compile("Termination condition (.+)")
  with open (fname, "r") as myfile:
    for l in myfile.readlines():
      m=regex.match(l)
      if m:
        TermReason=m.group(1)
#  print "Termination Reason = {}".format(TermReason)
  return TermReason
      

def ListLevs(dir, abort_on_empty=True):
  """In directory dir, look for subdirectories with names 
consistent with SpEC BBH simulations.  Return a list of
all the different Lev's present.  The returned list contains only
the digit(s) of Lev.

Also search over directories with names Condensed_LevN_{XX,XXX}
"""
  # catch level N as named group 'Lev', account for negative N
  regex = re.compile("^(Condensed_|)Lev(?P<Lev>-?[0-9]+)_[A-Z]{2,3}$")

  # apply regex to each directory inside dir
  matches = [regex.match(name) for name in os.listdir(dir)
             if os.path.isdir(os.path.join(dir,name))]

  # extract the integer lev for those directories that matched
  Levs = [int(m.group('Lev')) for m in matches if m]

  if not Levs and abort_on_empty:
    print("Levs={}".format(Levs))
    error("Did not find any directories matching "+regex.pattern)

  # uniquify, sort and return
  return sorted(list(set(Levs)))


################################################################
def SegmentListNoRingdown(dir, Lev, abort_on_empty):
  """Given directory 'dir' and integer level 'Lev',
determine all subdirectories that combine to the
evolution on Lev'Lev'. Does NOT descend into LevN_Ringdown/.

Search for LevN_{XX,XXX} (takes priority) and Condensed_LevN_{XX,XXX}.
"""
  if not os.path.isdir(dir):
    error("path {} is not a directory".format(dir))

  # Algorithm:  Find all matching subdirectory names, and
  # put them into a dictionary indexed by the final letters
  # 'AA', 'AB', 'ACD' etc.  Then sort keys and return entries
  segments = {}
  dirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir,d))]
  LevRE = re.compile("^(Condensed_|)Lev"+str(Lev)+"_([A-Z]{2,3})$")
  # Sort dirs so that "Condensed_*" get processed first
  for name in sorted(dirs):
    m = LevRE.match(name)
    if m:
      key = m.group(2)
      if key in segments:
        warning("Segment {} found twice ({}, {}), using the latter"
                .format(key, segments[key], name))
      segments[key] = name
  if abort_on_empty and not segments:
    error("Did not find any directories matching "+LevRE.pattern)

  # sort by letters 'AA', with presedence on number of chars
  def cmp(x,y):
    if len(x)<len(y): return -1
    if len(x)>len(y): return +1
    if x<y: return -1
    if x>y: return +1
    return 0
  sorted_segments = [segments[k] for k in sorted(segments.keys(),
                                                 key=cmp_to_key(cmp))]
  return sorted_segments


################################################################
def SegmentList(dir, Lev, abort_on_empty=True):
  """Given directory 'dir' and integer level 'Lev',
determine all subdirectories that combine to the
evolution on Lev'Lev'. Does descend into LevN_Ringdown/.
"""
  segs = SegmentListNoRingdown(dir, Lev, False)
  rd = "Lev{}_Ringdown".format(Lev)
  rd_dir = os.path.join(dir,rd)
  if os.path.isdir(rd_dir):
    # get ringdown-segments, and prepend LevN_Ringdown
    ringdown = [os.path.join(rd, p)
                for p in SegmentListNoRingdown(rd_dir, Lev, False)]
    segs.extend(ringdown)
    if abort_on_empty:
      if len(segs)==0:
        error("SegmentList.  Did not find any directories matching Lev='{}'".format(Lev))
  return segs
