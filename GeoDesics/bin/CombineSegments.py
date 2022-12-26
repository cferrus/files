#!/usr/bin/env python

"""
Combine segments from a standard SpEC evolution. Must be run in the directory
that contains the segments. Will look for matching files in Run directories
of both the inspiral and ringdown. Files will NOT be joined if the path
contains the string 'Checkpoints'.

e.g.: CombineSegments -o JoinedLev3 -L 3 -e dat h5 -f 'Constraint.*L2'
"""

from __future__ import print_function

import sys, os
import re
import time
import tempfile
sys.path.insert(0,os.path.realpath(__file__+'/../../Python'))
from Utils import System, error, warning, ReadH5
from ManageRuns import SegmentList

from UpdateH5DataVersion import *

#------------------------------------------------------------------------------

def IsCachedH5FileFormat(File):
  """Identify CachedH5Files by the use of specific dataset name extensions"""
  f = ReadH5(File)
  root_dsets = f.keys()
  known_exts = ["dir", "dat", "tdm", "txt", "ver"]
  # Every root dataset name must match at least one known extension
  for dset in root_dsets:
    if not any(re.search("\.{}$".format(ext), dset) for ext in known_exts):
      return False
  else:
    return True

def main():
  # Check that this version supports Argparse
  if sys.version_info <= (2,7):
      error("This script is only supported in Python version 2.7 or higher.")
  import argparse

  p = argparse.ArgumentParser(
    prog="CombineSegments",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=__doc__,
  )

  # Required:
  pr = p.add_argument_group('required arguments')
  pr.add_argument("-o", type=str, metavar="DIR", required=True, nargs=1,
                  help="Output directory for joined segment data. "
                  "Must not exist prior to execution.")
  pr.add_argument("-L", metavar="LEV", required=True, nargs=1,
                  help="Level to join. If e.g. L=3, segments are assumed "
                  "to match 'Lev3_[A-Z]{2,3}'.")

  # Optional:
  p.add_argument("-d", type=str, metavar="DIR", nargs=1, default=["Run"],
                 help="Directory to search under inside Lev* "
                 "(default: %(default)s). If empty, searches all directories "
                 "inside Lev*.")
  p.add_argument("-f", type=str, metavar="FILE", nargs='+', default=[".*"],
                 help="Files to join (default: %(default)s). May be a regex "
                "matching any part of relative file paths. Use single quotes "
                 "to avoid shell expansion. Files will only be joined if they "
                 "satisfy both -f and -e options.")
  p.add_argument("-e", type=str, metavar="EXT", nargs='+',
                 choices=["dat","out","h5","dump"],
                 default=["dat","out","h5","dump"],
                 help="Only join files with these extensions "
                 "(default: %(default)s).")
  p.add_argument("-x", type=str, metavar="FILE", nargs='+', default=[],
                 help="Files to exclude, even if they satisfy the -f and -e "
                 "options.  May be a regex like -f.")
  p.add_argument("-v", type=int, metavar="INT", default=1,
                 help="Verbosity level from 0 to 2 (default: %(default)s).")
  p.add_argument("--no-wipe", action="store_false", dest="Wipe", default=True,
                 help="If specified, do not use the flags that wipe "
                 "non-monotonic time steps. This makes the script faster.")
  p.add_argument("--search-first", action="store_true", dest="SearchFirst",
                 default=False,
                 help="If specified, search only the first segment for files "
                 "to join instead of all. This makes the script faster, but "
                 "may be incomplete if all files are not written in the first "
                 "segment.")
  p.add_argument("--use-existing", action="store_true", dest="UseExisting",
                 default=False,
                 help="If specified, will allow output directory to exist. "
                 "Will skip files that already exist in output directory.")
  p.add_argument("--skip-recent-segments", action="store_true",
                 dest='SkipRecentSegments', default=False,
                 help="If given, ignore segments that were modified within "
                 "the last 3 hours. This avoids errors for ongoing runs, "
                 "where the last segment may be written to while "
                 "CombineSegments runs.")
  p.add_argument("--gnuparallel", action="store_true", default=False,
                 help="If specified, use gnu parallel to execute N/2 "
                 "Join-commands in parallel, where N is the number of cores")
  p.add_argument("--poor_mans_cache", action="store_true", default=False,
                 help="If specified, 'cat > /dev/null' each .h5 file before "
                 "the JoinH5 command.  Harald finds this speeds up JoinH5 by "
                 "an order of magnitude on CITA 128GB memory servers")
  p.add_argument("--cachesize", type=float, default=1e9,
                 help="Size of h5 cache in bytes.(default: %(default)s).")

  opts = p.parse_args()

  # Give input arguments descriptive names
  Exts = opts.e
  DirToSearch = opts.d[0].strip()
  JoinDir = opts.o[0]
  Lev = opts.L[0]

  # Only look inside the Run directory
  CurDir = os.getcwd()
  dirs = SegmentList(CurDir, Lev)
  if opts.v > 0:
    print("CurDir={}".format(CurDir))
    print("dirs={}".format(dirs))
  if DirToSearch: dirs = [ os.path.join(i,DirToSearch) for i in dirs ]

  # Make sure the directories are readable
  no_read = []
  for dirname in dirs:
    if not os.access(dirname, os.R_OK):
      no_read.append(dirname)
  if no_read:
    error("Could not access the following directories:\n{}".format(no_read))

  # skip recent segments (if so asked) This option is intented to
  #   skip the latest segment when joining *ongoing* simulations,
  #   where the files in the latest segment may be written to.  For
  #   convenience, the option is phrased in terms of file-age to allow
  #   to *always* use this option, and get reasonable results, no
  #   matter whether the run is still ongoing or not (in the latter it
  #   will join all files).
  if opts.SkipRecentSegments:
    now = time.time()
    new_dirs = [d for d in dirs if now-os.path.getmtime(d)>3*3600]
    discarded=sorted(list(set(dirs)-set(new_dirs)))
    if len(discarded)>0:
      if opts.v>0:
        print("--skip-recent-segments:  skipping {}".format(discarded))
      dirs = new_dirs

  # Regexp to match files
  FileRE = [re.compile(f) for f in opts.f]
  # Regexp to exclude files
  ExcludeFileRE = [re.compile(f) for f in opts.x]
  # Regexp to match extensions
  ExtRE = {ext: re.compile(".*."+ext+"$") for ext in Exts}

  # Make sure the output directory can be made
  if not opts.UseExisting:
    try:
      os.makedirs(JoinDir)
    except OSError:
      error("Will not overwrite output directory. To continue, remove "+JoinDir
            +"\nor use the '--use-existing' flag to add data to an existing "
            +" output directory.")

  # Walk through the search segment(s) to get all the filenames
  Files = {}  # Using a dict avoids duplictes (necessary if len(SearchDirs)>1)
  SearchDirs = [dirs[0]] if opts.SearchFirst else dirs
  for SearchDir in SearchDirs:
    for froot, fdirs, ffiles in os.walk(SearchDir):
      # NewRoot is the root without the segment in the path
      NewRoot = froot.replace(SearchDir,'',1).lstrip("/")
      # Get all the files that need to be joined
      for File in ffiles:
        NewFile = os.path.join(NewRoot,File)
        # Check restrictions on file path
        # => file must NOT match 'Checkpoint'
        # => file must match at least one FileRE requirement
        # => file must NOT match any ExcludeFileRE requirement
        if re.search("Checkpoint",NewFile): continue
        if not any(RE.search(NewFile) for RE in FileRE): continue
        if any(RE.search(NewFile) for RE in ExcludeFileRE): continue
        # Make sure each file matches at least one ExtRE requirement
        # Save the filename and the associated extension
        for ext,RE in ExtRE.items():
          if RE.match(File):
            Files[NewFile] = ext
            break

  # Flags to wipe non-monotonic time steps (MUST HAVE SPACES!)
  WipeDat = "-w" if opts.Wipe else ""
  WipeH5  = "" if opts.Wipe else "-NoWipe"
  WipeH5old = "-l" if opts.Wipe else ""
  if opts.Wipe and ('out' in Files.values() or 'dump' in Files.values()):
     warning("Wiping non-monotonic time steps for dump and out files not "
             "implemented!")

  # Do the joining for each type of file
  BinDir = os.path.realpath(__file__+'/../../bin')
  Exec = {}
  Exec["dat"] = BinDir+"/JoinDatFiles --take_first_header "+WipeDat
  Exec["out"] = "cat"
  Exec["dump"] = "cat"
  Exec["h5"] = BinDir+"/JoinH5 "+WipeH5+ \
               " --cachesize {} ".format(opts.cachesize)
  Exec["h5old"] = BinDir+"/JoinH5Files "+WipeH5old

  # Join the segments
  AllCommands="" # only needed for --gnuparallel
  for File in sorted(Files):
    ext = Files[File]
    OutFile = os.path.join(JoinDir,File)

    # mkdir output directory, if needed
    OutDir = os.path.dirname(OutFile)
    if not os.path.exists(OutDir):
      os.mkdir(OutDir)

    # Skip output files that already exist
    if os.path.isfile(OutFile):
      warning("Skipping '%s' because it is already joined." %File)
      continue

    # Decide which segments to use in joining this file
    JoinList = []
    for Dir in dirs:
      NewFile = os.path.join(Dir,File)
      # Only add the file if it exists and is not empty
      if os.path.exists(NewFile) and os.path.getsize(NewFile)>0:
        JoinList.append(NewFile)

    try: # check if JoinList is a list of H5 files
      h5py.File(JoinList[0],'r')
      is_HDF5_file = True
    except IOError:
      is_HDF5_file = False
      pass

    if is_HDF5_file and not ensure_same_version(JoinList):
      error("Trying to combine files across different versions. "
            "Please use UpdateH5DataVersion.py to find and update files.")

    if len(JoinList) > 0:
      if ext == "h5":
        if not IsCachedH5FileFormat(JoinList[0]): ext = "h5old"
        cmd = Exec[ext]+" -o "+OutFile+" "+" ".join(JoinList)
      else:
        cmd = Exec[ext]+" "+" ".join(JoinList)+" > "+OutFile
      if opts.poor_mans_cache:
        cmd = "cat "+" ".join(JoinList)+" > /dev/null;" + cmd
    else:
      cmd = "echo 'Not joining "+File+". All files are empty.'"
    if opts.v > 0: print("*** Joining "+File)
    if opts.v > 1: print(cmd)
    if opts.gnuparallel:
      AllCommands=AllCommands+cmd+"\n"
    else:
      System(cmd)

  if opts.gnuparallel:
    if opts.v>0:
      print("opts.gnuparallel={}.  Calling gnu parallel".format(
        opts.gnuparallel))
    f=tempfile.NamedTemporaryFile(delete=False)
    f.write(AllCommands)
    f.close()
    delay=""
    if opts.poor_mans_cache:
      # stagger 'cat' commands. Harald hopes this reduces strain
      # on file system
      delay="--delay 5"
    cmd="parallel --jobs 50% {} :::: {}".format(delay, f.name)
    print("about to execute cmd='{}'".format(cmd))
    System(cmd)
    os.unlink(f.name)

if __name__ == "__main__":
  main()
