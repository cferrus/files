#!/usr/bin/env python
from __future__ import division, print_function

import sys

if sys.version_info[0] < 3:
  from ConfigParser import SafeConfigParser
elif sys.version_info[1] < 2:
  from configparser import SafeConfigParser
else:
  from configparser import ConfigParser as SafeConfigParser
import argparse
import os
import re
from math import ceil
from collections import OrderedDict
from Utils import call_perl, error, warning

def ConfigFile():
  """Returns the config file"""
  return 'MakeSubmit.input'
def SubmitFile():
  """Returns the submission file"""
  return 'Submit.sh'

# For case-sensitive config file
class CaseConfigParser(SafeConfigParser):
  optionxform = str
  def __init__(self, defaults=None):
    SafeConfigParser.__init__(self, defaults)

def _read_cfg_file(parser, File=None):
  if not File: File = ConfigFile()
  if not os.path.isfile(File):
    error("{} must exist in {}".format(File, os.getcwd()))
  parser.read(File)

def Query(field, value_type):
  """Return the value of the specified field in MakeSubmit.input"""
  # DEPRECATED_CODE: default for ForceCores should be removed on 2015/06/15.
  # DEPRECATED_CODE: default for Preempt should be removed on 2017/09/01.
  # The default is needed for old runs restarting with a new bin directory.
  cfg = CaseConfigParser({'ForceCores': 'True', 'Preempt': '0'})
  _read_cfg_file(cfg)
  if value_type is bool:
    return cfg.getboolean('main', field)
  else:
    value = cfg.get('main', field)
    return value_type(value)

def QueryAll():
  return {k: Query(k,v) for k,v in allowed_types.items()}

def Update(field_update_dict, do_read, force=False):
  """Updates any number of fields in MakeSubmit.input"""
  cfg = CaseConfigParser()
  if do_read:
    _read_cfg_file(cfg)
    # Don't allow CoresPerNode to be updated if ForceCoresPerNode=True
    if 'CoresPerNode' in field_update_dict.keys() and \
       Query('ForceCoresPerNode', bool):
      if not force:
        warning("Requested change in CoresPerNode, but ForceCoresPerNode=True."
                "\nWill not change CoresPerNode.")
        del field_update_dict['CoresPerNode']
      else:
        warning("Requested change in CoresPerNode will be made, even though"
                "\nForceCoresPerNode=True, because update was passed '-f'.")
    # Don't allow Cores to be updated if ForceCores=True
    if 'Cores' in field_update_dict.keys() and Query('ForceCores', bool):
      if not force:
        warning("Requested change in Cores, but ForceCores=True.\n"
                "Will not change Cores.")
        del field_update_dict['Cores']
      else:
        warning("Requested change in Cores will be made, even though"
                "\nForceCores=True, because update was passed '-f'.")
  else:
    cfg.add_section('main')
  for k,v in field_update_dict.items():
    cfg.set('main', k, str(v))
  with open(ConfigFile(), 'w') as f:
    f.write('# This input file is used to write {}.\n'.format(SubmitFile())+
            '# If you edit this file by hand, run `MakeSubmit.py update`.\n')
    cfg.write(f)

  Convert()

def Create(args):
  """Create new submission scripts"""

  # Get the cores per node
  # cpn is consistency checked with args.cores in Convert()
  cpn = args.cpn if args.cpn else suggested_cpn(args.cores)

  # Get the script from file
  if args.f:
    with open(args.f) as f:
      script = f.read()
  elif args.g:
    script = args.g
  else:
    script = 'EvolutionWrapper -v -a="__EmailAddresses__" ' \
             '-f="__TerminationInfoFile__"'

  # Ensure jobname does not start with a number (see ticket #882)
  jobname = args.J
  if re.match("^\d", args.J):
    warning("Jobname '{}' starts with a number. Prepending '_'.".format(args.J))
    jobname = "_{}".format(args.J)

  fields = OrderedDict({
    'Account': args.A,
    'Jobname': jobname,
    'Hours': args.H,
    'Queue': args.q,
    'Cores': args.cores,
    'ForceCores': args.force_cores,
    'CoresPerNode': cpn,
    'ForceCoresPerNode': args.force_cpn,
    'AllowUnusedCores': args.i,
    'Preempt': args.Preempt,
    'Script': script,
  })

  # Verify keys are only allowed ones
  assert set(fields.keys()) == set(allowed_types.keys()), \
    "Inconsistency in config fields"

  # Write the files
  Update(fields, do_read=False)

def Convert():
  """Converts MakeSubmit.input to Submit.sh"""
  host = _call_machines('GetHost()')
  print("Writing submission script for host '{}'".format(host))

  # Since conf fields can be changed manually, we must sanity check them
  fields  = QueryAll()
  cores   = fields['Cores']
  cpn     = fields['CoresPerNode']
  hours   = fields['Hours']
  preempt = fields['Preempt']

  if _call_machines('GetSharedNodeQueue()') == fields['Queue']:
    shared_node_queue=True
    if cores>=cpn:
      error("Can't run on {cores} cores in a shared node queue because it \n"
            "would take >= 1 node. host {host} has {cpn} cores per node!"\
            .format(**locals()))  
  else:
    shared_node_queue=False

  # Get the submission file template
  header = _call_machines("GetBatchHeader({},{}"
                          ")".format(cores if shared_node_queue else cpn,
                                     preempt))
  header += """
# DO NOT MANUALLY EDIT THIS FILE! See `MakeSubmit.py update -h`.
# This is for submitting a batch job on '{host}'.
umask 0022
. bin/this_machine.env || echo 'Not using env file.'
set -x
export PATH=$(pwd -P)/bin:$PATH\n\n""".format(**locals())

  # Derived quantities
  nodes  = int(ceil(cores/cpn))
  node_cores = int(nodes*cpn)
  unused_cores = node_cores - cores

  # Make sure CoresPerNode is compatable
  if not fields['AllowUnusedCores']:
    if host == "local" and nodes > 1:
      error("Can't run on {cores} cores, host {host} only has {cpn} cores!\n"
            "Either use {cpn} (or fewer) cores or use the -i option (in\n"
            "create mode) or the --AllowUnusedCores option (in update mode)."\
            .format(**locals()))
    if host != "local" and unused_cores > 0:
      error("You have requested {cores} cores. However, host {host} has\n"
          "{cpn} cores per node, so you cannot fill an integer number"
          " of nodes\nwith {cores} cores.  Either change the number of"
          " cores to a multiple\nof {cpn}, or use the -i (in create mode)"
          " or --AllowUnusedCores (in update mode) option,\n"
          "which will ignore this error message and allocate\n"
          "{nodes} nodes ({node_cores} cores), using only {cores} of them, so "
          "{unused_cores} cores will be wasted.".format(**locals()))

  # Make sure Hours is compatible
  max_hours = max_walltime()
  if hours > max_hours:
    error("You have requested {hours} hours of runtime,\n"
          "but the host '{host}' supports {max_hours} hour jobs at most."\
          .format(**locals()))

  # Convert hours to integer hours and minutes
  int_hours = int(hours)
  int_minutes = int((hours-int_hours)*60)

  int_inseconds = int(hours*3600.)

  def _replace_or_remove(text, old, new):
    """Replace 'old' with 'new', or delete lines with 'old' if 'new' is empty"""
    if new:
      return text.replace(old, new)
    else:
      new_lines = []
      for line in text.split('\n'):
        if re.search(old, line): continue
        new_lines.append(line)
      return '\n'.join(new_lines)

  # Replace any placeholders in the header
  header += fields['Script'] + '\n'
  header = header.replace("__hours__", str(int_hours))
  header = header.replace("__minutes__", str(int_minutes))
  header = header.replace("__maxwalltimeinseconds__", str(int_inseconds))
  header = header.replace("__nodes__", str(nodes))
  header = header.replace("__cores__",
                          str(cores if shared_node_queue else node_cores))
  header = header.replace("__SpEC__", fields['Jobname'])
  header = _replace_or_remove(header, "__Account__", fields['Account'])
  header = _replace_or_remove(header, "__Queue__", fields['Queue'])

  # Write the submission file
  with open(SubmitFile(), 'w') as f:
    f.write(header)

# Special types for argparse
def positive_int(x):
  if not int(x)>0 or not int(x)==float(x): raise
  return int(x)
def positive_float(x):
  if not float(x)>0: raise
  return float(x)
def str2bool(x):
  if x == "True": return True
  if x == "False": return False
  print("str2bool values must be either 'True' or 'False'")
  raise

def _call_machines(func, want_array=False):
  return call_perl('Machines', 'new Machines()->{}'.format(func), want_array)

def allowed_cpn():
  return [int(x) for x in _call_machines('GetAllowedCPN()', want_array=True)]
def max_walltime():
  return float(_call_machines('GetMaxRuntimeHours()'))
def default_queue():
  return _call_machines('GetDefaultQueue()')
def default_account():
  return _call_machines('GetDefaultAccount()')

def suggested_cpn(cores):
  cmd = "GetNodeSetup(Nmin=>{}, Nmax=>{}, N=>{})".format(cores,cores,cores)
  N, cpn, _ = call_perl('Machines', cmd, want_array=True)
  assert int(N) == cores, "Expected N={} to be {}".format(N, cores)
  return int(cpn)

#-----------------------------------------------------------------------------

if __name__ == "__main__":

  p = argparse.ArgumentParser()
  subparsers = p.add_subparsers(
    dest = 'subparser_name',
    help = 'Mode to run. Use `MakeSubmit.py <mode> -h` for more info.')

  # For all modes
  p.add_argument('-d', metavar='DIR', default='.',
    help = "Directory to run the desired mode in. (default: %(default)s)")

  # For 'create'
  p_create = subparsers.add_parser('create')
  pr = p_create.add_argument_group(title="required arguments")
  pr.add_argument('--cores', metavar='INT', type=positive_int, required=True,
    help = "Total number of cores. Use --force-cores to prevent this from"\
           "being changed in future segments.")
  pr.add_argument('-J', metavar='JOBNAME', required=True,
    help = "Jobname.")
  p_create_grp = p_create.add_mutually_exclusive_group(required=False)
  p_create_grp.add_argument('-f', metavar='FILE',
    help = "File containing script to append to the batch header. "\
           "Defaults to the standard EvolutionWrapper stub.")
  p_create_grp.add_argument('-g',
    help = "Literal text to append to the batch header. "\
           "Defaults to the standard EvolutionWrapper stub.")
  p_create.add_argument('-H', metavar='HOURS', type=positive_float,
    default = max_walltime(),
    help = "Walltime in hours (default from Machines.pm: %(default)s).")
  p_create.add_argument('-A', metavar='ACCOUNT',
    default = default_account(),
    help = "Account (default from Machines.pm: %(default)s).")
  p_create.add_argument('-q', metavar='QUEUE',
    default = default_queue(),
    help = "Queue (default from Machines.pm: %(default)s).")
  p_create.add_argument('--cpn', type=positive_int,
    choices = allowed_cpn(),
    help = "Cores per node. Use --force-cpn to prevent this from being "\
           "auto-overwritten in future segments (on variable CPN hosts). "
           "Default from Machines.pm. (choices: %(choices)s)")
  p_create.add_argument('--force-cpn', action='store_true',
    help = "Do not allow the CoresPerNode value to change.")
  p_create.add_argument('--force-cores', action='store_true',
    help = "Do not allow the Cores value to change.")
  p_create.add_argument('-i', action='store_true',
    help = "Allow --cores to not be a multiple of --cpn. However, jobs "\
           "will still request an integer number of nodes.")
  p_create.add_argument('--Preempt', type=int, default=0, metavar="INT",
    help = "Should this job be a preemptee (-1), preemptor (1) or normal (0)? "\
           "Will only have an effect if preemption is set up on your machine. "
           "(default: %(default)s)")

  # For 'update'
  allowed_types = {
    'Account': str,
    'Jobname': str,
    'Hours': positive_float,
    'Queue': str,
    'Cores': positive_int,
    'ForceCores': bool,
    'CoresPerNode': positive_int,
    'ForceCoresPerNode': bool,
    'AllowUnusedCores': bool,
    'Preempt': int,
    'Script': str,
  }
  p_update = subparsers.add_parser('update')
  pf = p_update.add_argument_group(title="field arguments",
    description = "These options correspond directly to the "\
                  "fields in {}.".format(ConfigFile()))
  for k,v in allowed_types.items():
    if v is bool: v = str2bool  # modify bool types here (see ticket #924)
    pf.add_argument('--{}'.format(k), type=v, metavar=v.__name__)
  p_update.add_argument('-f', action='store_true',
    help = "Update the fields regardless of any directives in the file "\
           "(such as ForceCores or ForceCoresPerNode)")

  # For 'query'
  p_query = subparsers.add_parser('query')
  p_query.add_argument('field', type=str, metavar='FIELD',
    choices = allowed_types.keys(),
    help = "Field to query from {} (choices: %(choices)s)".format(ConfigFile()))

  args = p.parse_args()

  os.chdir(args.d)

  if args.subparser_name == 'create':
    print(vars(args))
    Create(args)
  elif args.subparser_name == 'update':
    print(vars(args))
    field_update_dict = {k:v for k,v in vars(args).items() \
                         if k in allowed_types.keys() and v is not None}
    Update(field_update_dict, do_read=True, force=args.f)
  elif args.subparser_name == 'query':
    # This call must not print any diagnostics to stdout!
    print(Query(args.field, allowed_types[args.field]))
