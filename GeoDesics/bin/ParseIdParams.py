#!/usr/bin/env python
import os,re
from Utils import error

def ParseIdParams(path,file="ID_Params.perl",array=False):
  """Parses an ID_Params.perl file and returns a dictionary of {key: value} 
pairs. Arguments are:
  path                # path to file to parse
  file=ID_Params.perl # filename of file to parse
  array=False         # if True, output values in array context, not string
                      # context, where 
                      #   string context => 'value1, value 2'
                      #   array context  => ['value1','value2']
"""

  fullpath = os.path.join(path,file)
  try:
    text = open(fullpath, 'r').read()
  except IOError:
    error("ParseIdParams: Could not read file %s!" %fullpath)

  BaseRegexp = r'([^\s]+)\s*= *([^\n]*)'
  Regexp     = "^(?:@|\$)%s;\s*$" %BaseRegexp

  # Grab all (key,value) pairs and put them in a dictionary
  PairList = re.findall(Regexp, text, re.MULTILINE)
  PairDict = dict(PairList)
  for key in PairDict:
    # strip leading/trailing spaces, overall parens, and ALL quotes
    PairDict[key] = PairDict[key].strip(' ()').replace('"','').replace("'",'')
  if array:
    for key in PairDict:
      PairDict[key] = [i.strip() for i in PairDict[key].split(',')]

  return ParseIdParamsDict(PairDict)

# Returns a more useful error message when dict key does not exist
class ParseIdParamsDict(dict):
  def __getitem__(self, key):
    try:
      val = dict.__getitem__(self, key)
    except KeyError:
        error("Key '{}' does not exist in this ParseIdParams dict.\n"
              "Available keys are {}".format(key,self.keys()))
    return val
