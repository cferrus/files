#!/usr/bin/env python
import sys, os
from datetime import datetime
from Utils import System, call_perl

if __name__ == "__main__":
    StartTime = datetime.now()

EvDir = os.getcwd()
EvInputPath = EvDir

def call_machines(func):
    return call_perl("Machines", "new Machines()->%s" %(func),want_array=False)

def EvolveAfterID():
    print("Attempting to submit evolution.")
    ResubCmd = call_machines("GetResubCmd(('WorkDir'=>'%s','SubCmd'=>'./StartJob.sh'))" % EvDir)
    System(ResubCmd, verbose=True)

EvolveAfterID()
