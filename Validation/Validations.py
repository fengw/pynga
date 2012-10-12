#!/usr/bin/env python 

import os, sys
from ValidationUtils import * 

opt = sys.argv[1]
nga = sys.argv[2] 

wrk = '/Users/fengw/local/pylib/pynga/Validation'
OFilePth = wrk + '/NGAmodelsTestFiles'

if nga in ['BA','CB','CY','AS']:
    PFilePth = wrk + '/NGAmodelsPyNGA'
else: 
    PFilePth = wrk + '/DistancesPyNGA'

# initialize the class
V = ValidationUtils(wrk, OFilePth, PFilePth, nga)

# different options
if opt == 'CalcPyNGA':
    V.CalcNGA_Py()

if opt == 'PlotRMS':
    V.PlotRMS()

if opt == 'PlotDiff':
    V.PlotDiff()
