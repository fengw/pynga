#!/usr/bin/env python 

import os, sys
from ValidationUtils import * 

opt = sys.argv[1]
nga = sys.argv[2] 

wrk = '/Users/fengw/local/pylib/pynga/Validation'
OFilePth = wrk + '/NGAmodelsTestFiles'
PFilePth = wrk + '/NGAmodelsPyNGA'

# initialize the class
NGAsV = NGAsValidation(wrk, OFilePth, PFilePth, nga)

if opt == 'CalcPyNGA':
    NGAsV.CalcNGA_Py()

if opt == 'PlotRMS':
    NGAsV.PlotRMS()

if opt == 'PlotDiff':
    NGAsV.PlotDiff()
