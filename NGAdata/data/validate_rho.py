#!/usr/bin/env python
"""
Valiate Python compueted NGA and Matlab computed NGA
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt

import rpy2.robjects as robjects

import cst

from my_util.image import *

eps = sys.argv[1]   # total, inter, intra

epsdict = {'total':0, 'inter':1, 'intra':2}

NGAmodels = ['AS',]
clr = ['r',]
NGAmodels = ['CB',]
clr = ['r',]
NGAmodels = ['CB','BA','CY','AS']
clr = ['b','r','g','#808080']


wrk = '/Users/fengw/local/pylib/pynga'
data = os.path.join( wrk, 'data' )
plots = os.path.join( wrk, 'plots/Validation_NGAs' )
if not os.path.exists( plots ):
    os.mkdir( plots )

T=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 5, 7.5, 10];
periods = T[1:]

fig = plt.figure(1,(12,6))
ax = fig.add_subplot(111)
ax.set_title( 'NGA models cross-verification using correlation of SA and PGA' )
plotname = 'Validation_rho_%s.pdf'%eps
lines = []
for inga, nga in enumerate( NGAmodels ):
    
    # rho computed from Python
    metafile = data + '/meta_eps_rho_%s.py'%eps
    meta = cst.util.load( metafile )
    rho01 = meta.rho_PGA_SA

    # rho computed from Matlab
    meta2file = data + '/meta_Rhos_%s.txt'%nga
    meta2 = np.loadtxt( meta2file )
    line0 = ax.semilogx( periods, rho01[nga], color=clr[inga], linestyle='solid', lw=1.5 )
    lines.append(line0)
    plt.hold(True) 
    ax.semilogx( periods,  meta2[:,epsdict[eps]], color=clr[inga], linestyle='dashed', lw=3 )
    plt.hold(True)

    plt.grid(True)
    plt.grid( b = True, which='minor' )
    
lg = plt.legend( lines, NGAmodels, loc=0 )
lg.draw_frame(False)
ax.set_ylim([0.0,1.0])
ax.set_xlabel( 'periods (sec)' )
ax.set_ylabel( r'$\rho_{lnPGA,lnSA}$' )
fig.savefig( plots + '/' + plotname, savetype='pdf' )

