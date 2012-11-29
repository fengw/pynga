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

nga = sys.argv[1]

wrk = '/Users/fengw/local/pylib/pynga'
data = os.path.join( wrk, 'data' )
plots = os.path.join( wrk, 'plots/Validation_NGAs' )
if not os.path.exists( plots ):
    os.mkdir( plots )
plots = os.path.join( plots, 'Validation_%s'%nga )
if not os.path.exists( plots ):
    os.mkdir( plots )


T=[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 5, 7.5, 10];
Ndata = 1842

# compute R median value
meta2file = data + '/Rmeta_IMs_%s.txt'%nga
#robjects.r['source'](os.path.expanduser('~/work/Project/CyberShake_analysis/utils/nga_util/R_computeNGA.R'))
#robjects.r['R_computeNGA'](data+'/FlatFile_Subset_file.txt', T, nga, meta2file )

meta2 = np.loadtxt( meta2file )

meta = cst.util.load( data + '/meta_IMs_%s.py'%nga )
PGA_IMs = meta.PGA_IMs
SA_IMs = meta.SA_IMs

meta1 = np.loadtxt( data + '/meta_IMs_%s.txt'%nga )
Nsub = 6; 
keys = ['median','tau','sigma','epsilonT','eta','epsilon']
inds = [0,2,3,4,5,6]

T0 = float( sys.argv[2] )   # 0: PGA

ind = (np.array(T)==T0).nonzero()[0]
if len(ind) == 0:
    print 'T0 is not in T list'
    raise ValueError
else:
    ip = ind[0]

fig = init_fig( num =1, figsize=(14,10), dpi=100 )
fig.clf()
axs = init_subaxes( fig, subs=(2,3), basic=(0.7,0.7,0.7,0.7) )

if ip == 0:
    ind0 = 2
    for ie in xrange( 6 ):
	P_IMs = np.array( PGA_IMs[keys[ie]] )
	P_IMs0 = np.sort( P_IMs )
	M_IMs = meta1[:,ind0+inds[ie]]
	
	ax = fig.add_axes( axs[ie] )
	if ie == 0:
	    R_IMs = np.log(meta2[:,ip+1])
	    ax.plot( P_IMs, M_IMs, 'r.', P_IMs, R_IMs,'b.', P_IMs0, P_IMs0, 'k' )
	    print 'Record ID', 'Python Computed', 'Matlab Computed'
	    for irecord in xrange( Ndata ):
		if abs( P_IMs[irecord] - M_IMs[irecord] ) > 1.e-3:
		    print meta1[irecord,0], np.exp(P_IMs[irecord]), np.exp(M_IMs[irecord])

	    #ax.plot( P_IMs, M_IMs, 'r.', P_IMs0, P_IMs0, 'k' )
	else:
	    ax.plot( P_IMs, M_IMs, 'r.', P_IMs0, P_IMs0, 'k' )
	ax.set_title('%s'%keys[ie])

    fig.text( 0.2,0.97,'PGA IMs Validation, Python vs Matlab' )
    fig.savefig( plots + '/Validation_PGA_%s.pdf'%(nga), savetype='pdf' )

else:
    key = '%s'%'%5.3f'%T[ip]
    ind0 = 2 + 7*ip 
    for ie in xrange( 6 ):
	P_IMs = np.array( SA_IMs[key][keys[ie]] )
	P_IMs0 = np.sort( P_IMs )
	M_IMs = meta1[:,ind0+inds[ie]]
	ax = fig.add_axes( axs[ie] )
	if ie == 0:
	    R_IMs = np.log(meta2[:,ip+1])
	    ax.plot( P_IMs, M_IMs, 'r.', P_IMs,R_IMs,'b.', P_IMs0, P_IMs0, 'k' )
	    #ax.plot( P_IMs, M_IMs, 'r.', P_IMs0, P_IMs0, 'k' )
	    print 'Record ID', 'Python Computed', 'Matlab Computed'
	    for irecord in xrange( Ndata ):
		if abs( P_IMs[irecord] - M_IMs[irecord] ) > 1.e-3:
		    print meta1[irecord,0], np.exp(P_IMs[irecord]), np.exp(M_IMs[irecord])
	else:
	    ax.plot( P_IMs, M_IMs, 'r.', P_IMs0, P_IMs0, 'k' )
	ax.set_title('%s'%keys[ie])

    fig.text( 0.2,0.97,'SA at %s sec, IMs Validation, Python vs Matlab'%'%3.2f'%T[ip] )
    fig.savefig( plots + '/Validation_SA_%s_period%s.pdf'%(nga,'%3.2f'%T[ip]), savetype='pdf' )


