#!/usr/bin/env python 
"""
Verificiation between Python compute SC08 and Matlab's
As published on ES 2008 by Spudich and Chiou
"""
import os
import numpy as np 
import matplotlib.pyplot as plt 

from pynga.SC08 import * 

inpth = './Matlab_SC08/'

pltpth0 = './plots/'
pltpth = pltpth0 + '/VerificiationSC08' 
for f in [ pltpth0, pltpth, ]: 
    if not os.path.exists( f ):
	os.mkdir(f)

Nh = 10    # for this verificiation 
Ts = [2.0, 3.0, 5.0, 10.0] 

# parameters need in SC08
Mw = 7.45  # this is just used for magnitude taper
model_name = 'BA'
NewCoefs = None
cuteps = {'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[200,250],'m_cut':[5.6,6.0]}

# compure results for each hypocenter
fig = plt.figure(1)
pfmt = 'png'
for ih in xrange(Nh):
    fig.clf()

    SC08_outputs = inpth + 'hypo%s_BA6.txt'%(ih+1)
    outputs = np.loadtxt( SC08_outputs, delimiter=',', skiprows=1 ) 
    Rrup = outputs[:,3]
    Rfn = outputs[:,5]
    Rfp = outputs[:,6]
    s = outputs[:,7]
    h = outputs[:,8]
    ctildepr = outputs[:,9]
    IDP = outputs[:,10]

    # read in information
    for it in xrange( len(Ts) ):
	Ti = Ts[it] 
	SC08M = outputs[:,it+15] 

	SC08 = SC08_model(model_name+'08',cuteps=cuteps )
	kwds = {'NewCoefs':NewCoefs}

	mapfunc( SC08, Mw, Rrup, ctildepr, s, h, Rfn, Rfp, Ti,**kwds )
	SC08P = np.log(SC08.fD)
	if it == 0: 
	    IDPP = SC08.IDP 
	
	# plot 
	ax = fig.add_subplot(2,2,it+1)
	ax.plot( SC08P, SC08M, 'r.' )
    fig.savefig( pltpth + '/SC08_PythonMatlab_hypo%s_fD.%s'%(ih,pfmt), format=pfmt )
    
    fig.clf()

    ax = fig.add_subplot(111)
    ax.plot(IDPP, IDPP, 'r.')
    fig.savefig( pltpth + '/S08_PythonMatlab_hypo%s_IDP.%s'%(ih,pfmt), format=pfmt )

