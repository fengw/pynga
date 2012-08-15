#!/usr/bin/env python 
# Verify pynga using OpenSHA results

import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt 

from pynga import *
from pynga.utils import *

opt = sys.argv[1]   # Compute, Plot 

# OpenSHA path
wrk = os.path.realpath( os.getcwd() )
OpenSHA_pth = wrk + '/OpenSHA_NGA08/'

output_pth = wrk + '/outputs' 
plot_pth = wrk + '/plots' 
for f in [ output_pth, plot_pth, ]: 
    if not os.path.exists( f ):
	os.mkdir( f )
metafile = output_pth + '/NGAcptVer.py'

# NGA computing related
periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	   0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
NGAmodels = ['CB','BA','CY','AS']
NGAs = {'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)},\
	}, \

# source info (one fault)
sid = 79

if opt == 'Prepare':
    # you need have OpenSHA install on the local machine to do so 
    filen = OpenSHA_pth + '/nga_inputs'
    fid = open( filen, 'w' )
    for ip in xrange( len(periods) ):
        Tkey = '%.3f'%periods[ip]
        fid.write( '%s 0 %s\n'%(sid,Tkey) )
    fid.close() 
    os.chdir(OpenSHA_pth)
    os.system( 'nga_comparison_calc nga_inputs' )
    os.chdir(wrk)

if opt == 'Compute': 
    # read source info from metafile which has all information about the fault
    meta = load( wrk+ '/SrcInfo/meta_rup_35_3_%s.py'%sid )
    meta_rup = meta.rups_info['MR'] 
    Mws = meta_rup[0]
    rake,dip,Ztor,Zbom = meta_rup[3:]
    W = Zbom-Ztor

    rid = 0
    Mw = Mws[rid]

    # read site info and distance info with OpenSHA compute PSA values (median and std)
    # site would be the same for all periods
    OpenSHA_output = OpenSHA_pth + '/ERF%s_src%s_rup%s_SA%s.csv'%(35, sid, rid, '%2.1f'%(3.0))
    lines = open(OpenSHA_output, 'r').readlines()
    Vs30 = []; Z25 = []; Z10 = []
    Rjb = []; Rrup = []; Rx = []
    for il in range( 1, len(lines) ):
	spl = lines[il].strip().split(',')
	stanam = spl[1]
	#Vs30 (m/s), Z2.5 (km),     Z1.0 (m),     Rjb,     Rrup,    Rx (km)
	Vs30.append( float(spl[-18]) ) 
	Z25.append( float(spl[-17]) ) 
	Z10.append( float(spl[-15]) )
	Rjb.append( float(spl[-12]) )
	Rrup.append( float(spl[-11]) )
	Rx.append( float(spl[-9]) )
    
    ngaP = {}; ngaO = {}
    for inga, nga in enumerate( NGAmodels ):
	ngaP[nga] = {}; ngaO[nga] = {}
	index1 = 2*inga
	
	# all OpenSHA_NGA08 files to get correct file name
	filenames = glob.glob( OpenSHA_pth + '/ERF%s_src%s_rup%s_SA*.csv'%(35, sid, rid) ) 
	for file1 in filenames: 
	    fspl = file1.strip().split('/')[-1].split('_')[-1]
	    Ti = float(fspl[2:-4])
	    print Ti
	    Tkey = '%.3f'%Ti
	    ngaP[Tkey] = []; ngaO[Tkey] = []
	    median, std, tau, sigma = NGA08(nga, Mw, Rjb, Vs30, Ti, rake=rake,Mech=None, \
					    Rrup=Rrup, Rx=Rx, dip=dip,W=W,Ztor=Ztor,Z25=Z25,Z10=Z10,Fas=0,AB11=None,VsFlag=0)
	    ngaP[nga][Tkey] = list(np.log(median)), list(np.log(std))    # all in natural log

	    lines = open(file1, 'r').readlines()
	    tmpO_median = []; tmpO_sigmaT = []
	    for il in range( 1, len(lines) ):
		spl = lines[il].strip().split(',')
		tmpO_median.append(float(spl[index1-8]))
		tmpO_sigmaT.append(float(spl[index1-7]))   # all in natural log 
	    ngaO[nga][Tkey] = tmpO_median, tmpO_sigmaT

    # save ngaP and ngaO into metafile for further visual comparison 
    meta = dict( 
	    ngaP = ngaP, 
	    ngaO = ngaO, 
	    ) 
    save( metafile, meta, header='#NGA comparison\n' )

if opt == 'Plot':
    Ti = float(sys.argv[2])
    Tkey = '%.3f'%Ti 

    pfmt = 'png' 
    try: 
	meta = load( metafile )
    except: 
	print 'Compute NGAs and generate metafile for Plot option'
	raise ValueError
    ngaP = meta.ngaP 
    ngaO = meta.ngaO 
    
    xlab, ylab = 'Python', 'OpenSHA'
    nameS = 'Mean','sigmaT'
    clr = 'r', 'b'
    for i in xrange( len( nameS ) ):
	fig = plt.figure(i+1) 
	for inga, nga in enumerate( NGAmodels ):
	    ax = fig.add_subplot( 2,2, inga+1) 
	    ngaP1 = ngaP[nga][Tkey][i] 
	    ngaO1 = ngaO[nga][Tkey][i] 
	    ax.plot( ngaP1, ngaO1, clr[i]+'.' )
	    if inga == 0: 
		ax.set_title( nga+' '+nameS[i] )
	    else: 
		ax.set_title( nga )
            if inga in [0,2]: 
		ax.set_ylabel( 'OpenSHA' )
	    if inga in [2,3]:
		ax.set_xlabel( 'Python' )
	    plt.axis('equal')
	fig.savefig( plot_pth + '/NGAs_PSA%s_%s.%s'%(Tkey,nameS[i], pfmt), format=pfmt )
