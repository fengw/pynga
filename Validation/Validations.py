#!/usr/bin/env python 

import os, sys
from ValidationUtils import * 

import numpy as np 
import matplotlib.pyplot as plt 

#from my_util.functools import * 
from pynga.utils import * 


opt = sys.argv[1]
nga = sys.argv[2] 

wrk = '/Users/fengw/local/pylib/pynga/Validation'
if nga in ['BA','CB','CY','AS']:
    OFilePth = wrk + '/NGAmodelsTestFiles'
    PFilePth = wrk + '/NGAmodelsPyNGA'

    # initialize the class
    V = ValidationUtils(wrk, OFilePth, PFilePth, nga)

    # different options
    if opt == 'CalcPyNGA':
	V.CalcNGA_Py()

    if opt == 'PlotRMS':
	V.PlotRMS()

    if opt == 'PlotDiff':
	V.PlotDiff()

else: 

    OFilePth = wrk + '/DistancesTestFiles'
    inpth = OFilePth + '/inputs'
    outpth = OFilePth + '/outputs'
    plotpth0 = wrk + '/plots' 
    plotpth = plotpth0 + '/ValidationDistances' 
    for f in [plotpth0, plotpth,]: 
	if not os.path.exists(f):
	    os.mkdir(f) 

    sid = int(nga)   # test source which are in CyberShake
    metafile = outpth + '/CptDistSrc%s.py'%sid   
    
    DistKey = ['Rjb','Rrup','Rx'] 

    if opt == 'TestDist': 
	# Test distance calculation (analytical and discretized) 

	# Generate fault geometry (fault surface)
	FaultGeom = [] 
        
	SiteGeom = [-118.243,34.053,0.0] 
        Rjb0, Rrup0, Rx0 = [1.95,3.58,-1.93]   # computed by OpenSHA (Rjb,Rrup,Rx)
        
	SiteGeom = [-118.286, 34.0192,0.0]
	Rjb0, Rrup0, Rx0 = [7.16,7.76,-7.16]
	print 'OpenSHA calculation: '
	print Rjb0, Rrup0, Rx0 

	point0 = SiteGeom 
	
	# Method 1 (discretized)
	# read from srf file (Elysian Park, upper)
	SRFfile = inpth + '/158_0/158_0.txt.variation-s0000-h0000'
	srfFaultSurface = srfFaultSurfaceExtract( SRFfile ) 
	FaultGeom = srfFaultSurface['FaultGeom'] 
	Nrow,Ncol,Nelm =  FaultGeom.shape

	# method 1 (using srf file)
	Rjb1, Rrup1, Rx1 = DistanceToEvenlyGriddedSurface( SiteGeom, FaultGeom )
	print 'Using SRF file:' 
	print Rjb1, Rrup1, Rx1 

	# method 2 (simple fault)
	UpperSeisDepth = 3.0 
	LowerSeisDepth = 15.0 
	AveDip = 50.  
	FaultTrace0 = [[34.06834,-118.09977],[34.061222,-118.130356],[34.067455,-118.234761],[34.112583,-118.29677]]  # Elysian Park (Upper)
	
	Npoints = len(FaultTrace0)
	FaultTrace1 = []
	for ipoint in xrange( Npoints ):
	    FaultTrace1.append( [FaultTrace0[ipoint][1],FaultTrace0[ipoint][0], UpperSeisDepth] )
	
	daa = None; ddd = None 
	Rjb2, Rrup2, Rx2 = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1, UpperSeisDepth,LowerSeisDepth,AveDip, GridSpaceAlongStrike=daa, GridSpaceDownDip=ddd,Fast=True)
	print 'Faster way: '
	print Rjb2, Rrup2, Rx2 

        # evenlygrided surface
	daa = ddd = 1.0
	Rjb3, Rrup3, Rx3 = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1, UpperSeisDepth,LowerSeisDepth,AveDip, GridSpaceAlongStrike=daa, GridSpaceDownDip=ddd,Fast=True)
	print 'Slow way: '
	print Rjb3,Rrup3,Rx3



    if opt == 'CalcDist': 
	# Generate fault geometry (fault surface)
	FaultGeom = [] 
	if 0: 
	    #(from simple fault data like UCERF2) (Like OpenSHA for fault extension using Frankel and Stirling method)
	    
	    # first read from UCERF FM 2.0 to get the fault trace information and depth, average dip info 
	    # ...
	    pass 

	else: 
	    # read from srf file 
	    SRFfile = inpth + '/158_0/158_0.txt.variation-s0000-h0000'
	    srfFaultSurface = srfFaultSurfaceExtract( SRFfile ) 
	    FaultGeom = srfFaultSurface['FaultGeom'] 
	    Nrow,Ncol,Nelm =  FaultGeom.shape

	# read in all site info (locations)
	siteinfofile = inpth + '/cs_site_types.txt'
	lines = open( siteinfofile, 'r' ).readlines() 
	site_info = {}
	for il in xrange(1,len(lines)): 
	    spl = lines[il].strip().split() 
	    stanam = spl[0] 
	    site_info[stanam] = [float(spl[2]),float(spl[1]),0.0]   # lon, lat, dep

	# read OpenSHA computed distances (as references)
	OpenSHA_output = inpth + '/ERF35_src%s_rup0_SA3.0.csv'%(sid)
	lines = open(OpenSHA_output, 'r').readlines()
	sitesO = {}; SA = {}
	for il in range( 1, len(lines) ):
	    spl = lines[il].strip().split(',')
	    stanam = spl[1]
	    try: 
		# save all CyberShake sites 
		#  Rjb,     Rrup,    Rx
		sitesO[stanam] = [site_info[stanam][0], site_info[stanam][1], float(spl[-12]),float(spl[-11]),float(spl[-9])]
	    except: 
		continue

	# compute using pynga.utils DistanceToEvenlyGriddedSurface(SiteGeo, FaultGeo)
	sitesP = {}
	err = {}
	errR = {}
	siteErr = {}
	for ikey in xrange( len(DistKey) ): 
	    err[DistKey[ikey]] = []
	    errR[DistKey[ikey]] = []

	Nsta = len(sitesO.keys()) 
	for isite,SiteName in enumerate( sitesO.keys() ): 
	    try: 
		SiteGeo = site_info[SiteName]
	    except: 
		continue
	    Rjb, Rrup, Rx = DistanceToEvenlyGriddedSurface( SiteGeo, FaultGeom )
	    sitesP[SiteName] = [SiteGeo[0], SiteGeo[1], Rjb,Rrup,Rx]
	    
	    siteErr[SiteName] = [SiteGeo[0], SiteGeo[1]] 
	    for ikey in xrange( len(DistKey) ): 
		DistP = sitesP[SiteName][ikey+2]
		DistO = sitesO[SiteName][ikey+2]
		err0 = DistP-DistO
		
		err1 = err0/(DistO+0.00000001) * 100   # if DistO is very small, then this will be very large
		
		siteErr[SiteName].append( err0 )
		err[DistKey[ikey]].append( err0 )
		errR[DistKey[ikey]].append( err1 )  # Relative error (in percentage)

	meta = dict( 
		# fault information
		FaultGeom = FaultGeom.tolist(), 
		# distance calculated by different models
		SitesO = sitesO, 
		SitesP = sitesP, 
		# difference in distance calculation
		DistErr = err, 
		DistErrR = errR, 
		DistSiteErr = siteErr, 
		)
	save( metafile, meta, header='#Metafile For distance validation between OpenSHA and pynga\n' )

    if opt == 'PlotDist': 

	meta = load( metafile ) 
	FaultGeom = np.array( meta.FaultGeom ) 
	Nrow, Ncol, Nelm = FaultGeom.shape
	Fault = FaultGeom.reshape((Nrow*Ncol,3))
	
	err = meta.DistErr 
	errR = meta.DistErrR 
	siteErr = meta.DistSiteErr 
      
	SitesP = meta.SitesP 
	SitesO = meta.SitesO 

	#plot 
	pfmt = 'png' 

	# err (P-O)
	fig = plt.figure(1) 
	clr=['r','g','b']
	for ikey in xrange(len(DistKey) ):
	    key = DistKey[ikey] 
	    ax = fig.add_subplot( 2,2,ikey+1 )
	    ax.plot(err[key], clr[ikey]+'o')
	    ax.set_title(key)
	    if ikey == 0: 
		#ax.set_ylabel('Relative Error (%)')
		ax.set_ylabel('Error' )
	fig.savefig( plotpth + '/CptDistErrSrc%s.%s'%(sid,pfmt),format=pfmt )

        if 0:
	    # Relative err (P-O)/O
	    fig = plt.figure(2) 
	    clr=['r','g','b']
	    for ikey in xrange(len(DistKey) ):
		key = DistKey[ikey] 
		ax = fig.add_subplot( 2,2,ikey+1 )
		ax.plot(errR[key], clr[ikey]+'o')
		ax.set_title(key)
		if ikey == 0: 
		    #ax.set_ylabel('Relative Error (%)')
		    ax.set_ylabel('Error' )
	    fig.savefig( plotpth + '/CptDistErrRelativeSrc%s.%s'%(sid,pfmt),format=pfmt )

        if 1: 
	    # Scatter results (to see the distribution)
	    lons = []; lats = []
	    Rjb = []; Rrup = []; Rx = []
	    for isite, SiteName in enumerate( siteErr.keys() ): 
		lons.append( siteErr[SiteName][0] )
		lats.append( siteErr[SiteName][1] )
		Rjb.append( siteErr[SiteName][2] )
		Rrup.append( siteErr[SiteName][3] )
		Rx.append( siteErr[SiteName][4] )

	    Rjb = np.array(Rjb)
	    Rrup = np.array(Rrup)
	    Rx = np.array(Rx) 

	    fig = plt.figure(3) 
	    i = 1
	    ax = fig.add_subplot( 2,2,i )
	    ax.plot( Fault[:,0], Fault[:,1], 'k.' ) # fault surface projection
	    ax.set_title('fault surface projection')
	    plt.axis('equal')
	    for f in [Rjb, Rrup, Rx]: 
		ax = fig.add_subplot( 2,2,i+1 )
		#ax.plot( Fault[:,0], Fault[:,1], 'r.' ) # fault surface projection
		sc = ax.scatter( lons, lats, c=f, s=abs(f)*50, edgecolor='w' )
		ax.set_title(DistKey[i-1])
		#plt.axis('equal')
		ax.set_ylim([33,35.5])
		ax.set_xlim([-119.5,-116.5])
		fig.colorbar(sc)
		i += 1   
	    fig.savefig( plotpth + '/CptDistDiffScatterSrc%s.%s'%(sid,pfmt), format=pfmt )

	# P distance
	lons = []; lats = []
	Rjb = []; Rrup = []; Rx = []
	for isite, SiteName in enumerate( siteErr.keys() ): 
	    lons.append( SitesP[SiteName][0] )
	    lats.append( SitesP[SiteName][1] )
	    Rjb.append( SitesP[SiteName][2] )
	    Rrup.append( SitesP[SiteName][3] )
	    Rx.append( SitesP[SiteName][4] )

	Rjb = np.array(Rjb)
	Rrup = np.array(Rrup)
	Rx = np.array(Rx) 

	fig = plt.figure(4) 
	i = 1
	ax = fig.add_subplot( 2,2,i )
	ax.plot( Fault[:,0], Fault[:,1], 'k.' ) # fault surface projection
	ax.set_title('fault surface projection')
	plt.axis('equal')
	#ax.set_ylim([33,35.5])
	#ax.set_xlim([-119.5,-116.5])
	for f in [Rjb, Rrup, Rx]: 
	    ax = fig.add_subplot( 2,2,i+1 )
	    sc = ax.scatter( lons, lats, c=f, s=abs(f)*3, edgecolor='w' )
	    ax.set_title(DistKey[i-1])
	    ax.set_ylim([33,35.5])
	    ax.set_xlim([-119.5,-116.5])
	    #plt.axis('equal')
	    fig.colorbar(sc)
	    i += 1   
	fig.savefig( plotpth + '/CptDistPscatterSrc%s.%s'%(sid,pfmt), format=pfmt )

	# O distance
	lons = []; lats = []
	Rjb = []; Rrup = []; Rx = []
	for isite, SiteName in enumerate( siteErr.keys() ): 
	    lons.append( SitesO[SiteName][0] )
	    lats.append( SitesO[SiteName][1] )
	    Rjb.append( SitesO[SiteName][2] )
	    Rrup.append( SitesO[SiteName][3] )
	    Rx.append( SitesO[SiteName][4] )

	Rjb = np.array(Rjb)
	Rrup = np.array(Rrup)
	Rx = np.array(Rx) 

	fig = plt.figure(5) 
	i = 1
	ax = fig.add_subplot( 2,2,i )
	ax.plot( Fault[:,0], Fault[:,1], 'k.' ) # fault surface projection
	ax.set_title('fault surface projection')
	plt.axis('equal')
	#ax.set_ylim([33,35.5])
	#ax.set_xlim([-119.5,-116.5])
	for f in [Rjb, Rrup, Rx]: 
	    ax = fig.add_subplot( 2,2,i+1 )
	    sc = ax.scatter( lons, lats, c=f, s=abs(f)*3, edgecolor='w' )
	    ax.set_title(DistKey[i-1])
	    ax.set_ylim([33,35.5])
	    ax.set_xlim([-119.5,-116.5])
	    #plt.axis('equal')
	    fig.colorbar(sc)
	    i += 1   
	fig.savefig( plotpth + '/CptDistOscatterSrc%s.%s'%(sid,pfmt), format=pfmt )


