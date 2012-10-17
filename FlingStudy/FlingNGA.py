#!/usr/bin/env python
"""
FlingStudy to test the application of pynga
FlingStudy Fault trace model is similar to BBP *.src file 
"""
# required modules
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt

# pynga
from pynga import *
from pynga.utils import *

sid = sys.argv[1]   # Scenario ID (140,141,142)
vid = sys.argv[2]   # fault id (hypocenters)
opt = sys.argv[3]   # ComputeDist, ComputePSA, PlotPSA

dflag = sys.argv[4]  # proj or extend

# paths
wkd = os.getcwd() 

inpth0 = wkd + '/inputs/'
plotpth = wkd + '/plots/'
metapth = wkd + '/outputs/'
for f in [metapth, plotpth,]:
    if not os.path.exists( f ):
	os.mkdir( f )

inpth = inpth0 + 'Scenario%s'%sid
if not os.path.exists( plotpth ):
    os.mkdir( plotpth )

srcpth = inpth + '/SrcDesc1/'
stapth = inpth + '/StatInfo/'

# read in the station list and site information
stafiles = glob.glob( stapth + '*.ll' )
stafile = stafiles[0]   # only one
lines = open( stafile ).readlines()
rlon = []; rlat = []; rnam = []
for il in xrange( len( lines) ):
    spl = lines[il].strip().split()
    rlon.append( float( spl[0] ) )
    rlat.append( float( spl[1] ) )
    rnam.append( spl[2] )

# read from 1D model and get basin depth and Vs30
vfile = inpth0 + 'VelModel1D'
lines = open( vfile ).readlines()
layer = []; vs = []
for il in range( 3, len(lines) ):
    spl = lines[il].strip().split()
    layer.append( float( spl[0] ) )
    vs.append( float( spl[2] ) )
layer = np.array( layer )
vs = np.array( vs )
index = ( layer <= 30./1000.).nonzero()[0]
Vs30 = np.mean( vs[index] ) * 1000 # from km/s to m/s

# use this one as the Vs30
Vs30 = 865   # m/s (as shown in the velocity model could be obtained by travel time averaging?)

# read fault information and hypocenter locations
vpth = srcpth + '/v%s/'%( vid )
sfiles = glob.glob( vpth + '*.src' )
lines = open( sfiles[0] ).readlines() 
values = []
for il in xrange( len(lines) ):
    spl = lines[il].strip().split()
    values.append( float( spl[2] ) )

M, Fl, dfl, Fw, dfw, ztor, strike, rake, dip, lat0, lon0, hypoAS, hypoDD = values[:-3]

# write into metafile (for further visualization)
plotpth0 = plotpth + 'PSA/'
plotpth1 = plotpth0 + 'Scenario%s/'%sid
plotpth2 = plotpth1 + 'v%s/'%vid 

plotpth00 = plotpth + 'FaultGeom/'
plotpth10 = plotpth00 + 'Scenario%s/'%sid
plotpth20 = plotpth10 + 'v%s/'%vid 

metapth0 = metapth + 'PSA/'
metapth1 = metapth0 + 'Scenario%s/'%sid
metapth2 = metapth1 + 'v%s/'%vid
for f in [metapth0, metapth1, metapth2,\
	  plotpth00, plotpth10, plotpth20, \
	  plotpth0, plotpth1, plotpth2, ]:
    if not os.path.exists( f ):
	os.mkdir( f )
metafileD = metapth2 + 'Distance.Scenario%s.v%s_%s.py'%(sid,vid,dflag) 
metafilePSA = metapth2 + 'PSA.Scenario%s.v%s.py'%(sid,vid) 

NGA_models = ['CB','BA','CY','AS']
clrs = ['r','b','g','m']
syms = ['o','v','^','*']

if opt == 'ComputeDist': 
    # =================================
    # Get the Fault Geometry and compute Rjb, Rrup, and Rx distance for AS NGA computing
    IDs = sid, vid
    Dims = Fl, dfl, Fw, dfw, ztor
    Mech = strike, dip, rake
    HypoLoc = hypoAS, hypoDD
    origin = lon0, lat0

    if dflag == 'proj': 
	# projection property
	zone = 10
	rot = strike
	inverse = True  # from xy to ll
	ProjDict = {'zone':zone,'origin':origin,'rot':rot,'inverse':inverse}
	VisualDict={'SiteLoc':(rlon,rlat),'savetype':'pdf','plotpth':plotpth20} 
	VisualDict=None
	from FlingUtils import *
	FaultGeom = FaultGeom(IDs, Dims, Mech, HypoLoc, ProjDict,VisualDict=VisualDict) 
	SiteGeo = rlon, rlat   # all site Geometry
	start_time = HourMinSecToSec(BlockName=opt)
	DistDict = calc_distances( SiteGeo, Dims, Mech, ProjDict, Rrup=True, Rx=True )
	end_time = HourMinSecToSec()
	SecToHourMinSec( end_time-start_time,BlockName=opt )
   
    else: 
	start_time = HourMinSecToSec(BlockName=opt)
	
	FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, GridSpaceAlongStrike, GridSpaceDownDip = FaultTraceGen( origin, Dims, Mech ) 
	FaultGeom = SimpleFaultSurface(FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, GridSpaceAlongStrike=GridSpaceAlongStrike, GridSpaceDownDip=GridSpaceDownDip ) 

	DistDict = {}
	for key in ['Rjb','Rrup','Rx']: 
	    DistDict[key] = []

	Rjb = []; Rrup = []; Rx =[]
	for isite in xrange( len(rlon) ): 
	    SiteGeom = [rlon[isite], rlat[isite], 0.0]
	    # simple calculation  (short time)
	    Rjb0, Rrup0, Rx0 = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1,UpperSeisDepth,LowerSeisDepth,AveDip)
	    # discretized calculation (very long time)
	    #Rjb0, Rrup0, Rx0 = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1,UpperSeisDepth,LowerSeisDepth,AveDip,GridSpaceAlongStrike=GridSpaceAlongStrike,GridSpaceDownDip=GridSpaceDownDip)
	    DistDict['Rjb'].append(Rjb0) 
	    DistDict['Rrup'].append(Rrup0) 
	    DistDict['Rx'].append(Rx0) 

	end_time = HourMinSecToSec()
	SecToHourMinSec( end_time-start_time,BlockName=opt )


    meta = dict( 
	    distance=DistDict, 
	    FaultGeom=FaultGeom,
	    )
    save( metafileD, meta, header = '# Distance and FaultGeom MetaData \n' ) 


if opt == 'PlotDist':
    
    meta1 = load( metafileD )
    Rx = meta1.distance['Rx']
    Rrup = meta1.distance['Rrup']
    Rjb = meta1.distance['Rjb']
    
    if dflag == 'proj':
	azimuth = meta1.distance['azimuth'] 
	slon2dS, slat2dS = meta1.FaultGeom['FaultSurface'] 
	slon2dS = np.array(slon2dS) 
	slat2dS = np.array(slat2dS) 
	fig = plt.figure(1) 
	i = 1
	tls = ['azimuth','Rjb','Rrup','Rx']
	for f in [azimuth, Rjb, Rrup, Rx]:
	    f = np.array(f) 
	    ax = fig.add_subplot( 2,2,i )
	    img = ax.scatter(rlon,rlat,c=f,s=2*abs(f),edgecolor='w')
	    fig.colorbar(img)
	    ax.plot( [slon2dS[0,0], slon2dS[0,-1],slon2dS[0,-1], slon2dS[0,0], slon2dS[0,0]], \
		     [slat2dS[0,0], slat2dS[0,0],slat2dS[-1,0], slat2dS[-1,0], slat2dS[0,0]],'k' )
	    ax.set_title( tls[i-1] )
	    i = i + 1
	fig.savefig( plotpth + 'DistanceAzimuthScatter%s.pdf'%dflag, format='pdf' )   
    
    else: 

	F = np.array( meta1.FaultGeom )
	fig = plt.figure(1) 
	i = 1
	ax = fig.add_subplot( 2,2,1 )
	ax.plot( rlon, rlat, 'k^' ) 
	ax.plot( lon0, lat0, 'yo' )
	ax.plot( [F[0,0,0],F[0,-1,0],F[-1,-1,0],F[-1,0,0],F[0,0,0]],\
		 [F[0,0,1],F[0,-1,1],F[-1,-1,1],F[-1,0,1],F[0,0,1]],'b')
	ax.set_title( 'FaultSurface and Sites' )
	tls = ['Rjb','Rrup','Rx']
	for f in [Rjb, Rrup, Rx]:
	    f = np.array(f) 
	    ax = fig.add_subplot( 2,2,i+1 )
	    ax.plot( [F[0,0,0],F[0,-1,0],F[-1,-1,0],F[-1,0,0],F[0,0,0]],\
		     [F[0,0,1],F[0,-1,1],F[-1,-1,1],F[-1,0,1],F[0,0,1]],'k')
	    img = ax.scatter(rlon,rlat,c=f,s=2*abs(f),edgecolor='w')
	    fig.colorbar(img)
	    ax.set_title( tls[i-1] )
	    i = i + 1
	fig.savefig( plotpth + 'DistanceAzimuthScatter%s.pdf'%dflag, format='pdf' )   



if opt=='ComputePSA':

    # compute SA value using NGA models
    periods = [0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	       0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
    meta1 = load( metafileD )
    Rx = meta1.distance['Rx']
    Rrup = meta1.distance['Rrup']
    Rjb = meta1.distance['Rjb']

    start_time = HourMinSecToSec(BlockName=opt)
    print rake 

    psaP = {}
    for inga in xrange( len(NGA_models) ):
	nga = NGA_models[inga]
	print 'Compute PSA for %s'%nga
	psaP[nga] = {} 
	for ip in xrange( len(periods) ):
	    T = periods[ip]
	    if T == 0.0: 
		Tp = -1.0   # PGA for Python
                Tr = 0.0
            else: 
		Tp = T 
		Tr = T

	    key = '%s'%'%.2f'%T
	    median, std, tau, sigma = NGA08(nga, M, Rjb, Vs30, Tp, \
		    rake = rake, Mech=None, dip = dip, W = Fw, Ztor = ztor, Rx =Rx, Rrup = Rrup )
	    psaP[nga][key] = list( median )  
            
    end_time = HourMinSecToSec()
    SecToHourMinSec( end_time-start_time, BlockName='finish '+opt )
    
    meta = dict( 
	    distance=meta1.distance, 
	    psaP = psaP,
	    )

    save( metafilePSA, meta, header = '# PSA metadata\n' )

if opt == 'PlotPSA':
    # plot for given periods (psa vs Rjb)
    T0 = '%s'%'%.2f'%float(sys.argv[4])   # as float
    
    savetype = 'pdf'

    # read from metafile for test purpose
    meta = load( metafilePSA )
    Rjb = meta.distance['Rjb']
    Rx = meta.distance['Rx']  # problem of compute Rx (debug by plot azimuth and Rjb first, then and Rrup, Rx)
    Rrup = meta.distance['Rrup'] 
    azimuth = meta.distance['azimuth']

    fig = plt.figure(1)
    ax1 = fig.add_subplot( 221 )
    ax2 = fig.add_subplot( 222 )
    ax3 = fig.add_subplot( 223 )
    #ax4 = fig.add_subplot( 224 )
    lines = []
    for inga in xrange( len( NGA_models) ):
	nga = NGA_models[inga]
	psa0 = meta.psaP[nga][T0]
	line = ax1.loglog( Rjb, psa0, clrs[inga]+syms[inga] )
	lines.append( line )
	if nga != 'BA' and nga != 'CB':
	    ax2.plot( Rx, np.log(psa0), clrs[inga]+syms[inga] )
	if nga != 'BA':
	    ax3.loglog( Rrup, psa0, clrs[inga]+syms[inga] )
    
    ax1.set_xlabel( r'$R_{JB}$ (km)' )
    ax1.set_ylabel( r'5%-damped PSA (g)' )
    ax2.set_xlabel( r'$R_{x}$ (km)' )
    ax1.set_ylabel( r'5%-damped PSA (g)' )
    ax3.set_xlabel( r'$R_{rup}$ (km)' )
    ax3.set_ylabel( r'5%-damped PSA (g)' )

    ax1.set_title( 'PSA at T=%s sec'%T0 )
    ax1.legend( lines, NGA_models, loc=0 )
    plotn ='PSA.Scenario%s.v%s.T%s.NGAs.%s'%(sid,vid,T0,savetype)
    fig.savefig( plotpth2 + plotn, format=savetype )

    #plt.show()
