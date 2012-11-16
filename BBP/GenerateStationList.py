#!/usr/bin/env python 
"""
Generate Station List based on Rrup 
"""
import os,sys 
import numpy as np 

from pynga import * 
from pynga.utils import * 

R0 = float(sys.argv[1])   # in km (distance extend from the both ends of the rupture, used to generate the grid space)
GridSize = float(sys.argv[2])   # grid size to generate station list 
NstaSelected = int(sys.argv[3])   # number of station to choose 
VisualFlag = int(sys.argv[4])   # plot or not (for check the distance calculations)

# SRC file (specify the location of the source)
sfile = './inputs/bbp_validation_part_b_scenario_src_files'
lines = open( sfile ).readlines()
values = []
for il in xrange( len(lines) ):
    spl = lines[il].strip().split()
    values.append( float( spl[1] ) )
M, L, dl, W, dw, ztor, strike, rake, dip, lat0, lon0, hypoAS, hypoDD = values[:13]

origin = lon0, lat0    # center 
Dims = L, dl, W, dw, ztor
Mech = strike, dip, rake 
FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, dl,dw = FaultTraceGen( origin, Dims, Mech )
FaultTrace, FaultSeg, AveStrike = SimpleFaultSurface( FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip ) 

# generate grid of stations (for grid search)  in lon/lat 
strike = strike*np.pi/180.
loc0 = [lon0, lat0, 0.0]

vD = 0.0
dhD = GridSize   # in km (grid size) 
hDx = 60.   # 20~60 km
hDy = 0.5*L+R0+dhD*2    # along strike

az = strike + np.pi 
vector = [az,hDy,vD]
loc1 = EndLocation( loc0, vector )

az = strike + np.pi*3./2.
vector =[ az, hDx, vD ] 
loc11 = EndLocation( loc1, vector )  # new origin 

Nx = int(2*hDx/dhD + 1)
Ny = int(2*hDy/dhD + 1)

LocY = loc11
Loc2D = []
for iy in xrange( Ny ): 
    Loc2D.append( LocY )
    az = strike + np.pi/2.
    hD = dhD
    vector = [az, hD, vD]
    LocX = LocY 
    for ix in xrange( Nx-1 ): 
	LocX = EndLocation(LocX, vector )
	Loc2D.append( LocX )
    az = strike
    hD = dhD
    vector = [az, hD, vD]
    LocY = EndLocation(LocY, vector )

Nloc = len(Loc2D) 
LocS0 = np.array( Loc2D ) 
print 'Total station locations: ', Nloc 

Rrups = []; Rxs = []
for iloc in xrange( len(Loc2D) ): 
    SiteGeom = Loc2D[iloc] 
    Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1,UpperSeisDepth,LowerSeisDepth,AveDip)    # in km 
    Rrups.append( Rrup ) 
    Rxs.append( Rx ) 

if VisualFlag:
    import matplotlib.pyplot as plt 
    from mpl_toolkits.mplot3d import Axes3D
    vertsClosed = FaultTrace + [FaultTrace[0]]
    verts1 = np.array( vertsClosed )
    fig = plt.figure(1) 
    ax = fig.add_subplot( 111 ) 
    img = ax.scatter( LocS0[:,0], LocS0[:,1], c=np.array( Rrups ), s=np.array(Rrups), edgecolor='none' )
    ax.plot( verts1[:,0], verts1[:,1], 'k' )
    ax.set_xlabel('lon')
    ax.set_ylabel('lat') 
    ax.set_title( 'All station grid and Rrup with domain: %s*%s km^2, and grid size: %s km'%(2*hDx,2*hDy,dhD) ) 
    fig.colorbar(img)
    fig.savefig( './plots/AllStationRrup_GridSize%s.pdf'%('%.2f'%GridSize), format='pdf' )

    fig = plt.figure(2) 
    ax = fig.add_subplot( 111 ) 
    img = ax.scatter( LocS0[:,0], LocS0[:,1], c=np.array( Rxs ), s=abs(np.array(Rxs)), edgecolor='none' )
    ax.plot( verts1[:,0], verts1[:,1], 'k' )
    ax.set_xlabel('lon')
    ax.set_ylabel('lat') 
    ax.set_title( 'All station grid and Rx' ) 
    fig.colorbar(img)
    fig.savefig( './plots/AllStationRx_GridSize%s.pdf'%('%.2f'%GridSize), format='pdf' )


# grid search (allowing +- 0.5km range around Rrup0. Rrup0 is the cutoff distance to select stations from the grid)
for Rrup0 in [30, 50,]:
    Loc2D1 = []; Rrups1 = []
    for iloc in xrange( Nloc ):
	if Rrup0-dhD/2.<= Rrups[iloc] <= Rrup0+dhD/2. and Rxs[iloc] <= 0:
	    Loc2D1.append( Loc2D[iloc] ) 
	    Rrups1.append( Rrups[iloc] ) 
	else: 
	    continue
    Nloc1 = len(Loc2D1) 
    print 'Stations with distance Rrup=%s: '%(Rrup0), Nloc1 
    if Nloc1 < NstaSelected: 
	print 'Change your grid size smaller in order to do the subsampling'
	raise ValueError

    LocS1 = np.array( Loc2D1 ) 
    
    # Generate random integer (total number would be NstaSelected) between 0 and Nloc 
    LocS2 = []; i = 0
    index = []; Rrups2 = []
    while i < NstaSelected: 
	index0 = np.random.randint( 0, Nloc1-1 ) 
	if not index0 in index: 
	    LocS2.append( LocS1[index0] )
	    Rrups2.append( Rrups1[index0] )
	    index.append( index0 )
	    i = i + 1 
	else: 
	    continue
    LocS2 = np.array( LocS2 ) 
    
    # write into file
    fid = open( './outputs/bbp_validation_part_b_scenario_StationList_FW_RrupCutoff%s'%Rrup0, 'w' )
    fid.write('#lon lat Rrup\n')
    for ista in xrange( NstaSelected ): 
	fid.write( '%s %s %s\n'%(LocS2[ista,0], LocS2[ista,1],Rrups2[ista]) )
    fid.close() 
    


    if VisualFlag:

	fig = plt.figure(3) 
	fig.clf()
	ax = Axes3D(fig) 
	DepFactor = 1./6371*180./np.pi
	DepFactor = 1.0
	ax.plot( verts1[:,0], verts1[:,1], -verts1[:,2]*DepFactor, 'k--', label='Fault Surface' )
	ax.plot( verts1[:,0], verts1[:,1], 0.0, 'k', label='Fault Surface Projection')
	ax.plot( LocS0[:,0], LocS0[:,1], 0.0, 'g+', label='Station Grid' ) 
	ax.text( verts1[3,0], verts1[3,1], -verts1[3,2]*DepFactor, 'dip=%s'%('%.1f'%AveDip) ) 

	ax.set_xlabel('lon')
	ax.set_ylabel('lat')
	ax.set_zlabel('dep (km)') 
	ax.set_zlim3d([-110,0])

	ax.plot( LocS1[:,0], LocS1[:,1], 0.0, 'b^', label='Station Rrup=~%s km'%Rrup0)
	
	ax.plot( LocS2[:,0], LocS2[:,1], 0.0, 'r*', markersize=12, label='Selected %s Stations (Rrup=~%s km)'%(NstaSelected,Rrup0) ) 
	lg = ax.legend(loc=0)
	lg.draw_frame(False)
	ltext = plt.gca().get_legend().get_texts()
	plt.setp( ltext, fontsize = 8 )
	fig.savefig('./plots/StationSelection_GridSize%sNsta%sRrupAt%s.pdf'%('%.2f'%GridSize,NstaSelected,Rrup0),format='pdf')


