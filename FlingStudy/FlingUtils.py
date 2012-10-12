#!/usr/bin/env python
"""
Utilities for Fling
"""
# common modules
import numpy as np
import matplotlib.pyplot as plt

# extra modules needed
from pynga.utils import *


def FaultGeom(IDs, Dims, Mech, HypoLoc, ProjDict, VisualDict=None):
    """
    from Fling source files to get the fault geometry
    """
    ScenarioID, VID = IDs  # e.g. 142, 00
    Fl,dfl,Fw,dfw,ztor  = Dims   # fault length along strike, fault width down dip 
    strike, dip, rake = Mech   # fault mechanism 
    hypoAS, hypoDD = HypoLoc   # relative hypocenter loation (along strike, down dip)

    # UTM projection property
    lon0, lat0 = ProjDict['origin']   # projection origin
    kwds = ProjDict
    
    # Create fault mesh:
	# along strike direction: y-axis
	# along dip direction: x-axis
    fx = np.arange( 0, Fw+dfw, dfw )  # along dip (x direction)
    fy = np.arange( 0, Fl+dfl, dfl ) - Fl/2   # along strike (y direction)
    fx, fy = fx*1000, fy*1000
    fxx, fyy = np.meshgrid( fx, fy )
    sdep2d = fxx * np.sin( dip*np.pi/180. ) / 1000  # in km
    slon2d, slat2d = projection( fxx, fyy, **kwds )
    
    # surface projection (change in dip direction)
    fxS = fx * np.cos( dip*np.pi/180. )
    fxxS, fyy = np.meshgrid( fxS, fy )
    slon2dS, slat2dS = projection( fxxS, fyy, **kwds )

    # get the lon/lat location of hypocenter (slon, slat, hdep)
    hy = hypoAS*1000
    hx = hypoDD*1000
    slon, slat = projection( hx, hy, **kwds )
    hdep = hx*np.sin( dip*np.pi/180. ) / 1000. # (in km)
    
    # epicenter location
    hxS = hx*np.cos( dip*np.pi/180. )
    slonS, slatS = projection( hxS, hy, **kwds )

    # visual test:
    if VisualDict != None:
	print 'test plt'
	from mpl_toolkits.mplot3d import Axes3D

	rlon, rlat = VisualDict['SiteLoc']   # site locaitons for visual analysis
	savetype = VisualDict['savetype']
        
	plotpth = VisualDict['plotpth']

	# 1. plot surface projection
	fig = plt.figure(1) 
	ax = fig.add_subplot( 111 )
	ax.plot( rlon, rlat, 'k^' )
	ax.plot( slonS, slatS, 'r*', ms=10 )
	ax.plot( lon0, lat0, 'yo', ms=12 )
	ax.plot( [slon2dS[0,0], slon2dS[0,-1],slon2dS[0,-1], slon2dS[0,0], slon2dS[0,0]], \
		 [slat2dS[0,0], slat2dS[0,0],slat2dS[-1,0], slat2dS[-1,0], slat2dS[0,0]],'b' )
	ax.set_title( 'surface projection of fault surface (blue box)\nstation distribution (black triangles), yellew circle: origin' )
	ax.set_xlabel( 'longitude' )
	ax.set_ylabel( 'latitude' )

	plotn = 'SurfaceProjection.Scenario%s.FaultV%s.Station.Hypo.%s'%(ScenarioID, VID, savetype) 
	fig.savefig( plotpth + plotn, format=savetype )
        
	# 2. plot in 3D
	fig = plt.figure(2) 
	ax = Axes3D(fig) 
	
	ax.plot( rlon, rlat, np.zeros(len(rlon)), 'k^', ms=8 )
	ax.plot( [slon], [slat], [-hdep], 'r*', ms=12 )
	ax.plot( [lon0], [lat0], [0], 'yo', ms=12 )
	
	linec = ax.plot_wireframe( slon2d, slat2d, -sdep2d )
	linec.set_color('b')
	
	ax.set_zlim3d(-50, 0) 
	ax.set_xlabel( 'lon' )
	ax.set_ylabel( 'lat' )
	ax.set_zlabel( 'depth (km)' )
	
	plotn = 'ThreeDimension.Scenario%s.FaultV%s.Station.Hypo.%s'%(ScenarioID, VID, savetype) 
	fig.savefig( plotpth + plotn, format=savetype )
    
    OutputDict = {}
    OutputDict['FaultPlane'] = slon2d.tolist(), slat2d.tolist(), sdep2d.tolist()
    OutputDict['FaultSurface'] = slon2dS.tolist(), slat2dS.tolist()
    OutputDict['HypoLoc'] = slon,slat,hdep 
    OutputDict['EpicLoc'] = slonS, slatS 
    
    return OutputDict
