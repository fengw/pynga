#!/usr/bin/env python
"""
Class with methods for Broadband platform
"""
import os, sys
import numpy as np 

from pynga.utils import * 

# ===============
# General methods 
# ===============
def DegreeToRad(angle, inverse=False): 
    if inverse: 
	# rad to degree
	return angle * 180./np.pi 
    else: 
	# degree to rad
	return angle * np.pi / 180.


class BBP: 
    
    def __init__(self): 
	# initialize the class (path, basic parameters which will be used by self functions)
	pass 


    def srcFileRead(self, srcFile): 
	"""
	Read in *.src file and return 
	fault parameters for distance calculation
	"""
	# Read file
	lines = open( srcFile ).readlines()
	Values = []
	for il in xrange( len(lines) ):
	    spl = lines[il].strip().split()
	    Values.append( float( spl[2] ) )
	
	# Parameter list in *.src file 
	Magnitude, \
	FaultLength, Dlen, FaultWidth, Dwid, DepthToTop, \
	Strike, Rake, Dip, \
	Lat0, Lon0, \
	HypoAS, HypoDD, \
	DT, Seed, CornerFrequency = Values 

	Dims = FaultLength, Dlen, FaultWidth, Dwid, DepthToTop 
	Mech = Strike, Dip, Rake 
	OriginLoc = Lon0, Lat0, DepthTop
	
	if HypoAS >= 0: 
	    vector = [DegreeToRad(Strike), HypoAS, 0.0]
	    HypoLocAS = EndLocation( OriginLoc, vector ) 
	else: 
	    vector = [DegreeToRad(Strike)+np.pi, HypoAS, 0.0]
	    HypoLocAS = EndLocation( OriginLoc, vector ) 

        d = DegreeToRad( Dip )
	vector = [DegreeToRad(Strike)+np.pi/2., HypoDD*np.cos(d), HypoDD*np.sin(d)]
	HypoLoc =  EndLocation( HypoLocAS, vector )

	return Dims, Mech, OriginLoc[:2], HypoLoc
        

    def ReadSite(self): 
	pass

    def ComputeGMEPs(self): 
	pass 

    def CompareGMPEsandSimulation(self):
	pass
