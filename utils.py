#!/usr/bin/env python
"""
Utilities used in NGA classes
"""
import os
import numpy as np


def cpt_sqrt( a,b ):
    a = np.array( a )
    b = np.array( b )
    return np.sqrt( a**2 + b**2 )

def logline(x1,x2,y1,y2,x):
    # linear interpolation
    k = (y2-y1)/(x2-x1)
    C = y1-k*x1
    y = k*x+C
    return y

def GetKey(key):
    return '%3.2f'%(key)


# save and load (useful metadata operations
class namespace:
    """
    Namespace with object attributes initialized from a dict
    Turn Dictionary keys into object attributes
    d['KeyName'] -> d.KeyName
    d is a dictioanry
    """
    def __init__( self, d ):
	self.__dict__.update(d) 

# load dictionary saved in the python files
def load( f, d=None ):
    """ 
    load variables from Python source files
    Input:
        f: file object (fid) 
	   or filename with full path (string)
	d: dictionary ( None default )
    Output: 
        namespace(d) : updated dictionary
    """
    if type(f) is not file:
	f = open( os.path.expanduser(f) )   # get the file object
    if d is None:
	d = {}
    exec f in d
    return namespace( d ) 


# save dictionary into python file (used in metadata manipulation)
def save( fd, d, expand=None, keep=None, header=''):
    """
    Write variables from a dict into a Python source file.
    """
    if type( d ) is not dict:
        d = d.__dict__

    if expand is None:
        expand = []
    out = header
    for k in sorted( d ):
	if k not in expand and (keep is None or k in keep):
	    out += '%s = %r\n' % (k, d[k])

    for k in expand:
        print 'test expand'
	if k in d:
            if type( d[k] ) is tuple:
                out += k + ' = (\n'
                for item in d[k]:
                    out += '    %r,\n' % (item,)
                out += ')\n'
            elif type( d[k] ) is list:
                out += k + ' = [\n'
                for item in d[k]:
                    out += '    %r,\n' % (item,)
                out += ']\n'
            elif type( d[k] ) is dict:
                out += k + ' = {\n'
                for item in sorted( d[k] ):
                    out += '    %r: %r,\n' % (item, d[k][item])
                out += '}\n'
            else:
                sys.exit( 'Cannot expand %s type %s' % ( k, type( d[k] ) ) )
    if fd is not None:
        if type( fd ) is not file:
            fd = open( os.path.expanduser( fd ), 'w' )
        fd.write( out )
    return out


# this function will be used a lot (general map function as in Python)
def mapfunc(func,*args,**kwds):
    """
    Modified function map tool
    Account for the single argument and multiple arguments
    Account for the keywords input
    """
    # arguments
    na = len(args)    # number of arguments
    args0 = {}
    for ina in xrange( na ):
	key1 = '%s'%ina
	args0[key1] = {}
	try:
	    tmp = len(args[ina])
	    if tmp == 1:
		# [1,], 'a', ['AB',],{}
		key2 = '%s'%(tmp-1)
		args0[key1][key2] = args[ina][tmp-1]
	    else:
		if isinstance( args[ina], str ):
		    # 'AB' case
		    key2 = '%s'%0
		    args0[key1][key2] = args[ina]
		else:
		    # [1,2,...],['a','b',...],['AB','CD',...] case
		    for il in xrange( tmp ):
			key2 = '%s'%il
			args0[key1][key2] = args[ina][il]
	except:
	    # single number as input
	    key2 = '%s'%0
	    args0[key1][key2] = args[ina]

    del(args)

    # keywords
    keys = kwds.keys()
    nk = len( keys )
     
    if nk != 0:
	# Has keywords input
	kwds0 = {}
	for ink,key1 in enumerate( keys ):
	    kwds0[key1] = {}
	    try: 
		tmp = len( kwds[key1] )   # elements each keyword has
		if tmp == 1:
		    # [1,], 'a', ['AB',]
		    key2 = '%s'%(tmp-1)
		    kwds0[key1][key2] = kwds[key1][tmp-1]
		else:
		    if isinstance( kwds[key1], str ):
			# 'AB' case
			key2 = '%s'%0
			kwds0[key1][key2] = kwds[key1]
		    else:
			# [1,2,...],['a','b',...],['AB','CD',...] case
			for il in xrange( tmp ):
			    key2 = '%s'%il
			    kwds0[key1][key2] = kwds[key1][il]
	    except:
		# single number as input
		key2 = '%s'%0
		kwds0[key1][key2] = kwds[key1]
        del( kwds )

	# get the maximum list length
	nl = 0; nla = 0; nlk = 0
	for ina in xrange( na ):
	    key1 = '%s'%ina
	    nla0 = len( args0[key1].keys() )
	    if nla0 >= nla:
		nla = nla0
	for ink in xrange( nk ):
	    key1 = keys[ink]
	    nlk0 = len( kwds0[key1].keys() )
	    if nlk0 >= nlk:
		nlk = nlk0
	nl = max( nlk, nla )
	
	# check input args and kwds
	for ina in xrange( na ):
	    key1 = '%s'%ina
	    nl0 = len(args0[key1].keys())
	    if nl0 != 1 and nl0 < nl:
		print 'input argument length error!'
		raise ValueError
	    
	for ink in xrange( nk ):
	    key1 = keys[ink]
	    nl0k = len(kwds0[key1].keys())  # number of elements for each arguments (for map)
	    if nl0k != 1 and nl0k < nl:
		print 'input kwds element length error!'
		raise ValueError

	# map function
	value = []
	for il in xrange( nl ):
	    arg0 = []; kwd0 = {}
	    for ina in xrange( na ):
		key1 = '%s'%ina
		nl0 = len( args0[key1].keys() )
		if nl0 == 1:
		    element = args0[key1]['0']
		else:
		    key2 = '%s'%il
		    element = args0[key1][key2]
		arg0.append(element)
	    for ink in xrange( nk ):
		nlk = len(kwds0[keys[ink]])  # number of elements for each arguments (for map)
		key1 = keys[ink]
		if nlk == 1:
		    kwd0[key1] = kwds0[key1]['0']
		else:
		    key2 = '%s'%il
		    kwd0[key1] = kwds0[key1][key2]
	    value.append( func( *arg0, **kwd0 ) )
    
    else:
	# No keywords input (use the default of the original function)
	nl = 0 
	for ina in xrange( na ):
	    key1 = '%s'%ina
	    nl0 = len( args0[key1].keys() )
	    if nl0 >= nl:
		nl = nl0 

	# check input args 
	for ina in xrange( na ):
	    key1 = '%s'%ina
	    nl0 = len(args0[key1].keys())
	    if nl0 != 1 and nl0 < nl:
		print 'input argument length error!'
		raise ValueError
	    
	# map function
	value = []
	for il in xrange( nl ):
	    arg0 = []; kwd0 = {}
	    for ina in xrange( na ):
		key1 = '%s'%ina
		nl0 = len( args0[key1].keys() )
		if nl0 == 1:
		    element = args0[key1]['0']
		else:
		    key2 = '%s'%il
		    element = args0[key1][key2]
		arg0.append(element)

	    value.append( func( *arg0 ) )

    return value   # return is a list even taking single number for each input (attention)


# NGA database related 
def RakeBin(rakes):
    # rakes in degree, list
       
    # rake: [-180,180]
    # 0: strike-slip, [-180,-150], [-30,30], [150,180]
    # 1: normal, [-120,-60]
    # 2: reverse, [60,120]
    # 3: reverse-oblique, [30,60], [120, 150]
    # 4: Normal-oblique, [-150,-120], [-60, -30]
    # These rules come from NGA flatfile
    group = {}
    groupnames = {'U':['k','Unknown'],'SS':['r','Strike-Slip'],'NM':['g','Normal'],\
	          'RV':['b','Reverse'],'NO':['#808080','Normal-Oblique'],'RO':['m','Reverse-Oblique']}
    for ig,groupname in enumerate( groupnames.keys() ):
	group[groupname] = []

    for ir in xrange( len(rakes) ):
	rake = rakes[ir]
	
	if rake>180. or rake < -180. or rake == None:
	    group['U'].append( rake )

	if -180<= rake <= -150 or -30<=rake<=30 or 150 <= rake <= 180:
	    group['SS'].append( rake )
	
	if -120<=rake<=-60:
	    group['NM'].append( rake )
	
	if 60 <= rake <= 120:
	    group['RV'].append( rake )
	
	if 30<rake<60 or 120<rake<150:
	    group['RO'].append( rake )
	
	if -150<rake<-120 or -60 < rake < -30:
	    group['NO'].append( rake )
    
    return group, groupnames

def Vs30Bin(Vs30s):
    # Vs30s in m/s, list

    # A: => 1500
    # B: [760, 1500)
    # C: [360, 760)
    # D: [180, 360)
    # E: < 180
    # This rules come from NGA flatfile   
    group = {}
    groupnames = {'A':['k','Hard Rock'],'B':['r','Rock'],'C':['g','Dense Soil and Soft Rock'],\
	          'D':['b','Stiff Soil'],'E':['m','Soft Soil']}
    for ikey, key in enumerate( groupnames.keys() ):
	group[key] = []
    
    for iv in xrange( len( Vs30s ) ):
	Vs30 = Vs30s[iv] 
	if Vs30 >= 1500.:
	    group['A'].append( Vs30 )
	if 760. <= Vs30 < 1500.:
	    group['B'].append( Vs30 )
	if 360. <= Vs30 < 760.:
	    group['C'].append( Vs30 )
	if 180. <= Vs30 < 360.:
	    group['D'].append( Vs30 )
	if Vs30 < 180.:
	    group['E'].append( Vs30 )
    return group, groupnames

# =========================================================================================
# Functions to compute exploratory variables
# =========================================================================================
# 1. Source related 
def rake2ftype_BA(rake):
    if rake == None:
	ftype = 'U'
    if rake != None:
	if -30. <= rake <= 30. or 150. <= rake <= 180. or -180. <= rake <=-150.:
	    ftype = 'SS' # strike-slip
	elif 30. < rake < 150.:
	    ftype = 'RS' # reverse
	elif -150. <= rake <= -30.:
	    ftype = 'NS' # normal
	else:
	    print 'Wrong rake angle!'
	    raise ValueError
    return ftype

def rake2ftype_CB(rake):
    Frv = 0; Fnm = 0
    if 30 < rake < 150:
	Frv=1
    if -150 < rake < -30:
	Fnm=1
    return Frv, Fnm

def rake2ftype_CY(rake):
    Frv, Fnm = 0, 0
    if 30 <= rake <= 150:
	Frv = 1
    elif -120 <= rake <= -60:
	Fnm = 1
    return Frv, Fnm

def rake2ftype_AS(rake):
    Frv, Fnm = 0, 0
    if 30 <= rake <= 150:
	Frv = 1
    elif -120 <= rake <= -60:
	Fnm = 1
    return Frv, Fnm


def calc_dip( rake ):
    """
    Empirical determination of dip angle from the faulting style
    Input: 
        rake in degree (-180<=rake<=180)
    Output: 
        dip angle in degree
    """
    if abs(rake) > 180: 
	print 'rake angle should be within -180 and 180'
	raise ValueError
    
    if abs(rake)<=30 or abs(rake)>=150:
	dip = 90
    elif -150 < rake < -30: 
	dip = 50
    elif 30 < rake < 150:
	dip = 40
    return dip 

def calc_Zhypo(M,rake):
    """
    Compute Zhypo from empirical relations 
    When Ztor is unknown from input models
    """
    if M < 0:
	print 'Magnitude should be larger than 0'
	raise ValueError
    if abs(rake) > 180:
	print 'rake angle should be within -180 and 180'                     
	raise ValueError
    
    if abs(rake) < 30 or abs(rake) > 150:
	# strike-slip
	Zhypo = 5.63 + 0.68 * M 
    else:
	Zhypo = 11.24 - 0.2 * M
    return Zhypo

def calc_W(M,rake):
    """
    Compute fault width when not specified by input
    """
    if M < 0:
	print 'Magnitude should be larger than 0'
	raise ValueError
    
    # In R
    if abs(rake) > 180:
	print 'rake angle should be within -180 and 180'                  
	raise ValueError
    if abs(rake) < 30 or abs( rake ) > 150:
	W = 10 ** (-0.76+0.27*M)
    elif -150 <= rake <= -30:
	W = 10 ** (-1.14+0.35*M)
    elif 30 <= rake <= 150: 
	W = 10 ** (-1.61+0.41*M)
    
    # In Matlab
    #W = 10**(-1.01+0.32*M)

    return W

def calc_Ztor(W,dip,Zhypo):
    """
    Compute Ztor if not specified by input
        dip should be in degree
    """
    if dip <= 0 or dip > 90:
	print 'dip angle should be with in (0,90]'
	raise ValueError
    if W <= 0:
	print 'Fault width should be larger than 0'
	raise ValueError
    if Zhypo < 0:
	print 'Zhypo should be larger than 0'
	raise ValueError
    Ztor = max( Zhypo-0.6*W*np.sin( dip * np.pi/ 180 ), 0 )
    return Ztor


# 2. Source-Site related (distances)
def calc_Rx(Rjb, Ztor, W, dip, azimuth, Rrup=None):
    """
    Compute distance parameter Rx from other inputs
    """
    if Rjb < 0:
	print 'Joyer-Boore distance Rjb should be larger than 0'
	raise ValueError
    if Ztor < 0:
	print 'Ztor should be larger than 0'
	raise ValueError
    if W <= 0:
	print 'Fault width should be larger than 0'
	raise ValueError
    if dip<=0 or dip > 90:
	print 'dip angle should be (0,90]'
	raise ValueError
    if abs( azimuth ) > 180.0:
	print 'azimuth should be width in -180.0 and 180.0'
	raise ValueError

    d = dip * np.pi/180    # degree to radius 
    a = azimuth * np.pi/180
    if dip != 90:
	if azimuth > 0: 
	    if azimuth == 90: 
		if Rjb == 0: 
		    if Rrup != None:
			if Rrup < Ztor / np.cos(d):
			    Rx = np.sqrt( Rrup**2 - Ztor**2 )
			else:
			    Rx = Rrup / np.sin( d ) - Ztor / np.tan( d )
		    else:
			Rx = W * np.cos(d)/2.   
			# empirical relation  (Rrup is easy to compute)
			# assume that the site is located at the center of the surface projection of the rupture plane
		else: 
		    Rx = Rjb + W * np.cos(d) 
	    else: 
		if Rjb*abs(np.tan(a)) <= W*np.cos(d):
		    Rx = Rjb * abs( np.tan(a) )
		else:
		    Rx = Rjb * np.tan(a) * np.cos( a-np.arcsin(W*np.cos(d)*np.cos(a)/Rjb) )
	else: 
	    Rx = Rjb * np.sin(a)
    else:
	Rx = Rjb * np.sin(a)
    
    return Rx


def calc_Rrup( Rx, Ztor, W, dip, azimuth, Rjb=None ):
    """
    Compute the closest distance from site the the surface the fault
    """
    if Ztor < 0:
	print 'Ztor should be larger than 0'
	raise ValueError
    if W <= 0:
	print 'Fault width should be larger than 0'
	raise ValueError
    if dip<=0 or dip > 90:
	print 'dip angle should be (0,90]'
	raise ValueError
    if abs( azimuth ) > 180:
	print 'azimuth should be width in -180 and 180'
	raise ValueError
    
    d = dip * np.pi/180    # degree to radius 
    a = azimuth * np.pi/180
   
    if dip == 90 and Rjb != None:
	Rrup = np.sqrt( Rjb**2 + Ztor**2 )
        return Rrup

    if dip != 90:
	if Rx < Ztor * np.tan( d ):
	    Rrup1 = np.sqrt( Rx**2 + Ztor**2 )
	elif Rx >= Ztor*np.tan(d) and Rx <= Ztor*np.tan(d)+W/np.cos(d):
	    Rrup1 = Rx*np.sin(d) + Ztor*np.cos(d)
	elif Rx > Ztor*np.tan(d) + W/np.cos(d):
	    Rrup1 = np.sqrt( (Rx-W*np.cos(d))**2 + (Ztor+W*np.sin(d))**2 )
    elif dip == 90:
	Rrup1 = sqrt( Rx**2 + Ztor**2 )

    if azimuth == 90 or azimuth == -90:
	Ry = 0
    elif azimuth == 0 or azimuth == 180 or azimuth == -180:
	if Rjb == None:
	    print 'Rjb cannot be None in the case azimuth == 0 or azimuth == 180 or azimuth == -180'
	    raise ValueError
	else:
	    Ry = Rjb
    else:
	Ry = abs( Rx/np.tan(a) )
    Rrup = np.sqrt( Rrup1**2 + Ry**2 )
    return Rrup

# Give fault geometry (fault plane coordinates) and site locations
# compute Rjb, azimuth, and then Rrup
# from those, computing Rx by the function defined above

# using fault plane as the xy coordinate is the key 
def projection(x,y,**kwds):
    """ 
    Projection of lon/lat to UTM or reverse direction
    input:
    x,y ( lon/lat or x/y )
    kwds: zone, origin, rot, inverse
    
    output:
    x,y ( x/y or lon/lat )
    """

    import pyproj

    zone = kwds['zone']
    origin = kwds['origin']
    rot = kwds['rot']
    inverse = kwds['inverse']

    if origin == None:
	return

    # geometrical origin
    x0 = origin[0]; y0 = origin[1]
    
    rot = rot*np.pi/180.
    c,s = np.cos(rot),np.sin(rot)
    
    x = np.array(x,'f')
    y = np.array(y,'f')
    
    # you can use other projections (modify here)
    proj = pyproj.Proj(proj='utm',zone=zone,ellps='WGS84')
    
    if inverse:
	x0,y0 = proj(x0,y0,inverse=False)
	x,y = c*x-s*y, s*x+c*y
	x,y = x+x0,y+y0
	x,y = proj(x,y,inverse=True)
    
    else:

	x0,y0 = proj(x0,y0,inverse=False)
	x,y = proj(x,y,inverse=False)
	x,y = x-x0,y-y0
	x,y = x*c+y*s, -s*x+c*y

    return x,y


def calc_distances(SiteGeo, Dims, Mech, ProjDict, Rrup=False, Rx=False):
    """
    Compute Rjb, Rrup, Rx given fault geometry and site location
    For one site
    FaultGeo: 3*NfaultPoints (slon,slat,sdep)
    SiteGeo: (rlon,rlat)
    """
    # UTM projection property
    lon0, lat0 = ProjDict['origin']   # projection origin
    ProjDict['inverse'] = False  # from ll to xy
    kwds = ProjDict
    
    # sites and compute the azimuth
    rlon, rlat = SiteGeo 
    rx, ry = projection( rlon, rlat, **kwds )
    
    # Fault dimension and focal mechanism 
    Fl,dfl,Fw,dfw,ztor = Dims   # fault length along strike, fault width down dip 
    strike, dip, rake = Mech   # fault mechanism 
    
    # create fault mesh:
    ProjDict['inverse'] = True  # from xy to ll
    kwds = ProjDict
	# along strike direction: y-axis
	# along dip direction: x-axis
    fx = np.arange( 0, Fw+dfw, dfw )  # along dip (x direction)
    fy = np.arange( 0, Fl+dfl, dfl ) - Fl/2   # along strike (y direction)
    fx, fy = fx*1000, fy*1000
    fxx, fyy = np.meshgrid( fx, fy )
    fzz = fxx * np.sin( dip * np.pi/180. )   # in meter
    sdep2d = fzz / 1000  # in km
    slon2d, slat2d = projection( fxx, fyy, **kwds )

    # surface projection (change in dip direction)
    fxS = fx * np.cos( dip*np.pi/180. )
    fxxS, fyy = np.meshgrid( fxS, fy )
    slon2dS, slat2dS = projection( fxxS, fyy, **kwds )

    Nlat, Nlon = slon2d.shape
    Nloc = Nlat * Nlon 
    fxxS1d = fxxS.reshape( Nloc ) / 1000.
    fxx1d = fxx.reshape( Nloc ) /1000.
    fyy1d = fyy.reshape( Nloc ) /1000.
    fzz1d = fzz.reshape( Nloc ) /1000.

    # compute Rjb using fault and  site locations and get the azimuth of all sites for later use to compute Rx
    Nsta = len(rlon)
    Rjb = []; Rrup = [];  azimuth = []
    print 'Computing Rjb, and Rrup, and azimuth...'
    for ista in xrange( Nsta ):
	rx0 = rx[ista] / 1000.
	ry0 = ry[ista] / 1000.
	
	# Rjb
	if fxxS1d.min() <= rx0 <= fxxS1d.max() and fyy1d.min() <= ry0 <= fyy1d.max(): 
	    Rjb.append( 0 )
	else:
	    distS = []
	    for iloc in xrange( Nloc ):
		dx = fxxS1d[iloc] - rx0
		dy = fyy1d[iloc] - ry0 
		distS.append( np.sqrt( dx**2 + dy**2 ) )
	    Rjb.append( min(distS) )   
	
	# Rrup
	dist = []
	for iloc in xrange( Nloc ):
	    dx = fxx1d[iloc] - rx0 
	    dy = fyy1d[iloc] - ry0 
	    dist.append( np.sqrt( dx**2 + dy**2 +  fzz1d[iloc]**2 ) )

	Rrup.append( min(dist) )
	
	# compute azimuth (different from common used) 
	# refer to James's defination:
	# The angle between positive fault strike direction and line connecting a site to the closest point on
	# the surface projection of the TOP EDGE of rupture, which clockwise assumed positive !
	# fy: the surface projection of the top edge of rupture 
        
	# different cases
	fymin = fy.min()/1000.; fymax= fy.max()/1000. 
	if fymin <= ry0 <= fymax: 
	    azimuth0 = -np.pi/2. * (rx0<0.0) + np.pi/2 * (rx0>=0.0)
	
	if ry0 > fymax: 
	    dy = rx0 - 0.0
	    dy = ry0 - fymax
	    if rx0 > 0.0: 
		azimuth0 = np.arctan( dx/dy )
	    elif rx0 < 0.0: 
		azimuth0 = -np.arctan( -dx/dy )
	    elif rx0 == 0.0: 
		azimuth0 = 0.0
	if ry0 < fymin: 
	    dx = rx0 - 0.0 
	    dy = fymin - ry0 
	    if rx0 > 0.0: 
		azimuth0 = np.pi - np.arctan( dx / dy )
	    elif rx0 < 0.0: 
		azimuth0 = np.arctan( -dx/dy ) - np.pi
	    elif rx0 == 0.0: 
		azimuth0 = np.pi
        
        azimuth.append( azimuth0*180./np.pi )

    if 0:
	# test projection (basis for Rjb and azimuth calculation0
	import matplotlib.pyplot as plt 
	fig = plt.figure(1) 
	ax = fig.add_subplot( 221 )
	ax.plot(rx,ry,'k^')
	ax.plot( [fxS[0],fxS[-1],fxS[-1],fxS[0],fxS[0]],[fy[0],fy[0],fy[-1],fy[-1],fy[0]],'b' )
	ax = fig.add_subplot( 222 )
	ax.plot( azimuth, 'k.' )
	ax = fig.add_subplot( 223 )
	ax.plot( azimuth, Rjb, 'k.' )
	ax = fig.add_subplot( 224 )
	ax.plot( azimuth, Rrup, 'k.' )
	plt.show()
    
    # Compute Rx from Rjb and Rrup
    Rx = mapfunc( calc_Rx, Rjb, ztor, Fw, dip, azimuth, Rrup=Rrup )
    
    OutputDict = {}
    OutputDict['Rjb'] = Rjb
    OutputDict['Rx'] = Rx
    OutputDict['Rrup'] = Rrup 
    OutputDict['azimuth'] = azimuth

    return OutputDict


# 3. Site related ( site condition ) Given Vs30 (inferred or measured)
def calc_Z1(Vs30,NGAmodel):
    """
    Compute Z1.0 parameter for AS and CY model if not specified
    Vs30 in m/s
    return in km
    """
    if Vs30 < 0:
	print ' Vs30 should be larger than 0'
	raise ValueError

    if NGAmodel == 'AS':
	if Vs30 < 180.:
	    Z10 = 6.745
	elif 180 <= Vs30 <= 500.:
	    Z10 = 6.745 - 1.35*np.log( Vs30/180. )
	elif Vs30 > 500:
	    Z10 = 5.394 - 4.48 * np.log( Vs30/500. )
    elif NGAmodel == 'CY':
	Z10 = 28.5 - 3.82/8 * np.log( Vs30**8 + 378.7**8 )
    return np.exp( Z10 )   # in meter


def calc_Z25( Vs30=None, Z1=None, Z15=None ):
    """
    input Vs30 m/s; Z1 in m; Z15 in m
    """

    if Vs30 == None and Z1== None and Z15==None:
	print 'Either Vs30, Z1.0 or Z1.5 should be specified'
	raise ValueError
    if Vs30 != None and Vs30 < 0:
	print 'Vs30 should be larger than 0'
	raise ValueError
    if Z1 != None and Z1 < 0:
	print 'Z1 should be larger than 0'
	raise ValueError
    if Z15 != None and Z15 < 0:
	print 'Z15 should be larger than 0'
	raise ValueError
    if Z15 != None:	
	Z25 = 636 + 1.549*Z15
    elif Z1 != None:
	Z25 = 519 + 3.595*Z1
    elif Vs30 != None:
	Z1 = calc_Z1( Vs30, 'AS' )
	Z25 = 519 + 3.595*Z1
    return Z25 / 1000.   # in km


# 4. inter- and intra- event residuals calculation using standard deivation and 
# residuals between the original observations and NGA predicted median values
def GetIntraInterResiduals(residualT, EQID, sigmaT, tau, sigma, AS=None):
    """
    Compute Normalized Total, Inter- and Intra residuals
    Input:
        residualT: un-normalized total residual
	EQID: earthquake id for each record (used for group bin)
	sigmaT: standard deviation (total residual) for each record
	tau: standard deviation (inter-event) for each record
	sigma: standard deviation (intra-event) for each record
	AS: None (just use average method): not None (Use AS 1992 method to compute)
    Output:
        epsilonT: normalized total residual
	eta: normalized inter-event residual
	epsilon: normalized intra-event residual
    """

    residualT = np.array( residualT )
    sigmaT = np.array( sigmaT )
    tau = np.array( tau )
    sigma = np.array( sigma )

    # group bin
    Events = group_list( EQID )
    Neq = len(Events)
    ID = []
    for ieq in xrange( Neq ):
	ID.append( Events[ieq][2] )

    # Compute residuals
    eta_EQinter = []
    for ieq in xrange( Neq ):
	index = (np.array(RID)==ID[ieq]).nonzero()[0]
	N_EQ = len(index)
        
	if AS != None:
	    tau0 = np.mean(tau[index]); 
	    sigma0 = np.mean( sigma[index] )                                            
	    residual_EQinter0 = tau0**2 * sum( residualT[index] ) / (N_EQ*tau0**2 + sigma0**2)
	else:
	    residual_EQinter0 = np.mean( residualT[index] )
	for isub in index:
	    eta_EQinter.append( residual_EQinter0 / tau[isub] )   # for each record in the event group

    eta = np.array( eta_EQinter )
    residual_EQintra = residual_total - eta * tau
    epsilon = residual_EQintra / sigma
    epsilonT = residual_total / sigmaT

    return epsilonT, eta, epsilon


