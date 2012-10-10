#!/usr/bin/env python
"""
Utilities used in NGA classes
"""
import os, time
import numpy as np

# optional for test plots 
import matplotlib.pyplot as plt 

# ===================
# General Functions
# ===================
def cpt_sqrt( a,b ):
    a = np.array( a )
    b = np.array( b )
    return np.sqrt( a**2 + b**2 )


def RMScalc( V1,V2, Ratio=True ): 
    N = len(V1) 
    if len(V2) != N: 
	print 'length of array2 should be the same as array1'
	raise ValueError
    else:
	if 0 in V2: 
	    print 'elements in array2 should not be zeros'
	    raise ValueError
	if Ratio:
	    Error = (V1-V2)/V2
	else: 
	    Error = (V1-V2)
	RMS = np.sqrt( sum(Error**2) / N )
    return RMS 


def logline(x1,x2,y1,y2,x):
    # linear interpolation
    k = (y2-y1)/(x2-x1)
    C = y1-k*x1
    y = k*x+C
    return y


def GetKey(key):
    return '%.3f'%(key)


def HourMinSecToSec(BlockName=None):
    hour, min, sec = time.localtime()[3:6]
    sec1 = hour*60*60 + min*60 + sec
    if  BlockName != None:
	print '%s'%BlockName
    return sec1

def SecToHourMinSec(sec1,BlockName=None):
    hour = sec1//3600
    min = (sec1-hour*3600)//60
    sec = sec1-hour*3600-min*60
    if BlockName == None:
	BlockName = 'the above block'
    print 'Time cost of %s is %s hr %s min %s sec'%(BlockName,hour,min,sec)
    return hour,min,sec        


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


# geometrical projection (general)
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



# =========================================================================================
# NGA database related 
# =========================================================================================
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
# Functions to compute exploratory variables used in GMPEs
# =========================================================================================

# ==============
# Fault type
# ==============
def rake2ftype_BA(rake):
    if rake == None:
	ftype = 'U'
    if rake != None:
	if -30. <= rake <= 30. or 150. <= rake <= 180. or -180. <= rake <=-150.:
	    ftype = 'SS' # strike-slip
	elif 30. < rake < 150.:
	    ftype = 'RS' # reverse
	elif -150. <= rake <= -30.:
	    ftype = 'NM' # normal
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


# ===================
# Fault geometry
# ===================
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


# ======================================
# Fault-Site distances (Rjb, Rx, Rrup)
# ======================================
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


# Simplified fault geometry (one trace) (all of these geometrical relations satisfy in this simple case)
# so as the above two functions (simple fault geometry)
def calc_distances(SiteGeo, Dims, Mech, ProjDict, Rrup=False, Rx=False):
    """
    Compute Rjb, Rrup, Rx implicitly given fault geometry and site location (in lon/lat)
    The grid generation is needed to get the explicit fault geometry for further calculation
    For Fling study and broadband platform (*.src file)
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




# =======================================================================================
# Utilities for general distance calculations (in spherical coordinates and earth-flatten)
# =======================================================================================
R = 6371.   # Earth radius (km)
tol = 1e-10 
def LonLatToAngleDistance( loc1, loc2, \
	CalcRadius=True, \
	CalcDist=True, Fast=True, \
        CalcAzimuth=True, Azimuth0to2PI=False ): 
    """ 
    Convert lon/lat to radius (earth center) and/or azimuth and great circle Distances
    between points on the spherical surface
    Inputs: 
	Loc1: lon1, lat1 in degree, dep1 in km (Starting point in azimuth calculation)
	Loc2: lon2, lat2 in degree, dep2 in km 
    Outputs: 
	Radius: angle between two points (central angle)
	Azimuth1to2: azimuth (relative to north pole) from point 1 to point2
	horzDistance: distance between two points (great circle along spherical surface, km)
	vertDistance: distance between two points vertically
    Angles in radius
    # converted from OpenSHA jave to python (credit to OpenSHA developers)
    org.opensha.commons.geo.Location*
    and in the website: 
    http://en.wikipedia.org/wiki/Haversine_formula
    http://www.movable-type.co.uk/scripts/latlong.html
    """
    loc1 = np.array( loc1 ) 
    loc2 = np.array( loc2 ) 
    
    lon1, lat1 = loc1[:2] * np.pi/180. 
    lon2, lat2 = loc2[:2] * np.pi/180.
    
    # initialization
    Radius = 0.
    horzDistance = 0.
    verzDistance = 0. 
    Azimuth1to2 = 0. 

    if CalcRadius: 
	sinDlatBy2 = np.sin( (lat2-lat1)/2.0 )
	sinDlonBy2 = np.sin( (lon2-lon1)/2.0 )
	c = (sinDlatBy2**2) + np.cos(lat1)*np.cos(lat2) * sinDlonBy2**2 
	Radius = 2.0 * np.arctan2( np.sqrt(c), np.sqrt(1-c) )  # in rad (to keep angle in between pi and -pi) 
	# central angle from point1 to point2
    
    if CalcDist:
	if Fast: 
	    dlat = lat1 - lat2 
	    dlon = (lon1 - lon2) * np.cos( (lat1+lat2)/2.0 )
	    horzDistance = R * np.sqrt( dlat**2 + dlon**2 )
	else: 
	    sinDlatBy2 = np.sin( (lat2-lat1)/2.0 )
	    sinDlonBy2 = np.sin( (lon2-lon1)/2.0 )
	    c = (sinDlatBy2**2) + np.cos(lat1)*np.cos(lat2) * sinDlonBy2**2 
	    Radius = 2.0 * np.arctan2( np.sqrt(c), np.sqrt(1-c) ) 
	    horzDistance = R * Radius 
	
	dep1, dep2 = loc1[2], loc2[2] 
	verzDistance = dep1-dep2  
    
    if CalcAzimuth: 
	# calculate azimuth (p1 to p2 vector relative to p1 to north pole)  
	dlon = lon2-lon1 
	cosLat2 = np.cos(lat2) 
	y1 = np.sin(dlon) * cosLat2
	y2 = np.cos( lat1 )*np.sin(lat2) - np.sin(lat1)*cosLat2*np.cos(dlon) 
	Azimuth1to2 = np.arctan2( y1, y2 ) 
	if Azimuth0to2PI:
	    Azimuth1to2 = (Azimuth1to2+2*np.pi)%(2*np.pi) 
    
    return Radius, horzDistance, verzDistance, Azimuth1to2


# ptLineDist and ptSegDist are referred from java/awt/geom/Line2D.java.htm
def ptLineDist(x1,y1,x2,y2,px,py): 
    """
    Compute the point (px,py) to line (x1,y1)=>(x2,y2) allowing infinitely-extending of the line 
    Not used
    """ 
    # adjust vectors relative to point (x1,y1)
    x2 -= x1
    y2 -= y1 
    px -= x1 
    py -= y1 

    # 1. projection using dot production of adjusted vector (px,py) and (x2,y2)
    dotprod = ( px * x2 + py * y2  ) 
    projLenSq = dotprod * dotprod / (x2*x2+y2*y2)  # length of the vector (x1,y1)=>projected point of (px,py) on the line

    # 2. subtraction to get the closet distance (length of the vector)
    lenSq = px*px + py*py - projLenSq 
    if lenSq < 0: 
	# (px,py) is in the line specified by (x1,y1) and (x2,y2)
	lenSq = 0 
    return np.sqrt( lenSq )


def ptSegDist(x1,y1,x2,y2,px,py): 
    """
    Compute the point (px,py) to line (x1,y1)=>(x2,y2) without infinitely-extending of the line 
    Distance measured is the distance between the specified point and the closest point between 
    the specified end points (x1,y1) and (x2,y2)
    """ 
    # adjust vectors relative to point (x1,y1)
    x2 -= x1
    y2 -= y1 
    px -= x1 
    py -= y1 

    # 1. projection using dot production of adjusted vector (px,py) and (x2,y2)
    dotprod = ( px * x2 + py * y2  ) 
    if dotprod <= 0.0: 
	# (px,py) is on the side of (x1,y1) (projection of the point (px,py) is not on the segment)
	projLenSq = 0.0 
    else:
	# check the other side relationship
	px = x2-px
	py = y2-py 
        dotprod = px*x2 + py*y2 
	if dotprod <= 0.0: 
	    # (px,py) is on the side of (x2,y2) (projection of the point (px,py) is not on the segment)
	    projLenSq = 0.0 
	else: 
	    # point (px,py) 's projection is in between (x1,y1) and (x2,y2) 
	    # same as the ptLineDist function 
	    projLenSq = dotprod * dotprod / (x2*x2 + y2*y2)

    # 2. subtraction to get the closet distance (length of the vector)
    lenSq = px*px + py*py - projLenSq    # if projLenSq = 0.0, then the distance would be either original (px,py) to (x1,y1) or (x2,y2) 
    if lenSq < 0: 
	# (px,py) is in the line specified by (x1,y1) and (x2,y2)
	lenSq = 0.0
    return np.sqrt( lenSq  )


def distToLine(loc1, loc2, loc3, Fast=False):
    """
    Compute the shortest distance between a point (loc3) and a line (great circle) 
    that extends infinitely in both directiions. Depth is ignored.
    refer to: http://williams.best.vwh.net/avform.htm#XTE
    # this distance could be postive or negative
    """
    lon1, lat1 = loc1[:2] * np.pi/180. 
    lon2, lat2 = loc2[:2] * np.pi/180.
    lon3, lat3 = loc3[:2] * np.pi/180.
    
    if Fast: 
	# with earth-flatten approximation (faster) 
	# used in shorter distance (<=200km)
	lonScale = np.cos( 0.5*lat3 + 0.25*lat1 + 0.25*lat2 )   # earth-flatten approximation factor
	x2 = (lon2-lon1)*lonScale 
	y2 = lat2-lat1
	x3 = (lon3-lon1)*lonScale
	y3 = lat3 - lat1
	# x1=y1 = 0 
	Term1 = (x2*(-y3)-(-x3)*y2)/np.sqrt(x2**2+y2**2) 
	# originally, Term1 = 
	# abs( (x3-x1)*(y2-y1) - (y3-y1)*(x2-x1) ) / np.sqrt((x2-x1)**2+(y2-y1)**2.) 
	# for x1=y1=0, Term1 = abs( x3*y2 - y3*x2 ) / [] = abs( -y3*x2 - (-x3)*y2 ) / []
	# but here, Term1 has sign which indicates which side of point3 is located relative to the line vector (1to2)
	# +: right; -:left
	return Term1 * R 

    else:
	# orignial method to compute the distance from a point to a line (spherical trigonometry)
	# sin(A)/sin(a) = sin(B)/sin(b) = sin(C)/sin(c) 
	# A, B, and C: angle between surface circle 
	# a, b, and c: center angle of surface circle (great)
	a13,hD,vD,az13 = LonLatToAngleDistance(loc1,loc3,CalcRadius=True,CalcDist=False,CalcAzimuth=True)
	a12,hD,vD,az12 = LonLatToAngleDistance(loc1,loc2,CalcRadius=True,CalcDist=False,CalcAzimuth=True)
	Daz13az12 = az13-az12 
	xtd = np.arcsin( np.sin(a13)*np.sin(Daz13az12) ) 
	if abs(xtd) < tol: 
	    return 0.0   # point3 is on the line 
	else: 
	    return xtd * R     # xtd could >0 or <0 to identify the location of the point3 relative to the line
        # you could use this to compute Rx without extend your fault trace to infinite, but it takes time


def distToLineSeg( loc1,loc2,loc3, Fast=False ): 
    """
    Compute distance between point3 and line defined by point1 and point2 
    loc1, loc2, loc3 are list with 3 elements
    There are three cases: the projection point of loc3 on: 
    the loc1 side, on the loc2 side, in between loc1 and loc2
    """
    loc1 = np.array( loc1 ) 
    loc2 = np.array( loc2 ) 
    loc3 = np.array( loc3 )
    lon1, lat1 = loc1[:2] * np.pi/180. 
    lon2, lat2 = loc2[:2] * np.pi/180.
    lon3, lat3 = loc3[:2] * np.pi/180.
    
    if Fast: 
	lonScale = np.cos( 0.5*lat3 + 0.25*lat1 + 0.25*lat2 ) 
	x2 = (lon2-lon1)*lonScale 
	y2 = lat2-lat1
	x3 = (lon3-lon1)*lonScale
	y3 = lat3 - lat1
	Term1 = ptSegDist( 0, 0, x2, y2, x3, y3 )   # always positive
	return Term1 * R   

    else:
	# use cos(c) = cos(a)cos(b) + sin(a)sin(b)cos(C) to get
	a13, hD13, vD13, az13 = LonLatToAngleDistance(loc1,loc3,Fast=Fast)
	a12, hD12, vD12, az12 = LonLatToAngleDistance(loc1,loc2,Fast=Fast)
	Daz13az12 = az13-az12 

        # cross-track distance (in radius)
	xtd = np.arcsin( np.sin(a13)*np.sin(Daz13az12) ) 

	# along-track distance (in km )
	atd = np.arccos( np.cos(a13)/np.cos(xtd) ) * R 
	a23, hD23, vD23, az23 = LonLatToAngleDistance(loc2,loc3,CalcRadius=False,CalcDist=True,Fast=Fast,CalcAzimuth=False)

	# check if beyond p3 (should be p2?) (different from the original Rx definition?)
	if atd > hD12: 
	    print 'Beyond p2'
	    return hD23 
	
	# check if beyond p1 
	if np.cos(Daz13az12)<0:
	    print 'Beyond p1'
	    return hD13 

        # projection of the point is within the two points
	if abs(xtd) < tol: 
	    return 0.0   # point3 is on the line 
	else: 
	    return abs(xtd) * R   


def minDistToLineSeg( loc, segs, Fast=False ): 
    """
    Compute minimum distance between loc and a line made of segments
    Segments in line are contrained by two points loc1, loc2
    """
    Npoints = len(segs) 
    minDist = 1000.
    for iseg in range(1,Npoints): 
	p1 = segs[iseg-1] 
	p2 = segs[iseg]
	dist = abs( distToLineSeg(p1,p2,loc,Fast=Fast) ) 
	if dist <= minDist: 
	    minDist = dist 
    return minDist  


def EndLocation(loc1, vector): 
    """
    Given Vector and its starting point, find the end point of the vector, where
    vector has information (azimuth,horizontal distance and vertical distance) 
    # converted from OpenSHA jave to python (credit to OpenSHA developers)
    org.opensha.commons.geo.Location*
    Line extension (limited, not to infinite)
    """
    loc1 = np.array( loc1 )
    loc1[:2] = loc1[:2] * np.pi / 180.

    lon1, lat1, dep1 = loc1 
    az, DH, DV = vector    # az in radius 
    
    # refer to http://williams.best.vwh.net/avform.htm#LL 
    sinLat1 = np.sin(lat1) 
    cosLat1 = np.cos(lat1) 
    ad = DH / R 
    sinD = np.sin(ad) 
    cosD = np.cos(ad) 

    # compute the location information for the end point loc2 
    lat2 = np.arcsin( sinLat1*cosD + cosLat1*sinD*np.cos(az) ) 
    dlon = np.arctan2(np.sin(az)*sinD*cosLat1, cosD-sinLat1*np.sin(lat2)) 
    lon2 = lon1 + dlon
    dep2 = dep1 + DV

    lat2 = lat2*180./np.pi 
    lon2 = lon2*180./np.pi
    
    return [lon2,lat2,dep2] 


def CheckPointInPolygon(point, verts):
    """
    check whether a point is inside a polygon (general coordiate and convex or non-convex) 
    refer to: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    inputs:
	point: test point (x,y,[z]), z could be None
	verts: points that define the polygon shape (the last one and the first one are not the same!)
    return: True or False
    """
    Ndim = len(point) 
    verts = np.array( verts )
    dim = verts.shape[1] 
    if Ndim != dim or Ndim != 2 or dim != 2: 
	print 'point and shape should be defined with two coordinates'
	raise ValueError
    
    # test point
    testx = point[0]
    testy = point[1] 
    
    vertx = verts[:,0] 
    verty = verts[:,1] 
    nvert = len(vertx) 
    check = False
    j = nvert - 1 
    for i in xrange( nvert ): 
	c1 = verty[i]>testy
	c2 = verty[j]>testy 
	factor = (vertx[j]-vertx[i])*(testy-verty[i])/(verty[j]-verty[i]) + vertx[i] 
	if c1 != c2 and testx < factor: 
	    check = not check 
	j = i    # edge is defined from j to i

    return check


def getDistances(SiteGeo, FaultGeo, Fast = True):
    """
    Compute Rjb, Rrup, Rx explicitly given fault geometry and site location (in lon/lat)
	Rx just use the fault trace along-strike and azimuth between the two points at the two ends of the fault and 
	   site location.
    Inputs: 
        FaultGeo: faultGeo
		list faultGeo has dim: (Nrow,Ncol,3), fault1 surface discretization: (Nrow,Ncol)
		Nrow: down-dip direction grid points; Ncol: along-strike direction grid points
        AveDip: average dip (to determine the points to use)
	SiteGeo: siteGeo
		list siteGeo has elements: rlon,rlat,rdep to generally specify the location of site 
    Outputs: 
        Rjb, Rrup, Rx 
	They all have the following shape (site-based): 
    """
    FaultGeo = np.array( FaultGeo ) 
    SiteGeo = np.array( SiteGeo ) 
    loc1 = SiteGeo

    minRjb = 1000.    # in km 
    minRrup = 1000.   # in km 
    minRx = 1000.     # in km 
    
    Nrow, Ncol, Nelm = FaultGeo.shape
    
    surf = FaultGeo.reshape((Nrow*Ncol,Nelm))
    for loc2 in surf:
	alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2, CalcRadius=False, CalcDist=True, Fast=Fast, CalcAzimuth=False) 
	if hD <= minRjb: 
	    minRjb = hD
	totalDist = np.sqrt( hD**2 + vD**2 )
	if totalDist <= minRrup: 
	    minRrup = totalDist
    
    Rrup = minRrup 

    # check site is within the surface projection of the fault 
    verts = []
    irow = 0
    for icol in xrange( Ncol ): 
	verts.append( FaultGeo[irow,icol][:2].tolist() )
    icol = 0
    for irow in xrange( Nrow ): 
	verts.append( FaultGeo[irow,icol][:2].tolist() )
    irow = Nrow-1
    for icol in xrange( Ncol ): 
	verts.append( FaultGeo[irow,icol][:2].tolist() )
    icol = Ncol-1 
    for irow in xrange( Nrow ): 
	verts.append( FaultGeo[irow,icol][:2].tolist() )
    check = CheckPointInPolygon( SiteGeo[:2], verts )
    if check: 
	Rjb = 0.0
    else: 
	Rjb = minRjb 
    
    if 0: 
	DDtmp = 0
	for irow in xrange( Nrow-1 ): 
	    loc1 = FaultGeo[irow,0] 
	    loc2 = FaultGeo[irow+1,0] 
	    alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2, CalcRadius=False, CalcDist=True, Fast=Fast) 
	    DDtmp += DDtmp + hD 
	DownDipGridSpace = DDtmp / (Nrow-1)
	AStmp = 0
	for icol in xrange( Ncol-1 ): 
	    loc1 = FaultGeo[0,icol] 
	    loc2 = FaultGeo[0,icol+1] 
	    alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2, CalcRadius=False, CalcDist=True, Fast=True) 
	    AStmp += AStmp + hD 
	AlongStrikeGridSpace = AStmp / (Ncol-1) 
	AveGridSpace = (DownDipGridSpace+AlongStrikeGridSpace)/2.0
	if (minRjb < AveGridSpace): 
	    Rjb = 0.0
	else: 
	    Rjb = minRjb 

    # ============
    # compute Rx by extending fault trace and fault surface projection
    # ============
    FaultTrace = FaultGeo[0] # first row to get the fault trace for Rx calculation 
    Ntrace = FaultTrace.shape[0]    # = Ncol 
    ps = FaultTrace[0]
    pe = FaultTrace[-1] 
    
    Radius, hD, vD, Azimuth1to2 = LonLatToAngleDistance(ps, pe, CalcRadius=False, CalcDist=True, Fast=Fast, CalcAzimuth=True, Azimuth0to2PI=True ) 
    vector = [Azimuth1to2, hD, vD]
    
    # define the extended fault traces (to project faults surface projections on!)
    vector[1] = 1000  # km 
    vector[2] = 0
    
    # problem here !!!
    Loc1 = EndLocation( ps, vector ) 
    vector[0] = Azimuth1to2 + np.pi   # flip over trace dir 
    Loc2 = EndLocation( pe, vector ) 

    # this step is key to identify the hanging wall effect
    vector[0] = Azimuth1to2 + np.pi/2.   # downdip to get the other two points (last row fault trace)
    
    Loc3 = EndLocation( Loc1, vector ) 
    Loc4 = EndLocation( Loc2, vector ) 

    # now  Loc1, Loc2, Loc3, Loc4 give the extended fault surface projection (1000*1000 area)
    # this defines a polygon on the spherical surface and you need to check the site location
    # inside or outside the polygon to compute Rx
    verts = []; segs = []
    verts.append(Loc1) 
    segs.append(Loc1) 
    for iseg in xrange( Ntrace ): 
	verts.append( FaultTrace[Ntrace-1-iseg] ) 
	segs.append( FaultTrace[Ntrace-1-iseg] ) 
    verts.append( Loc2 ) 
    segs.append( Loc2 ) 
    
    verts.append( Loc4 ) 
    verts.append( Loc3 ) 
    
    if 0:
    # Test extended fault trace and fault surface projection (plot) 
	verts1 = np.array( verts ) 
	fig = plt.figure(10) 
	ax = fig.add_subplot( 111 ) 
	ax.plot( verts1[:,0], verts1[:,1], 'ro' )
	ax.plot( [ps[0],pe[0]], [ps[1],pe[1]], 'bx' )  # plot the initial points (where to start from)
	raw_input() 
    if 0: 
	print 'Surface projection verts: '
	print verts
	print 'Fault trace segments: '
	print segs 
	print 'Fault geometry: '
	print FaultGeo
	print 'Site location:' 
	print SiteGeo
    
    # check site is inside the polygon 
    check = CheckPointInPolygon( SiteGeo[:2], np.array(verts)[:,:2] )

    # compute site to fault trace min distance  
    distToExtendedTrace = minDistToLineSeg(SiteGeo, segs, Fast=Fast)

    if check or distToExtendedTrace == 0.0: 
	Rx = distToExtendedTrace   # hanging wall
    else: 
	Rx = - distToExtendedTrace   # foot wall 

    return Rjb, Rrup, Rx 



# ==========================
# Site-Specific Parameters
# ==========================
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


# ====================
# Applications (other)
# ====================
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


