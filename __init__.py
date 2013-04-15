#!/usr/bin/env python
# This is the main module
# when you import pynga
# what it does is to do the following statements

# Package content
import CB08
import BA08
import CY08
import AS08
import SC08
from utils import *


# Period list (available for each NGA models) 
# -1.0: PGA; -2.0: PGV
TsDict = {
	'BA': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,-1,-2],   
	'CB': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,-1,-2],    
	'CY': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,-1,-2],    
	'AS': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
	      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,-1,-2],    
	}


# Compute NGAs (You can use your coefficients and term flags)
def NGA08(model_name, Mw, Rjb, Vs30, period, epislon=0, NGAs=None, \
	  rake=None, Mech=3, Ftype=None, Fnm=None, Frv=None, \
	  dip=None, W=None, Ztor=None, Zhypo=None, Fas=0, \
	  Rrup=None, Rx=None, Fhw=None, azimuth=None, \
	  VsFlag=0, Z25=None, Z15=None, Z10=None, \
	  AS09=None, AB11=None, ArbCB=0 ):
    """
    Combined function to compute median and standard deviation
    
    Arguments (has to be specified)
    ----------
    model_name : choose NGA model you want to use (AS,BA,CB,CY) 
    Mw : moment magnitude 
    Rjb: Joyner-Boore distance in km
         defined as the shortest distance from a site to the surface projection of the rupture surface   
    Vs30: The average shear-wave velocity between 0 and 30-meters depth (site condition) in m/s
    period: period at which you want to use NGA 
            This function allow to use periods that are not in the available periods (refer to TsDict) 
    
    Keywords 
    --------
    [*] shows the default value

    # ================
    # General Keywords
    # ================
    epislon : deviation from the median value [0]
    NGAs : dictionary to select terms in NGA models and use updated coefficents 
              default: 
		 {'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
		  'BA':{'NewCoefs':None,'terms':(1,1,1)},\
		  'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
		  'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}}\
    
    # ===============
    # Source Keywords
    # ===============
    rake: rake angle (in degree) [None] 
          used to determine the fault type
    Mech: Used in BA model [3]
          (0:Strike-slip, 1:Normal, 2:Reverse, 3:Unknown 
    Ftype: fault type string [None]
          'SS': Strike-slip, 'NM': Normal, 'RV': Reverse, 'U': Unknown (unknown is only used in BA model)
    Fnm : 0: not a normal fault; 1: Normal [None]
          default: None
    Frv : 0: not a reverse fault; 1: reverse [None]
          default: None
    dip : dip angle of the fault plane [None]
          default: None 
    W : Rupture width (down-dip) [None]
    Ztor : depth to the top of rupture [None]
    Zhypo: depth to the hypocenter location [None] 
    Fas : Aftershock flag [None] 
          0: Mainshock; 1: Aftershock 
    
    # ================
    # Path Keywords
    # ================
    Rrup: Rupture distance in km [None]
          defined as the distance from a site the to the fault plane
	  For simple fault geometry, function calc_Rrup in utils.py can be used to compute Rrup, otherwise 
	  use DistanceToEvenlyGriddedSurface function in utils.py to compute given fault geometry and site location
    Rx :  horizontal distance between a site and fault trace, in km [None]
	  defined by extending the fault trace (or the top edge of the rupture) to infinity in both directions. 
	  For simple fault geometry, function calc_Rx in utils.py can be used to compute Rrup, otherwise, 
	  use DistanceToEvenlyGriddedSurface function in utils.py to compute given fault geometry and site location
    Fhw : hanging wall flag [None] 
          0: in footwall; 1: in hanging wall  
    azimuth: source-to-site azimuth [None]
           defined as the angle between the positive fault strike direction and the line connecting 
	   a site to the closet point on the surface projection of the top edge of rupture (clockwise) 
	   (used in simple fault geometry) 
    
    # =================
    # Site Keywords
    # =================
    VsFlag : Vs30 inferred or measured flag [0]
            0: inferred Vs30; 1: measured Vs30 
    Z25: basin depth to S wave velocity equal to 2.5 km/s [None], in km 
         Z25 could be estimated by using calc_Z25 function in utils.py given Vs30
    Z15: basin depth to S wave velocity equal to 1.5 km/s [None], in km 
         used to estimate Z2.5 when Z2.5 = None
    Z10: basin depth to S wave velocity equal to 1.0 km/s [None], in meter
         Z10 could be estimated by using calc_Z1 function in utils.py given Vs30

    # =================
    # Updated models 
    # =================
    AS09 : Abrahamson and Silva 2009 updated model (taper5 hanging wall effect) [None]
    AB11 : Atkinson and Boore 2011 updated model with correction term (after more small magnitude events recordings)
    
    # =================
    # Other Keywords 
    # =================
    ArbCB: Campbell and Bozorgnia 2008 model standard deviation [0] 
           0: output total standard deviation is for GMRotIpp intensity measures (rotation-independent)
	   1: output total standard deviation is for arbitrary horizontal component

    """
    
    if NGAs == None:
	NGAs={'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	      'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	      'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	      'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}}\

    dict1 = NGAs
    itmp = 0
    
    # check the input period
    if period > 10.0 or 0<period<0.01:
	print 'invalid period value (it should be within [0.01,10] for SA or == -1,-2 for PGA and PGV'
	raise ValueError

    if model_name == 'BA':
	
	BAnga = BA08.BA08_nga()
	kwds = {'rake':rake,'Mech':Mech,'Ftype':Ftype,'AB11':AB11, 'CoefTerms':dict1[model_name]}   # OpenSHA doesn't have this
	
	periods = np.array(BAnga.periods)
	for ip in xrange( len(periods) ):
	    if abs( period-periods[ip] ) < 0.0001:
		# period is within the periods list
		itmp = 1
		break

	if itmp == 1:
	    # compute median, std directly for the existing period in the period list of the NGA model
	    values = mapfunc( BAnga, Mw, Rjb, Vs30, period, **kwds )
	    values = np.array( values )

	if itmp == 0:
	    # do the interpolation for periods that is not in the period list of the NGA model
	    ind_low = (periods < period).nonzero()[0]
	    ind_high = (periods > period).nonzero()[0]

	    period_low = max( periods[ind_low] )
	    period_high = min( periods[ind_high] )
	    
	    values_low = np.array( mapfunc( BAnga, Mw, Rjb, Vs30, period_low, **kwds ) )
	    values_high = np.array( mapfunc( BAnga, Mw, Rjb, Vs30, period_high, **kwds ) )
	    
	    N1,N2 = np.array( values_low).shape
	    
	    values = np.zeros( (N1,N2) )
	    for icmp in xrange( N2 ):
		if icmp != 0:
		    # stardand values are in ln (g)
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), values_low[:,icmp], values_high[:,icmp], np.log(period) )
		else:
		    # median value is in g
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), np.log(values_low[:,icmp]), np.log(values_high[:,icmp]), np.log(period) )
                    values[:,icmp] = np.exp( values[:,icmp] )    # change the median into g unit (logline gives the result in ln(g))

    if model_name == 'CB':
	
	CBnga = CB08.CB08_nga()
	kwds = {'Ftype':Ftype,'Rrup':Rrup,'Ztor':Ztor,'dip':dip,'Z25':Z25,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'Z15':Z15,'Arb':ArbCB,'CoefTerms':dict1[model_name]}
	
	periods = np.array(CBnga.periods)
	for ip in xrange( len(periods) ):
	    if abs( period-periods[ip] ) < 0.0001:
		# period is within the periods list
		itmp = 1
		break

	if itmp == 1:
	    values = mapfunc( CBnga, Mw, Rjb, Vs30, period, rake, **kwds )
	    values = np.array( values )
        
	if itmp == 0:
	    # do the interpolation for periods that is not in the period list of the NGA model
	    ind_low =  (periods < period).nonzero()[0]
	    ind_high = (periods > period).nonzero()[0]

	    period_low = max( periods[ind_low] )
	    period_high = min( periods[ind_high] )
	    
	    values_low = np.array( mapfunc( CBnga, Mw, Rjb, Vs30, period_low, rake, **kwds ) )
	    values_high = np.array( mapfunc( CBnga, Mw, Rjb, Vs30, period_high, rake, **kwds ) )
	    
	    N1,N2 = np.array( values_low).shape
	    
	    values = np.zeros( (N1,N2) )
	    for icmp in xrange( N2 ):
		if icmp != 0:
		    # stardand values are in ln (g)
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), values_low[:,icmp], values_high[:,icmp], np.log(period) )
		else:
		    # median value is in g
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), np.log(values_low[:,icmp]), np.log(values_high[:,icmp]), np.log(period) )
                    values[:,icmp] = np.exp( values[:,icmp] )    # change the median into g unit (logline gives the result in ln(g))

    if model_name == 'CY':
	
	CYnga = CY08.CY08_nga()
	kwds = {'Ftype':Ftype,'Rrup':Rrup,'Rx':Rx,'Ztor':Ztor,'dip':dip,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'AS':Fas,'VsFlag':VsFlag,'CoefTerms':dict1[model_name]}
	
	periods = np.array(CYnga.periods)
	for ip in xrange( len(periods) ):
	    if abs( period-periods[ip] ) < 0.0001:
		# period is within the periods list
		itmp = 1
		break

	if itmp == 1:
	    values = mapfunc( CYnga, Mw, Rjb, Vs30, period, rake, **kwds )
	    values = np.array( values )
	
	if itmp == 0:
	    print 'Do the interpolation at period = %s for NGA model: %s'%('%.3f'%period, model_name)
	    # do the interpolation for periods that is not in the period list of the NGA model
	    ind_low =  (periods < period).nonzero()[0]
	    ind_high = (periods > period).nonzero()[0]

	    period_low = max( periods[ind_low] )
	    period_high = min( periods[ind_high] )
	    
	    values_low = np.array( mapfunc( CYnga, Mw, Rjb, Vs30, period_low, rake, **kwds ) )
	    values_high = np.array( mapfunc( CYnga, Mw, Rjb, Vs30, period_high, rake, **kwds ) )
	    
	    N1,N2 = np.array( values_low).shape
	    
	    values = np.zeros( (N1,N2) )
	    for icmp in xrange( N2 ):
		if icmp != 0:
		    # stardand values are in ln (g)
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), values_low[:,icmp], values_high[:,icmp], np.log(period) )
		else:
		    # median value is in g
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), np.log(values_low[:,icmp]), np.log(values_high[:,icmp]), np.log(period) )
                    values[:,icmp] = np.exp( values[:,icmp] )    # change the median into g unit (logline gives the result in ln(g))


    if model_name == 'AS':                                                                                                                
	ASnga = AS08.AS08_nga()
	kwds = {'Ftype':Ftype,'Rrup':Rrup,'Rx':Rx,'Ztor':Ztor,'dip':dip,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'Fas':Fas,'VsFlag':VsFlag, 'CoefTerms':dict1[model_name]}
	
	periods = np.array(ASnga.periods)
	for ip in xrange( len(periods) ):
	    if abs( period-periods[ip] ) < 0.0001:
		# period is within the periods list
		itmp = 1
		break
        
	if itmp == 1:
	    values = mapfunc( ASnga, Mw, Rjb, Vs30, period, rake, **kwds )  
	    values = np.array( values )
            #print period, values[22]
	    #raw_input()

	if itmp == 0:
	    print 'Do the interpolation at period = %s for NGA model: %s'%('%.3f'%period, model_name)

	    # do the interpolation for periods that is not in the period list of the NGA model
	    ind_low =  (periods < period).nonzero()[0]
	    ind_high = (periods > period).nonzero()[0]

	    period_low = max( periods[ind_low] )
	    period_high = min( periods[ind_high] )
	    
	    values_low = np.array( mapfunc( ASnga, Mw, Rjb, Vs30, period_low, rake, **kwds ) )
	    values_high = np.array( mapfunc( ASnga, Mw, Rjb, Vs30, period_high, rake, **kwds ) )
	    
	    N1,N2 = np.array( values_low).shape
	    
	    values = np.zeros( (N1,N2) )
	    for icmp in xrange( N2 ):
		if icmp != 0:
		    # stardand values are in ln (g)
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), values_low[:,icmp], values_high[:,icmp], np.log(period) )
		else:
		    # median value is in g
		    values[:,icmp] = logline( np.log(period_low), np.log(period_high), np.log(values_low[:,icmp]), np.log(values_high[:,icmp]), np.log(period) )
                    values[:,icmp] = np.exp( values[:,icmp] )    # change the median into g unit (logline gives the result in ln(g))

    # outputs
    NGAsigmaT = values[:,1]
    NGAtau = values[:,2]
    NGAsigma = values[:,3]
    
    if epislon: 
	NGAmedian = np.exp( np.log(values[:,0]) + epislon * NGAsigmaT )
    else: 
	NGAmedian = values[:,0]  

    # returned quantities are all in g, not in log(g), event for the standard deviations
    return NGAmedian, np.exp( NGAsigmaT ), np.exp( NGAtau ), np.exp( NGAsigma )      # all in g, include the standard deviation


def BA08Test(): 
    # to reproduce BA model (shown in Earthquake Spectra 2008) 

    import matplotlib.pyplot as plt

    NGAs={'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	  'BA':{'NewCoefs':None,'terms':(1,1,1)},\
	  'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},\
	  'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}}\
    
    # validation with BA
    nga = 'BA'
    Mws = [5,6,7,8]
    Vs30 = 760
    FT = 'U'
    Rjb = np.arange( 0.1, 100, 0.5 )
    T = 3.0 
    
    fig = plt.figure(1) 
    ax = fig.add_subplot( 111 )
    lines = []
    for Mw in Mws:
	median, std, tau, sigma = NGA08( nga, Mw, Rjb, Vs30, T, Mech=3, NGAs=NGAs )
	line = ax.loglog( Rjb, median*100 * 9.8 )
        lines.append( line )
    ax.legend( lines, ('M=5','M=6','M=7','M=8'), loc=0 )
    ax.set_title(r"$V_{S30}$ = 760 m/s, mech='U'")
    ax.set_xlabel( r'$R_{JB}$ (km)' )
    ax.set_ylabel( r'5%-damped PSA (cm/s)' )
    plt.show()



def NGAtest(nga): 
    # simple test comparing with file: ./Validation/NGAmodelsTestFiles/nga_Sa_v19a.xls
    M = 6.93 
    Ztor = 3 
    Ftype = 'RV'
    W = 3.85 
    dip = 70 
    Rrup = Rjb = Rx = 30 
    Fhw = 0 
    Vs30 = 760 
    Z10 = 0.024 * 1000   # in meter 
    Z25 = 2.974    # in km 
    VsFlag = 0 

    periods = TsDict[nga] 
    NT = len(periods) 
    Medians = []; SigmaTs = []
    for ip in xrange( NT ): 
	Ti = periods[ip] 
	median, std, tau, sigma = NGA08( nga, M, Rjb, Vs30, Ti, Ftype=Ftype, W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Rx=Rx,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
        Medians.append( median ) 
	SigmaTs.append( np.log(std) ) 
    output = np.c_[ np.array( periods),  np.array( Medians ), np.array( SigmaTs ) ] 
    pth = './tmp'
    if not os.path.exists( pth ): 
	os.mkdir(pth) 
    np.savetxt( pth + '/SimpleTest%s.txt'%nga, output ) 

    print output 


# ====================
# self_application
# ====================
if __name__ == '__main__':

    import sys 

    #BA08Test()
    NGAtest(sys.argv[1])


