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

# =======================
# NGA Models
# =======================
# Compute NGAs (You can use your coefficients and term flags)
def NGA08(model_name, Mw,Rjb,Vs30, period,rake=None,Mech=3,NGAs=None, \
	  Rrup=None,Rx=None,dip=None,W=None,Ztor=None,Zhypo=None,\
	  Z25=None,Z15=None,Z10=None,azimuth=None, \
	  Fhw=None,Fnm=None,Frs=None,Fas=0,AB11=None,AS09=None,VsFlag=0):
    """
    Compute NGA model
    Flags: 
	Fhw: hanning wall flag (default: None) 
	     0: no hanging wall effect
	     1: with hanging wall effect
	Fnm: normal fault
	Frv: reverse fault
	Fas: aftershock
	AB11: updated model by Atkinson and Boore 2011
	AS09: Abrahamson and Silva 2009 (dip function)
	VsFlag: default is 0 (inferred); 1: measured
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
	kwds = {'rake':rake,'Mech':Mech,'AB11':AB11, 'CoefTerms':dict1[model_name]}   # OpenSHA doesn't have this
	
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
	kwds = {'Rrup':Rrup,'Ztor':Ztor,'dip':dip,'Z25':Z25,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'Z15':Z15,'CoefTerms':dict1[model_name]}
	
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
	kwds = {'Rrup':Rrup,'Rx':Rx,'Ztor':Ztor,'dip':dip,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'AS':Fas,'VsFlag':VsFlag,'CoefTerms':dict1[model_name]}
	
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
	    print 'do the interpolation...'
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
	kwds = {'Rrup':Rrup,'Rx':Rx,'Ztor':Ztor,'dip':dip,'W':W,'Zhypo':Zhypo,'azimuth':azimuth,'Fhw':Fhw,'Z10':Z10,'Fas':Fas,'VsFlag':VsFlag, 'CoefTerms':dict1[model_name]}
	
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
	    print 'AS interpolating period'

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
    NGAmedian = values[:,0]
    NGAsigmaT = values[:,1]
    NGAtau = values[:,2]
    NGAsigma = values[:,3]
    
    # returned quantities are all in g, not in log(g), event for the standard deviations
    return NGAmedian, np.exp( NGAsigmaT ), np.exp( NGAtau ), np.exp( NGAsigma )      # all in g, include the standard deviation


# ====================
# self_application
# ====================
if __name__ == '__main__':
    
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


