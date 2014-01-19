#!/usr/bin/env python 
"""
Workflow using pynga for broadband platform
"""

# import pynga and its utilities
from pynga import * 
from pynga.utils import * 

import matplotlib.pyplot as plt 

import sys 
R0 = int(sys.argv[1]) 
GridSize = int( sys.argv[2] )
RrupCut = int(sys.argv[3])   # 30 or 50 

plotpth = '/Users/fengw/local/pylib/pynga/BBP/plots/' 

# ==============================
# 1. Read Event Info from SRC file 
# assume the SRC file has full-path name: sfile
# ==============================
sfile = './inputs/bbp_validation_part_b_scenario_src_files'
lines = open( sfile ).readlines()
values = []
for il in xrange( len(lines) ):
    spl = lines[il].strip().split()
    values.append( float( spl[1] ) )
M, L, dl, W, dw, ztor, strike, rake, dip, lat0, lon0, hypoAS, hypoDD = values[:13]

# information needed to get distance parameters
origin = lon0, lat0    # center 
Dims = L, dl, W, dw, ztor
Mech = strike, dip, rake 

# =====================
# 2. Read Station list and compute distance parameters 
# Vs30, Z1.0 and Z2.5 are also needed from Station List 
# You can give distance file or give directly the site location (geom)
# =====================
data = np.loadtxt( './outputs/bbp_validation_part_b_scenario_StationList_FW_DistExtend%sGridSize%sRrupCutoff%s'%(R0,GridSize,RrupCut), skiprows=1 )
SiteLon = data[:,0]
SiteLat = data[:,1] 

# for validate
SiteRrup = data[:,2]
SiteRx = data[:,3]
SiteRjb = data[:,4] 

if 0:
    SiteRrup0 = []; SiteRx0 = []; SiteRjb0 = []
    for isite in xrange( len(SiteLon) ): 
	SiteGeom = [SiteLon[isite], SiteLat[isite], 0.0]
	FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, das, ddd = FaultTraceGen( origin, Dims, Mech )
	Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1,UpperSeisDepth,LowerSeisDepth,AveDip)    # in km 
	SiteRjb0.append( Rjb )
	SiteRrup0.append( Rrup ) 
	SiteRx0.append( Rx ) 

    fig = plt.figure(1) 
    ax = fig.add_subplot( 111 ) 
    ax.plot( SiteRjb, SiteRjb0, 'ro', label = 'Rjb' ) 
    ax.plot( SiteRrup, SiteRrup0, 'bo', label = 'Rrup' ) 
    ax.plot( SiteRx, SiteRx0, 'go', label = 'Rx' ) 
    plt.show() 

# assuming you get the following site condition from station list:
# just for one, you could use a list of station parameters
Vs30 = 863.   # in m/s
Z10 = None    # in meter
Z25 = None     # in km 

# ======================
# 3. Read periods list 
# ======================
periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,-1,-2]    
nga_models = ['CB','BA','CY','AS']    # you can use 'CB','CY', and 'AS' as well, or loop over all four NGA models
clrs = ['r','g','b','k']

# =====================
# 4. Compute PSA for the given station 
# =====================
for ip in xrange( len(periods) ):
    period = periods[ip] 
    if period < 0:
	if period == -1: 
	    titlelab = 'PGA (g)' 
	else: 
	    titlelab = 'PGV (cm/s)' 
    else: 
	titlelab = 'SA %.3f sec'%period 

    fig = plt.figure(1) 
    fig.clf() 
    ax = fig.add_subplot(111)
    for inga, nga in enumerate( nga_models ):
	median, sigmaT, tau, sigma = NGA08( nga, M, SiteRjb, Vs30, period, rake=rake, dip=dip, W=W, Ztor=ztor, Rrup=SiteRrup, Rx=SiteRx, Z10=Z10, Z25=Z25 )   # this input list should be enough based on what BBP has 
	ax.plot( SiteRrup, median, clrs[inga]+'o', label=nga ) 
	ax.set_xlabel( 'Rrup (km)' ) 
	ax.set_ylabel( 'median (g)' )
    ax.legend(loc=0)
    ax.set_title( titlelab )
    fig.savefig( plotpth + '/Median_Rrup%s_%s.png'%(RrupCut,'%.3f'%period) ) 


