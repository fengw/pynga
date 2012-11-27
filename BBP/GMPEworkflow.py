#!/usr/bin/env python 
"""
Workflow using pynga for broadband platform
"""

# import pynga and its utilities
from pynga import * 
from pynga.utils import * 

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
# =====================
SiteGeom = [-118.0,33.5,0.0]    # location of site
FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, das, ddd = FaultTraceGen( origin, Dims, Mech )
Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,FaultTrace1,UpperSeisDepth,LowerSeisDepth,AveDip)    # in km 

# assuming you get the following site condition from station list:
# just for one, you could use a list of station parameters
Vs30 = 863.   # in m/s
Z10 = None    # in meter
Z25 = None     # in km 

# ======================
# 3. Read periods list 
# ======================
periods = [0.11,]  # in seconds (could be any values in between 0.01 to 10) 

# =====================
# 4. Compute PSA for the given station 
# =====================
nga_model = 'CY'    # you can use 'CB','CY', and 'AS' as well, or loop over all four NGA models
StationMedian = []  # for a given event and a given station, this list has period-dependent PSA values
for ip in xrange( len(periods) ):
    period = periods[ip] 
    median, sigmaT, tau, sigma = NGA08( nga_model, M, Rjb, Vs30, period, rake=rake, dip=dip, W=W, Ztor=ztor, Rrup=Rrup, Rx=Rx, Z10=Z10, Z25=Z25 )   # this input list should be enough based on what BBP has 
    StationMedian.append( median[0] )    # median in (g) for PSA and PGA

# You can save StationMedian into file (one event, one station) for later uses (compare with simulation results for the same event and the same station)


