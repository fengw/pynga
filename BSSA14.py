#!/usr/bin/env python

from utils import *

class BSSA14_nga:
    """
    Class for Boore et al. 2013 NGA-West2 model
    (validated with results from original NGA modelers)
    """
    def __init__(self):
        self.filepth = './NGA_west2'    # change this in macbook pro
        self.CoefFile = self.filepth + '/BSSA14.csv'
        self.Coefs = {}
        self.ReadModelCoefs()         
        
        # put some period independent coefs here 
        self.Rref = 1.0 
        self.Mref = 4.5 
        self.Vref = 760.
        self.f1 = 0
        self.f3 = 0.1 

        self.faults = ['unspecified','strike-slip','normal','reverse','U','NM','SS','RV']
        self.Dregions = ['GlobalCATW', 'ChinaTurkey', 'ItalyJapan']   # considered regions
        self.countries = ['California','Japan']    # use for basin effect correction given region (default is California)
        
    def ReadModelCoefs(self): 
        self.CoefKeys = open(self.CoefFile,'r').readlines()[1].strip().split(',')[1:]
        inputs = np.loadtxt(self.CoefFile,skiprows=2,delimiter=',')
        self.periods = inputs[:,0]
        coefs = inputs[:,1:]
        for i in xrange( len(self.periods) ):
            T1 = self.periods[i]
            Tkey = GetKey(T1)
            
            # periods list ( -2: PGV, -1: PGA ) (mapping between the NGA models accordingly, -1: PGV, 0: PGA)
            if Tkey == '-1.000':
                Tkey = '-2.000'    # PGV
                self.periods[i] = -2
            if Tkey == '0.000':
                Tkey = '-1.000'    # PGA
                self.periods[i] = -1
            
            # assign to Coefs
            self.Coefs[Tkey] = {}
            for ikey in xrange(len(self.CoefKeys)):
                key = self.CoefKeys[ikey]
                cmd = "self.Coefs['%s']['%s'] = coefs[%i,%i]"%(Tkey,key,i,ikey)
                exec(cmd)
    
    def __call__( self, M ,Rjb ,Vs30, T, rake, Dregion='GlobalCATW', country='California', Mech=3, Ftype=None, Z10=None, CoefTerms={'terms':(1,1,1),'NewCoefs':None}):
	"""
	Compute IM for single period
	required inputs:
	M, Rjb, Vs30, T
	rake: rake angle (degree), default is None (Unspecified fault type)
	or give Mech instead of rake
	Mech: 
	     0: strike
	     1: normal
	     2: reverse
	     else: unspecified (U=1) (Default)
	Ftype = 'U', or 'SS', or 'RV', or 'NM'
        """
	# ==================
	# Input variables
	# ==================
	self.M = float(M)	     # moment magnitude
	self.Rjb = float(Rjb)	     # Joyner-Boore distance (km)
        self.Vs30 = float( Vs30 )    # 30 meter averaged S wave velocity (m/s)
        self.Z10 = Z10 

        self.region = Dregion
        self.country = country
        if self.region not in self.Dregions:
            print '%s is not in %s'%(self.region, self.Dregions)
        if self.country not in self.countries: 
            print '%s is not in %s'%(self.country, self.countries) 
            
        terms = CoefTerms['terms']
	NewCoefs = CoefTerms['NewCoefs']

	if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
	# check inputs
	if self.M == None or self.M < 0:
	    print 'Moment magnitude must be a postive number'
	    raise ValueError
	if self.Rjb == None or self.Rjb < 0:
	    print 'Joyner-Boore distance must be a non-negative number'
	    raise ValueError
	if self.Vs30 == None or self.Vs30 < 0:
	    print 'Vs30 must be a positive number'
	    raise ValueError

	self.rake = rake
	self.Mech = Mech
	
	if rake == None and Mech == None and Ftype == None:
	    print 'either rake or (U,SS,NM,RV) should be provided'
	    raise ValueError
	else: 
	    if Ftype != None: 
		self.U = 1*(Ftype == 'U')
		self.SS = 1*(Ftype == 'SS')
		self.NM = 1*(Ftype == 'NM')
		self.RV = 1*(Ftype == 'RV')
            else: 
		if Mech != None and rake != None:
		    # giveng Mech and rake at the same time, use Mech, not rake
		    rake = None

		if rake != None and Mech == None:
		    # Get ftype from rake
		    self.rake = rake
		    self.ftype()
		
		if rake == None and Mech != None:
		    self.U = 1*(Mech==0)
		    self.SS = 1*(Mech==1)
		    self.NM = 1*(Mech==2)
		    self.RV = 1*(Mech==3)

	# modify the coefficients
	if NewCoefs != None:
	    # only update Coefs given by NewCoefs (at self.T)
	    Tkey = GetKey( self.T )
	    NewCoefKeys = NewCoefs.keys()
	    for key in NewCoefKeys:
		self.Coefs[Tkey][key] = NewCoefs[key]
	
	# ======================
	# begin to compute IM
	# ======================
        IM = self.compute_im(terms=terms)     # Median ground motion [PGV: cm/s; PGA: ln (g); SA ln (g)]
        sigmaT, tau, sigma = self.compute_std()     # standard deviation [same unit as median GM]
        
	return IM, sigmaT, tau, sigma
                         
    # ============================
    # Functions used in the class
    # they could also be output for 
    # further regression analysis
    # ============================
    def ftype(self):
	FT = rake2ftype_BA( self.rake )   # change in this version
	if FT not in self.fault:
	    print 'Invalid fault type!'
	    print 'It should be in one of the following list:'
	    print self.fault
	    raise ValueError
	else:
	    if FT == 'unspecified' or FT == 'U':
		self.U = 1
	    else:
		self.U = 0
	    if FT == 'strike-slip' or FT == 'SS':
		self.SS = 1
	    else:
		self.SS = 0
	    if FT == 'normal' or FT == 'NM':
		self.NM = 1
	    else:
		self.NM = 0
	    if FT == 'reverse' or FT == 'RV':
		self.RV = 1
	    else:
		self.RV = 0
        return FT   

    def moment_function(self, Tother=None):
	"""
	Magnitude-Moment scaling
	"""
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)
        for key in ['e0','e1','e2','e3','e4','e5','e6','Mh']:
            cmd = "%s = self.Coefs['%s']['%s']"%(key,Ti,key)
            exec(cmd) 

	faulting = e0*self.U + e1*self.SS + e2*self.NM + e3*self.RV
        if self.M <= Mh:
	    return faulting + e4*(self.M-Mh) + e5*(self.M-Mh)**2
	else:
	    return faulting + e6*(self.M-Mh)


    def distance_function(self,Tother=None):
	"""
	Distance function
	Geometrical spreading? (yes ~ ln(R))
	"""
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)
        
	h = self.Coefs[Ti]['h']
	c1 = self.Coefs[Ti]['c1']
	c2 = self.Coefs[Ti]['c2']
	c3 = self.Coefs[Ti]['c3']
        D_c3 = self.Coefs[Ti]['D_c3_%s'%self.region] 

	R = np.sqrt( self.Rjb**2 + h**2 )
        term = (c1+c2*(self.M-self.Mref))*np.log(R/self.Rref)+c3*(R-self.Rref) + D_c3*(R-self.Rref)
        return term


    def Vs30toZ1pt0(self, Vs30):
        if self.country == 'California': 
            return (-7.15/4.)*np.log((Vs30**4+570.94**4)/(1360**4+570.94**4))
        if self.country == 'Japan': 
            return (-5.23/2.)*np.log((Vs30**2+412.39**2)/(1360**2+412.39**2))

    def soil_function(self, Vs30=None, Tother=None):
	"""
	Site Amplification Function
	"""
	if Vs30 != None: 
	    self.Vs30 = Vs30 
	if Tother != None: 
	    Ti = GetKey( Tother ) 
	else: 
	    Ti = GetKey(self.T )

	# ===============
	# linear term
        # ===============
	c = self.Coefs[Ti]['c']
        Vc = self.Coefs[Ti]['Vc']
        if self.Vs30 <= Vc: 
            flin = c * np.log(self.Vs30/self.Vref)
        else: 
            flin = c * np.log(Vc/self.Vref)
            
        # =================
	# non-linear term
	# =================
        # 1. compute pga4nl, which is defined as the media PGA when Vs30=Vref=760 m/s
	Tpga = -1.0    # compute PGA
	pga4nl = np.exp( self.moment_function(Tother=Tpga) + self.distance_function(Tother=Tpga) )
	
        # 2. compute nonlinear site effect
        f4 = self.Coefs[Ti]['f4']
	f5 = self.Coefs[Ti]['f5']
        f2 = f4*(np.exp(f5*(min([self.Vs30,760])-360))-np.exp(f5*(760-360)))
        fnl = self.f1 + f2*np.log((pga4nl+self.f3)/self.f3)

        self.Fsite = flin + fnl 

        # consider the correction due to basin apparence (big changes!!!)		        
        if self.Z10 == None:
            dZ1pt0 = 0.0 
            fdZ1 = 0.0  # no consideration of basin effects
        else: 
            f6 = self.Coefs[Ti]['f6']
            f7 = self.Coefs[Ti]['f7']
            Z1pt0_Vs30 = self.Vs30toZ1pt0(self.Vs30)  
            dZ1pt0 = abs(Z1pt0_Vs30 - self.Z10)
            if Ti <0.65: 
                fdZ1 = 0.0 
            elif Ti>=0.65: 
                if dZ1pt0 <= f7/f6:
                    fdZ1 = f6*dZ1pt0 
                else: 
                    fdZ1 = f7 

	return self.Fsite + fdZ1 


    def compute_im(self,terms=(1,1,1)):
        """
	Compute IM based on functional form of BA08 model
	"""
	IM =  np.exp(terms[0]*self.moment_function()+
	             terms[1]*self.distance_function()+
		     terms[2]*self.soil_function())
        # note: for PGA and PSA, IM has unit (g) here 
        #       for PGV, IM has unit (cm/s) !  (after np.exp)
	return IM


    def compute_std(self):
        Ti = GetKey(self.T )
        for key in ['phi1','phi2','tau1','tau2','R1','R2','V1','V2','D_phi_R','D_phi_V']: 
            cmd = "%s = self.Coefs['%s']['%s']"%(key, Ti,key)
            exec(cmd) 

        # compute intra-event sigma
        if self.M <= 4.5: 
            phiM = phi1 
        elif 4.5 < self.M <= 5.5:
            phiM = phi1 + (phi2-phi1)*(self.M-4.5) 
        else: 
            phiM = phi2

        if self.Rjb <= R1: 
            phiMR = phiM 
        elif R1 < self.Rjb <= R2: 
            phiMR = phiM + D_phi_R*(np.log(self.Rjb/R1)/np.log(R2/R1))
        else: 
            phiMR = phiM + D_phi_R 

        if self.Vs30 >= V2:
            sigma = phiMR
        elif V1 <= self.Vs30 < V2:
            sigma = phiMR - D_phi_V*(np.log(V2*1.0/self.Vs30)/np.log(V2*1.0/V1))
        else: 
            sigma = phiMR - D_phi_V 

        # compute inter-event
        if self.M <= 4.5: 
            tau = tau1
        elif 4.5 < self.M <= 5.5:
            tau = tau1 + (tau2-tau1)*(self.M-4.5)
        else: 
            tau = tau2
        sigmaT = np.sqrt(sigma**2+tau**2)

	return sigmaT, tau, sigma


def BSSA14nga_test(T,CoefTerms):
    """
    Basic Test of running of model (how to use it)
    Test of some specific inputs (debug)
    """
    # input parameter list
    Rjb = 0
    Rjb = np.arange(1,200,5)
    Vs30 = 748.0,1200.,345.,160.
    Vs30 = 760.
    Mw = 4.

    rake=180.
    rake = 0.
    Ftype='SS'
    Mech = 1
    kwds = {'Mech':Mech,'Ftype':Ftype, 'Z10':None, 'Dregion':'GlobalCATW', 'country':'California', 'CoefTerms':CoefTerms}
    BSSAnga = BSSA14_nga()    # BA08nga instance
    if 1:
        values = mapfunc( BSSAnga, Mw, Rjb, Vs30, T, rake, **kwds )
        for ivalue in xrange( len(values) ):
            print Rjb[ivalue], values[ivalue]
    else: 
        # debug mode (show each term)
        IM, sigmaT, tau, sigma = BSSAnga(Mw,Rjb,Vs30,T,rake, **kwds)
        print IM, sigmaT, tau, sigma
    
    return BSSAnga 

if __name__ == '__main__':
    import sys 
    T = 0.3; NewCoefs = None     # pure one
    print 'BA SA at %s second'%('%3.2f'%T)
    CoefTerms={'terms':(1,1,1),'NewCoefs':NewCoefs}
    BSSAnga = BSSA14nga_test(T,CoefTerms)

    T = -1.0
    CoefTerms={'terms':(1,1,1),'NewCoefs':None}
    print 'BA PGA at %s second'%('%3.2f'%T)
    BSSAnga = BSSA14nga_test(T,CoefTerms)

