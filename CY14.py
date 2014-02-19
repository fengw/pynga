#!/usr/bin/env python

from utils import *

class CY14_nga:
    """
    Class for Chiou and Youngs 2014 NGA model
    """
    def __init__(self):
        
        self.filepth = './NGA_west2'    # change this in macbook pro
        self.CoefFile = self.filepth + '/CY14.csv'
        self.Coefs = {}
        self.ReadModelCoefs() 
        self.countries = ['California', 'Japan'] 

        # period independent parameters 
        self.c2 = 1.06
        self.c4 = -2.1
        self.c4a = -0.5
        self.cRB = 50 
        self.c8 = 0.2153
        self.c8a = 0.2695 
        

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
              
            self.Coefs[Tkey] = {}
            for ikey in xrange(len(self.CoefKeys)):
                key = self.CoefKeys[ikey]
                cmd = "self.Coefs['%s']['%s'] = coefs[%i,%i]"%(Tkey,key,i,ikey)
                exec(cmd)


    # call the function 
    def __call__(self,M,Rjb,Vs30,T,rake, Ftype = None, \
	         Rrup=None,Rx=None,dip=None,Ztor=None,Z10=None,\
	         W=None,Zhypo=None,azimuth=None,Fhw=None, D_DPP=0,\
		 AS=0, VsFlag=1, country='California', \
		 CoefTerms={'terms':(1,1,1,1,1,1,1),'NewCoefs':None} \
		 ):

        if T == -1:
            T = 0.01     # for CY model, PGA's coefficients share with SA(0.01)
        if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
        # required inputs
	self.M = M         # Moment Magnitude
	self.Rjb = Rjb     # Joyner-Boore distance (km)
	self.rake = rake   # rake angle
	self.Vs30 = Vs30   # site-condition (m/s)
	self.AS = AS       # Aftershock flag (0 or 1)  (depends on the earthquake itself)
	self.VsFlag = VsFlag # 0: inferred Vs30; 1: measured Vs30
        self.country = country 

	terms = CoefTerms['terms']
	NewCoefs = CoefTerms['NewCoefs']

        # Obtain optional parameters
	if Ftype != None:
	    self.Fnm = 1*(Ftype == 'NM')
	    self.Frv = 1*(Ftype == 'RV')
	else: 
	    if rake == None or rake < -180 or rake > 180.:
		print 'rake angle should be within [-180,180]'
		raise ValueError
	    else: 
		self.Frv, self.Fnm = rake2ftype_CY( self.rake )

	if W == None:
	    if self.rake == None: 
		print 'you should give either the fault width W or the rake angle'
		raise ValueError
	    else:
		W = calc_W(self.M,self.rake)
	else: 
	    self.W = W 
	
	if dip == None:
	    if self.rake == None: 
		print 'you should give either the fault dip angle or the rake angle'
		raise ValueError
	    else:
		self.dip = calc_dip( self.rake )
	else:
	    self.dip = dip
	
        if Zhypo == None: 
            self.Zhypo = calc_Zhypo(self.M,self.rake) 
        else: 
            self.Zhypo = Zhypo 

	if Ztor == None:
	    if Zhypo == None:
		if self.rake == None: 
		    print 'you should give either the Ztor or the rake angle'
		    raise ValueError
		else:
		    Zhypo = calc_Zhypo( self.M, self.rake )
	    self.Ztor = calc_Ztor( W, self.dip, Zhypo )
        else:
	    self.Ztor = Ztor
	
	if Fhw == None:
	    if azimuth == None and Rx == None:
		print 'either one of azimuth angle, Rx and Fhw has to be specified'
		raise ValueError

	    if azimuth != None:
		if 0 <= azimuth <= 180. and dip != 90.:
		    Fhw = 1
		else:
		    Fhw = 0
	    
	    elif Rx != None:
		if Rx >=0 and dip != 90.:
		    Fhw = 1
		else:
		    Fhw = 0
	    
	    if dip == 90:
		Fhw = 0
	
	if azimuth == None:
	    if Fhw == 1:
		azimuth = 50
	    else:
		azimuth = -50.
	
	if self.Rjb == 0:
	    azimuth = 90.
	    Fhw = 1
	self.Fhw = Fhw 

	# Compute Rx and Rrup
	if Rx == None: 
	    self.Rx = calc_Rx( self.Rjb, self.Ztor, W, self.dip, azimuth, Rrup )
	else:
	    self.Rx = Rx
	if Rrup == None:
	    self.Rrup = calc_Rrup( self.Rx, self.Ztor, W, self.dip, azimuth, self.Rjb )
        else:
	    self.Rrup = Rrup 

	# Z10 (empirical relationship depends on dataset used to obtain the relationship)
	if Z10 == None:
	    self.Z10 = calc_Z1(self.Vs30,'CY')   # in meter
	else:
	    self.Z10 = Z10   # Z10 should be in meter  (for CY14 model)
        
        # directivity parameter 
        if D_DPP == None: 
            # compute D_DPP  Chiou and Spudich 2013 
            pass 
        else: 
            self.D_DPP = D_DPP 

        # update coeficient
	if NewCoefs != None:
	    NewCoefKeys = NewCoefs.keys()
	    Tkey = GetKey(self.T)
	    for key in NewCoefKeys:
		self.Coefs[Tkey][key] = NewCoefs[key]
        
	IM = self.compute_im()    # in g
	sigma, tau, sigmaT = self.calc_sigma_tau()     # in ln(g)

	return IM, sigmaT, tau, sigma

    
    def flt_function(self):
	Ti = GetKey(self.T)

	c1 = self.Coefs[Ti]['c1']
	c1a = self.Coefs[Ti]['c1a']
	c1b = self.Coefs[Ti]['c1b']
	c1c = self.Coefs[Ti]['c1c']
	c1d = self.Coefs[Ti]['c1d']
	c7 = self.Coefs[Ti]['c7']
	c7b = self.Coefs[Ti]['c7b']
	c11 = self.Coefs[Ti]['c11']
	c11b = self.Coefs[Ti]['c11b']
        tmp = np.cosh(2*max([self.M-4.5,0]))
        term0 = c1 + (c1a+c1c/tmp)*self.Frv + (c1b+c1d/tmp)*self.Fnm    # faulting type
        
        MeanZtor = self.calc_MeanZtor()
        D_Ztor = abs(self.Ztor-MeanZtor) 
        term1 = (c7+c7b/tmp)*D_Ztor             # Ztor 

        term2 = (c11+c11b/tmp)*(np.cos(self.dip*np.pi/180.))**2      # Dip related 

        return term0 + term1 + term2
    

    def calc_MeanZtor(self, M=None, Frv=None, Fnm=None):
        if M == None: 
            M = self.M
        if Frv == None: 
            Frv = self.Frv
        if Fnm == None: 
            Fnm = self.Fnm 
        MeanZtor = Frv*(max([2.704-1.226*max(M-5.849,0),0]))**2 + Fnm*(max([2.673-1.136*max(M-4.970,0),0]))**2
        return MeanZtor


    def moment_function(self):
	Ti = GetKey(self.T)
	c3 = self.Coefs[Ti]['c3']
	cn = self.Coefs[Ti]['cn']
	cM = self.Coefs[Ti]['cM']
	term2 = self.c2*(self.M-6)+(self.c2-c3)/cn*np.log(1+np.exp(cn*(cM-self.M)))
	return term2

    def distance_function( self ):
	Ti = GetKey(self.T)
	c5 = self.Coefs[Ti]['c5']
	c6 = self.Coefs[Ti]['c6']
	cHM = self.Coefs[Ti]['cHM']
	cg1 = self.Coefs[Ti]['cg1']
	cg2 = self.Coefs[Ti]['cg2']
	cg3 = self.Coefs[Ti]['cg3']
	term3 = self.c4*np.log(self.Rrup+c5*np.cosh(c6*max(self.M-cHM,0)))
	term4 = (self.c4a-self.c4)*np.log(np.sqrt(self.Rrup**2+self.cRB**2))
	term5 = (cg1+cg2/np.cosh(max(self.M-cg3,0)))*self.Rrup
	return term3+term4+term5

    def directivity_function(self): 
        Ti = GetKey(self.T)

        c8b = self.Coefs[Ti]['c8b']
        d_taper = max([1-max([self.Rrup-40,0])/30.,0])
        m_taper = min([max([self.M-5.5,0])/0.8,1])
        term6 = self.c8 * d_taper * m_taper * np.exp(-self.c8a*(self.M-c8b)**2) * self.D_DPP
        return term6 

    def hw_function(self):
	Ti = GetKey(self.T)
	
	c9 = self.Coefs[Ti]['c9']
	c9a = self.Coefs[Ti]['c9a']
	c9b = self.Coefs[Ti]['c9b']

	d = self.dip*np.pi/180.
        term7 = c9*self.Fhw*np.cos(d)*(c9a+(1-c9a)*np.tanh(self.Rx/c9b))*(1-np.sqrt(self.Rjb**2+self.Ztor**2)/(self.Rrup+1)) 
        return term7    
    
    
    def lnYref(self):
	# assume the site effects and basin effects are zero here
	return self.moment_function() + self.distance_function() + self.flt_function() + self.directivity_function() + self.hw_function()

    def site_function( self ):
	Ti = GetKey(self.T)

	f1 = self.Coefs[Ti]['phi1']
	f2 = self.Coefs[Ti]['phi2']
	f3 = self.Coefs[Ti]['phi3']
	f4 = self.Coefs[Ti]['phi4']

	lnY_ref = self.lnYref()
	term8 = f1*min([np.log(self.Vs30/1130.),0])
	term9 = f2*( np.exp(f3*(min(self.Vs30,1130)-360)) - np.exp(f3*(1130-360)) )*np.log((np.exp(lnY_ref)+f4)/f4)
	return term8 + term9
   

    def calc_MeanZ10(self, Vs30=None, country='California'): 
        if Vs30 == None: 
            Vs30 = self.Vs30 
        if country == 'California': 
            MeanLnZ10 = -7.15/4. * np.log((Vs30**4+571.**2)/(1360.**4+571.**4))
        elif country == 'Japan': 
            MeanLnZ10 = -5.23/2. * np.log((Vs30**2+412.**2)/(1360.**2+412.**2))
        else: 
            # for other region, just use default California
            MeanLnZ10 = -7.15/4. * np.log((Vs30**4+571.**2)/(1360.**4+571.**4))
        return np.exp(MeanLnZ10)/1000.
    

    def basin_function(self,Z10=None,Tother=None):
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )
	
	if Z10 != None: 
	    self.Z10 = Z10 
        
        MeanZ10 = self.calc_MeanZ10(country=self.country)

        D_Z10 = abs(self.Z10-MeanZ10) 
	
        phi5 = self.Coefs[Ti]['phi5']
	phi6 = self.Coefs[Ti]['phi6']

	term10  = phi5 * (1-np.exp(-D_Z10/phi6))
	return term10


    def compute_im(self, terms=(1,1,1,1,1,1,1)):
	# use this one
	return  np.exp( 
		terms[0]*self.moment_function() + \
		terms[1]*self.flt_function() + \
		terms[2]*self.hw_function() + \
		terms[3]*self.distance_function() + \
		terms[4]*self.directivity_function() + \
                terms[5]*self.basin_function() + \
		terms[6]*self.site_function() )


    def calc_NL(self):
	
	Ti = GetKey( self.T )
	
	f2 = self.Coefs[Ti]['phi2']
	f3 = self.Coefs[Ti]['phi3']
	f4 = self.Coefs[Ti]['phi4']
	
	yref = np.exp( self.lnYref() )

	b = f2 * ( np.exp( f3*(min(self.Vs30,1130)-360) ) - np.exp(f3*(1130-360)) )    # Eqn10

	return b*yref / (yref+f4)


    def calc_sigma_tau(self):
	Ti = GetKey(self.T)
        if self.VsFlag == 0:
	    Finfer = 1
	    Fmeasure = 0
	else:
	    Finfer = 0
	    Fmeasure = 1
	
        sigma1 = self.Coefs[Ti]['sigma1']
        sigma2 = self.Coefs[Ti]['sigma2']
        sigma3 = self.Coefs[Ti]['sigma3']
        tau1 = self.Coefs[Ti]['tau1'] 
        tau2 = self.Coefs[Ti]['tau2']
	NL = self.calc_NL()

        tmp = min([max([self.M,5]),7.25])-5
        tau = tau1 + (tau2-tau1)/2.25 * tmp
        sigma = (sigma1 + (sigma2-sigma1)/2.25 * tmp) * np.sqrt(sigma3*Finfer+0.7*Fmeasure+(1+NL)**2)
        
	# correct tau
	tauNL = (1+NL)*tau
	sigmaT = np.sqrt( sigma**2 + tauNL**2 )
        
	#return (sigma, tau, sigmaT)
	return (sigma, tauNL, sigmaT)


def CY14nga_test(T,CoefTerms):
    """
    Test CY nga model
    """
    M = 7.75
    Rjb = 10.0 

    #Vs30 = 748.0,1200.0, 356., 160.
    Vs30 = 865.0
    rake = 90    # for specific rupture
    
    W = 20

    Rrup = 21.0
    Rx = 20.0
    Ztor= 0.64274240
    dip = 45
    Z10 = 1000.0
    Z10 = None
    AS = 0
    VsFlag = 0

    CYnga = CY14_nga()
    
    kwds= {'Ztor':Ztor,'dip':dip,'Rrup':Rrup,'Rx':Rx,'Z10':Z10,'AS':AS,'VsFlag':VsFlag,'CoefTerms':CoefTerms} 
    values = mapfunc( CYnga, M, Rjb, Vs30, T, rake, **kwds )
    print 'Median, SigmaT, Tau, Sigma'
    for i in xrange( len(values) ):
	print values[i]
    return CYnga


if __name__ == '__main__':

    T = 0.1; NewCoefs={'c1':-0.5747, 'c1a':0.1}
    T = 0.1; NewCoefs={'c1':-0.6747, 'c1a':0.1}
    T = 0.01; NewCoefs=None
    
    CoefTerms = {'terms':(1,1,1,1,1,1,1),'NewCoefs':NewCoefs}
    print 'CY SA at %s'%('%3.2f'%T)
    CYnga = CY14nga_test(T,CoefTerms)
    
    T = -1
    print 'CY PGA:'
    CYnga = CY14nga_test(T,CoefTerms)


