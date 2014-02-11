#!/usr/bin/env python 
"""
NGA GMPEs comparison and evaluation in details
08 and 14 models
careful about the conversion between GMRotI50 and GMRotD50
""" 

import os, sys
from pynga import * 
from pynga.utils import * 

import matplotlib.pyplot as plt 

opt = sys.argv[1]  # magnitude, distance, hangingwall, vs30, basin, 
wrkpth = '/Users/fengw/local/pylib/pynga/' 
os.chdir(wrkpth) 

pltpth00 = './NGAcomparison/plots/'
pltpth0 = pltpth00 + 'NGA08_14'
pltpth = pltpth0 + '/%s'%opt
for f in [pltpth00, pltpth0, pltpth, ]:        
    if not os.path.exists(f): 
        os.mkdir(f)

pfmt = 'png' 
figsize = (10,8)

# ground motion based: 
Ts = [-1, 0.3, 1.0, 3.0]
IMTs = ['PGA','SA0.3','SA1.0','SA3.0']
#IMTs = ['PGA (g)','SA0.3 (g)', 'SA1.0 (g)', 'SA3.0 (g)']

NGA08_models = ['AS','BA','CB','CY']
NGA14_models = ['ASK','BSSA','CB','CY']

# simple parameters (doesn't change very often)
Ftype = 'SS'; Mech = 1; dip = 90; rake = 0.0     # for this case Rjb could equal to Rrup
Ztor = 3; W = 10
VsFlag = 0 

# 0. Test period variations (for SA only)
if opt == 'period':
    IMT = 'SA'
    Ts = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    #Ts = np.array(Ts,'f')    

    Mws = [5, 6, 7]
    Rjb = Rrup = 20
    Vs30 = 760.
    Z10 = Z25 = None
    Fhw = 0

    for Mw in Mws: 
        fig = plt.figure(1,figsize)
        fig.clf()
        title = 'IMT: %s (g), Mw=%s, Rjb=%skm, Rrup=%skm, Vs30=%sm/s, Z1.0=%s, Z2.5=%s'%(IMT, '%.2f'%Mw, '%.2f'%Rjb,'%.2f'%Rrup,'%.2f'%Vs30,Z10, Z25)
        fig.text(0.1,0.95,title,fontsize=16)
        IMs1 = []
        stds1 = []
        for ip in xrange(len(Ts)):
            Ti = Ts[ip]   
            IMs = []
            stds = []
            for inga in xrange(4):
                ax = fig.add_subplot(2,2,inga+1)
                nga1 = NGA08_models[inga]                
                median1, std1, tau, sigma = NGA08( nga1, Mw, Rjb, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )                                  
                
                nga2 = NGA14_models[inga]
                median2, std2, tau, sigma = NGA14( nga2, Mw, Rjb, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
                IMs.append([median1, median2])
                stds.append([std1,std2])
                
            IMs1.append(IMs)
            stds1.append(stds)
        IMs1 = np.array(IMs1)
        stds1 = np.array(stds1)    
        for inga in xrange(4):
            nga1 = NGA08_models[inga]
            nga2 = NGA14_models[inga]
            ax = fig.add_subplot(2,2,inga+1)
            ax.loglog(Ts, IMs1[:,inga,0], 'b-', label=nga1+'08')
            ax.loglog(Ts, IMs1[:,inga,0]*stds1[:,inga,0],'b--')
            ax.loglog(Ts, IMs1[:,inga,0]/stds1[:,inga,0],'b--') 
            ax.loglog(Ts, IMs1[:,inga,1], 'r-', label=nga2+'14')
            ax.loglog(Ts, IMs1[:,inga,1]*stds1[:,inga,1],'r--')
            ax.loglog(Ts, IMs1[:,inga,1]/stds1[:,inga,1],'r--') 
            ax.legend(loc=0).draw_frame(False) 
            ax.grid(b=True,which='major')
            ax.grid(b=True,which='minor')  
        pltnam = '/IMT_%s_periods_Mw_%s.%s'%(IMT,'%.2f'%Mw, pfmt)
        fig.savefig(pltpth+pltnam,format=pfmt)   

# test 1. Magnitude test for given distance and Vs30    
if opt == 'magnitude':
    Mws = np.arange(3.5, 9.0, 0.5)
    Rjbs = [10, 50, 150, 200]
    Rrups = Rjbs 
    Fhw = 0    # in this case, you don't need to provide the Rx 
    Vs30 = 760
    Z10 = Z25 = None
    for ip in xrange(len(Ts)):
        Ti = Ts[ip]
        IMT = IMTs[ip]
        for idist in xrange(len(Rjbs)):
            fig = plt.figure(1,figsize)
            fig.clf()
            Rjb = Rjbs[idist]
            Rrup = Rrups[idist]
            title = 'IMT: %s (g), Rjb=%skm, Rrup=%skm, Vs30=%sm/s, Z1.0=%s, Z2.5=%s'%(IMT,'%.2f'%Rjb,'%.2f'%Rrup,'%.2f'%Vs30,Z10, Z25)
            pltnam = '/IMT_%s_Mws_Rjb_%s_Rrup_%s.%s'%(IMT,'%.2f'%Rjb,'%.2f'%Rrup,pfmt)            
            fig.text(0.1,0.95,title,fontsize=16)
            for inga in xrange(4):
                ax = fig.add_subplot(2,2,inga+1)
                nga1 = NGA08_models[inga]                
                median, std, tau, sigma = NGA08( nga1, Mws, Rjb, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )                      
                ax.semilogy(Mws, median,'b', label=nga1+'08')
                ax.semilogy(Mws, median*std, 'b--')
                ax.semilogy(Mws, median/std, 'b--')
                
                nga2 = NGA14_models[inga]
                median, std, tau, sigma = NGA14( nga2, Mws, Rjb, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
                if 0:                
                    median = np.log(median) 
                    std = np.log(std)
                    tau = np.log(tau)
                    sigma = np.log(sigma)
                    ax.plot(Mws, median,'r', label=nga1+'08')
                    ax.plot(Mws, median+std, 'r--')
                    ax.plot(Mws, median-std, 'r--')     
                else: 
                    ax.semilogy(Mws, median,'r',label=nga2+'14')
                    ax.semilogy(Mws, median*std, 'r--')
                    ax.semilogy(Mws, median/std, 'r--')
                    
                ax.legend(loc=0).draw_frame(False)
                ax.grid(b=True,which='major')
                ax.grid(b=True,which='minor')  
            fig.savefig(pltpth+pltnam,format=pfmt)
            
# 2. Distance attenuation test for given magnitudes
if opt == 'distance0':
    Mw = 4.0 
    Rjbs = np.arange(1,200, 5)
    Rrups = Rjbs
    Fhw = 0    # in this case, you don't need to provide the Rx 
    Vs30 = 760
    Z10 = Z25 = None
    Ti = 0.3 
    IMT = 'SA (g)'
    fig = plt.figure(1,figsize)
    fig.clf()
    title = 'IMT: %s (g), Mw=%s, Vs30=%sm/s, Z1.0=%s, Z2.5=%s'%(IMT,'%.2f'%Mw,'%.2f'%Vs30,Z10, Z25)
    pltnam = '/IMT_%s_Rs_Mw_%s.%s'%(IMT,'%.2f'%Mw, pfmt)            
    fig.text(0.1,0.95,title,fontsize=16)
    for inga in xrange(4):
	ax = fig.add_subplot(2,2,inga+1)
	nga1 = NGA08_models[inga]                
	median, std, tau, sigma = NGA08( nga1, Mw, Rjbs, Vs30, Ti, Mech=Mech, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrups,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag,AB11=1)
	if nga1 == 'BA':
	    Rs = Rjbs
	    for i in xrange(len(Rjbs)):
		print Rjbs[i], median[i], std[i], tau[i],sigma[i]
	else: 
	    Rs = Rrups

	ax.loglog(Rs, median,'b', label=nga1+'08')
	ax.loglog(Rs, median*std, 'b--')
	ax.loglog(Rs, median/std, 'b--')
	
	nga2 = NGA14_models[inga]
	median, std, tau, sigma = NGA14( nga2, Mw, Rjbs, Vs30, Ti, Mech=Mech, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrups,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
	if nga2 == 'BSSA':
	    Rs = Rjbs 
	else: 
	    Rs = Rrups
	ax.loglog(Rs, median,'r', label=nga2+'14')
	ax.loglog(Rs, median*std, 'r--')
	ax.loglog(Rs, median/std, 'r--')
	ax.legend(loc=0).draw_frame(False)
	ax.grid(b=True,which='major')
	ax.grid(b=True,which='minor')  
    fig.savefig(pltpth+pltnam,format=pfmt)



if opt == 'distance':
    Mws= [4, 5, 6, 7, 8]
    Rjbs = np.arange(1,200, 5)
    Rrups = Rjbs
    Fhw = 0    # in this case, you don't need to provide the Rx 
    Vs30 = 760
    Z10 = Z25 = None

    for ip in xrange(len(Ts)):
        Ti = Ts[ip]
        IMT = IMTs[ip]
        for imw in xrange(len(Mws)):
            fig = plt.figure(1,figsize)
            fig.clf()
            Mw = Mws[imw]
            title = 'IMT: %s (g), Mw=%s, Vs30=%sm/s, Z1.0=%s, Z2.5=%s'%(IMT,'%.2f'%Mw,'%.2f'%Vs30,Z10, Z25)
            pltnam = '/IMT_%s_Rs_Mw_%s.%s'%(IMT,'%.2f'%Mw, pfmt)            
            fig.text(0.1,0.95,title,fontsize=16)
            for inga in xrange(4):
                ax = fig.add_subplot(2,2,inga+1)
                nga1 = NGA08_models[inga]                
                median, std, tau, sigma = NGA08( nga1, Mw, Rjbs, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrups,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )                                         
                if nga1 == 'BA':
                    Rs = Rjbs
                else: 
                    Rs = Rrups

                ax.loglog(Rs, median,'b', label=nga1+'08')
                ax.loglog(Rs, median*std, 'b--')
                ax.loglog(Rs, median/std, 'b--')
                
                nga2 = NGA14_models[inga]
                median, std, tau, sigma = NGA14( nga2, Mw, Rjbs, Vs30, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrups,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
                if nga2 == 'BSSA':
                    Rs = Rjbs 
                else: 
                    Rs = Rrups
                ax.loglog(Rs, median,'r', label=nga2+'14')
                ax.loglog(Rs, median*std, 'r--')
                ax.loglog(Rs, median/std, 'r--')
                ax.legend(loc=0).draw_frame(False)
                ax.grid(b=True,which='major')
                ax.grid(b=True,which='minor')  
            fig.savefig(pltpth+pltnam,format=pfmt)

if opt == 'vs30':
    Rjbs = [10, 50, 150, 200]
    Rrups = Rjbs 
    Fhw = 0    # in this case, you don't need to provide the Rx 
    Vs30s = [180,230,500,760,800,1300.]
    Mw = 6.5
    Z10 = Z25 = None
    for ip in xrange(len(Ts)):
        Ti = Ts[ip]
        IMT = IMTs[ip]
        for idist in xrange(len(Rjbs)):
            fig = plt.figure(1,figsize)
            fig.clf()
            Rjb = Rjbs[idist]
            Rrup = Rrups[idist]
            title = 'IMT: %s (g), Rjb=%skm, Rrup=%skm, Mw=%s, Z1.0=%s, Z2.5=%s'%(IMT,'%.2f'%Rjb,'%.2f'%Rrup,Mw,Z10, Z25)
            pltnam = '/IMT_%s_Vs30s_Mw_%s_Rjb_%s_Rrup_%s.%s'%(IMT,'%.2f'%Mw, '%.2f'%Rjb,'%.2f'%Rrup,pfmt)            
            fig.text(0.1,0.95,title,fontsize=16)
            for inga in xrange(4):
                ax = fig.add_subplot(2,2,inga+1)
                nga1 = NGA08_models[inga]                
                median, std, tau, sigma = NGA08( nga1, Mw, Rjb, Vs30s, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )                      
                ax.semilogy(Vs30s, median,'b', label=nga1+'08')
                ax.semilogy(Vs30s, median*std, 'b--')
                ax.semilogy(Vs30s, median/std, 'b--')
                
                nga2 = NGA14_models[inga]
                median, std, tau, sigma = NGA14( nga2, Mw, Rjb, Vs30s, Ti, Ftype=Ftype, rake=rake,W=W,Ztor=Ztor,dip=dip,Rrup=Rrup,Fhw=Fhw,Z10=Z10,Z25=Z25,VsFlag=VsFlag )
                ax.semilogy(Vs30s, median,'r',label=nga2+'14')
                ax.semilogy(Vs30s, median*std, 'r--')
                ax.semilogy(Vs30s, median/std, 'r--')
                    
                ax.legend(loc=0).draw_frame(False)
                ax.grid(b=True,which='major')
                ax.grid(b=True,which='minor')  
            fig.savefig(pltpth+pltnam,format=pfmt)    
