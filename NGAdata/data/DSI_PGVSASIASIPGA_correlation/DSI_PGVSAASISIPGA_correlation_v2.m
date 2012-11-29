function DSI_PGVSAASISIPGA_correlation_v2
clc
%Brendon Bradley 10 Aug 2010

%Purpose: To look at the empirical correlations between various intensity
%measures, DSI, PGV and SI, ASI, PGA, SA
nboot=100;

%get the observational data
data=observedIMintensities;
EQID=data(:,2);
%spectral periods to consider
T=[0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2 3 4 5 7.5 10];

h_wait=waitbar(0,'Running..');
for i=1:length(data)
    waitbar(i/length(data));
    
    %get Metadata
    rec_num(i)=data(i,1);   EQ_num(i)=data(i,2);
    M=data(i,3);    faultprop.dip=data(i,4);  faultprop.AS=0;
    faultprop.rake=data(i,5);    mech=data(i,6);
    faultprop.Ztor=data(i,11); 
    %use Wells and Coppersmith 'All' to get rup width given Mw
    faultprop.W=10^(-1.01+0.32*M);
    if mech==0
        faultprop.faultstyle='strikeslip';
    elseif mech==1
        faultprop.faultstyle='normal';
    elseif mech==2
        faultprop.faultstyle='reverse';
    else
        faultprop.faultstyle='other';
    end
    
    Rjb=data(i,14); 
    Rrup=data(i,15);     siteprop.Rjb=Rjb; siteprop.Rx=Rjb;
    siteprop.V30=data(i,16); siteprop.V30measured=0; siteprop.Z1pt0=-1; siteprop.Zvs=-1;
    siteprop.orientation='average'; siteprop.g=981; %in cm/s2
    
    %Intensity measures
    
    %     ------------DSI----------------
    %get the predicted value of DSI using BA08
    IMR=@BooreAtkinson_2007_nga;
    [DSI_BA08(i),sigma_DSI_BA08(i,1:3)]=Bradleyetal_2011_DSI(M,Rjb,siteprop,faultprop,IMR);
    %get the predicted value using CY08
    IMR=@ChiouYoungs_2008_nga;
    [DSI_CY08(i),sigma_DSI_CY08(i,1:3)]=Bradleyetal_2011_DSI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using CB08
    IMR=@CampbellBozorgina_2007_nga;
    [DSI_CB08(i),sigma_DSI_CB08(i,1:3)]=Bradleyetal_2011_DSI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using AS08
    IMR=@AbrahamsonSilva_2008_nga;
    [DSI_AS08(i),sigma_DSI_AS08(i,1:3)]=Bradleyetal_2011_DSI(M,Rrup,siteprop,faultprop,IMR);
    %compute DSI for the GM record
    Sa_forDSI=data(i,36:39); T_forDSI=[2 3 4 5];
    Sd_forDSI=(T_forDSI/(2*pi)).^2.*Sa_forDSI*siteprop.g;
    DSI_observed(i)=trapz(T_forDSI,Sd_forDSI);    
    
    freq_min=data(i,17);
    T_max=1/freq_min;
    if T_max>=5
        %compute zlnDSI
        z_lnDSI_BA08(i)=(log(DSI_observed(i))-log(DSI_BA08(i)))/sigma_DSI_BA08(i,1);
        z_lnDSI_CY08(i)=(log(DSI_observed(i))-log(DSI_CY08(i)))/sigma_DSI_CY08(i,1);
        z_lnDSI_CB08(i)=(log(DSI_observed(i))-log(DSI_CB08(i)))/sigma_DSI_CB08(i,1);
        z_lnDSI_AS08(i)=(log(DSI_observed(i))-log(DSI_AS08(i)))/sigma_DSI_AS08(i,1);
    else
        z_lnDSI_BA08(i)=-100;
        z_lnDSI_CY08(i)=-100;
        z_lnDSI_CB08(i)=-100;
        z_lnDSI_AS08(i)=-100;
    end
    
%     ------------ASI----------------
    %get the predicted value of ASI using BA08
    IMR=@BooreAtkinson_2007_nga;
    [ASI_BA08(i),sigma_ASI_BA08(i,1:3)]=Bradleyetal_2008_ASI(M,Rjb,siteprop,faultprop,IMR);
    %get the predicted value using CY08
    IMR=@ChiouYoungs_2008_nga;
    [ASI_CY08(i),sigma_ASI_CY08(i,1:3)]=Bradleyetal_2008_ASI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using CB08
    IMR=@CampbellBozorgina_2007_nga;
    [ASI_CB08(i),sigma_ASI_CB08(i,1:3)]=Bradleyetal_2008_ASI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using AS08
    IMR=@AbrahamsonSilva_2008_nga;
    [ASI_AS08(i),sigma_ASI_AS08(i,1:3)]=Bradleyetal_2008_ASI(M,Rrup,siteprop,faultprop,IMR);
    %compute ASI for the GM record
    Sa_forASI=data(i,26:32); T_forASI=[0.1 0.15 0.2 0.25 0.3 0.4 0.5];
    ASI_observed(i)=trapz(T_forASI,Sa_forASI);
    %compute zlnASI
    z_lnASI_BA08(i)=(log(ASI_observed(i))-log(ASI_BA08(i)))/sigma_ASI_BA08(i,1);
    z_lnASI_CY08(i)=(log(ASI_observed(i))-log(ASI_CY08(i)))/sigma_ASI_CY08(i,1);
    z_lnASI_CB08(i)=(log(ASI_observed(i))-log(ASI_CB08(i)))/sigma_ASI_CB08(i,1);
    z_lnASI_AS08(i)=(log(ASI_observed(i))-log(ASI_AS08(i)))/sigma_ASI_AS08(i,1);
    
    %------------SI----------------
    %get the predicted value of SI using BA08
    IMR=@BooreAtkinson_2007_nga; 
    [SI_BA08(i),sigma_SI_BA08(i,1:3)]=Bradleyetal_2008_SI(M,Rjb,siteprop,faultprop,IMR);
    %get the predicted value using CY08
    IMR=@ChiouYoungs_2008_nga;
    [SI_CY08(i),sigma_SI_CY08(i,1:3)]=Bradleyetal_2008_SI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using CB08
    IMR=@CampbellBozorgina_2007_nga;
    [SI_CB08(i),sigma_SI_CB08(i,1:3)]=Bradleyetal_2008_SI(M,Rrup,siteprop,faultprop,IMR);
    %get the predicted value using AS08
    IMR=@AbrahamsonSilva_2008_nga;
    [SI_AS08(i),sigma_SI_AS08(i,1:3)]=Bradleyetal_2008_SI(M,Rrup,siteprop,faultprop,IMR);
    %compute SI for the GM record
    %get the Sa at T=2.5
    Sa_obs=data(i,20:41);
    Sa2pt5_obs = interp1(T,Sa_obs,2.5);
    Sa_forSI=[data(i,26:36) Sa2pt5_obs]; T_forSI=[0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2.0 2.5];
    Sv_forSI=(T_forSI/(2*pi)).*Sa_forSI*siteprop.g;
    SI_observed(i)=trapz(T_forSI,Sv_forSI);
    %compute zlnSI
    z_lnSI_BA08(i)=(log(SI_observed(i))-log(SI_BA08(i)))/sigma_SI_BA08(i,1);
    z_lnSI_CY08(i)=(log(SI_observed(i))-log(SI_CY08(i)))/sigma_SI_CY08(i,1);
    z_lnSI_CB08(i)=(log(SI_observed(i))-log(SI_CB08(i)))/sigma_SI_CB08(i,1);
    z_lnSI_AS08(i)=(log(SI_observed(i))-log(SI_AS08(i)))/sigma_SI_AS08(i,1);
    
%     ------------PGA----------------
    %get the predicted value of PGA using BA08
    siteprop.period=0;
    [PGA_BA08(i),sigma_PGA_BA08(i,1:3)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
    %get the predicted value using CY08
    [PGA_CY08(i),sigma_PGA_CY08(i,1:3)]=ChiouYoungs_2008_nga(M,Rrup,siteprop,faultprop);
    %get the predicted value using CB08
    [PGA_CB08(i),sigma_PGA_CB08(i,1:3)]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);
    %get the predicted value using AS08
    [PGA_AS08(i),sigma_PGA_AS08(i,1:3)]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
    %get PGA for the GM record
    PGA_observed(i)=data(i,18);
    %compute zlnPGA
    z_lnPGA_BA08(i)=(log(PGA_observed(i))-log(PGA_BA08(i)))/sigma_PGA_BA08(i,1);
    z_lnPGA_CY08(i)=(log(PGA_observed(i))-log(PGA_CY08(i)))/sigma_PGA_CY08(i,1);
    z_lnPGA_CB08(i)=(log(PGA_observed(i))-log(PGA_CB08(i)))/sigma_PGA_CB08(i,1);
    z_lnPGA_AS08(i)=(log(PGA_observed(i))-log(PGA_AS08(i)))/sigma_PGA_AS08(i,1);
    
%     ------------PGV----------------
%     get the predicted value of PGV using BA08
    siteprop.period=-1;
    [PGV_BA08(i),sigma_PGV_BA08(i,1:3)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
%     get the predicted value using CY08
    [PGV_CY08(i),sigma_PGV_CY08(i,1:3)]=ChiouYoungs_2008_nga(M,Rrup,siteprop,faultprop);
%     get the predicted value using CB08
    [PGV_CB08(i),sigma_PGV_CB08(i,1:3)]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);
%     get the predicted value using AS08
    [PGV_AS08(i),sigma_PGV_AS08(i,1:3)]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
%     get PGA for the GM record
    PGV_observed(i)=data(i,19);
%     compute zlnPGA
    z_lnPGV_BA08(i)=(log(PGV_observed(i))-log(PGV_BA08(i)))/sigma_PGV_BA08(i,1);
    z_lnPGV_CY08(i)=(log(PGV_observed(i))-log(PGV_CY08(i)))/sigma_PGV_CY08(i,1);
    z_lnPGV_CB08(i)=(log(PGV_observed(i))-log(PGV_CB08(i)))/sigma_PGV_CB08(i,1);
    z_lnPGV_AS08(i)=(log(PGV_observed(i))-log(PGV_AS08(i)))/sigma_PGV_AS08(i,1);
    
%     ------------SA----------------
    for j=1:length(T)
        siteprop.period=T(j);
%         get the predicted value of PGA using BA08
        [SA_BA08(i,j),sigma_SA_BA08(i,1:3,j)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
%         get the predicted value using CY08
        [SA_CY08(i,j),sigma_SA_CY08(i,1:3,j)]=ChiouYoungs_2008_nga(M,Rrup,siteprop,faultprop);
%         get the predicted value using CB08
        [SA_CB08(i,j),sigma_SA_CB08(i,1:3,j)]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);
%         get the predicted value using AS08
        [SA_AS08(i,j),sigma_SA_AS08(i,1:3,j)]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
%         get PGA for the GM record
        SA_observed(i,j)=data(i,19+j);
%         if the usable frequency is ok compute residual else set residual =
%         -100;
        freq_min=data(i,17);
        T_max=1/freq_min;
        if T(j)<T_max
%             compute zlnSA
            z_lnSA_BA08(i,j)=(log(SA_observed(i,j))-log(SA_BA08(i,j)))/sigma_SA_BA08(i,1,j);
            z_lnSA_CY08(i,j)=(log(SA_observed(i,j))-log(SA_CY08(i,j)))/sigma_SA_CY08(i,1,j);
            z_lnSA_CB08(i,j)=(log(SA_observed(i,j))-log(SA_CB08(i,j)))/sigma_SA_CB08(i,1,j);
            z_lnSA_AS08(i,j)=(log(SA_observed(i,j))-log(SA_AS08(i,j)))/sigma_SA_AS08(i,1,j);
        else
            z_lnSA_BA08(i,j)=-100;
            z_lnSA_CY08(i,j)=-100;
            z_lnSA_CB08(i,j)=-100;
            z_lnSA_AS08(i,j)=-100;
        end
    end
end

%now determine the intra and inter event terms for all of the IMs
alleta=0; %see get_interintraeventterms
k=0;
for j=1:length(z_lnDSI_BA08)  %will be same for both BA08 and CY08 so dont need to check twice
    if z_lnDSI_BA08(j)~=-100
        k=k+1;

        EQID_ok(k)=EQID(j);
        z_lnDSI_BA08_ok(k)=z_lnDSI_BA08(j);
        sigma_DSI_BA08_ok(k,1:3)=sigma_DSI_BA08(j,1:3);
        z_lnDSI_CY08_ok(k)=z_lnDSI_CY08(j);
        sigma_DSI_CY08_ok(k,1:3)=sigma_DSI_CY08(j,1:3);
        z_lnDSI_CB08_ok(k)=z_lnDSI_CB08(j);
        sigma_DSI_CB08_ok(k,1:3)=sigma_DSI_CB08(j,1:3);
        z_lnDSI_AS08_ok(k)=z_lnDSI_AS08(j);
        sigma_DSI_AS08_ok(k,1:3)=sigma_DSI_AS08(j,1:3);
        
        z_lnASI_BA08_ok(k)=z_lnASI_BA08(j);
        sigma_ASI_BA08_ok(k,1:3)=sigma_ASI_BA08(j,1:3);
        z_lnASI_CY08_ok(k)=z_lnASI_CY08(j);
        sigma_ASI_CY08_ok(k,1:3)=sigma_ASI_CY08(j,1:3);
        z_lnASI_CB08_ok(k)=z_lnASI_CB08(j);
        sigma_ASI_CB08_ok(k,1:3)=sigma_ASI_CB08(j,1:3);
        z_lnASI_AS08_ok(k)=z_lnASI_AS08(j);
        sigma_ASI_AS08_ok(k,1:3)=sigma_ASI_AS08(j,1:3);
        
        z_lnPGA_BA08_ok(k)=z_lnPGA_BA08(j);
        sigma_PGA_BA08_ok(k,1:3)=sigma_PGA_BA08(j,1:3);
        z_lnPGA_CY08_ok(k)=z_lnPGA_CY08(j);
        sigma_PGA_CY08_ok(k,1:3)=sigma_PGA_CY08(j,1:3);
        z_lnPGA_CB08_ok(k)=z_lnPGA_CB08(j);
        sigma_PGA_CB08_ok(k,1:3)=sigma_PGA_CB08(j,1:3);
        z_lnPGA_AS08_ok(k)=z_lnPGA_AS08(j);
        sigma_PGA_AS08_ok(k,1:3)=sigma_PGA_AS08(j,1:3);
            
        z_lnSI_BA08_ok(k)=z_lnSI_BA08(j);
        sigma_SI_BA08_ok(k,1:3)=sigma_SI_BA08(j,1:3);
        z_lnSI_CY08_ok(k)=z_lnSI_CY08(j);
        sigma_SI_CY08_ok(k,1:3)=sigma_SI_CY08(j,1:3);
        z_lnSI_CB08_ok(k)=z_lnSI_CB08(j);
        sigma_SI_CB08_ok(k,1:3)=sigma_SI_CB08(j,1:3);
        z_lnSI_AS08_ok(k)=z_lnSI_AS08(j);
        sigma_SI_AS08_ok(k,1:3)=sigma_SI_AS08(j,1:3);
            
        z_lnPGV_BA08_ok(k)=z_lnPGV_BA08(j);
        sigma_PGV_BA08_ok(k,1:3)=sigma_PGV_BA08(j,1:3);
        z_lnPGV_CY08_ok(k)=z_lnPGV_CY08(j);
        sigma_PGV_CY08_ok(k,1:3)=sigma_PGV_CY08(j,1:3);
        z_lnPGV_CB08_ok(k)=z_lnPGV_CB08(j);
        sigma_PGV_CB08_ok(k,1:3)=sigma_PGV_CB08(j,1:3);
        z_lnPGV_AS08_ok(k)=z_lnPGV_AS08(j);
        sigma_PGV_AS08_ok(k,1:3)=sigma_PGV_AS08(j,1:3); 
    end
end

%now get inter intra residuals
%     ------------DSI----------------
[eta_lnDSI_BA08,eps_lnDSI_BA08]=get_interintraeventterms(z_lnDSI_BA08_ok,EQID_ok,sigma_DSI_BA08_ok(:,1)',sigma_DSI_BA08_ok(:,2)',sigma_DSI_BA08_ok(:,3)',alleta);
[eta_lnDSI_CY08,eps_lnDSI_CY08]=get_interintraeventterms(z_lnDSI_CY08_ok,EQID_ok,sigma_DSI_CY08_ok(:,1)',sigma_DSI_CY08_ok(:,2)',sigma_DSI_CY08_ok(:,3)',alleta);
[eta_lnDSI_CB08,eps_lnDSI_CB08]=get_interintraeventterms(z_lnDSI_CB08_ok,EQID_ok,sigma_DSI_CB08_ok(:,1)',sigma_DSI_CB08_ok(:,2)',sigma_DSI_CB08_ok(:,3)',alleta);
[eta_lnDSI_AS08,eps_lnDSI_AS08]=get_interintraeventterms(z_lnDSI_AS08_ok,EQID_ok,sigma_DSI_AS08_ok(:,1)',sigma_DSI_AS08_ok(:,2)',sigma_DSI_AS08_ok(:,3)',alleta);
sigmaintraT_DSI_BA08=mean(sigma_DSI_BA08_ok(:,3)./sigma_DSI_BA08_ok(:,1));
sigmaintraT_DSI_CY08=mean(sigma_DSI_CY08_ok(:,3)./sigma_DSI_CY08_ok(:,1));
sigmaintraT_DSI_CB08=mean(sigma_DSI_CB08_ok(:,3)./sigma_DSI_CB08_ok(:,1));
sigmaintraT_DSI_AS08=mean(sigma_DSI_AS08_ok(:,3)./sigma_DSI_AS08_ok(:,1));
%     ------------ASI----------------
[eta_lnASI_BA08,eps_lnASI_BA08]=get_interintraeventterms(z_lnASI_BA08_ok,EQID_ok,sigma_ASI_BA08_ok(:,1)',sigma_ASI_BA08_ok(:,2)',sigma_ASI_BA08_ok(:,3)',alleta);
[eta_lnASI_CY08,eps_lnASI_CY08]=get_interintraeventterms(z_lnASI_CY08_ok,EQID_ok,sigma_ASI_CY08_ok(:,1)',sigma_ASI_CY08_ok(:,2)',sigma_ASI_CY08_ok(:,3)',alleta);
[eta_lnASI_CB08,eps_lnASI_CB08]=get_interintraeventterms(z_lnASI_CB08_ok,EQID_ok,sigma_ASI_CB08_ok(:,1)',sigma_ASI_CB08_ok(:,2)',sigma_ASI_CB08_ok(:,3)',alleta);
[eta_lnASI_AS08,eps_lnASI_AS08]=get_interintraeventterms(z_lnASI_AS08_ok,EQID_ok,sigma_ASI_AS08_ok(:,1)',sigma_ASI_AS08_ok(:,2)',sigma_ASI_AS08_ok(:,3)',alleta);
sigmaintraT_ASI_BA08=mean(sigma_ASI_BA08_ok(:,3)./sigma_ASI_BA08_ok(:,1));
sigmaintraT_ASI_CY08=mean(sigma_ASI_CY08_ok(:,3)./sigma_ASI_CY08_ok(:,1));
sigmaintraT_ASI_CB08=mean(sigma_ASI_CB08_ok(:,3)./sigma_ASI_CB08_ok(:,1));
sigmaintraT_ASI_AS08=mean(sigma_ASI_AS08_ok(:,3)./sigma_ASI_AS08_ok(:,1));
%     ------------SI-----------------
[eta_lnSI_BA08,eps_lnSI_BA08]=get_interintraeventterms(z_lnSI_BA08_ok,EQID_ok,sigma_SI_BA08_ok(:,1)',sigma_SI_BA08_ok(:,2)',sigma_SI_BA08_ok(:,3)',alleta);
[eta_lnSI_CY08,eps_lnSI_CY08]=get_interintraeventterms(z_lnSI_CY08_ok,EQID_ok,sigma_SI_CY08_ok(:,1)',sigma_SI_CY08_ok(:,2)',sigma_SI_CY08_ok(:,3)',alleta);
[eta_lnSI_CB08,eps_lnSI_CB08]=get_interintraeventterms(z_lnSI_CB08_ok,EQID_ok,sigma_SI_CB08_ok(:,1)',sigma_SI_CB08_ok(:,2)',sigma_SI_CB08_ok(:,3)',alleta);
[eta_lnSI_AS08,eps_lnSI_AS08]=get_interintraeventterms(z_lnSI_AS08_ok,EQID_ok,sigma_SI_AS08_ok(:,1)',sigma_SI_AS08_ok(:,2)',sigma_SI_AS08_ok(:,3)',alleta);
sigmaintraT_SI_BA08=mean(sigma_SI_BA08_ok(:,3)./sigma_SI_BA08_ok(:,1));
sigmaintraT_SI_CY08=mean(sigma_SI_CY08_ok(:,3)./sigma_SI_CY08_ok(:,1));
sigmaintraT_SI_CB08=mean(sigma_SI_CB08_ok(:,3)./sigma_SI_CB08_ok(:,1));
sigmaintraT_SI_AS08=mean(sigma_SI_AS08_ok(:,3)./sigma_SI_AS08_ok(:,1));
%     ------------PGA----------------
[eta_lnPGA_BA08,eps_lnPGA_BA08]=get_interintraeventterms(z_lnPGA_BA08_ok,EQID_ok,sigma_PGA_BA08_ok(:,1)',sigma_PGA_BA08_ok(:,2)',sigma_PGA_BA08_ok(:,3)',alleta);
[eta_lnPGA_CY08,eps_lnPGA_CY08]=get_interintraeventterms(z_lnPGA_CY08_ok,EQID_ok,sigma_PGA_CY08_ok(:,1)',sigma_PGA_CY08_ok(:,2)',sigma_PGA_CY08_ok(:,3)',alleta);
[eta_lnPGA_CB08,eps_lnPGA_CB08]=get_interintraeventterms(z_lnPGA_CB08_ok,EQID_ok,sigma_PGA_CB08_ok(:,1)',sigma_PGA_CB08_ok(:,2)',sigma_PGA_CB08_ok(:,3)',alleta);
[eta_lnPGA_AS08,eps_lnPGA_AS08]=get_interintraeventterms(z_lnPGA_AS08_ok,EQID_ok,sigma_PGA_AS08_ok(:,1)',sigma_PGA_AS08_ok(:,2)',sigma_PGA_AS08_ok(:,3)',alleta);
sigmaintraT_PGA_BA08=mean(sigma_PGA_BA08_ok(:,3)./sigma_PGA_BA08_ok(:,1));
sigmaintraT_PGA_CY08=mean(sigma_PGA_CY08_ok(:,3)./sigma_PGA_CY08_ok(:,1));
sigmaintraT_PGA_CB08=mean(sigma_PGA_CB08_ok(:,3)./sigma_PGA_CB08_ok(:,1));
sigmaintraT_PGA_AS08=mean(sigma_PGA_AS08_ok(:,3)./sigma_PGA_AS08_ok(:,1));
%     ------------PGV----------------
[eta_lnPGV_BA08,eps_lnPGV_BA08]=get_interintraeventterms(z_lnPGV_BA08_ok,EQID_ok,sigma_PGV_BA08_ok(:,1)',sigma_PGV_BA08_ok(:,2)',sigma_PGV_BA08_ok(:,3)',alleta);
[eta_lnPGV_CY08,eps_lnPGV_CY08]=get_interintraeventterms(z_lnPGV_CY08_ok,EQID_ok,sigma_PGV_CY08_ok(:,1)',sigma_PGV_CY08_ok(:,2)',sigma_PGV_CY08_ok(:,3)',alleta);
[eta_lnPGV_CB08,eps_lnPGV_CB08]=get_interintraeventterms(z_lnPGV_CB08_ok,EQID_ok,sigma_PGV_CB08_ok(:,1)',sigma_PGV_CB08_ok(:,2)',sigma_PGV_CB08_ok(:,3)',alleta);
[eta_lnPGV_AS08,eps_lnPGV_AS08]=get_interintraeventterms(z_lnPGV_AS08_ok,EQID_ok,sigma_PGV_AS08_ok(:,1)',sigma_PGV_AS08_ok(:,2)',sigma_PGV_AS08_ok(:,3)',alleta);
sigmaintraT_PGV_BA08=mean(sigma_PGV_BA08_ok(:,3)./sigma_PGV_BA08_ok(:,1));
sigmaintraT_PGV_CY08=mean(sigma_PGV_CY08_ok(:,3)./sigma_PGV_CY08_ok(:,1));
sigmaintraT_PGV_CB08=mean(sigma_PGV_CB08_ok(:,3)./sigma_PGV_CB08_ok(:,1));
sigmaintraT_PGV_AS08=mean(sigma_PGV_AS08_ok(:,3)./sigma_PGV_AS08_ok(:,1));
%     ------------SA----------------
%done in next section now 
%%%%%%%%%%%%%%%%%%%%%%end of computations %%%%%%%%%%%%%%%%%%%%


%1) Correlation between DSI and ASI
%----------------------------------
%intra
[rho_epsDSIASI_BA08] = bootstrp(nboot,@corr,eps_lnDSI_BA08,eps_lnASI_BA08);
[rho_epsDSIASI_CY08] = bootstrp(nboot,@corr,eps_lnDSI_CY08,eps_lnASI_CY08);
[rho_epsDSIASI_CB08] = bootstrp(nboot,@corr,eps_lnDSI_CB08,eps_lnASI_CB08);
[rho_epsDSIASI_AS08] = bootstrp(nboot,@corr,eps_lnDSI_AS08,eps_lnASI_AS08);
%inter
[rho_etaDSIASI_BA08] = bootstrp(nboot,@corr,eta_lnDSI_BA08,eta_lnASI_BA08);
[rho_etaDSIASI_CY08] = bootstrp(nboot,@corr,eta_lnDSI_CY08,eta_lnASI_CY08);
[rho_etaDSIASI_CB08] = bootstrp(nboot,@corr,eta_lnDSI_CB08,eta_lnASI_CB08);
[rho_etaDSIASI_AS08] = bootstrp(nboot,@corr,eta_lnDSI_AS08,eta_lnASI_AS08);
%total - incorrect computation
[rho_zDSIASI_BA08] = bootstrp(nboot,@corr,z_lnDSI_BA08_ok,z_lnASI_BA08_ok);
[rho_zDSIASI_CY08] = bootstrp(nboot,@corr,z_lnDSI_CY08_ok,z_lnASI_CY08_ok);
[rho_zDSIASI_CB08] = bootstrp(nboot,@corr,z_lnDSI_CB08_ok,z_lnASI_CB08_ok);
[rho_zDSIASI_AS08] = bootstrp(nboot,@corr,z_lnDSI_AS08_ok,z_lnASI_AS08_ok);
%total - correct computation
rho_DSIASI_BA08 = rho_etaDSIASI_BA08*sqrt((1-sigmaintraT_DSI_BA08^2)*(1-sigmaintraT_ASI_BA08^2))+rho_epsDSIASI_BA08*sigmaintraT_DSI_BA08*sigmaintraT_ASI_BA08;
rho_DSIASI_CY08 = rho_etaDSIASI_CY08*sqrt((1-sigmaintraT_DSI_CY08^2)*(1-sigmaintraT_ASI_CY08^2))+rho_epsDSIASI_CY08*sigmaintraT_DSI_CY08*sigmaintraT_ASI_CY08;
rho_DSIASI_CB08 = rho_etaDSIASI_CB08*sqrt((1-sigmaintraT_DSI_CB08^2)*(1-sigmaintraT_ASI_CB08^2))+rho_epsDSIASI_CB08*sigmaintraT_DSI_CB08*sigmaintraT_ASI_CB08;
rho_DSIASI_AS08 = rho_etaDSIASI_AS08*sqrt((1-sigmaintraT_DSI_AS08^2)*(1-sigmaintraT_ASI_AS08^2))+rho_epsDSIASI_AS08*sigmaintraT_DSI_AS08*sigmaintraT_ASI_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(1);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaDSIASI_BA08 rho_etaDSIASI_CY08 rho_etaDSIASI_CB08 rho_etaDSIASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnDSI,lnASI}'); grid on;
% subplot(222);
% boxplot([rho_epsDSIASI_BA08 rho_epsDSIASI_CY08 rho_epsDSIASI_CB08 rho_epsDSIASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnDSI,lnASI}'); grid on;
% subplot(223);
% boxplot([rho_zDSIASI_BA08 rho_zDSIASI_CY08 rho_zDSIASI_CB08 rho_zDSIASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnDSI,lnASI}'); grid on;
% subplot(224);
% boxplot([rho_DSIASI_BA08 rho_DSIASI_CY08 rho_DSIASI_CB08 rho_DSIASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnDSI,lnASI}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

%2) Correlation between DSI and PGA
%----------------------------------
%intra
[rho_epsDSIPGA_BA08] = bootstrp(nboot,@corr,eps_lnDSI_BA08,eps_lnPGA_BA08);
[rho_epsDSIPGA_CY08] = bootstrp(nboot,@corr,eps_lnDSI_CY08,eps_lnPGA_CY08);
[rho_epsDSIPGA_CB08] = bootstrp(nboot,@corr,eps_lnDSI_CB08,eps_lnPGA_CB08);
[rho_epsDSIPGA_AS08] = bootstrp(nboot,@corr,eps_lnDSI_AS08,eps_lnPGA_AS08);
%inter
[rho_etaDSIPGA_BA08] = bootstrp(nboot,@corr,eta_lnDSI_BA08,eta_lnPGA_BA08);
[rho_etaDSIPGA_CY08] = bootstrp(nboot,@corr,eta_lnDSI_CY08,eta_lnPGA_CY08);
[rho_etaDSIPGA_CB08] = bootstrp(nboot,@corr,eta_lnDSI_CB08,eta_lnPGA_CB08);
[rho_etaDSIPGA_AS08] = bootstrp(nboot,@corr,eta_lnDSI_AS08,eta_lnPGA_AS08);
%total - incorrect computation
[rho_zDSIPGA_BA08] = bootstrp(nboot,@corr,z_lnDSI_BA08_ok,z_lnPGA_BA08_ok);
[rho_zDSIPGA_CY08] = bootstrp(nboot,@corr,z_lnDSI_CY08_ok,z_lnPGA_CY08_ok);
[rho_zDSIPGA_CB08] = bootstrp(nboot,@corr,z_lnDSI_CB08_ok,z_lnPGA_CB08_ok);
[rho_zDSIPGA_AS08] = bootstrp(nboot,@corr,z_lnDSI_AS08_ok,z_lnPGA_AS08_ok);
%total - correct computation
rho_DSIPGA_BA08 = rho_etaDSIPGA_BA08*sqrt((1-sigmaintraT_DSI_BA08^2)*(1-sigmaintraT_PGA_BA08^2))+rho_epsDSIPGA_BA08*sigmaintraT_DSI_BA08*sigmaintraT_PGA_BA08;
rho_DSIPGA_CY08 = rho_etaDSIPGA_CY08*sqrt((1-sigmaintraT_DSI_CY08^2)*(1-sigmaintraT_PGA_CY08^2))+rho_epsDSIPGA_CY08*sigmaintraT_DSI_CY08*sigmaintraT_PGA_CY08;
rho_DSIPGA_CB08 = rho_etaDSIPGA_CB08*sqrt((1-sigmaintraT_DSI_CB08^2)*(1-sigmaintraT_PGA_CB08^2))+rho_epsDSIPGA_CB08*sigmaintraT_DSI_CB08*sigmaintraT_PGA_CB08;
rho_DSIPGA_AS08 = rho_etaDSIPGA_AS08*sqrt((1-sigmaintraT_DSI_AS08^2)*(1-sigmaintraT_PGA_AS08^2))+rho_epsDSIPGA_AS08*sigmaintraT_DSI_AS08*sigmaintraT_PGA_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(100);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaDSIPGA_BA08 rho_etaDSIPGA_CY08 rho_etaDSIPGA_CB08 rho_etaDSIPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnDSI,lnPGA}'); grid on;
% subplot(222);
% boxplot([rho_epsDSIPGA_BA08 rho_epsDSIPGA_CY08 rho_epsDSIPGA_CB08 rho_epsDSIPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnDSI,lnPGA}'); grid on;
% subplot(223);
% boxplot([rho_zDSIPGA_BA08 rho_zDSIPGA_CY08 rho_zDSIPGA_CB08 rho_zDSIPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnDSI,lnPGA}'); grid on;
% subplot(224);
% boxplot([rho_DSIPGA_BA08 rho_DSIPGA_CY08 rho_DSIPGA_CB08 rho_DSIPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnDSI,lnPGA}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

%2) Correlation between DSI with SI
%intra
[rho_epsDSISI_BA08] = bootstrp(nboot,@corr,eps_lnDSI_BA08,eps_lnSI_BA08);
[rho_epsDSISI_CY08] = bootstrp(nboot,@corr,eps_lnDSI_CY08,eps_lnSI_CY08);
[rho_epsDSISI_CB08] = bootstrp(nboot,@corr,eps_lnDSI_CB08,eps_lnSI_CB08);
[rho_epsDSISI_AS08] = bootstrp(nboot,@corr,eps_lnDSI_AS08,eps_lnSI_AS08);
%inter
[rho_etaDSISI_BA08] = bootstrp(nboot,@corr,eta_lnDSI_BA08,eta_lnSI_BA08);
[rho_etaDSISI_CY08] = bootstrp(nboot,@corr,eta_lnDSI_CY08,eta_lnSI_CY08);
[rho_etaDSISI_CB08] = bootstrp(nboot,@corr,eta_lnDSI_CB08,eta_lnSI_CB08);
[rho_etaDSISI_AS08] = bootstrp(nboot,@corr,eta_lnDSI_AS08,eta_lnSI_AS08);
%total - incorrect computation
[rho_zDSISI_BA08] = bootstrp(nboot,@corr,z_lnDSI_BA08_ok,z_lnSI_BA08_ok);
[rho_zDSISI_CY08] = bootstrp(nboot,@corr,z_lnDSI_CY08_ok,z_lnSI_CY08_ok);
[rho_zDSISI_CB08] = bootstrp(nboot,@corr,z_lnDSI_CB08_ok,z_lnSI_CB08_ok);
[rho_zDSISI_AS08] = bootstrp(nboot,@corr,z_lnDSI_AS08_ok,z_lnSI_AS08_ok);
%total - correct computation
rho_DSISI_BA08 = rho_etaDSISI_BA08*sqrt((1-sigmaintraT_DSI_BA08^2)*(1-sigmaintraT_SI_BA08^2))+rho_epsDSISI_BA08*sigmaintraT_DSI_BA08*sigmaintraT_SI_BA08;
rho_DSISI_CY08 = rho_etaDSISI_CY08*sqrt((1-sigmaintraT_DSI_CY08^2)*(1-sigmaintraT_SI_CY08^2))+rho_epsDSISI_CY08*sigmaintraT_DSI_CY08*sigmaintraT_SI_CY08;
rho_DSISI_CB08 = rho_etaDSISI_CB08*sqrt((1-sigmaintraT_DSI_CB08^2)*(1-sigmaintraT_SI_CB08^2))+rho_epsDSISI_CB08*sigmaintraT_DSI_CB08*sigmaintraT_SI_CB08;
rho_DSISI_AS08 = rho_etaDSISI_AS08*sqrt((1-sigmaintraT_DSI_AS08^2)*(1-sigmaintraT_SI_AS08^2))+rho_epsDSISI_AS08*sigmaintraT_DSI_AS08*sigmaintraT_SI_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(1);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaDSISI_BA08 rho_etaDSISI_CY08 rho_etaDSISI_CB08 rho_etaDSISI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnDSI,lnSI}'); grid on;
% subplot(222);
% boxplot([rho_epsDSISI_BA08 rho_epsDSISI_CY08 rho_epsDSISI_CB08 rho_epsDSISI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnDSI,lnSI}'); grid on;
% subplot(223);
% boxplot([rho_zDSISI_BA08 rho_zDSISI_CY08 rho_zDSISI_CB08 rho_zDSISI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnDSI,lnSI}'); grid on;
% subplot(224);
% boxplot([rho_DSISI_BA08 rho_DSISI_CY08 rho_DSISI_CB08 rho_DSISI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnDSI,lnSI}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

%1) Correlation between DSI and PGV
%--------------------------------------
%intra
[rho_epsDSIPGV_BA08] = bootstrp(nboot,@corr,eps_lnDSI_BA08,eps_lnPGV_BA08);
[rho_epsDSIPGV_CY08] = bootstrp(nboot,@corr,eps_lnDSI_CY08,eps_lnPGV_CY08);
[rho_epsDSIPGV_CB08] = bootstrp(nboot,@corr,eps_lnDSI_CB08,eps_lnPGV_CB08);
[rho_epsDSIPGV_AS08] = bootstrp(nboot,@corr,eps_lnDSI_AS08,eps_lnPGV_AS08);
%inter
[rho_etaDSIPGV_BA08] = bootstrp(nboot,@corr,eta_lnDSI_BA08,eta_lnPGV_BA08);
[rho_etaDSIPGV_CY08] = bootstrp(nboot,@corr,eta_lnDSI_CY08,eta_lnPGV_CY08);
[rho_etaDSIPGV_CB08] = bootstrp(nboot,@corr,eta_lnDSI_CB08,eta_lnPGV_CB08);
[rho_etaDSIPGV_AS08] = bootstrp(nboot,@corr,eta_lnDSI_AS08,eta_lnPGV_AS08);
%total - incorrect computation
[rho_zDSIPGV_BA08] = bootstrp(nboot,@corr,z_lnDSI_BA08_ok,z_lnPGV_BA08_ok);
[rho_zDSIPGV_CY08] = bootstrp(nboot,@corr,z_lnDSI_CY08_ok,z_lnPGV_CY08_ok);
[rho_zDSIPGV_CB08] = bootstrp(nboot,@corr,z_lnDSI_CB08_ok,z_lnPGV_CB08_ok);
[rho_zDSIPGV_AS08] = bootstrp(nboot,@corr,z_lnDSI_AS08_ok,z_lnPGV_AS08_ok);
%total - correct computation
rho_DSIPGV_BA08 = rho_etaDSIPGV_BA08*sqrt((1-sigmaintraT_DSI_BA08^2)*(1-sigmaintraT_PGV_BA08^2))+rho_epsDSIPGV_BA08*sigmaintraT_DSI_BA08*sigmaintraT_PGV_BA08;
rho_DSIPGV_CY08 = rho_etaDSIPGV_CY08*sqrt((1-sigmaintraT_DSI_CY08^2)*(1-sigmaintraT_PGV_CY08^2))+rho_epsDSIPGV_CY08*sigmaintraT_DSI_CY08*sigmaintraT_PGV_CY08;
rho_DSIPGV_CB08 = rho_etaDSIPGV_CB08*sqrt((1-sigmaintraT_DSI_CB08^2)*(1-sigmaintraT_PGV_CB08^2))+rho_epsDSIPGV_CB08*sigmaintraT_DSI_CB08*sigmaintraT_PGV_CB08;
rho_DSIPGV_AS08 = rho_etaDSIPGV_AS08*sqrt((1-sigmaintraT_DSI_AS08^2)*(1-sigmaintraT_PGV_AS08^2))+rho_epsDSIPGV_AS08*sigmaintraT_DSI_AS08*sigmaintraT_PGV_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(1);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaDSIPGV_BA08 rho_etaDSIPGV_CY08 rho_etaDSIPGV_CB08 rho_etaDSIPGV_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnDSI,lnPGV}'); grid on;
% subplot(222);
% boxplot([rho_epsDSIPGV_BA08 rho_epsDSIPGV_CY08 rho_epsDSIPGV_CB08 rho_epsDSIPGV_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnDSI,lnPGV}'); grid on;
% subplot(223);
% boxplot([rho_zDSIPGV_BA08 rho_zDSIPGV_CY08 rho_zDSIPGV_CB08 rho_zDSIPGV_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnDSI,lnPGV}'); grid on;
% subplot(224);
% boxplot([rho_DSIPGV_BA08 rho_DSIPGV_CY08 rho_DSIPGV_CB08 rho_DSIPGV_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnDSI,lnPGV}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

clear EQID_ok;
clear z_lnDSI_BA08_ok z_lnDSI_CY08_ok z_lnDSI_CB08_ok z_lnDSI_AS08_ok;
clear sigma_DSI_BA08_ok sigma_DSI_CY08_ok sigma_DSI_CB08_ok sigma_DSI_AS08_ok;
%4) Correlation between DSI with SA(T)
for i=1:length(T)
    %need to check if each residual is not -100 indicating above fmin
    k=0;
    for j=1:length(z_lnSA_BA08)  %will be same for both BA08 and CY08 so dont need to check twice
        if (z_lnSA_BA08(j,i)~=-100)&(z_lnDSI_BA08(j)~=-100)
            k=k+1;
            
            EQID_ok(k)=EQID(j);
            z_lnSA_BA08_ok(k)=z_lnSA_BA08(j,i);
            sigma_SA_BA08_ok(k,1:3)=sigma_SA_BA08(j,1:3,i);
            z_lnSA_CY08_ok(k)=z_lnSA_CY08(j,i);
            sigma_SA_CY08_ok(k,1:3)=sigma_SA_CY08(j,1:3,i);
            z_lnSA_CB08_ok(k)=z_lnSA_CB08(j,i);
            sigma_SA_CB08_ok(k,1:3)=sigma_SA_CB08(j,1:3,i);
            z_lnSA_AS08_ok(k)=z_lnSA_AS08(j,i);
            sigma_SA_AS08_ok(k,1:3)=sigma_SA_AS08(j,1:3,i);
            
            z_lnDSI_BA08_ok(k)=z_lnDSI_BA08(j);
            sigma_DSI_BA08_ok(k,1:3)=sigma_DSI_BA08(j,1:3);
            z_lnDSI_CY08_ok(k)=z_lnDSI_CY08(j);
            sigma_DSI_CY08_ok(k,1:3)=sigma_DSI_CY08(j,1:3);
            z_lnDSI_CB08_ok(k)=z_lnDSI_CB08(j);
            sigma_DSI_CB08_ok(k,1:3)=sigma_DSI_CB08(j,1:3);
            z_lnDSI_AS08_ok(k)=z_lnDSI_AS08(j);
            sigma_DSI_AS08_ok(k,1:3)=sigma_DSI_AS08(j,1:3);
        end
    end
    
    [eta_lnSA_BA08_ok,eps_lnSA_BA08_ok]=get_interintraeventterms(z_lnSA_BA08_ok,EQID_ok,sigma_SA_BA08_ok(:,1)',sigma_SA_BA08_ok(:,2)',sigma_SA_BA08_ok(:,3)',alleta);
    [eta_lnSA_CY08_ok,eps_lnSA_CY08_ok]=get_interintraeventterms(z_lnSA_CY08_ok,EQID_ok,sigma_SA_CY08_ok(:,1)',sigma_SA_CY08_ok(:,2)',sigma_SA_CY08_ok(:,3)',alleta);
    [eta_lnSA_CB08_ok,eps_lnSA_CB08_ok]=get_interintraeventterms(z_lnSA_CB08_ok,EQID_ok,sigma_SA_CB08_ok(:,1)',sigma_SA_CB08_ok(:,2)',sigma_SA_CB08_ok(:,3)',alleta);
    [eta_lnSA_AS08_ok,eps_lnSA_AS08_ok]=get_interintraeventterms(z_lnSA_AS08_ok,EQID_ok,sigma_SA_AS08_ok(:,1)',sigma_SA_AS08_ok(:,2)',sigma_SA_AS08_ok(:,3)',alleta);
    sigmaintraT_SA_BA08_ok=mean(sigma_SA_BA08_ok(:,3)./sigma_SA_BA08_ok(:,1));
    sigmaintraT_SA_CY08_ok=mean(sigma_SA_CY08_ok(:,3)./sigma_SA_CY08_ok(:,1));
    sigmaintraT_SA_CB08_ok=mean(sigma_SA_CB08_ok(:,3)./sigma_SA_CB08_ok(:,1));
    sigmaintraT_SA_AS08_ok=mean(sigma_SA_AS08_ok(:,3)./sigma_SA_AS08_ok(:,1));
    [eta_lnDSI_BA08_ok,eps_lnDSI_BA08_ok]=get_interintraeventterms(z_lnDSI_BA08_ok,EQID_ok,sigma_DSI_BA08_ok(:,1)',sigma_DSI_BA08_ok(:,2)',sigma_DSI_BA08_ok(:,3)',alleta);
    [eta_lnDSI_CY08_ok,eps_lnDSI_CY08_ok]=get_interintraeventterms(z_lnDSI_CY08_ok,EQID_ok,sigma_DSI_CY08_ok(:,1)',sigma_DSI_CY08_ok(:,2)',sigma_DSI_CY08_ok(:,3)',alleta);
    [eta_lnDSI_CB08_ok,eps_lnDSI_CB08_ok]=get_interintraeventterms(z_lnDSI_CB08_ok,EQID_ok,sigma_DSI_CB08_ok(:,1)',sigma_DSI_CB08_ok(:,2)',sigma_DSI_CB08_ok(:,3)',alleta);
    [eta_lnDSI_AS08_ok,eps_lnDSI_AS08_ok]=get_interintraeventterms(z_lnDSI_AS08_ok,EQID_ok,sigma_DSI_AS08_ok(:,1)',sigma_DSI_AS08_ok(:,2)',sigma_DSI_AS08_ok(:,3)',alleta);
    sigmaintraT_DSI_BA08_ok=mean(sigma_DSI_BA08_ok(:,3)./sigma_DSI_BA08_ok(:,1));
    sigmaintraT_DSI_CY08_ok=mean(sigma_DSI_CY08_ok(:,3)./sigma_DSI_CY08_ok(:,1));
    sigmaintraT_DSI_CB08_ok=mean(sigma_DSI_CB08_ok(:,3)./sigma_DSI_CB08_ok(:,1));
    sigmaintraT_DSI_AS08_ok=mean(sigma_DSI_AS08_ok(:,3)./sigma_DSI_AS08_ok(:,1));
            
    %intra
    [rho_epsDSISA_BA08(:,i)] = bootstrp(nboot,@corr,eps_lnDSI_BA08_ok,eps_lnSA_BA08_ok);
    [rho_epsDSISA_CY08(:,i)] = bootstrp(nboot,@corr,eps_lnDSI_CY08_ok,eps_lnSA_CY08_ok);
    [rho_epsDSISA_CB08(:,i)] = bootstrp(nboot,@corr,eps_lnDSI_CB08_ok,eps_lnSA_CB08_ok);
    [rho_epsDSISA_AS08(:,i)] = bootstrp(nboot,@corr,eps_lnDSI_AS08_ok,eps_lnSA_AS08_ok);
    %inter
    [rho_etaDSISA_BA08(:,i)] = bootstrp(nboot,@corr,eta_lnDSI_BA08_ok,eta_lnSA_BA08_ok);
    [rho_etaDSISA_CY08(:,i)] = bootstrp(nboot,@corr,eta_lnDSI_CY08_ok,eta_lnSA_CY08_ok);
    [rho_etaDSISA_CB08(:,i)] = bootstrp(nboot,@corr,eta_lnDSI_CB08_ok,eta_lnSA_CB08_ok);
    [rho_etaDSISA_AS08(:,i)] = bootstrp(nboot,@corr,eta_lnDSI_AS08_ok,eta_lnSA_AS08_ok);
    %total - incorrect computation
    [rho_zDSISA_BA08(:,i)] = bootstrp(nboot,@corr,z_lnDSI_BA08_ok,z_lnSA_BA08_ok);
    [rho_zDSISA_CY08(:,i)] = bootstrp(nboot,@corr,z_lnDSI_CY08_ok,z_lnSA_CY08_ok);
    [rho_zDSISA_CB08(:,i)] = bootstrp(nboot,@corr,z_lnDSI_CB08_ok,z_lnSA_CB08_ok);
    [rho_zDSISA_AS08(:,i)] = bootstrp(nboot,@corr,z_lnDSI_AS08_ok,z_lnSA_AS08_ok);
    %total - correct computation
    rho_DSISA_BA08(:,i) = rho_etaDSISA_BA08(:,i)*sqrt((1-sigmaintraT_DSI_BA08_ok^2)*(1-sigmaintraT_SA_BA08_ok^2))+rho_epsDSISA_BA08(:,i)*sigmaintraT_DSI_BA08_ok*sigmaintraT_SA_BA08_ok;
    rho_DSISA_CY08(:,i) = rho_etaDSISA_CY08(:,i)*sqrt((1-sigmaintraT_DSI_CY08_ok^2)*(1-sigmaintraT_SA_CY08_ok^2))+rho_epsDSISA_CY08(:,i)*sigmaintraT_DSI_CY08_ok*sigmaintraT_SA_CY08_ok;
    rho_DSISA_CB08(:,i) = rho_etaDSISA_CB08(:,i)*sqrt((1-sigmaintraT_DSI_CB08_ok^2)*(1-sigmaintraT_SA_CB08_ok^2))+rho_epsDSISA_CB08(:,i)*sigmaintraT_DSI_CB08_ok*sigmaintraT_SA_CB08_ok;
    rho_DSISA_AS08(:,i) = rho_etaDSISA_AS08(:,i)*sqrt((1-sigmaintraT_DSI_AS08_ok^2)*(1-sigmaintraT_SA_AS08_ok^2))+rho_epsDSISA_AS08(:,i)*sigmaintraT_DSI_AS08_ok*sigmaintraT_SA_AS08_ok;
    
    %comparison of correlations (inter, intra, z, total)
%     Tnum=10;
%     if i==Tnum
%         fig1=figure(1);
%         axes('Parent',gcf,'FontSize',16);
%         subplot(221);
%         boxplot([rho_etaDSISA_BA08(:,Tnum) rho_etaDSISA_CY08(:,Tnum) rho_etaDSISA_CB08(:,Tnum) rho_etaDSISA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
%         h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
%         ylabel('Corr coeff, \rho_{eta,lnDSI,lnSA}'); grid on;
%         subplot(222);
%         boxplot([rho_epsDSISA_BA08(:,Tnum) rho_epsDSISA_CY08(:,Tnum) rho_epsDSISA_CB08(:,Tnum) rho_epsDSISA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
%         h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
%         ylabel('Corr coeff, \rho_{eps,lnDSI,lnSA}'); grid on;
%         subplot(223);
%         boxplot([rho_zDSISA_BA08(:,Tnum) rho_zDSISA_CY08(:,Tnum) rho_zDSISA_CB08(:,Tnum) rho_zDSISA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
%         h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
%         ylabel('Corr coeff, \rho_{z,lnDSI,lnSA}'); grid on;
%         subplot(224);
%         boxplot([rho_DSISA_BA08(:,Tnum) rho_DSISA_CY08(:,Tnum) rho_DSISA_CB08(:,Tnum) rho_DSISA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
%         h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
%         ylabel('Corr coeff, \rho_{lnDSI,lnSA}'); grid on;
%         set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
%         a=1;
%     end    
end


%Transformation of rho to fisher z values
z_rho_DSIASI(:,1)=FisherZTransformation(rho_DSIASI_BA08);
z_rho_DSIASI(:,2)=FisherZTransformation(rho_DSIASI_CY08);
z_rho_DSIASI(:,3)=FisherZTransformation(rho_DSIASI_CB08);
z_rho_DSIASI(:,4)=FisherZTransformation(rho_DSIASI_AS08);
% 
z_rho_DSIPGA(:,1)=FisherZTransformation(rho_DSIPGA_BA08);
z_rho_DSIPGA(:,2)=FisherZTransformation(rho_DSIPGA_CY08);
z_rho_DSIPGA(:,3)=FisherZTransformation(rho_DSIPGA_CB08);
z_rho_DSIPGA(:,4)=FisherZTransformation(rho_DSIPGA_AS08);
% 
z_rho_DSISI(:,1)=FisherZTransformation(rho_DSISI_BA08);
z_rho_DSISI(:,2)=FisherZTransformation(rho_DSISI_CY08);
z_rho_DSISI(:,3)=FisherZTransformation(rho_DSISI_CB08);
z_rho_DSISI(:,4)=FisherZTransformation(rho_DSISI_AS08);

z_rho_DSIPGV(:,1)=FisherZTransformation(rho_DSIPGV_BA08);
z_rho_DSIPGV(:,2)=FisherZTransformation(rho_DSIPGV_CY08);
z_rho_DSIPGV(:,3)=FisherZTransformation(rho_DSIPGV_CB08);
z_rho_DSIPGV(:,4)=FisherZTransformation(rho_DSIPGV_AS08);
% 
z_rho_DSIASI_all = [z_rho_DSIASI(:,1)' z_rho_DSIASI(:,2)' z_rho_DSIASI(:,3)' z_rho_DSIASI(:,4)']';
z_rho_DSIPGA_all = [z_rho_DSIPGA(:,1)' z_rho_DSIPGA(:,2)' z_rho_DSIPGA(:,3)' z_rho_DSIPGA(:,4)']';
z_rho_DSISI_all = [z_rho_DSISI(:,1)' z_rho_DSISI(:,2)' z_rho_DSISI(:,3)' z_rho_DSISI(:,4)']';
z_rho_DSIPGV_all = [z_rho_DSIPGV(:,1)' z_rho_DSIPGV(:,2)' z_rho_DSIPGV(:,3)' z_rho_DSIPGV(:,4)']';

for i=1:length(T)
    z_rho_DSISA(:,i,1)=FisherZTransformation(rho_DSISA_BA08(:,i));
    z_rho_DSISA(:,i,2)=FisherZTransformation(rho_DSISA_CY08(:,i));
    z_rho_DSISA(:,i,3)=FisherZTransformation(rho_DSISA_CB08(:,i));
    z_rho_DSISA(:,i,4)=FisherZTransformation(rho_DSISA_AS08(:,i));
    
    z_rho_DSISA_all(:,i) = [z_rho_DSISA(:,i,1)' z_rho_DSISA(:,i,2)' z_rho_DSISA(:,i,3)' z_rho_DSISA(:,i,4)']';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputting of results
median_rho_DSIASI=InvFisherZTransformation(mean(z_rho_DSIASI,1));
median_rho_DSIPGA=InvFisherZTransformation(mean(z_rho_DSIPGA,1));
median_rho_DSISI=InvFisherZTransformation(mean(z_rho_DSISI,1));
median_rho_DSIPGV=InvFisherZTransformation(mean(z_rho_DSIPGV,1));

median_rho_DSIASI_all=InvFisherZTransformation(mean(z_rho_DSIASI_all));
median_rho_DSIPGA_all=InvFisherZTransformation(mean(z_rho_DSIPGA_all));
median_rho_DSISI_all=InvFisherZTransformation(mean(z_rho_DSISI_all));
median_rho_DSIPGV_all=InvFisherZTransformation(mean(z_rho_DSIPGV_all));

std_z_DSIASI=std(z_rho_DSIASI_all);
std_z_DSIPGA=std(z_rho_DSIPGA_all);
std_z_DSISI=std(z_rho_DSISI_all);
std_z_DSIPGV=std(z_rho_DSIPGV_all);

for i=1:length(T)
    median_rho_DSISA(i,:)=InvFisherZTransformation(mean(z_rho_DSISA(:,i,:),1));
    median_rho_DSISA_all(i)=InvFisherZTransformation(mean(z_rho_DSISA_all(:,i),1));
    std_z_DSISA(i)=std(z_rho_DSISA_all(:,i));
end



fprintf('Results of correlation analyses \n');
fprintf('\n');
fprintf('-----------------------  Median correlations   ---------------------- | ----------Overall   -------\n');
fprintf('             BA08   |     CY08      |     CY08      |     AS08      |   median     |  Z std   \n');
fprintf('DSI,ASI    %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_DSIASI,median_rho_DSIASI_all,std_z_DSIASI);
fprintf('DSI,PGA    %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_DSIPGA,median_rho_DSIPGA_all,std_z_DSIPGA);
fprintf('DSI,SI     %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_DSISI,median_rho_DSISI_all,std_z_DSISI);
fprintf('DSI,PGV    %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_DSIPGV,median_rho_DSIPGV_all,std_z_DSIPGV);
fprintf('\n');
fprintf('              median DSI - SA correlation       |                          |  \n');
fprintf('T      |    BA08  |    CY08 |    CB08 |    AS08 |   median   |  Z variance |  \n');
for i=1:length(T)
fprintf('%6.3f     %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   \n',...
                                                        T(i),median_rho_DSISA(i,:),median_rho_DSISA_all(i),std_z_DSISA(i));
end
fprintf('------------------------------------- \n');

%plotting
%DSI with ASI,SI, PGA, PGV correlations
fig1=figure(1);
axes('Parent',gcf,'FontSize',16);
subplot(221);
boxplot([rho_DSIASI_BA08 rho_DSIASI_CY08 rho_DSIASI_CB08 rho_DSIASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnDSI,lnASI}'); grid on;
subplot(222);
boxplot([rho_DSIPGA_BA08 rho_DSIPGA_CY08 rho_DSIPGA_CB08 rho_DSIPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnDSI,lnPGA}'); grid on;
subplot(223);
boxplot([rho_DSISI_BA08 rho_DSISI_CY08 rho_DSISI_CB08 rho_DSISI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnDSI,lnSI}'); grid on;
subplot(224);
boxplot([rho_DSIPGV_BA08 rho_DSIPGV_CY08 rho_DSIPGV_CB08 rho_DSIPGV_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnDSI,lnPGV}'); grid on;
set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.85 0.85]);
% 
% %DSI correlation with SA
fig2=figure(2);
axes('Parent',gcf,'FontSize',16);
semilogx(T,mean(rho_DSISA_BA08,1),'LineWidth',3,'LineStyle','-','Color',[1 0 0]); hold on;
semilogx(T,prctile(rho_DSISA_BA08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[1 0 0]);
semilogx(T,mean(rho_DSISA_CY08,1),'LineWidth',3,'LineStyle','-','Color',[0 1 0]);
semilogx(T,prctile(rho_DSISA_CY08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0 1 0]);
semilogx(T,mean(rho_DSISA_CB08,1),'LineWidth',3,'LineStyle','-','Color',[0 0 1]);
semilogx(T,prctile(rho_DSISA_CB08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0 0 1]);
semilogx(T,mean(rho_DSISA_AS08,1),'LineWidth',3,'LineStyle','-','Color',[0.5 0.5 0.5]);
semilogx(T,prctile(rho_DSISA_AS08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5]);
ylabel('Corr coeff, \rho_{lnDSI,lnSA}'); xlabel('Period, T (s)'); grid on;
set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.5 0.4]);

%show adequacy of z having a normal distribution
z_rho_ij=z_rho_DSIASI_all; %which is considered

cdf=0.01:0.01:0.99;
X=norminv(cdf,mean(z_rho_ij),std(z_rho_ij));

alpha=0.1;
N=length(z_rho_ij);
emp_cdf=1/(N+1):1/(N+1):N/(N+1);
[H,P,KSSTAT,CV] = kstest(sort(z_rho_ij),[sort(z_rho_ij) emp_cdf'],alpha); 

fig4=figure(4);
axes('Parent',gcf,'FontSize',16);
plot(X,cdf,'LineWidth',3,'LineStyle','-','Color',[1 0 0]); hold on; 
plot(X,cdf+CV,'Color',[1 0 0],'LineWidth',2,'LineStyle','--'); 
plot(X,cdf-CV,'Color',[1 0 0],'LineWidth',2,'LineStyle','--');
plot(sort(z_rho_ij),emp_cdf,'LineWidth',3,'LineStyle','-','Color',[0 0 1]);
xlabel('Fisher Z value'); ylabel('Cumulative probability');
ylim([0 1]);
grid on;

close(h_wait)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z]=FisherZTransformation(rho)
%Computes the Fisher Z Transformation of the correlation

z=0.5*log((1+rho)./(1-rho));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho]=InvFisherZTransformation(z)
%Computes the Inverse Fisher Z Transformation of the correlation

rho=(exp(2*z)-1)./(1+exp(2*z));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
