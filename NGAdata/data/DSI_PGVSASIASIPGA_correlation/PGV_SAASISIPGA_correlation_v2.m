function PGV_SAASISIPGA_correlation
clc
%Brendon Bradley 10 Aug 2010

%Purpose: To look at the empirical correlations between various intensity
%measures, PGV and SI, ASI, PGA, SA
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
    siteprop.orientation='average';
    
    %Intensity measures
    
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
    IMR=@BooreAtkinson_2007_nga; siteprop.g=981; %in cm/s2
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
%     ------------ASI----------------
[eta_lnASI_BA08,eps_lnASI_BA08]=get_interintraeventterms(z_lnASI_BA08,EQID,sigma_ASI_BA08(:,1)',sigma_ASI_BA08(:,2)',sigma_ASI_BA08(:,3)',alleta);
[eta_lnASI_CY08,eps_lnASI_CY08]=get_interintraeventterms(z_lnASI_CY08,EQID,sigma_ASI_CY08(:,1)',sigma_ASI_CY08(:,2)',sigma_ASI_CY08(:,3)',alleta);
[eta_lnASI_CB08,eps_lnASI_CB08]=get_interintraeventterms(z_lnASI_CB08,EQID,sigma_ASI_CB08(:,1)',sigma_ASI_CB08(:,2)',sigma_ASI_CB08(:,3)',alleta);
[eta_lnASI_AS08,eps_lnASI_AS08]=get_interintraeventterms(z_lnASI_AS08,EQID,sigma_ASI_AS08(:,1)',sigma_ASI_AS08(:,2)',sigma_ASI_AS08(:,3)',alleta);
sigmaintraT_ASI_BA08=mean(sigma_ASI_BA08(:,3)./sigma_ASI_BA08(:,1));
sigmaintraT_ASI_CY08=mean(sigma_ASI_CY08(:,3)./sigma_ASI_CY08(:,1));
sigmaintraT_ASI_CB08=mean(sigma_ASI_CB08(:,3)./sigma_ASI_CB08(:,1));
sigmaintraT_ASI_AS08=mean(sigma_ASI_AS08(:,3)./sigma_ASI_AS08(:,1));
%     ------------SI-----------------
[eta_lnSI_BA08,eps_lnSI_BA08]=get_interintraeventterms(z_lnSI_BA08,EQID,sigma_SI_BA08(:,1)',sigma_SI_BA08(:,2)',sigma_SI_BA08(:,3)',alleta);
[eta_lnSI_CY08,eps_lnSI_CY08]=get_interintraeventterms(z_lnSI_CY08,EQID,sigma_SI_CY08(:,1)',sigma_SI_CY08(:,2)',sigma_SI_CY08(:,3)',alleta);
[eta_lnSI_CB08,eps_lnSI_CB08]=get_interintraeventterms(z_lnSI_CB08,EQID,sigma_SI_CB08(:,1)',sigma_SI_CB08(:,2)',sigma_SI_CB08(:,3)',alleta);
[eta_lnSI_AS08,eps_lnSI_AS08]=get_interintraeventterms(z_lnSI_AS08,EQID,sigma_SI_AS08(:,1)',sigma_SI_AS08(:,2)',sigma_SI_AS08(:,3)',alleta);
sigmaintraT_SI_BA08=mean(sigma_SI_BA08(:,3)./sigma_SI_BA08(:,1));
sigmaintraT_SI_CY08=mean(sigma_SI_CY08(:,3)./sigma_SI_CY08(:,1));
sigmaintraT_SI_CB08=mean(sigma_SI_CB08(:,3)./sigma_SI_CB08(:,1));
sigmaintraT_SI_AS08=mean(sigma_SI_AS08(:,3)./sigma_SI_AS08(:,1));
%     ------------PGA----------------
[eta_lnPGA_BA08,eps_lnPGA_BA08]=get_interintraeventterms(z_lnPGA_BA08,EQID,sigma_PGA_BA08(:,1)',sigma_PGA_BA08(:,2)',sigma_PGA_BA08(:,3)',alleta);
[eta_lnPGA_CY08,eps_lnPGA_CY08]=get_interintraeventterms(z_lnPGA_CY08,EQID,sigma_PGA_CY08(:,1)',sigma_PGA_CY08(:,2)',sigma_PGA_CY08(:,3)',alleta);
[eta_lnPGA_CB08,eps_lnPGA_CB08]=get_interintraeventterms(z_lnPGA_CB08,EQID,sigma_PGA_CB08(:,1)',sigma_PGA_CB08(:,2)',sigma_PGA_CB08(:,3)',alleta);
[eta_lnPGA_AS08,eps_lnPGA_AS08]=get_interintraeventterms(z_lnPGA_AS08,EQID,sigma_PGA_AS08(:,1)',sigma_PGA_AS08(:,2)',sigma_PGA_AS08(:,3)',alleta);
sigmaintraT_PGA_BA08=mean(sigma_PGA_BA08(:,3)./sigma_PGA_BA08(:,1));
sigmaintraT_PGA_CY08=mean(sigma_PGA_CY08(:,3)./sigma_PGA_CY08(:,1));
sigmaintraT_PGA_CB08=mean(sigma_PGA_CB08(:,3)./sigma_PGA_CB08(:,1));
sigmaintraT_PGA_AS08=mean(sigma_PGA_AS08(:,3)./sigma_PGA_AS08(:,1));
%     ------------PGV----------------
[eta_lnPGV_BA08,eps_lnPGV_BA08]=get_interintraeventterms(z_lnPGV_BA08,EQID,sigma_PGV_BA08(:,1)',sigma_PGV_BA08(:,2)',sigma_PGV_BA08(:,3)',alleta);
[eta_lnPGV_CY08,eps_lnPGV_CY08]=get_interintraeventterms(z_lnPGV_CY08,EQID,sigma_PGV_CY08(:,1)',sigma_PGV_CY08(:,2)',sigma_PGV_CY08(:,3)',alleta);
[eta_lnPGV_CB08,eps_lnPGV_CB08]=get_interintraeventterms(z_lnPGV_CB08,EQID,sigma_PGV_CB08(:,1)',sigma_PGV_CB08(:,2)',sigma_PGV_CB08(:,3)',alleta);
[eta_lnPGV_AS08,eps_lnPGV_AS08]=get_interintraeventterms(z_lnPGV_AS08,EQID,sigma_PGV_AS08(:,1)',sigma_PGV_AS08(:,2)',sigma_PGV_AS08(:,3)',alleta);
sigmaintraT_PGV_BA08=mean(sigma_PGV_BA08(:,3)./sigma_PGV_BA08(:,1));
sigmaintraT_PGV_CY08=mean(sigma_PGV_CY08(:,3)./sigma_PGV_CY08(:,1));
sigmaintraT_PGV_CB08=mean(sigma_PGV_CB08(:,3)./sigma_PGV_CB08(:,1));
sigmaintraT_PGV_AS08=mean(sigma_PGV_AS08(:,3)./sigma_PGV_AS08(:,1));
%     ------------SA----------------
%done in next section now 
% for j=1:length(T)
%     [eta_lnSA_BA08(:,j),eps_lnSA_BA08(:,j)]=get_interintraeventterms(z_lnSA_BA08(:,j)',EQID,sigma_SA_BA08(:,1,j)',sigma_SA_BA08(:,2,j)',sigma_SA_BA08(:,3,j)',alleta);
%     [eta_lnSA_CY08(:,j),eps_lnSA_CY08(:,j)]=get_interintraeventterms(z_lnSA_CY08(:,j)',EQID,sigma_SA_CY08(:,1,j)',sigma_SA_CY08(:,2,j)',sigma_SA_CY08(:,3,j)',alleta);
%     [eta_lnSA_CB08(:,j),eps_lnSA_CB08(:,j)]=get_interintraeventterms(z_lnSA_CB08(:,j)',EQID,sigma_SA_CB08(:,1,j)',sigma_SA_CB08(:,2,j)',sigma_SA_CB08(:,3,j)',alleta);
%     [eta_lnSA_AS08(:,j),eps_lnSA_AS08(:,j)]=get_interintraeventterms(z_lnSA_AS08(:,j)',EQID,sigma_SA_AS08(:,1,j)',sigma_SA_AS08(:,2,j)',sigma_SA_AS08(:,3,j)',alleta);
%     sigmaintraT_SA_BA08(j)=mean(sigma_SA_BA08(:,3,j)./sigma_SA_BA08(:,1,j));
%     sigmaintraT_SA_CY08(j)=mean(sigma_SA_CY08(:,3,j)./sigma_SA_CY08(:,1,j));
%     sigmaintraT_SA_CB08(j)=mean(sigma_SA_CB08(:,3,j)./sigma_SA_CB08(:,1,j));
%     sigmaintraT_SA_AS08(j)=mean(sigma_SA_AS08(:,3,j)./sigma_SA_AS08(:,1,j));
% end
    
% %temp plots to check PGV calculations
% figure(1);
% ecdf(eps_lnPGA_BA08); hold on;
% ecdf(eps_lnPGA_CY08); 
% ecdf(eps_lnPGA_CB08);
% plot(norminv([0.01:0.01:0.99],0,1),[0.01:0.01:0.99],'r-');

%%%%%%%%%%%%%%%%%%%%%%end of computations %%%%%%%%%%%%%%%%%%%%

%processing of residuals
%1) Correlation between PGV and ASI
%----------------------------------
%intra
[rho_epsPGVASI_BA08] = bootstrp(nboot,@corr,eps_lnPGV_BA08,eps_lnASI_BA08);
[rho_epsPGVASI_CY08] = bootstrp(nboot,@corr,eps_lnPGV_CY08,eps_lnASI_CY08);
[rho_epsPGVASI_CB08] = bootstrp(nboot,@corr,eps_lnPGV_CB08,eps_lnASI_CB08);
[rho_epsPGVASI_AS08] = bootstrp(nboot,@corr,eps_lnPGV_AS08,eps_lnASI_AS08);
%inter
[rho_etaPGVASI_BA08] = bootstrp(nboot,@corr,eta_lnPGV_BA08,eta_lnASI_BA08);
[rho_etaPGVASI_CY08] = bootstrp(nboot,@corr,eta_lnPGV_CY08,eta_lnASI_CY08);
[rho_etaPGVASI_CB08] = bootstrp(nboot,@corr,eta_lnPGV_CB08,eta_lnASI_CB08);
[rho_etaPGVASI_AS08] = bootstrp(nboot,@corr,eta_lnPGV_AS08,eta_lnASI_AS08);
%total - incorrect computation
[rho_zPGVASI_BA08] = bootstrp(nboot,@corr,z_lnPGV_BA08,z_lnASI_BA08);
[rho_zPGVASI_CY08] = bootstrp(nboot,@corr,z_lnPGV_CY08,z_lnASI_CY08);
[rho_zPGVASI_CB08] = bootstrp(nboot,@corr,z_lnPGV_CB08,z_lnASI_CB08);
[rho_zPGVASI_AS08] = bootstrp(nboot,@corr,z_lnPGV_AS08,z_lnASI_AS08);
%total - correct computation
rho_PGVASI_BA08 = rho_etaPGVASI_BA08*sqrt((1-sigmaintraT_PGV_BA08^2)*(1-sigmaintraT_ASI_BA08^2))+rho_epsPGVASI_BA08*sigmaintraT_PGV_BA08*sigmaintraT_ASI_BA08;
rho_PGVASI_CY08 = rho_etaPGVASI_CY08*sqrt((1-sigmaintraT_PGV_CY08^2)*(1-sigmaintraT_ASI_CY08^2))+rho_epsPGVASI_CY08*sigmaintraT_PGV_CY08*sigmaintraT_ASI_CY08;
rho_PGVASI_CB08 = rho_etaPGVASI_CB08*sqrt((1-sigmaintraT_PGV_CB08^2)*(1-sigmaintraT_ASI_CB08^2))+rho_epsPGVASI_CB08*sigmaintraT_PGV_CB08*sigmaintraT_ASI_CB08;
rho_PGVASI_AS08 = rho_etaPGVASI_AS08*sqrt((1-sigmaintraT_PGV_AS08^2)*(1-sigmaintraT_ASI_AS08^2))+rho_epsPGVASI_AS08*sigmaintraT_PGV_AS08*sigmaintraT_ASI_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(1);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaPGVASI_BA08 rho_etaPGVASI_CY08 rho_etaPGVASI_CB08 rho_etaPGVASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnPGV,lnASI}'); grid on;
% subplot(222);
% boxplot([rho_epsPGVASI_BA08 rho_epsPGVASI_CY08 rho_epsPGVASI_CB08 rho_epsPGVASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnPGV,lnASI}'); grid on;
% subplot(223);
% boxplot([rho_zPGVASI_BA08 rho_zPGVASI_CY08 rho_zPGVASI_CB08 rho_zPGVASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnPGV,lnASI}'); grid on;
% subplot(224);
% boxplot([rho_PGVASI_BA08 rho_PGVASI_CY08 rho_PGVASI_CB08 rho_PGVASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnPGV,lnASI}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

%2) Correlation between PGV and PGA
%----------------------------------
%intra
[rho_epsPGVPGA_BA08] = bootstrp(nboot,@corr,eps_lnPGV_BA08,eps_lnPGA_BA08);
[rho_epsPGVPGA_CY08] = bootstrp(nboot,@corr,eps_lnPGV_CY08,eps_lnPGA_CY08);
[rho_epsPGVPGA_CB08] = bootstrp(nboot,@corr,eps_lnPGV_CB08,eps_lnPGA_CB08);
[rho_epsPGVPGA_AS08] = bootstrp(nboot,@corr,eps_lnPGV_AS08,eps_lnPGA_AS08);
%inter
[rho_etaPGVPGA_BA08] = bootstrp(nboot,@corr,eta_lnPGV_BA08,eta_lnPGA_BA08);
[rho_etaPGVPGA_CY08] = bootstrp(nboot,@corr,eta_lnPGV_CY08,eta_lnPGA_CY08);
[rho_etaPGVPGA_CB08] = bootstrp(nboot,@corr,eta_lnPGV_CB08,eta_lnPGA_CB08);
[rho_etaPGVPGA_AS08] = bootstrp(nboot,@corr,eta_lnPGV_AS08,eta_lnPGA_AS08);
%total - incorrect computation
[rho_zPGVPGA_BA08] = bootstrp(nboot,@corr,z_lnPGV_BA08,z_lnPGA_BA08);
[rho_zPGVPGA_CY08] = bootstrp(nboot,@corr,z_lnPGV_CY08,z_lnPGA_CY08);
[rho_zPGVPGA_CB08] = bootstrp(nboot,@corr,z_lnPGV_CB08,z_lnPGA_CB08);
[rho_zPGVPGA_AS08] = bootstrp(nboot,@corr,z_lnPGV_AS08,z_lnPGA_AS08);
%total - correct computation
rho_PGVPGA_BA08 = rho_etaPGVPGA_BA08*sqrt((1-sigmaintraT_PGV_BA08^2)*(1-sigmaintraT_PGA_BA08^2))+rho_epsPGVPGA_BA08*sigmaintraT_PGV_BA08*sigmaintraT_PGA_BA08;
rho_PGVPGA_CY08 = rho_etaPGVPGA_CY08*sqrt((1-sigmaintraT_PGV_CY08^2)*(1-sigmaintraT_PGA_CY08^2))+rho_epsPGVPGA_CY08*sigmaintraT_PGV_CY08*sigmaintraT_PGA_CY08;
rho_PGVPGA_CB08 = rho_etaPGVPGA_CB08*sqrt((1-sigmaintraT_PGV_CB08^2)*(1-sigmaintraT_PGA_CB08^2))+rho_epsPGVPGA_CB08*sigmaintraT_PGV_CB08*sigmaintraT_PGA_CB08;
rho_PGVPGA_AS08 = rho_etaPGVPGA_AS08*sqrt((1-sigmaintraT_PGV_AS08^2)*(1-sigmaintraT_PGA_AS08^2))+rho_epsPGVPGA_AS08*sigmaintraT_PGV_AS08*sigmaintraT_PGA_AS08;
%comparison of correlations (inter, intra, z, total)
fig1=figure(100);
axes('Parent',gcf,'FontSize',16);
subplot(221);
boxplot([rho_etaPGVPGA_BA08 rho_etaPGVPGA_CY08 rho_etaPGVPGA_CB08 rho_etaPGVPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{eta,lnPGV,lnPGA}'); grid on;
subplot(222);
boxplot([rho_epsPGVPGA_BA08 rho_epsPGVPGA_CY08 rho_epsPGVPGA_CB08 rho_epsPGVPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{eps,lnPGV,lnPGA}'); grid on;
subplot(223);
boxplot([rho_zPGVPGA_BA08 rho_zPGVPGA_CY08 rho_zPGVPGA_CB08 rho_zPGVPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{z,lnPGV,lnPGA}'); grid on;
subplot(224);
boxplot([rho_PGVPGA_BA08 rho_PGVPGA_CY08 rho_PGVPGA_CB08 rho_PGVPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnPGV,lnPGA}'); grid on;
set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
a=1;

%2) Correlation between PGV with SI
%intra
[rho_epsPGVSI_BA08] = bootstrp(nboot,@corr,eps_lnPGV_BA08,eps_lnSI_BA08);
[rho_epsPGVSI_CY08] = bootstrp(nboot,@corr,eps_lnPGV_CY08,eps_lnSI_CY08);
[rho_epsPGVSI_CB08] = bootstrp(nboot,@corr,eps_lnPGV_CB08,eps_lnSI_CB08);
[rho_epsPGVSI_AS08] = bootstrp(nboot,@corr,eps_lnPGV_AS08,eps_lnSI_AS08);
%inter
[rho_etaPGVSI_BA08] = bootstrp(nboot,@corr,eta_lnPGV_BA08,eta_lnSI_BA08);
[rho_etaPGVSI_CY08] = bootstrp(nboot,@corr,eta_lnPGV_CY08,eta_lnSI_CY08);
[rho_etaPGVSI_CB08] = bootstrp(nboot,@corr,eta_lnPGV_CB08,eta_lnSI_CB08);
[rho_etaPGVSI_AS08] = bootstrp(nboot,@corr,eta_lnPGV_AS08,eta_lnSI_AS08);
%total - incorrect computation
[rho_zPGVSI_BA08] = bootstrp(nboot,@corr,z_lnPGV_BA08,z_lnSI_BA08);
[rho_zPGVSI_CY08] = bootstrp(nboot,@corr,z_lnPGV_CY08,z_lnSI_CY08);
[rho_zPGVSI_CB08] = bootstrp(nboot,@corr,z_lnPGV_CB08,z_lnSI_CB08);
[rho_zPGVSI_AS08] = bootstrp(nboot,@corr,z_lnPGV_AS08,z_lnSI_AS08);
%total - correct computation
rho_PGVSI_BA08 = rho_etaPGVSI_BA08*sqrt((1-sigmaintraT_PGV_BA08^2)*(1-sigmaintraT_SI_BA08^2))+rho_epsPGVSI_BA08*sigmaintraT_PGV_BA08*sigmaintraT_SI_BA08;
rho_PGVSI_CY08 = rho_etaPGVSI_CY08*sqrt((1-sigmaintraT_PGV_CY08^2)*(1-sigmaintraT_SI_CY08^2))+rho_epsPGVSI_CY08*sigmaintraT_PGV_CY08*sigmaintraT_SI_CY08;
rho_PGVSI_CB08 = rho_etaPGVSI_CB08*sqrt((1-sigmaintraT_PGV_CB08^2)*(1-sigmaintraT_SI_CB08^2))+rho_epsPGVSI_CB08*sigmaintraT_PGV_CB08*sigmaintraT_SI_CB08;
rho_PGVSI_AS08 = rho_etaPGVSI_AS08*sqrt((1-sigmaintraT_PGV_AS08^2)*(1-sigmaintraT_SI_AS08^2))+rho_epsPGVSI_AS08*sigmaintraT_PGV_AS08*sigmaintraT_SI_AS08;
%comparison of correlations (inter, intra, z, total)
% fig1=figure(1);
% axes('Parent',gcf,'FontSize',16);
% subplot(221);
% boxplot([rho_etaPGVSI_BA08 rho_etaPGVSI_CY08 rho_etaPGVSI_CB08 rho_etaPGVSI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eta,lnPGV,lnSI}'); grid on;
% subplot(222);
% boxplot([rho_epsPGVSI_BA08 rho_epsPGVSI_CY08 rho_epsPGVSI_CB08 rho_epsPGVSI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{eps,lnPGV,lnSI}'); grid on;
% subplot(223);
% boxplot([rho_zPGVSI_BA08 rho_zPGVSI_CY08 rho_zPGVSI_CB08 rho_zPGVSI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{z,lnPGV,lnSI}'); grid on;
% subplot(224);
% boxplot([rho_PGVSI_BA08 rho_PGVSI_CY08 rho_PGVSI_CB08 rho_PGVSI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
% h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
% ylabel('Corr coeff, \rho_{lnPGV,lnSI}'); grid on;
% set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
% a=1;

%4) Correlation between PGV with SA(T)
for i=1:length(T)
    %need to check if each residual is not -100 indicating above fmin
    k=0;
    for j=1:length(z_lnSA_BA08)  %will be same for both BA08 and CY08 so dont need to check twice
        if z_lnSA_BA08(j,i)~=-100
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
    
    [eta_lnSA_BA08_ok,eps_lnSA_BA08_ok]=get_interintraeventterms(z_lnSA_BA08_ok,EQID_ok,sigma_SA_BA08_ok(:,1)',sigma_SA_BA08_ok(:,2)',sigma_SA_BA08_ok(:,3)',alleta);
    [eta_lnSA_CY08_ok,eps_lnSA_CY08_ok]=get_interintraeventterms(z_lnSA_CY08_ok,EQID_ok,sigma_SA_CY08_ok(:,1)',sigma_SA_CY08_ok(:,2)',sigma_SA_CY08_ok(:,3)',alleta);
    [eta_lnSA_CB08_ok,eps_lnSA_CB08_ok]=get_interintraeventterms(z_lnSA_CB08_ok,EQID_ok,sigma_SA_CB08_ok(:,1)',sigma_SA_CB08_ok(:,2)',sigma_SA_CB08_ok(:,3)',alleta);
    [eta_lnSA_AS08_ok,eps_lnSA_AS08_ok]=get_interintraeventterms(z_lnSA_AS08_ok,EQID_ok,sigma_SA_AS08_ok(:,1)',sigma_SA_AS08_ok(:,2)',sigma_SA_AS08_ok(:,3)',alleta);
    sigmaintraT_SA_BA08_ok=mean(sigma_SA_BA08_ok(:,3)./sigma_SA_BA08_ok(:,1));
    sigmaintraT_SA_CY08_ok=mean(sigma_SA_CY08_ok(:,3)./sigma_SA_CY08_ok(:,1));
    sigmaintraT_SA_CB08_ok=mean(sigma_SA_CB08_ok(:,3)./sigma_SA_CB08_ok(:,1));
    sigmaintraT_SA_AS08_ok=mean(sigma_SA_AS08_ok(:,3)./sigma_SA_AS08_ok(:,1));
    [eta_lnPGV_BA08_ok,eps_lnPGV_BA08_ok]=get_interintraeventterms(z_lnPGV_BA08_ok,EQID_ok,sigma_PGV_BA08_ok(:,1)',sigma_PGV_BA08_ok(:,2)',sigma_PGV_BA08_ok(:,3)',alleta);
    [eta_lnPGV_CY08_ok,eps_lnPGV_CY08_ok]=get_interintraeventterms(z_lnPGV_CY08_ok,EQID_ok,sigma_PGV_CY08_ok(:,1)',sigma_PGV_CY08_ok(:,2)',sigma_PGV_CY08_ok(:,3)',alleta);
    [eta_lnPGV_CB08_ok,eps_lnPGV_CB08_ok]=get_interintraeventterms(z_lnPGV_CB08_ok,EQID_ok,sigma_PGV_CB08_ok(:,1)',sigma_PGV_CB08_ok(:,2)',sigma_PGV_CB08_ok(:,3)',alleta);
    [eta_lnPGV_AS08_ok,eps_lnPGV_AS08_ok]=get_interintraeventterms(z_lnPGV_AS08_ok,EQID_ok,sigma_PGV_AS08_ok(:,1)',sigma_PGV_AS08_ok(:,2)',sigma_PGV_AS08_ok(:,3)',alleta);
    sigmaintraT_PGV_BA08_ok=mean(sigma_PGV_BA08_ok(:,3)./sigma_PGV_BA08_ok(:,1));
    sigmaintraT_PGV_CY08_ok=mean(sigma_PGV_CY08_ok(:,3)./sigma_PGV_CY08_ok(:,1));
    sigmaintraT_PGV_CB08_ok=mean(sigma_PGV_CB08_ok(:,3)./sigma_PGV_CB08_ok(:,1));
    sigmaintraT_PGV_AS08_ok=mean(sigma_PGV_AS08_ok(:,3)./sigma_PGV_AS08_ok(:,1));
            
    %intra
    [rho_epsPGVSA_BA08(:,i)] = bootstrp(nboot,@corr,eps_lnPGV_BA08_ok,eps_lnSA_BA08_ok);
    [rho_epsPGVSA_CY08(:,i)] = bootstrp(nboot,@corr,eps_lnPGV_CY08_ok,eps_lnSA_CY08_ok);
    [rho_epsPGVSA_CB08(:,i)] = bootstrp(nboot,@corr,eps_lnPGV_CB08_ok,eps_lnSA_CB08_ok);
    [rho_epsPGVSA_AS08(:,i)] = bootstrp(nboot,@corr,eps_lnPGV_AS08_ok,eps_lnSA_AS08_ok);
    %inter
    [rho_etaPGVSA_BA08(:,i)] = bootstrp(nboot,@corr,eta_lnPGV_BA08_ok,eta_lnSA_BA08_ok);
    [rho_etaPGVSA_CY08(:,i)] = bootstrp(nboot,@corr,eta_lnPGV_CY08_ok,eta_lnSA_CY08_ok);
    [rho_etaPGVSA_CB08(:,i)] = bootstrp(nboot,@corr,eta_lnPGV_CB08_ok,eta_lnSA_CB08_ok);
    [rho_etaPGVSA_AS08(:,i)] = bootstrp(nboot,@corr,eta_lnPGV_AS08_ok,eta_lnSA_AS08_ok);
    %total - incorrect computation
    [rho_zPGVSA_BA08(:,i)] = bootstrp(nboot,@corr,z_lnPGV_BA08_ok,z_lnSA_BA08_ok);
    [rho_zPGVSA_CY08(:,i)] = bootstrp(nboot,@corr,z_lnPGV_CY08_ok,z_lnSA_CY08_ok);
    [rho_zPGVSA_CB08(:,i)] = bootstrp(nboot,@corr,z_lnPGV_CB08_ok,z_lnSA_CB08_ok);
    [rho_zPGVSA_AS08(:,i)] = bootstrp(nboot,@corr,z_lnPGV_AS08_ok,z_lnSA_AS08_ok);
    %total - correct computation
    rho_PGVSA_BA08(:,i) = rho_etaPGVSA_BA08(:,i)*sqrt((1-sigmaintraT_PGV_BA08_ok^2)*(1-sigmaintraT_SA_BA08_ok^2))+rho_epsPGVSA_BA08(:,i)*sigmaintraT_PGV_BA08_ok*sigmaintraT_SA_BA08_ok;
    rho_PGVSA_CY08(:,i) = rho_etaPGVSA_CY08(:,i)*sqrt((1-sigmaintraT_PGV_CY08_ok^2)*(1-sigmaintraT_SA_CY08_ok^2))+rho_epsPGVSA_CY08(:,i)*sigmaintraT_PGV_CY08_ok*sigmaintraT_SA_CY08_ok;
    rho_PGVSA_CB08(:,i) = rho_etaPGVSA_CB08(:,i)*sqrt((1-sigmaintraT_PGV_CB08_ok^2)*(1-sigmaintraT_SA_CB08_ok^2))+rho_epsPGVSA_CB08(:,i)*sigmaintraT_PGV_CB08_ok*sigmaintraT_SA_CB08_ok;
    rho_PGVSA_AS08(:,i) = rho_etaPGVSA_AS08(:,i)*sqrt((1-sigmaintraT_PGV_AS08_ok^2)*(1-sigmaintraT_SA_AS08_ok^2))+rho_epsPGVSA_AS08(:,i)*sigmaintraT_PGV_AS08_ok*sigmaintraT_SA_AS08_ok;
    
    %comparison of correlations (inter, intra, z, total)
    Tnum=10;
    if i==Tnum
        fig1=figure(1);
        axes('Parent',gcf,'FontSize',16);
        subplot(221);
        boxplot([rho_etaPGVSA_BA08(:,Tnum) rho_etaPGVSA_CY08(:,Tnum) rho_etaPGVSA_CB08(:,Tnum) rho_etaPGVSA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
        h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
        ylabel('Corr coeff, \rho_{eta,lnPGV,lnSA}'); grid on;
        subplot(222);
        boxplot([rho_epsPGVSA_BA08(:,Tnum) rho_epsPGVSA_CY08(:,Tnum) rho_epsPGVSA_CB08(:,Tnum) rho_epsPGVSA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
        h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
        ylabel('Corr coeff, \rho_{eps,lnPGV,lnSA}'); grid on;
        subplot(223);
        boxplot([rho_zPGVSA_BA08(:,Tnum) rho_zPGVSA_CY08(:,Tnum) rho_zPGVSA_CB08(:,Tnum) rho_zPGVSA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
        h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
        ylabel('Corr coeff, \rho_{z,lnPGV,lnSA}'); grid on;
        subplot(224);
        boxplot([rho_PGVSA_BA08(:,Tnum) rho_PGVSA_CY08(:,Tnum) rho_PGVSA_CB08(:,Tnum) rho_PGVSA_AS08(:,Tnum)],'labels',{'BA08';'CY08';'CB08';'AS08'});
        h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
        ylabel('Corr coeff, \rho_{lnPGV,lnSA}'); grid on;
        set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.9]);
        a=1;
    end
    
    clear z_lnSA_BA08_ok z_lnSA_CY08_ok z_lnSA_CB08_ok z_lnSA_AS08_ok;
    clear z_lnPGV_BA08_ok z_lnPGV_CY08_ok z_lnPGV_CB08_ok z_lnPGV_AS08_ok;
    clear sigma_SA_BA08_ok sigma_SA_CY08_ok sigma_SA_CB08_ok sigma_SA_AS08_ok;
    clear sigma_PGV_BA08_ok sigma_PGV_CY08_ok sigma_PGV_CB08_ok sigma_PGV_AS08_ok;
    clear eps_lnSA_BA08_ok eps_lnSA_CY08_ok eps_lnSA_CB08_ok eps_lnSA_AS08_ok;
    clear eps_lnPGV_BA08_ok eps_lnPGV_CY08_ok eps_lnPGV_CB08_ok eps_lnPGV_AS08_ok;
    clear eta_lnSA_BA08_ok eta_lnSA_CY08_ok eta_lnSA_CB08_ok eta_lnSA_AS08_ok;
    clear eta_lnPGV_BA08_ok eta_lnPGV_CY08_ok eta_lnPGV_CB08_ok eta_lnPGV_AS08_ok;
    clear z_lnSA_BA08_ok z_lnSA_CY08_ok z_lnSA_CB08_ok z_lnSA_AS08_ok;
    clear z_lnPGV_BA08_ok z_lnPGV_CY08_ok z_lnPGV_CB08_ok z_lnPGV_AS08_ok;
    
end


%Transformation of rho to fisher z values
z_rho_PGVASI(:,1)=FisherZTransformation(rho_PGVASI_BA08);
z_rho_PGVASI(:,2)=FisherZTransformation(rho_PGVASI_CY08);
z_rho_PGVASI(:,3)=FisherZTransformation(rho_PGVASI_CB08);
z_rho_PGVASI(:,4)=FisherZTransformation(rho_PGVASI_AS08);
% 
z_rho_PGVPGA(:,1)=FisherZTransformation(rho_PGVPGA_BA08);
z_rho_PGVPGA(:,2)=FisherZTransformation(rho_PGVPGA_CY08);
z_rho_PGVPGA(:,3)=FisherZTransformation(rho_PGVPGA_CB08);
z_rho_PGVPGA(:,4)=FisherZTransformation(rho_PGVPGA_AS08);
% 
z_rho_PGVSI(:,1)=FisherZTransformation(rho_PGVSI_BA08);
z_rho_PGVSI(:,2)=FisherZTransformation(rho_PGVSI_CY08);
z_rho_PGVSI(:,3)=FisherZTransformation(rho_PGVSI_CB08);
z_rho_PGVSI(:,4)=FisherZTransformation(rho_PGVSI_AS08);
% 
z_rho_PGVASI_all = [z_rho_PGVASI(:,1)' z_rho_PGVASI(:,2)' z_rho_PGVASI(:,3)' z_rho_PGVASI(:,4)']';
z_rho_PGVPGA_all = [z_rho_PGVPGA(:,1)' z_rho_PGVPGA(:,2)' z_rho_PGVPGA(:,3)' z_rho_PGVPGA(:,4)']';
z_rho_PGVSI_all = [z_rho_PGVSI(:,1)' z_rho_PGVSI(:,2)' z_rho_PGVSI(:,3)' z_rho_PGVSI(:,4)']';

for i=1:length(T)
    z_rho_PGVSA(:,i,1)=FisherZTransformation(rho_PGVSA_BA08(:,i));
    z_rho_PGVSA(:,i,2)=FisherZTransformation(rho_PGVSA_CY08(:,i));
    z_rho_PGVSA(:,i,3)=FisherZTransformation(rho_PGVSA_CB08(:,i));
    z_rho_PGVSA(:,i,4)=FisherZTransformation(rho_PGVSA_AS08(:,i));
    
    z_rho_PGVSA_all(:,i) = [z_rho_PGVSA(:,i,1)' z_rho_PGVSA(:,i,2)' z_rho_PGVSA(:,i,3)' z_rho_PGVSA(:,i,4)']';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputting of results
median_rho_PGVASI=InvFisherZTransformation(mean(z_rho_PGVASI,1));
median_rho_PGVPGA=InvFisherZTransformation(mean(z_rho_PGVPGA,1));
median_rho_PGVSI=InvFisherZTransformation(mean(z_rho_PGVSI,1));

median_rho_PGVASI_all=InvFisherZTransformation(mean(z_rho_PGVASI_all));
median_rho_PGVPGA_all=InvFisherZTransformation(mean(z_rho_PGVPGA_all));
median_rho_PGVSI_all=InvFisherZTransformation(mean(z_rho_PGVSI_all));

std_z_PGVASI=std(z_rho_PGVASI_all);
std_z_PGVPGA=std(z_rho_PGVPGA_all);
std_z_PGVSI=std(z_rho_PGVSI_all);

for i=1:length(T)
    median_rho_PGVSA(i,:)=InvFisherZTransformation(mean(z_rho_PGVSA(:,i,:),1));
    median_rho_PGVSA_all(i)=InvFisherZTransformation(mean(z_rho_PGVSA_all(:,i),1));
    std_z_PGVSA(i)=std(z_rho_PGVSA_all(:,i));
end



fprintf('Results of correlation analyses \n');
fprintf('\n');
fprintf('-----------------------  Median correlations   ---------------------- | ----------Overall   -------\n');
fprintf('             BA08   |     CY08      |     CY08      |     AS08      |   median     |  Z std   \n');
fprintf('PGV,ASI    %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_PGVASI,median_rho_PGVASI_all,std_z_PGVASI);
fprintf('PGV,PGA    %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_PGVPGA,median_rho_PGVPGA_all,std_z_PGVPGA);
fprintf('PGV,SI     %6.3f       %6.3f       %6.3f      %6.3f      |     %6.3f       %6.3f      \n',median_rho_PGVSI,median_rho_PGVSI_all,std_z_PGVSI);
fprintf('\n');
fprintf('              median PGV - SA correlation       |                          |  \n');
fprintf('T      |    BA08  |    CY08 |    CB08 |    AS08 |   median   |  Z variance |  \n');
for i=1:length(T)
fprintf('%6.3f     %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   \n',...
                                                        T(i),median_rho_PGVSA(i,:),median_rho_PGVSA_all(i),std_z_PGVSA(i));
end
fprintf('------------------------------------- \n');

%plotting
%PGV with ASI,SI, PGA correlations
fig1=figure(1);
axes('Parent',gcf,'FontSize',16);
subplot(131);
boxplot([rho_PGVASI_BA08 rho_PGVASI_CY08 rho_PGVASI_CB08 rho_PGVASI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnPGV,lnASI}'); grid on;
subplot(132);
boxplot([rho_PGVPGA_BA08 rho_PGVPGA_CY08 rho_PGVPGA_CB08 rho_PGVPGA_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnPGV,lnPGA}'); grid on;
subplot(133);
boxplot([rho_PGVSI_BA08 rho_PGVSI_CY08 rho_PGVSI_CB08 rho_PGVSI_AS08],'labels',{'BA08';'CY08';'CB08';'AS08'});
h = findobj(gca,'Type','Line'); for i=1:length(h); set(h(i),'LineWidth',2); end;
ylabel('Corr coeff, \rho_{lnPGV,lnSI}'); grid on;
set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.9 0.35]);
% 
% %PGV correlation with SA
fig2=figure(2);
axes('Parent',gcf,'FontSize',16);
semilogx(T,mean(rho_PGVSA_BA08,1),'LineWidth',3,'LineStyle','-','Color',[1 0 0]); hold on;
semilogx(T,prctile(rho_PGVSA_BA08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[1 0 0]);
semilogx(T,mean(rho_PGVSA_CY08,1),'LineWidth',3,'LineStyle','-','Color',[0 1 0]);
semilogx(T,prctile(rho_PGVSA_CY08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0 1 0]);
semilogx(T,mean(rho_PGVSA_CB08,1),'LineWidth',3,'LineStyle','-','Color',[0 0 1]);
semilogx(T,prctile(rho_PGVSA_CB08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0 0 1]);
semilogx(T,mean(rho_PGVSA_AS08,1),'LineWidth',3,'LineStyle','-','Color',[0.5 0.5 0.5]);
semilogx(T,prctile(rho_PGVSA_AS08,[5 95],1)','LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5]);
ylabel('Corr coeff, \rho_{lnPGVSI,lnSA}'); xlabel('Period, T (s)'); grid on;
set(gcf,'units','normalized'); set(gcf,'Position',[0.05 0.05 0.5 0.4]);

%show adequacy of z having a normal distribution
z_rho_ij=z_rho_PGVASI_all; %which is considered

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
