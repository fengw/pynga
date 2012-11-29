function Ngm_Neq_functionperiod
clc
%Brendon Bradley 10 Aug 2010

%Purpose: Determine the number of earthquakes and ground motions used as a function of period 
%in the empirical correlation equations
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
    
%     ------------SA----------------
    for j=1:length(T)
        siteprop.period=T(j);
%         get the predicted value of PGA using BA08
        [SA_BA08(i,j),sigma_SA_BA08(i,1:3,j)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
%         get PGA for the GM record
        SA_observed(i,j)=data(i,19+j);
%         if the usable frequency is ok compute residual else set residual =
%         -100;
        freq_min=data(i,17);
        T_max=1/freq_min;
        if T(j)<T_max
%             compute zlnSA
            z_lnSA_BA08(i,j)=(log(SA_observed(i,j))-log(SA_BA08(i,j)))/sigma_SA_BA08(i,1,j);
        else
            z_lnSA_BA08(i,j)=-100;
        end
    end
end

%4) Correlation between DSI with SA(T)
for i=1:length(T)
    %need to check if each residual is not -100 indicating above fmin
    k=0;
    for j=1:length(z_lnSA_BA08)  %will be same for both BA08 and CY08 so dont need to check twice
        if (z_lnSA_BA08(j,i)~=-100)
            k=k+1;
            
            EQID_ok(k)=EQID(j);
            z_lnSA_BA08_ok(k)=z_lnSA_BA08(j,i);
            sigma_SA_BA08_ok(k,1:3)=sigma_SA_BA08(j,1:3,i);
            
        end
    end
    alleta=0;
    [eta_lnSA_BA08_ok,eps_lnSA_BA08_ok]=get_interintraeventterms(z_lnSA_BA08_ok,EQID_ok,sigma_SA_BA08_ok(:,1)',sigma_SA_BA08_ok(:,2)',sigma_SA_BA08_ok(:,3)',alleta);
    Ngm(i)=length(eps_lnSA_BA08_ok);
    Neq(i)=length(eta_lnSA_BA08_ok);
    clear eta_lnSA_BA08_ok eps_lnSA_BA08_ok;
    clear EQID_ok z_lnSA_BA08_ok sigma_SA_BA08_ok
end

fig1=figure(1);
axes1=axes('Parent',fig1,'FontSize',14);
loglog(T,Ngm,'LineWidth',3,'Color',[0 0 0],'LineStyle','-'); hold on;
loglog(T,Neq,'LineWidth',3,'Color',[0.5 0.5 0.5],'LineStyle','--');
xlabel('Period, T (s)'); legend('N_{gm}','N_{eq}','Location','SouthWest'); 
xlim([0.01 10]); ylim([10 2000]);
grid on
