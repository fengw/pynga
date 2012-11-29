function input_VSI_directvsBradleymethod

%for comparison of the Danciu and Tselentis attenuation relations for
%computing VSI directly, and via the spectral attentuation relations
format long
%ERF = earthquake rupture forecast 
ERF=@ERF_onefault_spatialIMinvestigation; 
tint=1; %time interval over which to compute the solution

%IMR = intensity measure relationship
IMR{1}=@DanciuTselentis_2007_Sa;
IM=0.0:10:300; %%range of IM values to compute solution for

%site properties - soil type etc
siteprop.soiltype='rock'

faultprop.faultstyle='normal';
siteprop.g=1;
% 

    M=5.5; r=1:5:200;

    for i=1:length(r)
        R=r(i)
    %     [SI(i),sigma_SI]=Bradleyetal_2008_SI(M,R(i),siteprop,faultprop,IMR(2))
        [SI(i),sigma_SI(i)]=DanciuTselentis_2007_SI(M,R,siteprop,faultprop);
        SI(i)=SI(i)*2.4;
        [SIa(i),sigma_SIa(i)]=Bradleyetal_2008_SI(M,R,siteprop,faultprop,@DanciuTselentis_2007_Sa);
    end

    loglog(r,SI,'-r',r,SIa,'-b'); legend('Danicu Tselentis','via Sa attenuation')
    hold on
    loglog(r,SI.*exp(sigma_SI),'--r',r,SI.*exp(-sigma_SI),'--r')
    loglog(r,SIa.*exp(sigma_SIa),'--b',r,SIa.*exp(-sigma_SIa),'--b')
    
    figure(2)
    plot(r,sigma_SI,'--b',r,sigma_SIa,'-g'); legend('DT','this study')









