function input_ASI_dispersionmagnitude

%for comparison of the ground motion dispersion in ASI as compared to that
%of spectral ordinates
clc
format long
%ERF = earthquake rupture forecast 
ERF=@ERF_onefault_spatialIMinvestigation; 
tint=1; %time interval over which to compute the solution

%IMR = intensity measure relationship
% IMR{1}=@BooreAtkinson_2007_nga;    DanciuTselentis_2007_Sa
IMR{2}=@Bradleyetal_2008_ASI;
IM=0.0:0.01:1.5; %%range of IM values to compute solution for

%site properties - soil type etc
siteprop.soiltype='rock';

faultprop.faultstyle='normal';
siteprop.g=981;
siteprop.V30=300;
% 
att_type=1;     %=1 BA; =2 DT; 3=BA but with R varying ; 4=EB with M varying

M=6.5; R=10;
T=0.1:0.1:2.5;
    
if att_type==1
    for i=1:length(T)
        siteprop.period=T(i);
        [SA(i),sigma_SA(i)]=BooreAtkinson_2007_nga(M,R,siteprop,faultprop);
    end
    [ASI,sigma_ASI]=Bradleyetal_2008_ASI(M,R,siteprop,faultprop,@BooreAtkinson_2007_nga);
    [SI,sigma_SI]=Bradleyetal_2008_SI(M,R,siteprop,faultprop,@BooreAtkinson_2007_nga);
    siteprop.period=-1;
    [PGV,sigma_PGV]=BooreAtkinson_2007_nga(M,R,siteprop,faultprop);
elseif att_type==2
    for i=1:length(T)
        siteprop.period=T(i);
        [SA(i),sigma_SA(i)]=DanciuTselentis_2007_SA(M,R,siteprop,faultprop);
    end
    [SI,sigma_SI]=Bradleyetal_2008_SI(M,R,siteprop,faultprop,@DanciuTselentis_2007_SA);
    [SI_actual,sigma_SI_actual]=DanciuTselentis_2007_SI(M,R,siteprop,faultprop);
elseif att_type==3
    M=6.5; R=5:5:75;
    for i=1:length(R)
        siteprop.period=-1;  %PGV
        [PGV(i),sigma_PGV(i)]=BooreAtkinson_2007_nga(M,R(i),siteprop,faultprop);
        [SI(i),sigma_SI(i)]=Bradleyetal_2008_SI(M,R(i),siteprop,faultprop,@BooreAtkinson_2007_nga);
        [ASI(i),sigma_ASI(i)]=Bradleyetal_2008_ASI(M,R(i),siteprop,faultprop,@BooreAtkinson_2007_nga);
    end
    minmax_sigmaSI=[min(sigma_SI) max(sigma_SI)]
else 
    M=5:0.2:8; R=30;
    for i=1:length(M)
        siteprop.period=-1;  %PGV
        [PGV(i),sigma_PGV(i)]=BooreAtkinson_2007_nga(M(i),R,siteprop,faultprop);
        [SI(i),sigma_SI(i)]=Bradleyetal_2008_SI(M(i),R,siteprop,faultprop,@BooreAtkinson_2007_nga);
    end
    minmax_sigmaSI=[min(sigma_SI) max(sigma_SI)]
end

outplot=2;
if outplot==1
    %mean
%     title('median')
    loglog(T,SA,'-or',[min(T) max(T)],[SI SI],'-b'); legend('Danicu Tselentis SA','SI')
else
    %std
%     title('dispersion')
    if att_type==1
        semilogx(T,sigma_SA,'-or',[min(T) max(T)],[sigma_SI sigma_SI],'-b',[min(T) max(T)],[sigma_ASI sigma_ASI],'--r',[min(T) max(T)],[sigma_PGV sigma_PGV],'-g'); legend('SA','SI','ASI','PGV')
    elseif att_type==2
        semilogx(T,sigma_SA,'-or',[min(T) max(T)],[sigma_SI sigma_SI],'-b',[min(T) max(T)],[sigma_SI_actual sigma_SI_actual],'-g'); legend('SA','SI','SI_actual')
    elseif att_type==3
        plot(R,sigma_SI,'-or',R,sigma_PGV); legend('SI','PGV')
    else
        plot(M,sigma_SI,'-or',M,sigma_PGV); legend('SI','PGV')
    end

end









