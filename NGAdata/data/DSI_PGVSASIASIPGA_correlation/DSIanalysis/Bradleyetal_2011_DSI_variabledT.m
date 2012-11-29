function [DSI,sigma_DSI]=Bradleyetal_2011_DSI_variabledT(M,R,siteprop,faultprop,IMR)
%Brendon Bradley   1 Sept 2010

%Provides the attenuation relation for disp spectrum intensity based
%upon spectral acceleration relationships.

%Def:  DSI = integral[Sd(T,5%)] from T=0.5-5.0s

%Basis of approach is that DSI attenuation can be obtained from attenuation
%relation from SA using any general SA model.

%Provides the attenuation relation for disp spectrum intensity defined as:

%               /-T=5.0s
%      DSI =    |       Sd(T,5%)dT     
%              -/T=2.0s

%note that care should be taken in the value of 'g' used (use g=9.81m/s2 to
%give DSI in m.s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input Variables:
% M             = Moment magnitude (Mw)
% R             = Source-to-site distance (km) 
% siteprop      = properties of site (soil etc) 
%                 siteprop.g    -'acceleration of gravity (typ 9.81m/s2)
%                 (others required for the second IMR selected)
% faultprop     = properties of fault (strikeslip etc) (not required in the VSI computation here
%                 ,but input as for the second IMR)
% IMR           = handle of the Intensity Measure Relationship (attenuation
%                 relation) used for getting the spectral acceleration attenuation

%Output Variables:
% DSI           = median DSI 
% sigma_ASI     = lognormal standard deviation in DSI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first determine elastic spectral accelerations at various period ranges
g=siteprop.g;
dT_type=siteprop.dTtype;
if dT_type==0 %linear
    Nt=siteprop.NTi;
    T=2.0:(5-2.0)/(Nt-1):5;
elseif dT_type==1 %logarithmic
    Nt=siteprop.NTi;
    T=exp(log(2.0):(log(5)-log(2.0))/(Nt-1):log(5));
end

for i=1:length(T)
    siteprop.period=T(i);
    [SA(i),sigma_SA(i,1:3)]=feval(IMR,M,R,siteprop,faultprop);
end

%now convert to spectral acceleration mean and stadnard deviation (normal)
for i=1:length(T)
    SA_mean(i)=SA(i)*exp(0.5*sigma_SA(i,1)^2);
    Sd(i)=(T(i)/(2*pi)).^2.*SA_mean(i)*g;
    %convert lognormal standard deviation in SA to normal 
    sigma_SA_normal(i)=SA_mean(i)*sqrt(exp(sigma_SA(i)^2)-1);
    sigma_Sd_normal(i)=(T(i)/(2*pi)).^2.*sigma_SA_normal(i)*g;
end

%allocate integration weights for trapz rule 
weight(1)=(T(2)-T(1))/2;
for i=2:length(T)-1
    weight(i)=(T(i+1)-T(i-1))/2;
end
weight(Nt)=(T(Nt)-T(Nt-1))/2;

%integration discretely as a summation
%mean
DSI_mean=0;
for i=1:length(T)
    DSI_mean=DSI_mean+weight(i)*Sd(i);    
end

%stdev
var_DSI_normal=0;
for i=1:length(T)
    for j=1:i
        if i==j
            var_DSI_normal=var_DSI_normal+weight(i)^2*sigma_Sd_normal(i)^2;
        else
            %compute correlation
            [rholn]=SA_correlation(T(i),T(j));  %log correlation
            rhon=(exp(rholn*sigma_SA(i)*sigma_SA(j))-1)/sqrt((exp(sigma_SA(i)^2)-1)*(exp(sigma_SA(j)^2)-1));    %normal correlation
            var_DSI_normal=var_DSI_normal+2*rhon*weight(i)*weight(j)*sigma_Sd_normal(i)*sigma_Sd_normal(j);
        end
    end
end
sigma_DSI_normal=sqrt(var_DSI_normal);
%convert to median and lognormal standard deviation
sigma_DSI=sqrt(log((sigma_DSI_normal/DSI_mean)^2+1));
DSI=DSI_mean*exp(-0.5*sigma_DSI^2);
            
%end of attenuation relationship
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho]=SA_correlation(T1,T2);

% Created by Jack Baker, 2/28/07 (updated 6/25/2007)
% Compute the correlation of epsilons for the NGA ground motion models
%
% The function is strictly emperical, fitted over the range the range 0.01s <= T1, T2 <= 10s
%
% Documentation is provided in the following document:
% Baker, J.W. and Jayaram, N., "Correlation of spectral acceleration values from NGA ground 
% motion models," Earthquake Spectra, (in review).

% INPUT
%
%   T1, T2      = The two periods of interest. The periods may be equal,
%                 and there is no restriction on which one is larger.
%
% INPUT
%
%   rho         = The predicted correlation coefficient

T_min = min(T1, T2);
T_max = max(T1, T2);

C1 = (1-cos(pi/2 - log(T_max/max(T_min, 0.109)) * 0.366 ));
if T_max < 0.2
    C2 = 1 - 0.105*(1 - 1./(1+exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099);
end
if T_max < 0.109
    C3 = C2;
else
    C3 = C1;
end
C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi*(T_min)/(0.109)));

if T_max <= 0.109
    rho = C2;
elseif T_min > 0.109
    rho = C1;
elseif T_max < 0.2
    rho = min(C2, C4);
else
    rho = C4;
end
%end of script
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
