function [SI,sigma_SI]=Bradleyetal_2008_SI(M,R,siteprop,faultprop,IMR)
%Brendon Bradley   6 April 2008

%Provides the attenuation relation for velocity spectrum intensity based
%upon spectral acceleration relationships.

%Def: SI = VSI = integral[PSV(T,5%)] from T=0.1-2.5s

%Basis of approach is that SI attenuation can be obtained from attenuation
%relation from SA using any general SA model.

%Provides the attenuation relation for spectrum intensity defined as:

%               /-T=2.5s
%      SI =     |       PSV(T,5%)dT     
%              -/T=0.1s

%note that care should be taken in the value of 'g' used.  If the spectral
%acceleration attenuation relation used gives Sa in % of g then using
%g=9.81m/s2 will give SI in units of m.s/s (i.e. m).  g=981.0 will give SI
%in units cm.s/s (common)
%if Sa is given in units of cm/s2 then using 'g=1' will give SI in cm.s/s

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
% SI           = median SI 
% sigma_SI     = lognormal standard deviation in SI
                  %sigma_SI(1) = total sigma
                  %sigma_SI(2) = interevent sigma
                  %sigma_SI(3) = intraevent sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Coefficients
g=siteprop.g;   %acceleration of gravity

%first determine elastic spectral accelerations at various period ranges
dT=0.1;
T=0.1:dT:2.5;

for i=1:length(T)
    siteprop.period=T(i);
    [SA(i),sigma_SA(i,1:3)]=feval(IMR,M,R,siteprop,faultprop);
end

%now convert to spectral velocity
for i=1:length(T)
    SA_mean=SA(i)*exp(0.5*sigma_SA(i,1)^2);
    SA_mean_intra=SA(i)*exp(0.5*sigma_SA(i,3)^2); 
    SV(i)=(T(i)/(2*pi))*SA_mean*g;
    SV_intra(i)=(T(i)/(2*pi))*SA_mean_intra*g;
    %convert lognormal standard deviation in SA to normal 
    sigma_SA_normal=SA_mean*sqrt(exp(sigma_SA(i,1)^2)-1);
    sigma_SAintra_normal=SA_mean_intra*sqrt(exp(sigma_SA(i,3)^2)-1);
    sigma_SV_normal(i)=(T(i)/(2*pi))*sigma_SA_normal*g;
    sigma_SVintra_normal(i)=(T(i)/(2*pi))*sigma_SAintra_normal*g;
end

%allocate integration weights for trapz rule (dT/2 for first and last, dT otherwise)
weight=dT*ones(1,length(T)); 
weight(1)=weight(1)/2; weight(length(T))=weight(length(T))/2;

%integration discretely as a summation
%mean
SI_mean=0; 
SI_mean_intra=0;
for i=1:length(T)
    SI_mean=SI_mean+weight(i)*SV(i);    
    SI_mean_intra=SI_mean_intra+weight(i)*SV_intra(i);   
end

%stdev
var_SI_normal=0;
var_SIintra_normal=0;
for i=1:length(T)
    for j=1:i
        if i==j
            var_SI_normal=var_SI_normal+weight(i)^2*sigma_SV_normal(i)^2;
            var_SIintra_normal=var_SIintra_normal+weight(i)^2*sigma_SVintra_normal(i)^2;
        else
            %compute correlation
            [rholn]=SA_correlation(T(i),T(j));  %log correlation
            rhon=(exp(rholn*sigma_SA(i,1)*sigma_SA(j,1))-1)/sqrt((exp(sigma_SA(i,1)^2)-1)*(exp(sigma_SA(j,1)^2)-1));    %normal correlation
            rhonintra=(exp(rholn*sigma_SA(i,3)*sigma_SA(j,3))-1)/sqrt((exp(sigma_SA(i,3)^2)-1)*(exp(sigma_SA(j,3)^2)-1));    %normal correlation
            var_SI_normal=var_SI_normal+2*rhon*weight(i)*weight(j)*sigma_SV_normal(i)*sigma_SV_normal(j);
            var_SIintra_normal=var_SIintra_normal+2*rhonintra*weight(i)*weight(j)*sigma_SVintra_normal(i)*sigma_SVintra_normal(j);
        end
    end
end
sigma_SI_normal=sqrt(var_SI_normal);
sigma_SIintra_normal=sqrt(var_SIintra_normal);
%convert to median and lognormal standard deviation
sigma_SI(1)=sqrt(log((sigma_SI_normal/SI_mean)^2+1));   %total sigma
sigma_SI(3)=sqrt(log((sigma_SIintra_normal/SI_mean_intra)^2+1));   %intra-event sigma
sigma_SI(2)=sqrt(sigma_SI(1)^2-sigma_SI(3)^2); %inter event sigma
SI=SI_mean*exp(-0.5*sigma_SI(1)^2);
            
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
