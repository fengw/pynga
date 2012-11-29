function [ASI,sigma_ASI]=Bradleyetal_2008_ASI(M,R,siteprop,faultprop,IMR)
%Brendon Bradley   6 April 2008

%Provides the attenuation relation for accleration spectrum intensity based
%upon spectral acceleration relationships (gives both total, inter, and intra event residuals).

%Def:  ASI = integral[PSA(T,5%)] from T=0.1-0.5s

%Basis of approach is that ASI attenuation can be obtained from attenuation
%relation from SA using any general SA model.

%Provides the attenuation relation for acceleration spectrum intensity defined as:

%               /-T=0.5s
%      ASI =    |       PSA(T,5%)dT     
%              -/T=0.1s

%note that care should be taken in the value of 'g' used.  If the spectral
%acceleration attenuation relation used gives Sa in % of g then 
%ASI in units of g.s (if g=9.81m/s2 then ASI is m.s/s2)

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
% ASI           = median ASI 
% sigma_ASI     = lognormal standard deviation in ASI 
                  %sigma_ASI(1) = total sigma
                  %sigma_ASI(2) = interevent sigma
                  %sigma_ASI(3) = intraevent sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first determine elastic spectral accelerations at various period ranges
dT=0.1;
T=0.1:dT:0.5;

for i=1:length(T)
    siteprop.period=T(i);
    [SA(i),sigma_SA(i,1:3)]=feval(IMR,M,R,siteprop,faultprop);
end

%now convert to spectral acceleration mean and stadnard deviation (normal)
for i=1:length(T)
    SA_mean(i)=SA(i)*exp(0.5*sigma_SA(i,1)^2);
    SA_mean_intra(i)=SA(i)*exp(0.5*sigma_SA(i,3)^2); 
    %convert lognormal standard deviation in SA to normal 
    sigma_SA_normal(i)=SA_mean(i)*sqrt(exp(sigma_SA(i,1)^2)-1);
    sigma_SAintra_normal(i)=SA_mean_intra(i)*sqrt(exp(sigma_SA(i,3)^2)-1);
end

%allocate integration weights for trapz rule (dT/2 for first and last, dT otherwise)
weight=dT*ones(1,length(T)); 
weight(1)=weight(1)/2; weight(length(T))=weight(length(T))/2;

%integration discretely as a summation
%mean
ASI_mean=0;
ASI_mean_intra=0;
for i=1:length(T)
    ASI_mean=ASI_mean+weight(i)*SA_mean(i);    
    ASI_mean_intra=ASI_mean_intra+weight(i)*SA_mean_intra(i); 
end

%stdev
var_ASI_normal=0;
var_ASIintra_normal=0;
for i=1:length(T)
    for j=1:i
        if i==j
            var_ASI_normal=var_ASI_normal+weight(i)^2*sigma_SA_normal(i)^2;
            var_ASIintra_normal=var_ASIintra_normal+weight(i)^2*sigma_SAintra_normal(i)^2;
        else
            %compute correlation
            [rholn]=SA_correlation(T(i),T(j));  %log correlation
            rhon=(exp(rholn*sigma_SA(i,1)*sigma_SA(j,1))-1)/sqrt((exp(sigma_SA(i,1)^2)-1)*(exp(sigma_SA(j,1)^2)-1));    %normal correlation
            rhonintra=(exp(rholn*sigma_SA(i,3)*sigma_SA(j,3))-1)/sqrt((exp(sigma_SA(i,3)^2)-1)*(exp(sigma_SA(j,3)^2)-1));    %normal correlation
            var_ASI_normal=var_ASI_normal+2*rhon*weight(i)*weight(j)*sigma_SA_normal(i)*sigma_SA_normal(j);
            var_ASIintra_normal=var_ASIintra_normal+2*rhonintra*weight(i)*weight(j)*sigma_SAintra_normal(i)*sigma_SAintra_normal(j);
        end
    end
end
sigma_ASI_normal=sqrt(var_ASI_normal);
sigma_ASIintra_normal=sqrt(var_ASIintra_normal);
%convert to median and lognormal standard deviation
sigma_ASI(1)=sqrt(log((sigma_ASI_normal/ASI_mean)^2+1));        %total sigma
sigma_ASI(3)=sqrt(log((sigma_ASIintra_normal/ASI_mean_intra)^2+1));   %intra-event sigma
sigma_ASI(2)=sqrt(sigma_ASI(1)^2-sigma_ASI(3)^2); %inter-event sigma
ASI(1)=ASI_mean*exp(-0.5*sigma_ASI(1)^2); 
            
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
