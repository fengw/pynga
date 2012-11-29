function [DSI,sigma_DSI]=Bradleyetal_2011_DSI(M,R,siteprop,faultprop,IMR)
%Brendon Bradley   31 Aug 2011

%Provides the attenuation relation for displacement spectrum intensity based
%upon spectral acceleration relationships.

%Def: DSI = integral[Sd(T,5%)] from T=0.5-5s

%Basis of approach is that DSI attenuation can be obtained from attenuation
%relation from SA using any general SA model.

%Provides the attenuation relation for spectrum intensity defined as:

%               /-T=5.0s
%      DSI =    |       Sd(T,5%)dT     
%              -/T=2.0s

%note that care should be taken in the value of 'g' used.  If the spectral
%acceleration attenuation relation used gives Sa in % of g then using
%g=9.81m/s2 will give DSI in units of m.s/s (i.e. m) (common)
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
% DSI          = median DSI 
% sigma_DSI     = lognormal standard deviation in DSI
                  %sigma_DSI(1) = total sigma
                  %sigma_DSI(2) = interevent sigma
                  %sigma_DSI(3) = intraevent sigma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Coefficients
g=siteprop.g;   %acceleration of gravity

%first determine elastic spectral accelerations at various period ranges
Nt=21;
Tmin=2; Tmax=5;
T=exp(log(Tmin):(log(Tmax)-log(Tmin))/(Nt-1):log(Tmax));

for i=1:length(T)
    siteprop.period=T(i);
    [SA(i),sigma_SA(i,1:3)]=feval(IMR,M,R,siteprop,faultprop);
end

%now convert to spectral velocity
for i=1:length(T)
    SA_mean=SA(i)*exp(0.5*sigma_SA(i,1)^2);
    SA_mean_intra=SA(i)*exp(0.5*sigma_SA(i,3)^2); 
    Sd(i)=(T(i)/(2*pi)).^2.*SA_mean*g;
    Sd_intra(i)=(T(i)/(2*pi)).^2.*SA_mean_intra*g;
    %convert lognormal standard deviation in SA to normal 
    sigma_SA_normal=SA_mean*sqrt(exp(sigma_SA(i,1)^2)-1);
    sigma_SAintra_normal=SA_mean_intra*sqrt(exp(sigma_SA(i,3)^2)-1);
    sigma_Sd_normal(i)=(T(i)/(2*pi)).^2.*sigma_SA_normal*g;
    sigma_Sdintra_normal(i)=(T(i)/(2*pi)).^2.*sigma_SAintra_normal*g;
end

%allocate integration weights for trapz rule (dT/2 for first and last, dT otherwise)
weight(1)=(T(2)-T(1))/2;
for i=2:length(T)-1
    weight(i)=(T(i+1)-T(i-1))/2;
end
weight(Nt)=(T(Nt)-T(Nt-1))/2;

%integration discretely as a summation
%mean
DSI_mean=0; 
DSI_mean_intra=0;
for i=1:length(T)
    DSI_mean=DSI_mean+weight(i)*Sd(i);    
    DSI_mean_intra=DSI_mean_intra+weight(i)*Sd_intra(i);   
end

%stdev
var_DSI_normal=0;
var_DSIintra_normal=0;
for i=1:length(T)
    for j=1:i
        if i==j
            var_DSI_normal=var_DSI_normal+weight(i)^2*sigma_Sd_normal(i)^2;
            var_DSIintra_normal=var_DSIintra_normal+weight(i)^2*sigma_Sdintra_normal(i)^2;
        else
            %compute correlation
            [rholn]=SA_correlation(T(i),T(j));  %log correlation
            rhon=(exp(rholn*sigma_SA(i,1)*sigma_SA(j,1))-1)/sqrt((exp(sigma_SA(i,1)^2)-1)*(exp(sigma_SA(j,1)^2)-1));    %normal correlation
            rhonintra=(exp(rholn*sigma_SA(i,3)*sigma_SA(j,3))-1)/sqrt((exp(sigma_SA(i,3)^2)-1)*(exp(sigma_SA(j,3)^2)-1));    %normal correlation
            var_DSI_normal=var_DSI_normal+2*rhon*weight(i)*weight(j)*sigma_Sd_normal(i)*sigma_Sd_normal(j);
            var_DSIintra_normal=var_DSIintra_normal+2*rhonintra*weight(i)*weight(j)*sigma_Sdintra_normal(i)*sigma_Sdintra_normal(j);
        end
    end
end
sigma_DSI_normal=sqrt(var_DSI_normal);
sigma_DSIintra_normal=sqrt(var_DSIintra_normal);
%convert to median and lognormal standard deviation
sigma_DSI(1)=sqrt(log((sigma_DSI_normal/DSI_mean)^2+1));   %total sigma
sigma_DSI(3)=sqrt(log((sigma_DSIintra_normal/DSI_mean_intra)^2+1));   %intra-event sigma
sigma_DSI(2)=sqrt(sigma_DSI(1)^2-sigma_DSI(3)^2); %inter event sigma
DSI=DSI_mean*exp(-0.5*sigma_DSI(1)^2);
            
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
