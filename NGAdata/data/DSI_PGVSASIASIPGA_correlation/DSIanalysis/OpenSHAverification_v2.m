function [ASI,sigma_ASI]=Bradleyetal_2008_ASI_variabledT(M,R,siteprop,faultprop,IMR)
%Brendon Bradley   6 April 2008
%           modified 28 Aug 2009 - allows both linear and logarithmic
%           spacing of Ti
clc
%Provides the attenuation relation for accleration spectrum intensity based
%upon spectral acceleration relationships.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first determine elastic spectral accelerations at various period ranges
% dT_type=siteprop.dTtype;
% if dT_type==0 %linear
%     Nt=siteprop.NTi;
%     T=0.1:(0.5-0.1)/(Nt-1):0.5;
% elseif dT_type==1 %logarithmic
%     Nt=siteprop.NTi;
%     T=exp(log(0.1):(log(0.5)-log(0.1))/(Nt-1):log(0.5));
% end
% 
% for i=1:length(T)
%     siteprop.period=T(i);
%     [SA(i),sigma_SA(i)]=feval(IMR,M,R,siteprop,faultprop);
% end

OpenSHAdata=[2.0 0.03786425629767887 0.7
2.2427067839402772 0.03354166892106412 0.698587593122147
2.5148668593658705 0.029712548562039105 0.6971751862442941
2.820054483103209 0.026320560975339913 0.6957627793664412
3.16227766016838 0.023381378520994686 0.6955493591315574
3.5460307805812055 0.020839019063761986 0.6967437645441866
3.976353643835253 0.01857310145977488 0.6979381699568159
4.458897646187482 0.017710638577324307 0.720388868612844
5.000000000000001 0.01695063199460118 0.744];
T=OpenSHAdata(:,1); SA=OpenSHAdata(:,2); sigma_SA=OpenSHAdata(:,3);

%now convert to spectral acceleration mean and stadnard deviation (normal)
for i=1:length(T)
    SA_mean(i)=SA(i)*exp(0.5*sigma_SA(i)^2);
    %convert lognormal standard deviation in SA to normal 
    sigma_SA_normal(i)=SA_mean(i)*sqrt(exp(sigma_SA(i)^2)-1);
%     fprintf('%10.5f %10.5f %10.5f \n',T(i),SA_mean(i),sigma_SA_normal(i));
end

%allocate integration weights for trapz rule 
weight(1)=(T(2)-T(1))/2;
for i=2:length(T)-1
    weight(i)=(T(i+1)-T(i-1))/2;
end
weight(length(T))=(T(length(T))-T(length(T)-1))/2;

%integration discretely as a summation
%mean
ASI_mean=0;
for i=1:length(T)
    ASI_mean=ASI_mean+weight(i)*SA_mean(i);
    fprintf('%10.5f %10.5f %10.5f \n',T(i),weight(i),ASI_mean);
end

%stdev
var_ASI_normal=0;
for i=1:length(T)
    for j=1:i
        if i==j
            var_ASI_normal=var_ASI_normal+weight(i)^2*sigma_SA_normal(i)^2;
        else
            %compute correlation
            [rholn]=SA_correlation(T(i),T(j));  %log correlation
            rhon=(exp(rholn*sigma_SA(i)*sigma_SA(j))-1)/sqrt((exp(sigma_SA(i)^2)-1)*(exp(sigma_SA(j)^2)-1));    %normal correlation
            var_ASI_normal=var_ASI_normal+2*rhon*weight(i)*weight(j)*sigma_SA_normal(i)*sigma_SA_normal(j);
        end
    end
end
sigma_ASI_normal=sqrt(var_ASI_normal);
%convert to median and lognormal standard deviation
sigma_ASI=sqrt(log((sigma_ASI_normal/ASI_mean)^2+1));
ASI=ASI_mean*exp(-0.5*sigma_ASI^2);
fprintf('ASI median = %10.5f \n',ASI);
fprintf('sigma ASI = %10.5f \n',sigma_ASI);
            
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
