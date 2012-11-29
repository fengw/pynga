function [SA_CMSmedian,SA_CMSbeta]=CMSspectra_approximate(M,R,epsilon,targetSa,targetperiod,periodrange,IMR,siteprop,faultprop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Purpose: To compute the approximate CMS spectra by using the mean
%M,R,epsilon values from PSHA deaggregation

%required M-files in same folder
%   SA_correlation.m
%Reference: Baker, JW and Cornell, CA, 2006.  Spectral acceleration, record
%selection and epsilon.  Earthquake Engineering and Structural Dynamics.
%35: 1077-1095.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input data
M=6.5;
R=50;
eps=1.5;
targetSa=0.68;   %in g
targetperiod=0.5; %seconds

siteprop.soiltype='rock';
faultprop.faultstyle='normal';
siteprop.g=981;
siteprop.V30=300;

periodrange=0:0.05:4;
IMR=@BooreAtkinson_2007_nga;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get predicted SA
siteprop.period=targetperiod;
[SApredicted,sigma_SApredicted]=feval(IMR,M,R,siteprop,faultprop);

%adjust epsilon value to get target Sa
epsilonmod=log(targetSa/SApredicted)/sigma_SApredicted;

%get unconditional spectra and correlation at range of periods
for i=1:length(periodrange)
    siteprop.period=periodrange(i);
    [SA(i),sigma_SA(i)]=feval(IMR,M,R,siteprop,faultprop);
    [rho(i)]=SA_correlation(periodrange(i),targetperiod);
    %get CMS
    SA_CMSmedian(i)=exp(log(SA(i))+sigma_SA(i)*rho(i)*epsilonmod);
    SA_CMSbeta(i)=sigma_SA(i)*sqrt(1-rho(i)^2);
end

%84th and 16th percentiles CMS
SA_CMS16=exp(log(SA_CMSmedian)-SA_CMSbeta);
SA_CMS84=exp(log(SA_CMSmedian)+SA_CMSbeta);

%plotting
figure(1)
plot(periodrange,SA,'-r');
hold on
plot(periodrange,SA_CMSmedian,'b-');
plot(periodrange,SA_CMS16,'b--');
plot(periodrange,SA_CMS84,'b--');
plot(periodrange,exp(log(SA)+epsilonmod*sigma_SA),'--r');





