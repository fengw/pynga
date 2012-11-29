function [SA_CMSmedian,SA_CMSbeta]=CMSspectra_approximate_SI(M,R,epsilon,targetASI,periodrange,IMR,siteprop,faultprop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Purpose: To compute the approximate CMS spectra fixed at SI
%by using the mean M,R,epsilon values from PSHA deaggregation

%required M-files in same folder
%   SA_correlation.m

%Reference: 
%Baker, JW and Cornell, CA, 2006.  Spectral acceleration, record
%selection and epsilon.  Earthquake Engineering and Structural Dynamics.
%35: 1077-1095.
%Baker, JW.  2005.  Vector-valued Intensity Measures for probabilistic
%seismic demand analysis, PhD Thesis.  Stanford University, CA.
%Bradley, BA. 2009.  Ground motion prediction equation for spectrum intensity
%from spectral acceleration equations.  Bulletin of the Seismological Society of
%America.  (in press).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input data
M=6.5;
R=15;
eps=1;
targetSI=150;   %in cm.s/s


siteprop.soiltype='rock';
faultprop.faultstyle='normal';
siteprop.g=981;
siteprop.V30=300;

periodrange=0:0.02:5;
IMR{1}=@Bradleyetal_2008_SI;
IMR{2}=@BooreAtkinson_2007_nga;
g=siteprop.g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%integration data for SI
dT=0.1;
T=0.1:dT:2.5;
%allocate integration weights for trapz rule (dT/2 for first and last, dT otherwise)
weight=dT*ones(1,length(T)); 
weight(1)=weight(1)/2; weight(length(T))=weight(length(T))/2;

[SIpredicted,sigma_SIpredicted]=feval(IMR{1},M,R,siteprop,faultprop,IMR{2});
meanSI=SIpredicted*exp(0.5*sigma_SIpredicted^2);
sigma_SInormal=meanSI*sqrt(exp(sigma_SIpredicted^2)-1);

%adjust epsilon value to get target Sa
epsilonmod=log(targetSI/SIpredicted)/sigma_SIpredicted;

%get unconditional spectra and correlation at range of periods
for i=1:length(periodrange)
    siteprop.period=periodrange(i);
    [SA(i),sigma_SA(i)]=feval(IMR{2},M,R,siteprop,faultprop);
    rho_sum=0.0;
    for j=1:length(T)
        %get spectral acceleration at each of the T values of SI and
        %convert from lognormal to normal moments
        siteprop.period=T(j);
        [SA_SI,sigma_SA_SI]=feval(IMR{2},M,R,siteprop,faultprop);
%         plot(T(j),SA_ASI,'o'); hold on
        SA_SImean=SA_SI*exp(0.5*sigma_SA_SI^2);
        sigma_SA_SInormal=SA_SImean*sqrt(exp(sigma_SA_SI^2)-1);
        %compute LN correlation and convert to normal
        [rhoij(j)]=SA_correlation(periodrange(i),T(j));
        rhoijnormal=(exp(rhoij(j)*sigma_SA(i)*sigma_SA_SI)-1)/sqrt((exp(sigma_SA(i)^2)-1)*(exp(sigma_SA_SI^2)-1));
%         plot(T(j),rhoijnormal,'o'); hold on
        omega(j)=2*pi/T(j);
        %now add to summation
        rho_sum=rho_sum+weight(j)/omega(j)*g*rhoijnormal*sigma_SA_SInormal;
    end
    rhoSI_SAnormal=rho_sum/(sigma_SInormal);
    %now convert to LN correlation
    rhoSI_SA(i)=log(1+rhoSI_SAnormal*sqrt((exp(sigma_SA(i).^2)-1)*(exp(sigma_SIpredicted.^2)-1)))/(sigma_SA(i)*sigma_SIpredicted);
    
    %get CMS
    SA_CMSmedian(i)=exp(log(SA(i))+sigma_SA(i)*rhoSI_SA(i)*epsilonmod);
    SA_CMSbeta(i)=sigma_SA(i)*sqrt(1-rhoSI_SA(i)^2);
end
% plot(periodrange,SA);
% plot(periodrange,sigma_SA);
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

figure(2)
plot(periodrange,rhoSI_SA)



