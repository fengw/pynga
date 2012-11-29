function [SA_CMSmedian,SA_CMSbeta]=CMSspectra_approximate_ASI(M,R,epsilon,targetASI,periodrange,IMR,siteprop,faultprop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Purpose: To compute the approximate CMS spectra fixed at ASI
%by using the mean M,R,epsilon values from PSHA deaggregation

%required M-files in same folder
%   SA_correlation.m

%Reference: 
%Baker, JW and Cornell, CA, 2006.  Spectral acceleration, record
%selection and epsilon.  Earthquake Engineering and Structural Dynamics.
%35: 1077-1095.
%Baker, JW.  2005.  Vector-valued Intensity Measures for probabilistic
%seismic demand analysis, PhD Thesis.  Stanford University, CA.
%Bradley, BA. 2009.  Ground motion prediction and selection using
%acceleration spectrum intensity.  Bulletin of the Seismological Society of
%America.  (In preparation).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input data
M=6.5;
R=50;
eps=1.5;
targetASI=0.3;   %in g


siteprop.soiltype='rock';
faultprop.faultstyle='normal';
siteprop.g=981;
siteprop.V30=300;

periodrange=0:0.02:4;
IMR{1}=@Bradleyetal_2008_ASI;
IMR{2}=@BooreAtkinson_2007_nga;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%integration data for ASI
dT=0.05;
T=0.1:dT:0.5;
%allocate integration weights for trapz rule (dT/2 for first and last, dT otherwise)
weight=dT*ones(1,length(T)); 
weight(1)=weight(1)/2; weight(length(T))=weight(length(T))/2;

[ASIpredicted,sigma_ASIpredicted]=feval(IMR{1},M,R,siteprop,faultprop,IMR{2});
meanASI=ASIpredicted*exp(0.5*sigma_ASIpredicted^2);
sigma_ASInormal=meanASI*sqrt(exp(sigma_ASIpredicted^2)-1);

%adjust epsilon value to get target Sa
epsilonmod=log(targetASI/ASIpredicted)/sigma_ASIpredicted;
% epsilonmod=eps; %temp solution

%get unconditional spectra and correlation at range of periods
for i=1:length(periodrange)
    siteprop.period=periodrange(i);
    [SA(i),sigma_SA(i)]=feval(IMR{2},M,R,siteprop,faultprop);
    rho_sum=0.0;
    for j=1:length(T)
        %get spectral acceleration at each of the T values of ASI and
        %convert from lognormal to normal moments
        siteprop.period=T(j);
        [SA_ASI,sigma_SA_ASI]=feval(IMR{2},M,R,siteprop,faultprop);
%         plot(T(j),SA_ASI,'o'); hold on
        SA_ASImean=SA_ASI*exp(0.5*sigma_SA_ASI^2);
        sigma_SA_ASInormal=SA_ASImean*sqrt(exp(sigma_SA_ASI^2)-1);
        %compute LN correlation and convert to normal
        [rhoij(j)]=SA_correlation(periodrange(i),T(j));
        rhoijnormal=(exp(rhoij(j)*sigma_SA(i)*sigma_SA_ASI)-1)/sqrt((exp(sigma_SA(i)^2)-1)*(exp(sigma_SA_ASI^2)-1));
%         plot(T(j),rhoijnormal,'o'); hold on
        %now add to summation
        rho_sum=rho_sum+weight(j)*rhoijnormal*sigma_SA_ASInormal;
    end
    rhoASI_SAnormal=rho_sum/(sigma_ASInormal);
    %now convert to LN correlation
    rhoASI_SA(i)=log(1+rhoASI_SAnormal*sqrt((exp(sigma_SA(i).^2)-1)*(exp(sigma_ASIpredicted.^2)-1)))/(sigma_SA(i)*sigma_ASIpredicted);
    
    %get CMS
    SA_CMSmedian(i)=exp(log(SA(i))+sigma_SA(i)*rhoASI_SA(i)*epsilonmod);
    SA_CMSbeta(i)=sigma_SA(i)*sqrt(1-rhoASI_SA(i)^2);
end
% plot(periodrange,SA);
% plot(periodrange,sigma_SA);
%84th and 16th percentiles CMS
SA_CMS16=exp(log(SA_CMSmedian)-SA_CMSbeta);
SA_CMS84=exp(log(SA_CMSmedian)+SA_CMSbeta);

%plotting
figure(1)
loglog(periodrange,SA,'-r');
hold on
loglog(periodrange,SA_CMSmedian,'b-');
loglog(periodrange,SA_CMS16,'b--');
loglog(periodrange,exp(log(SA)+epsilonmod*sigma_SA),'--r');
loglog(periodrange,SA_CMS84,'b--');
loglog(periodrange,exp(log(SA)+epsilonmod*sigma_SA),'--r');

figure(2)
semilogx(periodrange,rhoASI_SA)



