function [H,P,KSSTAT,CV]=KStest_lognormality(alpha,inputdata,outplot,moments)

%purpose: to use the Kolmogrov-Smirnov goodness-of-fit test to illustrate
%the applicability of the lognormal distribution at fitting some data

%input: inputdata - one column array of raw data
%       alpha     - the confidence level to perform the KS test
%       outplot   - an indicator for plotting the KS test result
%       moments   - (optional) a 2x1 array giving the mean and std of the
%                   lognormal distribution to use in the KS test - if not supplied then
%                   the sample mean and std are used
%output:H         - the null hypotesis result (=1 reject null hypothesis)
%       P         - p-value for the significance test
%       KSSTAT    - the KS test statistic
%       CV        - cut-off value for determining if KSSTAT is significant
%       see KSTEST for further details on the output variables

format long

if nargin<4
    mu=mean(inputdata);
    stdev=std(inputdata,0,2);
    stdlnx=sqrt(log((stdev/mu)^2+1));
    mulnx=log(mu)-0.5*stdlnx^2;
else
    mulnx=moments(1);
    stdlnx=moments(2);
end

datafit=0.9*min(inputdata):(1.1*max(inputdata)-0.9*min(inputdata))/100:1.1*max(inputdata);
CDF=logncdf(datafit,mulnx,stdlnx);

[H,P,KSSTAT,CV] = kstest(inputdata',[datafit' CDF'],alpha);
if outplot==1
    figure(1); plot(datafit,CDF','-r','LineWidth',2); hold on
    plot(datafit,[(CDF+CV)' (CDF-CV)'],'--r','LineWidth',1.5);
    [h,stats]=cdfplot(inputdata);
    axis([0 max(datafit) 0 1]);
end




