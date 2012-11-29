function [ind,M,Pm,R,Pr,faultprop]=ERF_onefault_spatialIMinvestigation(tint,faultnum)

%Brendon Bradley   6 April 2008

%gives the ERF for the "onefault_spatialIMinvestigation" scenario

%This scenario consists of 1 fault at a distance 30km from the site which
%has a magnitude recurrence relationship of Pm=4.0-0.8M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input Variables:
% tint          = the time interval over which to compute the solution
% faultnum      = the fault number to set the rupture data for

%Output Variables:
% ind           = an indicator variable. =1 if fault exists, =0 otherwise
% M             = an array of magnitude values used in the computation
% Pm            = the probabilities of the various magnitudes in M
% R             = an array of distance values of the fault from the site
% Pr            = the probabilities of the various distances in R
% faultstyle    = the style of faulting 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if faultnum==1
    %first (and only fault)
    ind=faultnum;
    R=[30];
    Pr=[1.0];
    
    delM=0.1;
    M=[4.0:delM:8.0]';         %minimum M=4.0; Max=8.0; increment=0.1
    lambdaM=exp(4.0-0.8*M);
    Pm=(1-exp(-exp(4.0-0.8.*(M-delM/2))*tint))-(1-exp(-exp(4.0-0.8.*(M+delM/2))*tint));

    faultprop.faultstyle='strikeslip';
    faultprop.Rrup=20;
    faultprop.Ztor=15;
    faultprop.rake=15;
    faultprop.dip=40;
else
    ind=0;    M=0;    Pm=0;     R=0;    Pr=0;   faultprop.faultstyle='';
end

end %end of ERF: onefault_spatialIMinvestigation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%