function [PIM,IMdeagg]=PSHAtool(ERF,tint,IMR,siteprop,IM);
%PSHA tool version 2 6 april 2008.  
%Author: Brendon Bradley      
%
%Carries out the PSHA methodology

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input Variables:
% ERF           = the function handle for the earthquake rupture forecast
% tint          = time interval to compute the solution for
% IMR           = the Intensity Measure relationship
% siteprop      = site properties (soil type etc)
% IM            = the IM values to compute the solution for

%Output Variables:
% PIM           = the probability of the IM values 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
faultnum=1; nfeval=0;

while faultnum~=0      %iterate over each fault
    %CORE PSHA CODE
    
    %get data from ERF
    [faultnum,M,Pm,R,Pr,faultprop]=feval(ERF,tint,faultnum);
    if faultnum~=0
        for i=1:length(IM)
            im=IM(i);
            Pim(i,faultnum)=0;
            for k=1:length(M)
                for l=1:length(R)
                    nfeval=nfeval+1;
                    if length(IMR)==1
                        [medianIM,sigmaIM]=feval(IMR{1},M(k),R(l),siteprop,faultprop);
                    elseif length(IMR)==2    %allows for attenuation relations comprised of multiple relationships
                        [medianIM,sigmaIM]=feval(IMR{1},M(k),R(l),siteprop,faultprop,IMR{2});
                    end
                    epsilon=(log(im)-log(medianIM))/sigmaIM;
                    dPim=(1-normcdf(epsilon,0,1))*Pm(k)*Pr(l);
                    Pim(i,faultnum)=Pim(i,faultnum)+dPim;
                    IMdeagg(nfeval,1:4,i)=[M(k) R(l) epsilon dPim];
                end
            end
        end
    else
        break
    end
    %next fault
    faultnum=faultnum+1;
end
%END OF CORE CODE
PIM=1.0-prod((1-Pim),2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%