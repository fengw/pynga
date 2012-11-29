function [eta,epsilon]=get_interintraeventterms(residual_totalnormalised,EQ_num,sigmaT,sigma_inter,sigma_intra)
    
%Purpose: given the total residuals (normalised) obtain the inter and intra
%event terms
% 20 March 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input variables
%-----------
%residual_totalnormalised -- the total (normalised) residuals for each ground
%                       motion considered
%EQ_num            -- the EQID numbers for each GM considered
%sigmaT            -- the total logn stdev
%sigma_inter       -- the inter event logn stdev
%sigma_intra       -- the intra event logn stdev

%output variables
%-----------
%eta               -- a vector of the inter event terms for each GM
%epsilon           -- a vector of the intra event terms for each GM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
N_EQ=0;    %initialise
EQ_num_array=0;
Ngm=length(residual_totalnormalised);

residual_total=residual_totalnormalised.*sigmaT; %make non-normalised

for i=1:Ngm
    found_EQid=0;
    for k=1:N_EQ
        if EQ_num_array(k,1)==EQ_num(i);
            EQ_num_array(k,2)=EQ_num_array(k,2)+1;  %number of events
            EQ_num_array(k,3)=EQ_num_array(k,3)+residual_total(i);
            found_EQid=1;
            break
        end
    end
    
    if found_EQid==0
        N_EQ=N_EQ+1;
        EQ_num_array(N_EQ,1)=EQ_num(i);
        EQ_num_array(N_EQ,2)=1;  %number of events
        EQ_num_array(N_EQ,3)=residual_total(i);
    end
    
end

%now compute the inter and intra event residuals
for i=1:N_EQ
    N_GM=EQ_num_array(i,2); sum_residuals=EQ_num_array(i,3);
    residual_EQinter(i)=sigma_inter^2*sum_residuals/(N_GM*sigma_inter^2+sigma_intra^2);
    A=1/(1+(sigma_intra/sigma_inter)^2/N_GM)*sum_residuals/N_GM;
end

for i=1:Ngm
    %work out the EQID and corresponding inter event residual
    for j=1:N_EQ
        if EQ_num_array(j,1)==EQ_num(i);
            break
        end
    end
    eta(i)=residual_EQinter(j)/sigma_inter;
    epsilon(i)=(residual_total(i)-eta(i)*sigma_inter)/sigma_intra;
end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%