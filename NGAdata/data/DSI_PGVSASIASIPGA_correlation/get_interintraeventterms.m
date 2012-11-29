function [eta,epsilon]=get_interintraeventterms(residual_totalnormalised,EQ_num,sigmaT,sigma_inter,sigma_intra,alleta)
    
%Purpose: given the total residuals (normalised) obtain the inter and intra
%event terms
% 20 March 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input variables
%-----------
%residual_totalnormalised -- the total (normalised) residuals for each ground
%                       motion considered
%EQ_num            -- the EQID numbers for each GM considered
%sigmaT            -- the total logn stdev (for each GM)
%sigma_inter       -- the inter event logn stdev (for each GM)
%sigma_intra       -- the intra event logn stdev (for each GM)
%,alleta           -- whether to output [eta] as a vector of length Ngm,
%                     which will have the eta values for every GM (=1), or to
%                     output eta of length N_EQ (=0)

%output variables
%-----------
%eta               -- a vector of the inter event terms for each GM
%epsilon           -- a vector of the intra event terms for each GM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
N_EQ=0;    %initialise
EQ_num_array=0;
sigma_inter_sum=0;
sigma_intra_sum=0;
Ngm=length(residual_totalnormalised);

residual_total=residual_totalnormalised.*sigmaT; %make non-normalised

for i=1:Ngm
    found_EQid=0;
    for k=1:N_EQ
        if EQ_num_array(k,1)==EQ_num(i);
            EQ_num_array(k,2)=EQ_num_array(k,2)+1;  %number of events
            EQ_num_array(k,3)=EQ_num_array(k,3)+residual_total(i);
            sigma_inter_sum(k)=sigma_inter_sum(k)+sigma_inter(i);
            sigma_intra_sum(k)=sigma_intra_sum(k)+sigma_intra(i);
            found_EQid=1;
            break
        end
    end
    
    if found_EQid==0
        N_EQ=N_EQ+1;
        EQ_num_array(N_EQ,1)=EQ_num(i);
        EQ_num_array(N_EQ,2)=1;  %number of events
        EQ_num_array(N_EQ,3)=residual_total(i);
        sigma_inter_sum(N_EQ)=sigma_inter(i);
        sigma_intra_sum(N_EQ)=sigma_intra(i);
    end
    
end

%now compute the inter and intra event residuals
for i=1:N_EQ
    N_GM=EQ_num_array(i,2); sum_residuals=EQ_num_array(i,3);
    sigma_inter_avg(i)=sigma_inter_sum(i)/N_GM;
    sigma_intra_avg(i)=sigma_intra_sum(i)/N_GM;
    residual_EQinter(i)=sigma_inter_avg(i)^2*sum_residuals/(N_GM*sigma_inter_avg(i)^2+sigma_intra_avg(i)^2);
    eta_EQ(i)=residual_EQinter(i)/sigma_inter_avg(i);
%     A=1/(1+(sigma_intra_avg(i)/sigma_inter_avg(i))^2/N_GM)*sum_residuals/N_GM;
end

for i=1:Ngm
    %work out the EQID and corresponding inter event residual
    for j=1:N_EQ
        if EQ_num_array(j,1)==EQ_num(i);
            break
        end
    end
%     sigma_inter_avg_gm(i)=sigma_inter_avg(j); sigma_intra_avg_gm(i)=sigma_intra_avg(j);
    eta(i)=residual_EQinter(j)/sigma_inter(i);
    epsilon(i)=(residual_total(i)-eta(i)*sigma_inter(i))/sigma_intra(i);
end

if alleta==0
    clear eta
    eta = eta_EQ;
end

%debug plots
% figure(1)
% gmid=1:Ngm;
% subplot(211);
% plot(gmid,sigma_inter,'or'); hold on;
% plot(gmid,sigma_intra,'or');
% plot(gmid,sigma_inter_avg_gm,'xb'); 
% plot(gmid,sigma_intra_avg_gm,'xb');
% xlabel('gmid'); ylabel('std');
% ylim([0.2 0.7]);
% 
% subplot(212)
% plot(gmid,sigma_inter-sigma_inter_avg_gm,'og'); hold on;
% plot(gmid,sigma_intra-sigma_intra_avg_gm,'og');
% xlabel('gmid'); ylabel('std');
% ylim([-0.1 0.1]);

% figure(2)
% gmid=1:Ngm;
% plot(gmid,sigma_intra./sigmaT,'or'); hold on;
% plot([1 Ngm],[1 1]*mean(sigma_intra./sigmaT),'-');
% plot([1 Ngm],[1 1]*prctile(sigma_intra./sigmaT,16),'--');
% plot([1 Ngm],[1 1]*prctile(sigma_intra./sigmaT,84),'--');
% a=1;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
