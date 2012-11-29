function DSI_PGVSASIASIPGA_correlation_fit
clc
fclose all
%Brendon BRadley 1 July 2010
%TO fit the correlations parameterically

%               rho_{lnDSI,lnSA}       
% T  |        BA08  |    CY08 |    CB08  |    AS08 |    mean rho   |   stddev z    
rho_computed=[ 0.010      0.417    0.386    0.362    0.400    0.391    0.051   
 0.020      0.405    0.375    0.355    0.387    0.381    0.049   
 0.030      0.387    0.349    0.328    0.370    0.359    0.053   
 0.040      0.363    0.318    0.309    0.341    0.333    0.051   
 0.050      0.331    0.291    0.284    0.319    0.306    0.052   
 0.075      0.306    0.259    0.233    0.279    0.269    0.057   
 0.100      0.283    0.255    0.236    0.270    0.261    0.050   
 0.150      0.288    0.246    0.239    0.283    0.264    0.052   
 0.200      0.329    0.286    0.264    0.324    0.301    0.056   
 0.250      0.357    0.314    0.300    0.355    0.332    0.052   
 0.300      0.387    0.354    0.332    0.384    0.364    0.050   
 0.400      0.435    0.415    0.386    0.435    0.418    0.046   
 0.500      0.469    0.452    0.435    0.476    0.458    0.045   
 0.750      0.565    0.547    0.527    0.551    0.548    0.042   
 1.000      0.642    0.638    0.624    0.631    0.634    0.033   
 1.500      0.779    0.763    0.756    0.747    0.761    0.044   
 2.000      0.877    0.865    0.861    0.855    0.864    0.049   
 3.000      0.963    0.961    0.959    0.959    0.961    0.036   
 4.000      0.980    0.980    0.978    0.980    0.979    0.034   
 5.000      0.945    0.941    0.938    0.940    0.941    0.042   
 7.500      0.862    0.838    0.848    0.850    0.850    0.059   
10.000      0.821    0.821    0.828    0.823    0.823    0.062           
];

%rho DSI with  i) ASI; ii) PGA; iii) SI; iv) PGV
%                BA08  |    CY08 |    CB08  |    AS08 |    mean rho   |   std z   |
rho_computed2=[0.399        0.370        0.339       0.396       0.376        0.053      
               0.423        0.388        0.364       0.404       0.395        0.049      
               0.800        0.778        0.776       0.773       0.782        0.045      
               0.816        0.772        0.804       0.807       0.800        0.062      
];   

% -------------------------------------
Traw=rho_computed(:,1); rho_DSI_SA=rho_computed(:,2:5); rho_DSI_SA_mean=rho_computed(:,6);
sigma_z_DSI_SA=rho_computed(:,7);

rho_DSI_ASI=rho_computed2(1,1:4); rho_DSI_ASI_mean=rho_computed2(1,5);  sigma_z_DSI_ASI=rho_computed2(1,6);
rho_DSI_PGA=rho_computed2(2,1:4); rho_DSI_PGA_mean=rho_computed2(2,5);  sigma_z_DSI_PGA=rho_computed2(2,6);
rho_DSI_SI=rho_computed2(3,1:4); rho_DSI_SI_mean=rho_computed2(3,5);  sigma_z_DSI_SI=rho_computed2(3,6);
rho_DSI_PGV=rho_computed2(4,1:4); rho_DSI_PGV_mean=rho_computed2(4,5);  sigma_z_DSI_PGV=rho_computed2(4,6);


%compute +- one sigma correlations
[z_DSI_SA_mean]=FisherZTransformation(rho_DSI_SA_mean);
rho_DSI_SA_84=InvFisherZTransformation(z_DSI_SA_mean+sigma_z_DSI_SA);
rho_DSI_SA_16=InvFisherZTransformation(z_DSI_SA_mean-sigma_z_DSI_SA);

T=Traw;
T=0.01:0.01:10;
%parametric fit
% DSI,SA
for i=1:length(T)
%     a1=0.90; b1=0.8; c1=0.042; d1=2.2; e1=0.075; f1=1;
%     a2=0.8; b2=0.955;  c2=0.15; d2=2.2; e2=0.3;  f2=1;
%     a3=1.; b3=0.22; c3=0.85; d3=0.9; f3=1.2;
    %a=value for x'=0; b=value for x'=1; c=x' value at midpoint between a
    %and b; e=controls slope; f=controls log scaling
    a=[0.39 0.19 0.98];
    b=[0.265 1.20 0.82];
    c=[0.04 1.2 6.1];
    d=[1.8 0.6 3];
    e=[0.15 3.4 10];
    %sigma params
    as=0.05;
    x1=0.6; x2=1;
    y1=as; y2=0.13;
    m=(y2-y1)/(x2-x1);
    cs=as-x1*m;
    for j=1:3
        if T(i)<=e(j)
            %mean
            rho_DSI_SA_fit(i,1)=(a(j)+b(j))/2-(a(j)-b(j))/2*tanh(d(j)*log((T(i)/c(j))));
            %sigma
            [z_DSI_SA_fit]=FisherZTransformation(rho_DSI_SA_fit(i,1));
            sigma_DSI_SA_fit(i)=as; %max(as,m*rho_DSI_SA_fit(i,1)+cs);
            %16th and 84th percentiles
            rho_DSI_SA_fit(i,2)=InvFisherZTransformation(z_DSI_SA_fit+sigma_DSI_SA_fit(i));
            rho_DSI_SA_fit(i,3)=InvFisherZTransformation(z_DSI_SA_fit-sigma_DSI_SA_fit(i));
            break
        end
    end
end
%output values at the boundaries of the piece-wise predictions
fprintf('Continuity at the piece-wise boundaries of the DSI, SA model \n');
fprintf('T_e      rho_left     rho_right   \n');
for i=1:length(e)-1
    rho_left=(a(i)+b(i))/2-(a(i)-b(i))/2*tanh(d(i)*log((e(i)/c(i)))); 
    rho_right=(a(i+1)+b(i+1))/2-(a(i+1)-b(i+1))/2*tanh(d(i+1)*log((e(i)/c(i+1)))); 
    fprintf('%6.3f  %8.3f %8.3f \n',e(i),rho_left,rho_right);
end
fprintf('\n');

outtype=1;    %  = 1 : DSI, SA
              %  = 4 : Examine sigma as a function of mean correlation
              
if outtype==1

%     fig2=figure(2); 
%     axes('Parent',gcf,'FontSize',16);
%     semilogx(T,rho_DSI_SA(:,1),'LineWidth',3,'LineStyle','-','Color',[1 0 0]); hold on;
%     semilogx(T,rho_DSI_SA(:,2),'LineWidth',3,'LineStyle','-','Color',[0 1 0]);
%     semilogx(T,rho_DSI_SA(:,3),'LineWidth',3,'LineStyle','-','Color',[0 0 1]);
%     semilogx(T,rho_DSI_SA(:,4),'LineWidth',3,'LineStyle','-','Color',[0.5 0.5 0.5]);
%     semilogx(T,rho_DSI_SA_mean,'LineWidth',3,'LineStyle','-','Color',[0 0 0]);
%     semilogx(T,rho_DSI_SA_16,'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
%     semilogx(T,rho_DSI_SA_84,'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
%     ylabel('Corr coeff, \rho_{lnDSI,lnSA}'); xlabel('Period, T (s)'); grid on; 

    fig3=figure(3); %compare mean and parametric soln
    axes('Parent',gcf,'FontSize',16);
    % subplot(211);
    semilogx(Traw,rho_DSI_SA_mean,'LineWidth',3,'LineStyle','-','Color',[0 0 0],'Marker','o'); hold on;
    semilogx(T,rho_DSI_SA_fit(:,1),'-r','LineWidth',3); 
    
    semilogx(Traw,rho_DSI_SA_16,'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    semilogx(Traw,rho_DSI_SA_84,'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    
    semilogx(T,rho_DSI_SA_fit(:,2:3),'--r','LineWidth',3); 
    grid on;
    set(gcf,'units','normalized'); set(gcf,'Position',[0.1 0.5 0.6 0.4]);
    xlabel('Period, T (s)'); ylabel('Corr coeff, \rho_{lnDSI,lnSA}');
    xlim([0.01 10]); set(gcf,'units','normalized'); set(gcf,'Position',[0.5 0.5 0.5 0.4]);
    legend('Empirical','Parametric','Location','SouthEast');
    ylim([0 1]);
%     set(gcf,'units','normalized'); set(gcf,'Position',[0.5 0.05 0.45 0.85]);

    fig4=figure(4);
    axes('Parent',gcf,'FontSize',16);
    % subplot(212);
    semilogx(Traw,sigma_z_DSI_SA,'LineWidth',3,'LineStyle','-','Color',[0 0 0],'Marker','o'); hold on;
    semilogx(T,sigma_DSI_SA_fit,'-r','LineWidth',3); grid on;
    set(gcf,'units','normalized'); set(gcf,'Position',[0.1 0.05 0.6 0.4]);
    xlabel('Period, T (s)'); ylabel('Std. dev., \sigma_{z_{lnDSI,lnSA}}');
    xlim([0.01 10]); set(gcf,'units','normalized'); set(gcf,'Position',[0.5 0.5 0.5 0.4]);
    legend('Empirical','Parametric');
    
elseif outtype==4
    %examine uncertainty in correlations as a function of mean correlation
    fig1=figure(1)
    axes('Parent',gcf,'FontSize',16);
    plot(rho_DSI_SA_mean,sigma_z_DSI_SA,'ro','LineWidth',3); hold on;
    plot(rho_DSI_SI_mean,sigma_z_DSI_SI,'gsq','LineWidth',4);
    plot(rho_DSI_ASI_mean,sigma_z_DSI_ASI,'b^','LineWidth',4);
    plot(rho_DSI_PGA_mean,sigma_z_DSI_PGA,'k*','LineWidth',8);
    grid on;
    legend('DSI, SA','DSI, SI','DSI, ASI','DSI, PGA','Location','NorthWest');
    xlabel('Corr. coef., \rho_{lnIM_i,lnIM_j}'); ylabel('Std dev. z, \sigma_z');
    
    %model used for SI, ASI etc correlation
    x1=0.3; y1=0.08;
    x2=0.9; y2=0.045;
      %sigma_z=a*rho+b form
    a=(y2-y1)/(x2-x1)
    b=y1+a*(0-x1)
    plot([x1 x2],[y1 y2],'LineWidth',3,'Color',[0 0 0],'LineStyle','--');
    
    %model used for this DSI - SA correlation
    rho=0.2:0.02:0.97;
    sigz=as*ones(length(rho),1); %max(as,m*rho+cs);
    plot(rho,sigz,'LineWidth',3,'Color',[0 0 0],'LineStyle','-');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z]=FisherZTransformation(rho)
%Computes the Fisher Z Transformation of the correlation

z=0.5*log((1+rho)./(1-rho));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho]=InvFisherZTransformation(z)
%Computes the Inverse Fisher Z Transformation of the correlation

rho=(exp(2*z)-1)./(1+exp(2*z));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

