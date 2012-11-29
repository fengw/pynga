function dTsensitivity

%to investigate the sensitivity of the solution to the time step used
clc
format long
%ERF = earthquake rupture forecast 
ERF=@ERF_onefault_spatialIMinvestigation; 
tint=1; %time interval over which to compute the solution

%IMR = intensity measure relationship
% IMR{1}=@BooreAtkinson_2007_nga;    DanciuTselentis_2007_Sa
IMR{2}=@Bradleyetal_2011_DSI;
IM=0.0:10:300; %%range of IM values to compute solution for

%site properties - soil type etc
% siteprop.soiltype='rock';

faultprop.faultstyle='normal';
siteprop.g=981;
siteprop.V30=760;
% 

    M=[7.0 6.5]; R=[30 20];
    Nt=[3 5 9 21 41 81];
    for j=1:length(M)
        for i=1:length(Nt)
            for k=1:2
                siteprop.dTtype=k-1;
                siteprop.NTi=Nt(i);
                [DSI(i,j,k),sig]=Bradleyetal_2011_DSI_variabledT(M(j),R(j),siteprop,faultprop,@BooreAtkinson_2007_nga);
                sigma_DSI(i,j,k)=sig(1);
            end
        end
    end
    

outplot=2;
if outplot==1
    %mean
    fig1=figure(1);
    axes('Parent',fig1,'FontSize',14); 
    semilogx(Nt,DSI(:,1,1),'Marker','o','LineStyle','-','LineWidth',3,'Color',[0 0.498 0]);  hold on;
    semilogx(Nt,DSI(:,2,1),'Marker','o','LineStyle','-','LineWidth',3,'Color',[1 0 0]); 
    semilogx(Nt,DSI(:,1,2),'Marker','o','LineStyle','--','LineWidth',3,'Color',[0 0.498 0]);
    semilogx(Nt,DSI(:,2,2),'Marker','o','LineStyle','--','LineWidth',3,'Color',[1 0 0]);
    for j=1:length(M)
        semilogx([min(Nt) max(Nt)],[DSI(length(Nt),j,1) DSI(length(Nt),j,1)],'LineStyle','--','Color',[0.502 0.502 0.502])
    end
    legend('M=7.0, R=30 km','M=6.5, R=20 km','Location','SouthEast')
    xlabel('Number of integration points, N_{pts}');
    ylabel('Median disp. spectrum intensity,  {\itDSI}_{50} (cm.s)');
%     ylim([0.08 0.12]);
    set(fig1,'Units','normalized'); set(fig1,'Position',[0.3 0.2 0.4 0.5]);
    grid on;
else
    %std
    fig1=figure(1);
    axes('Parent',fig1,'FontSize',14); 
    semilogx(Nt,sigma_DSI(:,1,1),'Marker','o','LineStyle','-','LineWidth',3,'Color',[0 0.498 0]);  hold on;
    semilogx(Nt,sigma_DSI(:,2,1),'Marker','o','LineStyle','-','LineWidth',3,'Color',[1 0 0]); 
    semilogx(Nt,sigma_DSI(:,1,2),'Marker','o','LineStyle','--','LineWidth',3,'Color',[0 0.498 0]);
    semilogx(Nt,sigma_DSI(:,2,2),'Marker','o','LineStyle','--','LineWidth',3,'Color',[1 0 0]);
    
    for j=1:length(M)
        semilogx([min(Nt) max(Nt)],[sigma_DSI(length(Nt),j,1) sigma_DSI(length(Nt),j,1)],'--g')
    end
    legend('M=7.0, R=30 km','M=6.5, R=20 km','Location','SouthEast')
    xlabel('Number of integration points, N_{pts}');
    ylabel('Sigma disp. spectrum intensity,  \sigma_{lnDSI} '); grid on;
%     ylim([0.53 0.58]);
    set(fig1,'Units','normalized'); set(fig1,'Position',[0.3 0.2 0.4 0.5]);
end









