function DSI_Sa_intraAndTotalStandardDeviation

%Purpose: To compare the intra event and total standard deviations of ASI
%and Sa values

%define arbitrary scenario details
M=6.5;
Rjb=30;
siteprop.V30=400;
period=0.01:0.01:10;
period2=2.0:0.01:5; period3=0.1:0.01:0.5; period4=0.1:0.01:2.5;
faultprop.faultstyle='strikeslip';
outputplot=4;

%consider total and interevent standard deviations
for i=1:length(period)
    siteprop.period=period(i);
    [SA(i),sigma_SA(i,1:3)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
end

if outputplot==1
    figure(1); hold on;
    plot(period,SA,'LineWidth',2,'LineStyle','-','Color',[0 0 1])
    plot(period,SA.*exp(sigma_SA(:,1)'),'LineWidth',1,'LineStyle','-','Color',[0 0 1])
    plot(period,SA.*exp(sigma_SA(:,2)'),'LineWidth',1,'LineStyle','-','Color',[1 0 0])
    plot(period,SA.*exp(-sigma_SA(:,1)'),'LineWidth',1,'LineStyle','-','Color',[0 0 1])
    plot(period,SA.*exp(-sigma_SA(:,2)'),'LineWidth',1,'LineStyle','-','Color',[1 0 0])
    legend('Median','Median','Median \pm total std','Median \pm inter std'); xlabel('Period'); ylabel('Sa (g)');
elseif outputplot==2
    figure(2); hold on;
    plot(period,sigma_SA(:,1),'LineWidth',2,'LineStyle','-','Color',[0 0 1])
    plot(period,sigma_SA(:,2)','LineWidth',2,'LineStyle','--','Color',[1 0 0]);
    legend('Total stdev','Inter stdev'); xlabel('Period'); ylabel('sigma');
end

%now look at how the intra and total stdev effects the ASI computations
R=Rjb; siteprop.g=9.81;
[DSI,sigma_DSI]=Bradleyetal_2011_DSI(M,R,siteprop,faultprop,@BooreAtkinson_2007_nga);
[ASI,sigma_ASI]=Bradleyetal_2008_ASI(M,R,siteprop,faultprop,@BooreAtkinson_2007_nga);
siteprop.g=981; 
[SI,sigma_SI]=Bradleyetal_2008_SI(M,R,siteprop,faultprop,@BooreAtkinson_2007_nga);

if outputplot==3
    figure(3); hold on;
    plot(1,[DSI*exp(-sigma_DSI(1)) DSI DSI*exp(sigma_DSI(1))],'LineWidth',3,'Marker','o','Color',[0 0 1],'MarkerSize',10);
    plot(1,[DSI*exp(-sigma_DSI(2)) DSI DSI*exp(sigma_DSI(2))],'LineWidth',3,'Marker','sq','Color',[1 0 0],'MarkerSize',10);
elseif outputplot==4
    h=figure(4); 
    axes('FontSize',14,'Parent',h); 
    semilogx(period,sigma_SA(:,1),'LineWidth',2,'LineStyle','-','Color',[1 0 0]); hold on; box on;
    semilogx(period,sigma_SA(:,3)','LineWidth',2,'LineStyle','--','Color',[1 0 0]);
    semilogx(period2,sigma_DSI(1)*ones(1,length(period2)),'LineWidth',3,'LineStyle','-','Color',[0.749 0.749 0]);
    semilogx(period2,sigma_DSI(3)*ones(1,length(period2)),'LineWidth',3,'LineStyle','--','Color',[0.749 0.749 0]);
    semilogx(period3,sigma_ASI(1)*ones(1,length(period3)),'LineWidth',3,'LineStyle','-','Color',[0 1 0]);
    semilogx(period3,sigma_ASI(3)*ones(1,length(period3)),'LineWidth',3,'LineStyle','--','Color',[0 1 0]);
    semilogx(period4,sigma_SI(1)*ones(1,length(period4)),'LineWidth',3,'LineStyle','-','Color',[0 0 1]);
    semilogx(period4,sigma_SI(3)*ones(1,length(period4)),'LineWidth',3,'LineStyle','--','Color',[0 0 1]);
    legend('Sa (total)','Sa (intra)','DSI (total)','DSI (intra)','ASI (total)','ASI (intra)','SI (total)','SI (intra)'); xlabel('Period (s)','FontSize',12); ylabel('Dispersion, \sigma_{lnX}','FontSize',12);
end

