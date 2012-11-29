function PulsePeriodForDSI

%Purpose: To examine the variation of pulse period with magnitude given by
%some empirical models - then compared with the definition of DSI

M=5:0.1:8;

%ShahiBaker2010
a=[  -5.73 -6.37 -2.2*log(10)];
b=[   0.99 1.03 0.4*log(10)];
% sig=[ 0.56 0.57 0];

for j=1:length(a)
    for i=1:length(M)
        mulnTp_M(i,j) = a(j) + b(j)*M(i);
    end
end

fig1=figure(1);
axes('Parent',fig1,'FontSize',16);

%now plot DSI range from 2.0 - 5.0
x=[5.5 8]';
y=[2 2;5-2 5-2]';
% y = [y1; (y2-y1)]'; 
ha = area(x, y); hold on;
set(ha(1), 'FaceColor', 'none') % this makes the bottom area invisible
set(ha, 'LineStyle', 'none')
set(ha(2), 'FaceColor', [0.7 0.7 0.7])
set(gca, 'Layer', 'top')
set(gca,'YScale','log');
grid on

color=[1 0 0;
       0 1 0;
       0 0 1];
for i=1:length(a)
    semilogy(M,exp(mulnTp_M(:,i)),'LineWidth',3,'LineStyle','-','Color',color(i,:)); hold on;
%     semilogy(M,exp(mulnTp_M(:,i)+sig(i)),'LineWidth',2,'LineStyle','--','Color',color(i,:));
%     semilogy(M,exp(mulnTp_M(:,i)-sig(i)),'LineWidth',2,'LineStyle','--','Color',color(i,:));
end
xlim([5 8]); ylim([0.2 10]); grid on;
ylabel('Pulse period, T_p (s)');  xlabel('Magnitude, M');



