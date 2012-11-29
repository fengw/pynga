function run_SaSpatialCorrelation_JayaramBaker
%Purpose: To test the above relationship

T=0:0.1:10;
h=1;
casei=[1 1 1 2 2 2];
bound=[0 1 -1 0 1 -1];

fig1=figure(1); axes('Parent',fig1,'FontSize',14); hold on;

for j=1:6
    for i=1:length(T)
        [rho(i,j)]=SaSpatialCorrelation_JayaramBaker(h,T(i),casei(j),bound(j))
    end
    b(:,j)=-3*h./log(rho(:,j));
    plot(T,b(:,j),'LineWidth',3);
end
%now invert to get the range




xlabel('Period, T (s)'); ylabel('Range, b (km)');
ylim([0 100]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%