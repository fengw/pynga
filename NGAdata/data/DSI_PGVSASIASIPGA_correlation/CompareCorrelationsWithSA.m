function CompareCorrelationsWithSA

%Purpose: To compare the correlations of ASI, SI, and DSI with SA as a
%function of SA period

T=[0.01	0.02	0.03	0.04	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10];

for i=1:length(T)
    a=[0.39 0.19 0.98];
    b=[0.265 1.20 0.82];
    c=[0.04 1.2 6.1];
    d=[1.8 0.6 3];
    e=[0.15 3.4 10];
    for j=1:3
        if T(i)<=e(j)
            %mean
            rho_DSI_SA_fit(i,1)=(a(j)+b(j))/2-(a(j)-b(j))/2*tanh(d(j)*log((T(i)/c(j))));
            break
        end
    end
end

%parametric fit
% ASI,SA
for i=1:length(T)
    %a=value for x'=0; b=value for x'=1; c=x' value at midpoint between a
    %and b; e=controls slope; f=controls log scaling
    a_1=[0.927 0.823 1.05];
    b_1=[0.823 0.962 0.29];
    c_1=[0.04 0.14 0.8];
    d_1=[1.8 2.2 1.0];
    e_1=[0.075 0.3 10];
    for j=1:3
        if T(i)<=e_1(j)
            %mean
            rho_ASI_SA_fit(i,1)=(a_1(j)+b_1(j))/2-(a_1(j)-b_1(j))/2*tanh(d_1(j)*log((T(i)/c_1(j)))); 
            break
        end
    end
end
% SI,SA
for i=1:length(T)
    a_2=[0.6 0.38 0.95];
    b_2=[0.38 0.94 0.68];
    c_2=[0.045 0.33 3.10];
    d_2=[1.5 1.4 1.6];
    e_2=[0.1 1.4 10];
    for j=1:3
        if T(i)<=e_2(j)
            %mean
            rho_SI_SA_fit(i,1)=(a_2(j)+b_2(j))/2-(a_2(j)-b_2(j))/2*tanh(d_2(j)*log((T(i)/c_2(j)))); 
            break
        end
    end
end

%plot
fig1=figure(1)
axes('Parent',gcf,'FontSize',16);
semilogx(T,rho_ASI_SA_fit,'LineWidth',3,'LineStyle','-','Color',[1 0 0]); hold on;
semilogx(T,rho_SI_SA_fit,'LineWidth',3,'LineStyle','--','Color',[0 1 0]);
semilogx(T,rho_DSI_SA_fit,'LineWidth',3,'LineStyle',':','Color',[0 0 1]);
legend('\rho_{lnASI,lnSa}','\rho_{lnSI,lnSa}','\rho_{lnDSI,lnSa}');
xlabel('Period, T (s)','FontSize',16); ylabel('Corr. Coeff., \rho','FontSize',16);
grid on;
