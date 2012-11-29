function [SA,sigma_SA]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop)
%Brendon Bradley   24 Aug 2010

%Provides the attenuation relation for Sa in units of g 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input Variables:
% M             = Moment magnitude (Mw)
% R             = Source-to-site distance (km) (Rrup distance)
% siteprop      = properties of site (soil etc)
%                 siteprop.V30   -'(any real variable)' shear wave
%                                   velocity(m/s)
%                 siteprop.V30measured - yes =1 (i.e. from Vs tests); no =
%                   0 (i.e. estimated from geology)
%                 siteprop.Rx -'closest horiz distance coseismic rupture measured perpendicular to strike (km)
%                           strike from surface projection of updip edge of the fault rupture (+ve in downdip dir) (km)
%                 siteprop.Rjb -distancehorizontal dist measured
%                                   perpendicular to fault
%                 siteprop.period -'(-1),(0),(real variable)' period of vibration =-1->PGV; =0->PGA; >0->SA
%                 siteprop.Z1pt0 -'depth to the 1.0km/s shear wave velocity horizon (optional, uses default relationship otherwise)

% faultprop     = properties of fault (strikeslip etc)
%                 faultprop.Ztor -'depth to top of coseismic rupture (km)
%                 faultprop.rake -'rake angle in degrees
%                 faultprop.dip -'avg dip angle in degrees
%                 faultprop.AS - aftershock indicator (=0 no; =1 yes)
%                 faultprop.W - down dip rupture width (km)

%Output Variables:
% SA           = median SA  (or PGA or PGV)
% sigma_SA     = lognormal standard deviation of SA
%                %sigma_SA(1) = total std
                 %sigma_SA(2) = interevent std
                 %sigma_SA(3) = intraevent std

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coefficients
% 	// coefficients (index 22 is PGA and 23 is PGV):

%period independent constants
c1 = 6.75;
c4 = 4.5;
a3 = 0.265;
a4 = -0.231;
a5 = -0.398;
n = 1.18;
c = 1.88;
c2 = 50;

period = [     -1       0    0.01    0.02    0.03    0.04     0.05    0.075     0.1    0.15     0.2    0.25     0.3     0.4     0.5    0.75       1     1.5       2       3       4       5     7.5      10];
Vlin =   [  400.0   865.1	865.1	865.1   907.8   994.5	1053.5	 1085.7	 1032.5	  877.6	  748.2	  654.3	  587.1	    503	  456.6	  410.5	    400	    400	    400	    400	    400	    400	    400	    400];
b =      [ -1.955  -1.186  -1.186  -1.219  -1.273  -1.308	-1.346	 -1.471	 -1.624	 -1.931	 -2.188	 -2.381	 -2.518	 -2.657	 -2.669	 -2.401	 -1.955	 -1.025	 -0.299	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000];
a1 =     [ 5.7578  0.8040  0.8110  0.8550  0.9620  1.0370	1.1330	 1.3750	 1.5630	 1.7160	 1.6870	 1.6460	 1.6010	 1.5110	 1.3970	 1.1370	 0.9150	 0.5100	 0.1920	 -0.280	 -0.639	 -0.936	 -1.527	 -1.993];
a2 =     [-0.9046 -0.9679 -0.9679 -0.9774 -1.0024 -1.0289  -1.0508	-1.0810	-1.0833	-1.0357	-0.9700	-0.9202	-0.8974	-0.8677	-0.8475	-0.8206	-0.8088	-0.7995	-0.7960	-0.7960	-0.7960	-0.7960	-0.7960	-0.7960];
a8 =     [-0.1200 -0.0372 -0.0372 -0.0372 -0.0372 -0.0315  -0.0271	-0.0191	-0.0166	-0.0254	-0.0396	-0.0539	-0.0656	-0.0807	-0.0924	-0.1137	-0.1289	-0.1534	-0.1708	-0.1954	-0.2128	-0.2263	-0.2509	-0.2683];
a10=     [ 1.5390  0.9445  0.9445  0.9834  1.0471  1.0884	1.1333	 1.2808	 1.4613	 1.8071	 2.0773	 2.2794	 2.4201  2.5510	 2.5395	 2.1493	 1.5705	 0.3991	-0.6072	-0.9600	-0.9600	-0.9208	-0.7700	-0.6630];
a12 =    [ 0.0800  0.0000  0.0000  0.0000  0.0000  0.0000	0.0000   0.0000	 0.0000	 0.0181	 0.0309	 0.0409	 0.0491	 0.0619	 0.0719	 0.0800	 0.0800	 0.0800	 0.0800	 0.0800	 0.0800	 0.0800	 0.0800	 0.0800];
a13 =    [-0.0600 -0.0600 -0.0600 -0.0600 -0.0600 -0.0600  -0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600	-0.0600 -0.0600	-0.0600	-0.0600];
a14 =    [ 0.7000  1.0800  1.0800  1.0800  1.1331  1.1708	1.2000	 1.2000	 1.2000	 1.1683	 1.1274	 1.0956	 1.0697	 1.0288	 0.9971	 0.9395	 0.8985	 0.8409	 0.8000	 0.4793	 0.2518	 0.0754	 0.0000	 0.0000];
a15 =    [  -0.39 -0.3500 -0.3500 -0.3500 -0.3500 -0.3500  -0.3500	-0.3500	-0.3500	-0.3500	-0.3500	-0.3500	-0.3500	-0.3500	-0.3191	-0.2629	-0.2230	-0.1668	-0.1270	-0.0708	-0.0309	 0.0000  0.0000	 0.0000];
a16 =    [   0.63  0.9000  0.9000  0.9000  0.9000  0.9000	0.9000	 0.9000	 0.9000	 0.9000	 0.9000	 0.9000	 0.9000	 0.8423	 0.7458	 0.5704	 0.4460	 0.2707	 0.1463	-0.0291	-0.1535	-0.2500	-0.2500	-0.2500];
a18 =    [   0.00 -0.0067 -0.0067 -0.0067 -0.0067 -0.0067  -0.0076	-0.0093	-0.0093	-0.0093	-0.0083	-0.0069	-0.0057	-0.0039	-0.0025	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000	 0.0000];
s1_e=    [   0.59    0.59	 0.59	 0.59	0.605	0.615	 0.623	   0.63	   0.63	   0.63	   0.63	   0.63	   0.63	   0.63	   0.63	   0.63	   0.63	  0.615	  0.604	  0.589	  0.578	   0.57	  0.611	   0.64];
s2_e=    [   0.47    0.47	 0.47	 0.47	0.478	0.483	 0.488	  0.495	  0.501	  0.509	  0.514	  0.518	  0.522	  0.527	  0.532	  0.539	  0.545	  0.552	  0.558	  0.565	   0.57	  0.587	  0.618	   0.64];
s1_m=    [  0.576   0.576	0.576	0.576	0.591	0.602	  0.61	  0.617	  0.617	  0.616	  0.614	  0.612	  0.611	  0.608	  0.606	  0.602	  0.594	  0.566	  0.544	  0.527	  0.515	   0.51	  0.572	  0.612];
s2_m=    [  0.453   0.453	0.453	0.453	0.461	0.466	 0.471	  0.479	  0.485	  0.491	  0.495	  0.497	  0.499	  0.501	  0.504	  0.506	  0.503	  0.497	  0.491	    0.5	  0.505	  0.529	  0.579	  0.612];
s3 =     [   0.42    0.47	 0.42	 0.42	0.462	0.492	 0.515	   0.55	   0.55	   0.55	   0.52	  0.497	  0.479	  0.449	  0.426	  0.385	   0.35	   0.35	   0.35	   0.35    0.35	   0.35	   0.35	   0.35];
s4 =     [    0.3     0.3	  0.3	  0.3	0.305	0.309	 0.312	  0.317	  0.321	  0.326	  0.329	  0.332	  0.335	  0.338	  0.341   0.346	   0.35	   0.35	   0.35	   0.35	   0.35	   0.35	   0.35	   0.35];
rho_T_PGA=[  0.74       1	    1	    1	0.991	0.982	 0.973	  0.952	  0.929	  0.896	  0.874	  0.856	  0.841	  0.818	  0.783	   0.68	  0.607	  0.504	  0.431	  0.328	  0.255	    0.2	    0.2	    0.2];
    
%%%%%%%%%%%%%%%%%%%  
Rjb=siteprop.Rjb;
Rx=siteprop.Rx;
Vs30=siteprop.V30;

if siteprop.Z1pt0<0
    %use Z10 = Z10bar
    if Vs30 < 180
        Z10 = exp(6.75);
    elseif Vs30 <500
        Z10 = exp(6.745 - 1.35*log(Vs30/180));
    else
        Z10 = exp(5.394 - 4.48*log(Vs30/500));
    end
else
    Z10=siteprop.Z1pt0; %depth to 1.0km/s Vs horizon
end

delta=faultprop.dip; %dip in degrees
lambda=faultprop.rake; %rake in degrees
Fas = faultprop.AS;
Ztor=faultprop.Ztor;
W = faultprop.W;


frv = lambda >= 30 & lambda <= 150; % frv: 1 for lambda between 30 and 150, 0 otherwise
fnm = lambda >= -120 & lambda <= -60; % fnm: 1 for lambda between -120 and -60, 0 otherwise
HW = Rx>=0;

T=siteprop.period;

%compute PGA reference on vs=1100 m/s rock
if T~=0 | Vs30~=1100;
    siteprop.period=0;
    siteprop.V30=1100;
    [PGA1100,sigma_PGA1100]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
    siteprop.period=T;
    siteprop.V30=Vs30;
end


% interpolate between periods if neccesary
if (length(find(abs((period-T))<0.0001))==0)
    T_low=max(period(find(period<T)));
    T_hi=min(period(find(period>T)));
    
    siteprop.period=T_low;
    [SA_low,sigma_SA_low]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
    siteprop.period=T_hi;
    [SA_high,sigma_SA_high]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
    siteprop.period=T;
    
    if T_low>eps %log interpolation
        x=[log(T_low) log(T_hi)];
        Y_sa=[log(SA_low) log(SA_high)];
        SA_sigma=[sigma_SA_low' sigma_SA_high'];
        SA=exp(interp1(x,Y_sa,log(T)));
        for i=1:3
            sigma_SA(i) = interp1(x,SA_sigma(i,:),log(T));
        end
    else    %inear interpolation
        x=[T_low T_hi];
        Y_sa=[SA_low SA_high];
        SA_sigma=[sigma_SA_low' sigma_SA_high'];
        SA=interp1(x,Y_sa,T);
        for i=1:3
            sigma_SA(i) = interp1(x,SA_sigma(i,:),T);
        end
    end
    
else
    i = find(abs((period - T)) < 0.0001); % Identify the period
    
    %calculate terms in the median computation
    
    %magnitude and distance scaling
    R=sqrt(Rrup^2+c4^2);
    if M<=c1
        f1 =  a1(i) + a4*(M - c1) + a8(i)*(8.5-M)^2 + (a2(i) + a3*(M-c1))*log(R);
    else
        f1 =  a1(i) + a5*(M - c1) + a8(i)*(8.5-M)^2 + (a2(i) + a3*(M-c1))*log(R);
    end
    
    %f5 : site response scaling
    if T==-1 %PGV
        V1=862;
    elseif T<=0.5
        V1=1500;
    elseif T<=1
        V1=exp(8.0-0.795*log(T/0.21));
    elseif T<2
        V1=exp(6.76-0.297*log(T));
    elseif T>=2
        V1=700;
    end
    
    if Vs30<V1
        Vs30star = Vs30;
    else
        Vs30star = V1;
    end
    
    if Vs30<Vlin(i)
        f5 = a10(i)*log(Vs30star/Vlin(i)) - b(i)*log(PGA1100+c) + b(i)*log(PGA1100+c*(Vs30star/Vlin(i))^n);
    else
        f5 = (a10(i)+b(i)*n)*log(Vs30star/Vlin(i));
    end
    
    %f4 : hanging wall term
    T1 = max(1-Rjb/30,0);
    if Rx > W*cos(delta) | delta==90
        T2 = 1;
    else
        T2 = 0.5 + Rx/(2*W*cos(delta*pi/180));
    end
    
    if Rx>=Ztor
        T3 = 1;
    else
        T3 = Rx/Ztor;
    end
    
    T4 = max(min(M-6,1),0);
    T5 = min(1, 1-(delta-70)/20);

    f4 = a14(i)*T1*T2*T3*T4*T5;
    
    %f6 : depth to top of rutpure
    f6 = a16(i)*min(1, Ztor/10);
    
    %f8 : large distance scaling
    T6 = min(1,max(0.5*(6.5-M)+0.5,0.5));
    f8 = a18(i)*max(0,(Rrup-100))*T6;
    
    %f10: soil depth model
    if Vs30 < 180
        Z10bar = exp(6.75);
    elseif Vs30 <500
        Z10bar = exp(6.745 - 1.35*log(Vs30/180));
    else
        Z10bar = exp(5.394 - 4.48*log(Vs30/500));
    end
    
    if T < 0.35 | Vs30 > 1000
        e2 = 0;
    elseif T <= 2
        e2 = -0.25*log(Vs30/1000)*log(T/0.35);
    else
        e2 = -0.25*log(Vs30/1000)*log(2/0.35);
    end
    
    if Vs30>=1000
        a21 = 0;
    elseif (a10(i)+b(i)*n)*log(Vs30star/min(V1,1000))+e2*log((Z10+c2)/(Z10bar+c2)) < 0
        a21 = -(a10(i)+b(i)*n*log(Vs30star/min(V1,1000))/log((Z10+c2)/(Z10bar+c2)));
    else
        a21 = e2;
    end
    
    if T==-1
        a22 = 0;
    else
        a22 = 0.0625*max(T-2,0);
    end
    
    if Z10<200
        f10=a21*log((Z10 + c2)/(Z10bar + c2));
    else
        f10 = a21*log((Z10 + c2)/(Z10bar + c2)) + a22*log(Z10/200);
    end
    
    %constant displacement corner period
    TD = 10^(-1.25+0.3*M);

    % Compute median
    %%%%%Need to account for TD
    if T <= TD | (Vs30==1100 & Z10==0)
        Sa = exp(f1 + a12(i)*frv + a13(i)*fnm + a15(i)*Fas + f5 + HW*f4 + f6 +f8 + f10);
    else
       %get the closest period to TD
        Tdminus = period(length(find(period<TD)));
        siteprop.period=Tdminus;
        siteprop.V30 = 1100;
        siteprop.Z1pt0 = 0;
        [SaTdminus,sigma_SaTdminus]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
        
        Tdplus = period(length(find(period<TD))+1);
        siteprop.period=Tdplus;
        siteprop.V30 = 1100;
        siteprop.Z1pt0 = 0;
        [SaTdplus,sigma_SaTdplus]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop);
        
        SaTd=exp( log(SaTdminus) + log(SaTdplus/SaTdminus)*log(TD/Tdminus)/log(Tdplus/Tdminus) );
        siteprop.period=T;
        siteprop.V30 = Vs30;
        siteprop.Z1pt0 = Z10;
        
        %calc f5,1100
        if 1100<V1
            Vs301100star = 1100;
        else
            Vs301100star = V1;
        end
    
        if 1100<Vlin(i)
            f51100 = a10(i)*log(Vs301100star/Vlin(i)) - b(i)*log(PGA1100+c) + b(i)*log(PGA1100+c*(Vs301100star/Vlin(i))^n);
        else
            f51100 = (a10(i)+b(i)*n)*log(Vs301100star/Vlin(i));
        end
    
        %const disp model with f51100 subtracted
        Sa = exp(log(SaTd*(TD/T)^2) -f51100  + f5 + f10);
    end
    

    % Compute standard deviation
    if siteprop.V30measured==1
        s1 = s1_m(i);
        s2 = s2_m(i);
        s1_PGA = s1_m(2);
        s2_PGA = s2_m(2);
    else
        s1 = s1_e(i);
        s2 = s2_e(i);
        s1_PGA = s1_e(2);
        s2_PGA = s2_e(2);
    end
    
    if Vs30 >= Vlin(i)
        dlnAmp_dlnPGA = 0;
    else
        dlnAmp_dlnPGA=b(i)*PGA1100*(-1.0/(PGA1100+c)   +1.0/(PGA1100+ c*(Vs30/Vlin(i))^n)); %Note that error in AS08 Eq Spectra publications here
    end
    
    %intra-event
    sigma_amp=0.3;
    sigma0 = min(s1, max(s2, s1+(s2-s1)/2*(M-5) ) );
    sigma0_PGA = min(s1_PGA, max(s2_PGA, s1_PGA+(s2_PGA-s1_PGA)/2*(M-5) ) );
    sigmaB = sqrt(sigma0^2 - sigma_amp^2);
    sigmaB_PGA = sqrt(sigma0_PGA^2 - sigma_amp^2);
    sigma = sqrt( sigmaB^2 + sigma_amp^2 + dlnAmp_dlnPGA^2*sigmaB_PGA^2 + 2*dlnAmp_dlnPGA*sigmaB*sigmaB_PGA*rho_T_PGA(i) );  %Note that error in AS08 Eq Spectra publications here
    
    %inter-event
    tau0 = min(s3(i), max(s4(i), s3(i)+(s4(i)-s3(i))/2*(M-5) ) );
    tau0_PGA = min(s3(2), max(s4(2), s3(2)+(s4(2)-s3(2))/2*(M-5) ) );
    tauB = tau0;
    tauB_PGA = tau0_PGA;
    tau = sqrt( tau0^2 + dlnAmp_dlnPGA^2*tau0_PGA^2 + 2*dlnAmp_dlnPGA*tauB*tauB_PGA*rho_T_PGA(i) );
    
    %outputs
    SA=Sa;
    sigma_SA(1)=sqrt(tau^2+sigma^2);
    sigma_SA(2)=tau;  %%%%%%%%%%%%%%%%%%%%%%%
    sigma_SA(3)=sigma; %%%%%%%%%%%%%%%%%%%%%%%
 

end

