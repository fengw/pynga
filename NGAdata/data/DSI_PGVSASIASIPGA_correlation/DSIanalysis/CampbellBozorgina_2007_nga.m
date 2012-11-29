function [SA,sigma_SA]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop)
%Brendon Bradley   6 April 2008

%Provides the attenuation relation for peak horizontal acceleration, PHA
%in units of g 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input Variables:
% M             = Moment magnitude (Mw)
% Rrup          = 'closest distance coseismic rupture (km)
% siteprop      = properties of site (soil etc)
%                 siteprop.Rjb -Source-to-site distance (km) (Joyner Boore
%                 distance)
%                 siteprop.V30   -'(any real variable)' shear wave velocity(m/s)
%                 siteprop.period -'(-1,-10),(0),(real variable)' period of
%                                   vibration =-1->PGV; -10->PGD; =0->PGA; >0->SA
%                 siteprop.Zvs -'depth to the 2.5km/s shear wave velocity
%                                  horizon
%                 siteprop.orientation -'random'
%                                      -'average'

% faultprop     = properties of fault (strikeslip etc)
%                 faultprop.Ztor -'depth to top of coseismic rupture (km)
%                 faultprop.rake -'rake angle in degrees
%                 faultprop.dip -'avg dip angle in degrees

%Output Variables:
% SA           = median SA  (or PGA or PGV)
% sigma_SA     = lognormal standard deviation in SA
                 %sigma_SA(1) = total std
                 %sigma_SA(2) = interevent std
                 %sigma_SA(3) = intraevent std

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coefficients
period = [0.01	  0.02	  0.03	  0.05	 0.075	   0.1	  0.15	   0.2	  0.25	   0.3	   0.4	   0.5	  0.75	     1	   1.5	     2	     3	      4	      5	    7.5	     10      0      -1     -10];
c0 = [  -1.715	 -1.68	-1.552	-1.209	-0.657	-0.314	-0.133	-0.486	 -0.89	-1.171	-1.466	-2.569	-4.844	-6.406	-8.692	-9.701 -10.556	-11.212	-11.684	-12.505	-13.087	-1.715	 0.954	 -5.27];
c1 = [     0.5	   0.5	   0.5	   0.5	   0.5	   0.5	   0.5	   0.5	   0.5	   0.5	   0.5	 0.656	 0.972	 1.196	 1.513	   1.6	   1.6	    1.6	    1.6	    1.6	    1.6	   0.5	 0.696	   1.6];
c2 = [   -0.53	 -0.53	 -0.53	 -0.53	 -0.53	 -0.53	 -0.53	-0.446	-0.362	-0.294	-0.186	-0.304	-0.578	-0.772	-1.046	-0.978  -0.638	 -0.316	  -0.07	  -0.07	  -0.07	 -0.53	-0.309	 -0.07];
c3 = [  -0.262	-0.262	-0.262	-0.267	-0.302	-0.324	-0.339	-0.398	-0.458	-0.511	-0.592	-0.536	-0.406	-0.314	-0.185	-0.236  -0.491	  -0.77	 -0.986	 -0.656	 -0.422	-0.262	-0.019  	 0];
c4 = [  -2.118	-2.123	-2.145	-2.199	-2.277	-2.318	-2.309	 -2.22	-2.146	-2.095	-2.066	-2.041	    -2	    -2	    -2	    -2	    -2	     -2	     -2	     -2	     -2	-2.118	-2.016	    -2];
c5 = [    0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	  0.17	   0.17	   0.17	   0.17	   0.17	  0.17	  0.17	  0.17];
c6 = [     5.6	   5.6	   5.6	  5.74	  7.09	  8.05	  8.79	   7.6	  6.58	  6.04	   5.3	  4.73	     4	     4	     4	     4	     4	      4	      4	      4	      4	   5.6	     4	     4];
c7 = [    0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	  0.28	 0.255	 0.161	 0.094	     0	      0	      0	      0	      0	  0.28	 0.245	     0];
c8 = [   -0.12	 -0.12	 -0.12	 -0.12	 -0.12	-0.099	-0.048	-0.012	     0	     0	     0	     0	     0	     0	     0	     0	     0	      0	      0	      0	      0	 -0.12	     0	     0];
c9 = [    0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	  0.49	 0.371	 0.154	      0	      0	      0	      0	  0.49	 0.358	     0];
c10 = [  1.058	 1.102	 1.174	 1.272	 1.438	 1.604	 1.928	 2.194	 2.351	  2.46	 2.587	 2.544	 2.133	 1.571	 0.406	-0.456	 -0.82	  -0.82	  -0.82	  -0.82	  -0.82	 1.058	 1.694	 -0.82];
c11 = [   0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	  0.04	 0.077	  0.15	 0.253	   0.3	   0.3	    0.3	    0.3	    0.3	    0.3	  0.04	 0.092	   0.3];
c12 = [   0.61	  0.61	  0.61	  0.61	  0.61	  0.61	  0.61	  0.61	  0.61	  0.61	  0.61	 0.883	     1	     1	     1	     1	     1	      1	      1	      1	      1	  0.61	     1	     1];
k1 = [     865	   865	   908	  1054	  1086	  1032	   878	   748	   654	   587	   503	   457	   410	   400	   400	   400	   400	    400	    400	    400	    400	   865	   400	   400];
k2 = [  -1.186	-1.219	-1.273	-1.346	-1.471	-1.624	-1.931	-2.188	-2.381	-2.518	-2.657	-2.669	-2.401	-1.955	-1.025	-0.299	     0	      0	      0	      0	      0	-1.186	-1.955	     0];
k3 = [   1.839	  1.84	 1.841	 1.843	 1.845	 1.847	 1.852	 1.856	 1.861	 1.865	 1.874	 1.883	 1.906	 1.929	 1.974	 2.019	  2.11	    2.2	  2.291	  2.517	  2.744	 1.839	 1.929	 2.744];
slny = [ 0.478	  0.48	 0.489	  0.51	  0.52	 0.531	 0.532	 0.534	 0.534	 0.544	 0.541	  0.55	 0.568	 0.568	 0.564	 0.571	 0.558	  0.576	  0.601	  0.628	  0.667	 0.478	 0.484	 0.667];
tlny = [ 0.219	 0.219	 0.235	 0.258	 0.292	 0.286	  0.28	 0.249	  0.24	 0.215	 0.217	 0.214	 0.227	 0.255	 0.296	 0.296	 0.326	  0.297	  0.359	  0.428	  0.485	 0.219	 0.203	 0.485];
sigmac =[0.166	 0.166	 0.165	 0.162	 0.158	  0.17	  0.18	 0.186	 0.191	 0.198	 0.206	 0.208	 0.221	 0.225	 0.222	 0.226	 0.229	  0.237	  0.237	  0.271	   0.29	 0.166	  0.19	  0.29];
rhos = [     1	 0.999	 0.989	 0.963	 0.922	 0.898	  0.89	 0.871	 0.852	 0.831	 0.785	 0.735	 0.628	 0.534	 0.411	 0.331	 0.289	  0.261	    0.2	  0.174	  0.174	     1	 0.691	 0.174];
rhot = [     1	 0.994	 0.979	 0.927	  0.88	 0.871	 0.885	 0.913	 0.873	 0.848	 0.756	 0.631	 0.442	 0.29	  0.29	  0.29	  0.29	   0.29	   0.29	   0.29	   0.29	     1	 0.538	  0.29];

c=1.88;
n=1.18;

T=siteprop.period;
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
    i=find(period==T); 
    
    Rjb=siteprop.Rjb;    %coseismic rupture distance
    Ztor=faultprop.Ztor;    %depth to top of coseismic rupture surface
    lambda=faultprop.rake;  %rake angle in deg
    delta=faultprop.dip;    %dip angle in deg
    
          
    V30=siteprop.V30;       %avg shear wave velocity to 30m depth
    
    if siteprop.Zvs<0       %distance to 2.5km/s shear wave horizon
        Z10=exp(28.5-3.82/8*log(V30^8+378.7^8)); %CY08 estimate
        Zvs = 0.519 + 3.595 * (Z10/1000); %CB08 estimate of Zvs - note that the CY08 Z10 is in m not km
    else
        Zvs=siteprop.Zvs;  
    end
    
    
    % Magnitude dependence
    if M<=5.5
        fmag=c0(i)+c1(i)*M;
    else if M<=6.5
            fmag=c0(i)+c1(i)*M+c2(i)*(M-5.5); 
        else    
            fmag=c0(i)+c1(i)*M+c2(i)*(M-5.5)+c3(i)*(M-6.5); 
        end
    end 
    
    % Distance dependence
    fdis=(c4(i)+c5(i)*M)*log(sqrt(Rrup^2+c6(i)^2));

    % Style of faulting
    if Ztor<1
        ffltz=Ztor;
    else    
        ffltz=1;
    end 

    Frv=(lambda>30&lambda<150);
    Fnm=(lambda>-150&lambda<-30);

    fflt=c7(i)*Frv*ffltz+c8(i)*Fnm;

    % Hanging-wall effects
    if Rjb==0
        fhngr=1;
    else    
        if Ztor<1
            fhngr=(max(Rrup,sqrt(Rjb^2+1))-Rjb)/max(Rrup,sqrt(Rjb^2+1));
        else    
            fhngr=(Rrup-Rjb)/Rrup;
        end
    end 

    if M<=6
        fhngm=0;
    else
        if M<6.5
            fhngm=2*(M-6);
        else    
            fhngm=1;
        end
    end 

    fhngz=((20-Ztor)/20)*(Ztor>=0&Ztor<20);
    fhngdelta=(delta<=70)+((90-delta)/20)*(delta>70);
    fhng=c9(i)*fhngr*fhngm*fhngz*fhngdelta;

    % Site conditions
    if V30<k1(i)
        siteprop.V30=1100; siteprop.period=0;
        A1100=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);
        siteprop.V30=V30; siteprop.period=T;
        fsite=c10(i)*log(V30/k1(i))+k2(i)*(log(A1100+c*(V30/k1(i))^n)-log(A1100+c));
    else    
        fsite=(c10(i)+k2(i)*n)*log(V30/k1(i));
    end 

    % Sediment effects
    if Zvs<1
        fsed=c11(i)*(Zvs-1);
    else if Zvs<=3
            fsed=0;
        else    
            fsed=c12(i)*k3(i)*exp(-0.75)*(1-exp(-0.25*(Zvs-3)));
        end
    end 
    
    % Median value
    SA=exp(fmag+fdis+fflt+fhng+fsite+fsed);
    
    % Standard deviation computations
    if (V30<k1(i))
        alpha1=k2(i)*A1100*((A1100+c*(V30/k1(i))^n)^(-1)-(A1100+c)^(-1));
    else
        alpha1=0;
    end
    
    slnyB=sqrt(slny(i)^2-0.3^2);
    slnyA=sqrt(slny(find(period==0))^2-0.3^2);
    sigma_sq=slny(i)^2+alpha1^2*slny(find(period==0))^2+2*alpha1*rhos(i)*slnyB*slnyA;
    tau_sq=tlny(i)^2;
    
    
    if length(regexp(siteprop.orientation,'random','match'))~=0  %random/arbitrary component
        sigma_SA(1)=sqrt(sigma_sq+tau_sq+sigmac(i)^2);
        sigma_SA(2)=sqrt(tau_sq);
        sigma_SA(3)=sqrt(sigma_sq+sigmac(i)^2);
    else
        sigma_SA(1)=sqrt(sigma_sq+tau_sq);
        sigma_SA(2)=sqrt(tau_sq);
        sigma_SA(3)=sqrt(sigma_sq);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%