

nga_model = input( 'Please input NGA model (CB,BA,CY or AS):\n' );

outpth = '/Users/fengw/local/pylib/pynga/data';

clc

%get the observational data
data=observedIMintensities;
EQID=data(:,2);
%spectral periods to consider
T=[0.01 0.02 0.03 0.04 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1.0 1.5 2 3 4 5 7.5 10];

h_wait=waitbar(0,'Running..');
for i=1:length(data)
    
    waitbar(i/length(data));
    
    %get Metadata
    rec_num(i)=data(i,1);   EQ_num(i)=data(i,2);
    M=data(i,3);    faultprop.dip=data(i,4);  faultprop.AS=0;
    faultprop.rake=data(i,5);    mech=data(i,6);
    faultprop.Ztor=data(i,11);
    %use Wells and Coppersmith 'All' to get rup width given Mw
    faultprop.W=10^(-1.01+0.32*M);
    if mech==0
        faultprop.faultstyle='strikeslip';
    elseif mech==1
        faultprop.faultstyle='normal';
    elseif mech==2
        faultprop.faultstyle='reverse';
    else
        faultprop.faultstyle='other';
    end
    
    Rjb=data(i,14); 
    Rrup=data(i,15);     siteprop.Rjb=Rjb; siteprop.Rx=Rjb;
    siteprop.V30=data(i,16); siteprop.V30measured=0; siteprop.Z1pt0=-1; siteprop.Zvs=-1;
    siteprop.orientation='average'; siteprop.g=981;
    
    %Intensity measures
    meta(i,1) = data(i,1);
    meta(i,2) = EQID(i);    
  
    %     ------------PGA----------------
    siteprop.period=0;
    switch nga_model
        case 'BA'
            %get the predicted value of DSI using BA08
            [PGAmedian(i),PGAstd(i,1:3)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
        case 'CY'
            %get the predicted value using CY08
            [PGAmedian(i),PGAstd(i,1:3)]=ChiouYoungs_2008_nga(M,Rrup,siteprop,faultprop);

        case 'CB'
            %get the predicted value using CB08
            [PGAmedian(i),PGAstd(i,1:3)]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);

        case 'AS'
            %get the predicted value using AS08
            if i == 209
                disp([M Rjb siteprop.V30])
                [PGAmedian(i),PGAstd(i,1:3)]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop, 1);
                ['PGA ' num2str(PGAmedian(i))]
                ['PGA std ' num2str(PGAstd(i,1:3))]
                break
            end
                
    end
    PGA_observed(i)=data(i,18);
    PGAresidual_total(i) = (log(PGA_observed(i))-log(PGAmedian(i)))/PGAstd(i,1);  % normalized residual (total)

%     ------------SA----------------
    for j=1:length(T)
        siteprop.period=T(j);
        switch nga_model
            case 'BA'
                %get the predicted value of DSI using BA08
                [SAmedian(i,j),SAstd(i,j,1:3)]=BooreAtkinson_2007_nga(M,Rjb,siteprop,faultprop);
            case 'CY'
                %get the predicted value using CY08
                [SAmedian(i,j),SAstd(i,j,1:3)]=ChiouYoungs_2008_nga(M,Rrup,siteprop,faultprop);

            case 'CB'
                %get the predicted value using CB08
                [SAmedian(i,j),SAstd(i,j,1:3)]=CampbellBozorgina_2007_nga(M,Rrup,siteprop,faultprop);

            case 'AS'
                %get the predicted value using AS08
                [SAmedian(i,j),SAstd(i,j,1:3)]=AbrahamsonSilva_2008_nga(M,Rrup,siteprop,faultprop,0);
        end
        SA_observed(i,j)=data(i,19+j); 
        SAresidual_total(i,j) = (log(SA_observed(i,j))-log(SAmedian(i,j)))/SAstd(i,j,1);
    end

end

alleta=1;
[PGAeta, PGAepsilon] = get_interintraeventterms(PGAresidual_total, EQID, PGAstd(:,1)', PGAstd(:,2)', PGAstd(:,3)', alleta);
meta(1:length(data),3:9) = [log(PGAmedian)' PGAstd(:,1) PGAstd(:,2) PGAstd(:,3) PGAresidual_total' PGAeta' PGAepsilon'];    % all in log
Nsub = 7;
for j = 1:length(T)    
    [SAeta, SAepsilon] = get_interintraeventterms(SAresidual_total(:,j)', EQID, SAstd(:,j,1)', SAstd(:,j,2)', SAstd(:,j,3)', alleta);
    start = 3+7*j; 
    end0 =  9+7*j;
    meta(:,start:end0) = [log(SAmedian(:,j)) SAstd(:,j,1), SAstd(:,j,2), SAstd(:,j,3) SAresidual_total(:,j) SAeta' SAepsilon'];
    
    tmp = corrcoef( PGAresidual_total, SAresidual_total(:,j) );    
    rho01(1,j) = tmp(1,2);
    tmp = corrcoef( PGAeta, SAeta );    
    rho01(2,j) = tmp(1,2);
    tmp = corrcoef( PGAepsilon, SAepsilon );    
    rho01(3,j) = tmp(1,2);
end

% Save to file as txt
filepth = [ outpth '/meta_IMs_' nga_model '.txt' ];
dlmwrite( filepth, meta, 'delimiter', ' ' );

filepth1 = [outpth '/meta_Rhos_' nga_model '.txt' ];
dlmwrite( filepth1, rho01', 'delimiter', ' ' )
if 0
    for isub =1:3
        subplot( 3,1,isub) 
        semilogx( T, rho01(isub,:),'rx' )
        grid on;
        ylim([0,1.0])
    end
end

close(h_wait)

