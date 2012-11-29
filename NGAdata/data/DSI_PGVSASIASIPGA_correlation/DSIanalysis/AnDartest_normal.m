function [AD2a,AD2acrit,P,alpha] = AnDartest_normal(x,alpha)
%ANDARTEST Anderson-Darling test for assessing normality of a sample data.
% The Anderson-Darling test (Anderson and Darling, 1952) is used to test if 
% a sample of data comes from a specific distribution. It is a modification
% of the Kolmogorov-Smirnov (K-S) test and gives more weight to the tails
% than the K-S test. The K-S test is distribution free in the sense that the
% critical values do not depend on the specific distribution being tested.
% The Anderson-Darling test makes use of the specific distribution in calculating
% critical values. This has the advantage of allowing a more sensitive test
% and the disadvantage that critical values must be calculated for each
% distribution.
% The Anderson-Darling test is only available for a few specific distributions.
% The test is calculated as: 
%              
%        AD2 = integral{[F_o(x)-F_t(x)]^2/[F_t(x)(1-F_t(x)0]}dF_t(x)
%
%        AD2a = AD2*a
%
% Note that for a given distribution, the Anderson-Darling statistic may be
% multiplied by a constant, a (which usually depends on the sample size, n).
% These constants are given in the various papers by Stephens (1974, 1977a,
% 1977b, 1979, 1986). This is what should be compared against the critical 
% values. Also, be aware that different constants (and therefore critical 
% values) have been published. You just need to be aware of what constant 
% was used for a given set of critical values (the needed constant is typically
% given with the critical values). 
% The critical values for the Anderson-Darling test are dependent on the 
% specific distribution that is being tested. Tabulated values and formulas
% have been published for a few specific distributions (normal, lognormal, 
% exponential, Weibull, logistic, extreme value type 1). The test is a one-sided
% test and the hypothesis that the distribution is of a specific form is 
% rejected if the test statistic, AD2a, is greater than the critical value. 
% Here we develop the m-file for detecting departure from normality. It is 
% one of the most powerful statistics for test this.
%
% Syntax: function AnDartest(X) 
%      
%     Input:
%          x - data vector
%      alpha - significance level (default = 0.05)
%
%     Output:
%       AD2a - anderson darling test statistic
%   AD2acrit - critical AD2a value for alpha confidence
%          P - probability of emperical distrbution following normal dist
%      alpha - alpha value (may have been modified internally, see below)

switch nargin
    case{2}
        if isempty(x) == false && isempty(alpha) == false
            if (alpha <= 0 || alpha >= 1)
                fprintf('Warning: Significance level error; must be 0 < alpha < 1 \n');
                return;
            end
        end
    case{1}
        alpha = 0.05;
    otherwise 
        error('Requires at least one input argument.');
end

n = length(x);

if n < 7,
    disp('Sample size must be greater than 7.');
    return,
else
    x = x(:);
    x = sort(x);
    fx = normcdf(x,mean(x),std(x));
    i = 1:n;
    
    S = sum((((2*i)-1)/n)*(log(fx)+log(1-fx(n+1-i))));
    AD2 = -n-S;
    
    AD2a = AD2*(1 + 0.75/n + 2.25/n^2); %correction factor for small sample sizes: case normal
    
    if (AD2a >= 0.00 && AD2a < 0.200);
        P = 1 - exp(-13.436 + 101.14*AD2a - 223.73*AD2a^2);
    elseif (AD2a >= 0.200 && AD2a < 0.340);
        P = 1 - exp(-8.318 + 42.796*AD2a - 59.938*AD2a^2);
    elseif (AD2a >= 0.340 && AD2a < 0.600);
        P = exp(0.9177 - 4.279*AD2a - 1.38*AD2a^2);
    else (AD2a >= 0.600 && AD2a <= 13);
        P = exp(1.2937 - 5.709*AD2a + 0.0186*AD2a^2);
    end
    
    siglevel=[0.2 0.1 0.05 0.025 0.01 0.005];
    aalpha = [0.5091 0.6305 0.7514 0.8728 1.0348 1.1578];
    b0=[-0.756 -0.75 -0.795 -0.881 -1.013 -1.063];
    b1=[-0.39 -0.8 -0.89 -0.94 -0.93 -1.34];
    
    if (length(find(abs((siglevel-alpha))<0.001))==0)
        i=min(period(find(period<T)));  %select lower period option
        alpha=siglevel(i);
    else
        i=find(abs((siglevel-alpha)) < 0.001); % Identify the sign level
    end
    
    calpha=aalpha(i)*(1+b0(i)/n+b1(i)/n^2);
    AD2acrit=calpha;
end

