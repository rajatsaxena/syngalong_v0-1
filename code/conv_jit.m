function [m,s,p] = conv_jit(t,n,isgausswin)

if nargin<3 || ~isgausswin,
    % allow for non-integer n
    jit_filt = ones(n-(mod(n+1,2)),1);
    if (n-sum(jit_filt))>0
        jit_filt = [(n-sum(jit_filt))/2; jit_filt; (n-sum(jit_filt))/2];
    end
else
    jit_filt = normpdf(linspace(-length(t),length(t),2*length(t)+1),0,n);
    jit_filt = jit_filt/max(jit_filt);
end

tmpn = conv(t,jit_filt,'same');
basn = conv(ones(size(t)),jit_filt,'same');
s = sqrt(tmpn*1./basn.*(1-1./basn));
m=tmpn./basn;

if mean(t)>100 % use normal approximation
    p = normcdf(t,ceil(tmpn)./basn,sqrt(ceil(tmpn)./basn.*(1-1./basn)));
else
    p = min(binocdf(t,ceil(tmpn),1./basn),binocdf(t,ceil(tmpn),1./basn,'upper')+binopdf(t,ceil(tmpn),1./basn));
end
