
function [f,df,lam,bout,stat] = composite_model_loss_binom(bin,y,N,cov_construct,penalty,params)

X = cov_construct(bin);
X(~isfinite(X))=0;
% X = X(:,var(X)>0);

if nargin>5
    [bout,f,stat]=glmfit(X(params.mask,:),[y(params.mask) repmat(N,sum(params.mask),1)],'binomial','constant','off');
else
    [bout,f,stat]=glmfit(X,[y repmat(N,size(y))],'binomial','constant','off');
end
if isIllConditioned(decomposition(X))
%     keyboard
    f=Inf;
    df=[];
end

if nargin>4
    f = f+penalty(bin);
end
lam = 1./(1+exp(-X*bout));
df = [];