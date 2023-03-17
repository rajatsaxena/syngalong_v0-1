% Fits a composite model of a correlogram of the form...
%   yhat = logistic(X*b1 + f(c)*b2)
% in the typical case X is a smooth basis and f(c) is a model of fast
% structure (e.g. putative synapses)
%
% Input...
%   ij                  [1x2] indices for the pair of neurons
%   Tlist               cell-array of event times - correlogram is based on Tlist{ij(1)} vs Tlist{ij(2)}
%   params              struct specifying how to compute correlogram and how to fit
%   model_name          {'i>j','j>i','recip','null'}
%   pre_spk_id          (optional) a logical array of which events in Tlist{ij(1)} to include
%   feat_params_in      (optional) inputs "c" to use directly instead of optimizing with f(c)
%
% Output...
%   fit                 struct with results of fit
%
% Notes...
%  model_name='i>j' does a basic convolutional synapse
%  the form can be set using params.syn_type={'alpha','gamma','double_exp'}
%  everything is fit using a binomial observation model
%  X is a cubic B-spline basis with params.nknots equally-spaced knots
%  when fitting feat_params we use nrestarts random restarts to avoid local minima
%  note that composite_model_constructor includes penalties on c by default

function fit = composite_model_fit(ij,Tlist,params,model_name,pre_spk_id,feat_params_in)

if nargin<5
    pre_spk_id = ones(size(Tlist{ij(1)}))>0;
end
if nargin<6
    feat_params_in=[];
end

% Cross correlogram
[c,~] = spk_xcorr(Tlist{ij(1)}(pre_spk_id),Tlist{ij(2)},params.interval(1),params.interval(2),params.num_bins);
N = sum(pre_spk_id);

% Auto correlograms (for convolutional models)
a1 = spk_xcorr(Tlist{ij(1)}(pre_spk_id),Tlist{ij(1)},2*params.interval(1)+params.dt/2,2*params.interval(2)-params.dt/2,2*(params.num_bins-1)+1);
a2 = spk_xcorr(Tlist{ij(2)},Tlist{ij(2)},2*params.interval(1)+params.dt/2,2*params.interval(2)-params.dt/2,2*(params.num_bins-1)+1);


% Get the model structure...
if params.nknots<2
    params.nknots=-1;
end
[mdl,atoms] = composite_model_constructor(model_name,params,a1,a2);


if ~isfield(mdl,'penalty') || isempty(mdl.penalty)
    mdl.penalty=@(p) 0;
end

options=[];
options.MaxIter=100;
options.MaxFunEvals=10000;
options.method='lbfgs';
options.progTol = 10e-6;
options.Display = 'off';
options.numDiff = 2;

if ~isempty(feat_params_in)
    % Use the nonlinear fit parameters that are give...
    fit.feat_params=feat_params_in;
elseif strcmp(lower(model_name),'null') || strcmp(lower(model_name),'smooth') 
    fit.feat_params=[];
else
    % Optimize the parameters for the nonlinear fit...
    fopt = Inf; ftrace=[]; eftrace=[];
    nrestarts=40;
    for restart=1:nrestarts
        p = log(gamrnd(2,0.002,mdl.num_free_params,1));
        [b,ftrace(restart)] =  minFunc(@composite_model_loss_binom,p,options,c,N,mdl.cov_construct,mdl.penalty,params);
        if ftrace(restart)<fopt && ftrace(restart)>0
            fit.feat_params=b;
            fopt=ftrace(restart);
        end
    end
end

[fit.f,~,fit.yhat,fit.glm_params,fit.stat]=composite_model_loss_binom(fit.feat_params,c,N,mdl.cov_construct,@(p)0,params);
mdl.X = mdl.cov_construct(fit.feat_params);
if strcmp(lower(model_name),'i>j')
    if ~params.asym
        fit.lamsmoo = 1./(1+exp(-mdl.X(:,1:(end-1))*fit.glm_params(1:(end-1))));
    else
        fit.lamsmoo = 1./(1+exp(-mdl.X(:,1:(end-2))*fit.glm_params(1:(end-2))));
    end
end

% Estimate efficacy
if strcmp(lower(model_name),'i>j')
    tid = atoms.syn1(exp(fit.feat_params))>.01;
    fit.eff=sum(fit.yhat(tid)-fit.lamsmoo(tid))*N/sum(pre_spk_id);
    
    bsamp = mvnrnd(fit.glm_params,fit.stat.covb,100);
    yhatsamp = 1./(1+exp(-mdl.X*bsamp'));
    smoosamp = 1./(1+exp(-mdl.X(:,1:(params.nknots+2))*bsamp(:,1:(params.nknots+2))'));
    effsamp = sum(yhatsamp(tid,:)-smoosamp(tid,:),1)*N/sum(pre_spk_id);
    fit.effse=std(effsamp);
end

% preserve model detail
fit.mdl=mdl;
fit.atoms=atoms;

% plot
bar(params.t,c/N,1,'EdgeColor','flat')
hold on
plot(params.t,fit.yhat,'r','LineWidth',3)
if strcmp(lower(model_name),'i>j')
    plot(params.t,fit.lamsmoo)
end
hold off
title(mdl.name)
drawnow