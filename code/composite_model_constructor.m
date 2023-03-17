
function [mdl,atoms] = composite_model_constructor(model_name,params,a1,a2)

tau_penalty = 100;
mdl.tau_penalty = tau_penalty;

at = linspace(2*params.interval(1)+params.dt/2,2*params.interval(2)-params.dt/2,2*(params.num_bins-1)+1)';
asym=(at<0);
asym2=(at>0);

% Basis for slow fluctuation
x = linspace(0,1,params.num_bins);
nknots = params.nknots;
if nknots>=1
    weights = ones(nknots+1,1);
    knots = linspace(-2/nknots,1+2/nknots,nknots+5);
    s = fastBSpline(knots,weights);
    X = s.getBasis(x);
    X = orth([ones(size(X,1),1) X]);
elseif nknots==-1
    X = params.X;
else
    X = [ones(length(x),1)];
end

tp=params.t;
tn=-params.t;
if ~isfield(params,'syn_type')
    params.syn_type='alpha';
end
if ~isfield(params,'asym')
    params.asym=false;
end
if ~isfield(params,'fixed_center')
    params.fixed_center=false;
end

switch lower(params.syn_type)
    case 'double_exp'
        atoms.syn1 = @(p) (tp>p(1)).*(-exp(-(tp-p(1))/p(2))+exp(-(tp-p(1))/p(3)));
        atoms.syn2 = @(p) (tn>p(1)).*(-exp(-(tn-p(1))/p(2))+exp(-(tn-p(1))/p(3)));
        atoms.syn_params = 3;
    case 'gamma'
        if ~isfield(params,'gammak')
            params.gammak=2;
        end
        atoms.syn1 = @(p) (tp>p(1)).*gampdf(tp-p(1),params.gammak,p(2))/gampdf((params.gammak-1)*p(2),params.gammak,p(2));
        atoms.syn2 = @(p) (tn>p(1)).*gampdf(tn-p(1),params.gammak,p(2))/gampdf((params.gammak-1)*p(2),params.gammak,p(2));
        atoms.syn_params = 2;
    otherwise % alpha function
        atoms.syn1 = @(p) (tp>p(1)).*((tp-p(1))/p(2).*exp(1-(tp-p(1))/p(2)));
        atoms.syn2 = @(p) (tn>p(1)).*((tn-p(1))/p(2).*exp(1-(tn-p(1))/p(2)));
        atoms.syn_params = 2;
end

if params.fixed_center
    atoms.central_gauss = @(p) exp(-tp.^2./p(1).^2);
    atoms.central_params = 1;
    atoms.central_canbeneg = [];
else
    atoms.central_gauss = @(p) exp(-(tp-p(1)).^2./p(2).^2);
    atoms.central_params = 2;
    atoms.central_canbeneg = [1];
end

switch lower(model_name)
    case 'null'
        mdl.name = 'null';
        mdl.cov_construct = @(p) [ones(length(x),1)];
    case 'smooth'
        mdl.name = 'smooth';
        mdl.cov_construct = @(p) X;
    case 'i>j'
        mdl.name = 'i>j';
        if params.asym
            mdl.cov_construct = @(p) [X conv(atoms.syn1(exp(p)),a1/max(a1),'same')' conv(atoms.syn1(exp(p)),a1/max(a1).*asym,'same')'];
        else
            mdl.cov_construct = @(p) [X conv(atoms.syn1(exp(p)),a1/max(a1),'same')'];
        end
        mdl.num_free_params = atoms.syn_params;
        mdl.penalty = @(p) tau_penalty*exp(p(2));
        mdl.param_transform = @(p) exp(p);
    case 'j>i'
        mdl.name = 'j>i';
        if params.asym
            mdl.cov_construct = @(p) [X conv(atoms.syn2(exp(p)),a2/max(a2),'same')' conv(atoms.syn2(exp(p)),a2/max(a2).*asym2,'same')'];
        else
            mdl.cov_construct = @(p) [X conv(atoms.syn2(exp(p)),a2/max(a2),'same')'];
        end
        mdl.num_free_params = atoms.syn_params;
        mdl.penalty = @(p) tau_penalty*exp(p(2));
    case 'recip'
        mdl.name= 'i>j and j>i';
        if params.asym
            mdl.cov_construct = @(p) [X conv(atoms.syn1(exp(p(1:2))),a1/max(a1),'same')' conv(atoms.syn2(exp(p(3:4))),a2/max(a2),'same')' conv(atoms.syn1(exp(p)),a1/max(a1).*asym,'same')' conv(atoms.syn2(exp(p)),a2/max(a2).*asym2,'same')'];
        else
            mdl.cov_construct = @(p) [X conv(atoms.syn1(exp(p(1:2))),a1/max(a1),'same')' flipud(conv(fliplr(atoms.syn2(exp(p(3:4)))),a2/max(a2),'same')')];
        end
        mdl.num_free_params = 2*atoms.syn_params;
        mdl.penalty = @(p) tau_penalty*(exp(p(2))+exp(p(4)));
    otherwise
        error(['Unknown model_name ' model_name])
end