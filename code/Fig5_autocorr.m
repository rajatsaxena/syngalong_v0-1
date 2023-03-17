
%%

load('../data/abi_example_units.mat')

% select pair from examples
pair = 1;
ij = [pair*2-1 pair*2];

params.interval_centers = [-.02 .02];
params.num_bins = 201;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 2;
params.asym=true;

figure(2)
clf
fitout = composite_model_fit(ij,Tlist,params,'i>j');


%%

figure(3)
clf
isis1 = [NaN; diff(Tlist{ij(1)})];
isis2 = [diff(Tlist{ij(1)}); NaN];
bins = double(isis1>0.004 & isis2>0.004)+1;

generalizeToSplit(ij,Tlist,params,fitout,'i>j',bins,true)
