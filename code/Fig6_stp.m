
load('../data/abi_example_units.mat')

% select pair from examples
pair = 5;
ij = [pair*2-1 pair*2];

params=[];
% correlogram parameters
params.interval_centers = [-.02 .02];
params.num_bins = 101;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

% fit parameters
params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.0001)=0;
params.mask = params.mask>0;
params.nknots = 4;
params.asym=false;

figure(1); clf
[C,bins,n] = split_xcorr(ij,Tlist,params,'isi',5,true);

%% 
[C,bins,n] = split_xcorr(ij,Tlist,params,'isi',20,true);

fitout = composite_model_fit(ij,Tlist,params,'i>j');
fitout.eff
clf
isis = [NaN; diff(Tlist{ij(1)})];
xdat=[];ydat=[];
for i=1:size(C,1)
    pre_spk_id = bins==i;
    fitout = composite_model_fit(ij,Tlist,params,'i>j',pre_spk_id,fitout.feat_params);
    xdat(i,1) = mean(isis(pre_spk_id));
    ydat(i,1) = fitout.eff;
    ydat(i,2) = fitout.effse;
end
%%

% [A U tau_syn tau_fac tau_dep f]
sigm = @(p) 1./(1+exp(-p));
syn = @(d,x) sigm(d(1))*exp(-x/exp(d(2)));
fac = @(d,x) sigm(d(1))*exp(-x/exp(d(2)))+sigm(d(1))*(1-sigm(d(1))*exp(-x/exp(d(2))));
dep = @(d,x) 1-sigm(d(1))*exp(-x/exp(d(2)));
full_mod = @(p,x) p(1)*( syn(p(2:3),x) + fac(p([2 4]),x).*dep(p([2 5]),x) );
[beta,R,J] = nlinfit(xdat(:,1),ydat(:,1),full_mod,[1 0 log(0.01) log(20) log(0.5)]);

% Basic depression model (without facilitation)
% [A U tau_syn tau_dep]
% full_mod = @(p,x) p(1)*( syn(p(2:3),x) + dep(p([2 4]),x) );
% [beta,R,J] = nlinfit(xdat(:,1),ydat(:,1),full_mod,[1 0 log(0.01) log(0.5)]);

% Random restarts
minmse=Inf;
A0_pair = [NaN NaN NaN 1 1 -1];
for i=1:500
    p0 = [A0_pair(pair) 0.1 log(0.001) log(0.5) log(0.5)]+randn(1,5); % for Fig 5a
    [betap,Rp,Jp,Cp,mse] = nlinfit(xdat(:,1),ydat(:,1),full_mod,p0);
    if mse<minmse
        beta=betap;
        R=Rp;
        J=Jp;
        Covb=Cp;
    end
end

figure(2)
errorbar(xdat,ydat(:,1),ydat(:,2),'o','MarkerFaceColor','auto','CapSize',0)
set(gca,'XScale','log')
xl0=log10(xlim());
x0 = logspace(xl0(1),xl0(2),256);
[ypred,delta] = nlpredci(full_mod,x0,beta,R,'Jacobian',J);
hold on
plot(x0,ypred,x0,ypred+delta,'k',x0,ypred-delta,'k')
hold off
box off; set(gca,'TickDir','out')
xlabel('ISI [s]')
ylabel('Efficacy')
