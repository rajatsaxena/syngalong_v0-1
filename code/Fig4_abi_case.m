

dat_dir = 'D:\allen_neuropixel\';

prefix = '715093703'
load([dat_dir prefix '_spikes.mat'])
tab = readtable([dat_dir prefix '_units.csv']);
epoch = readtable([dat_dir prefix '_epoch.csv']);
running_speed = readtable([dat_dir prefix '_running_speed.csv']);

Tlist = cellfun(@(x) (x+(rand(size(x))-.5)*0.0001)',Tlist,'UniformOutput',false)';
data.ids = tab.snr>2 & cellfun(@length,Tlist)>1000;
Tlist=Tlist(data.ids);
tab_sub = tab(data.ids,:);

params=[];
params.trange = [-0.025 0.025];
params.tn = 102;
params.t = linspace(params.trange(1),params.trange(2),params.tn);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;
params.jit_min_t = 0.0011;
params.ssi_min_t = 0.00049;
params.ssi_max = 0.4;
params.alph = 0.1;
params.jit_conv_window = 2;
params.mask=abs(params.t)>0.0000001;

data.anatomy = tab.ecephys_structure_acronym;
data.xorig = tab.probe_horizontal_position;
data.y = tab.probe_vertical_position;
data.probe = tab.probe_id;
[~,~,data.probeu]=unique(data.probe);
data.x = data.xorig+data.probeu*500-35;
%%

result = detection_pipeline(Tlist,params);

%%

n = cellfun(@(x) length(x),Tlist);
p = n./max(cellfun(@max,Tlist))/100;

ep = p.*(1-p).*(norminv(0.975)-norminv(.2))./(sqrt(n.*p.*(1-p))-(p-.5)*norminv(.2));
en = p.*(1-p).*(norminv(0.025)-norminv(.8))./(sqrt(n.*p.*(1-p))-(p-.5)*norminv(.8));
histogram(ep)
hold on
histogram(en)
hold off


%%

sidx = find(result.detected_cnx(:));
sidx = sidx(kidx);
figure(1)

offset=101;
plotPairs(result,params,sidx,offset)

for i=0:19
    subplot(5,4,i+1)
    set(gca,'XTickLabel','')
    ylabel([num2str(tab_sub{ijlist(i+offset,1),'unit_id'}) '>' num2str(tab_sub{ijlist(i+offset,2),'unit_id'})])
    title([tab_sub{ijlist(i+offset,1),'ecephys_structure_acronym'}{1} '>' (tab_sub{ijlist(i+offset,2),'ecephys_structure_acronym'}{1})])
end


%%


[r,c] = find(result.detected_cnx);

params=[];
params.interval_centers = [-.02 .02];
params.num_bins = 201;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 6;
params.asym=false;

figure(2)
clf
clear fitout
for i=1:length(r)
    fprintf('%03i...\n',i)
    fitout(i) = composite_model_fit([r(i) c(i)],Tlist,params,'i>j');
    drawnow
end

%% screen

% pseudo-Rsquared
rsq=[];
for i=1:length(fitout)
    N = length(Tlist{r(i)});
    y = fitout(i).yhat(params.mask)*N+fitout(i).stat.resid;
    yhat = fitout(i).yhat(params.mask);
    llm = sum(y.*log(yhat)+(N-y).*log(1-yhat));
    ll0 = sum(y.*log(mean(y)/N)+(N-y).*log(1-mean(y)/N));
    lls = sum(y.*log(y./N)+(N-y).*log(1-y./N));
    rsq(i)=(llm-ll0)./(lls-ll0);
end

% latency and duration [ms] and efficacy
lat = exp(arrayfun(@(x) x.feat_params(1),fitout))*1000;
dur = exp(arrayfun(@(x) x.feat_params(2),fitout))*1000;
eff = arrayfun(@(x) x.eff,fitout);
ex = arrayfun(@(x) x.eff>0,fitout);
in = arrayfun(@(x) x.eff<0,fitout);

% select subset of candidate pairs
kidx = (rsq>.5) & (dur<2) & (lat<6) & (lat>.1) & (abs(eff)<.5);
ijlist = [r(kidx) c(kidx)];

%% assess stability across stimuli (Fig 4a)

tid = epoch.duration>300;
time_blocks = [epoch.start_time(tid) epoch.stop_time(tid)];
[u,~,ic] = unique(epoch{tid,'stimulus_name'});

figure(2)
clf

params=[];
params.interval_centers = [-.01 .01];
params.num_bins = 101;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 3;
params.asym=false;

% fig 4a
% example 1 -- [33 38] (ij=7, CA1>CA1)
% example 2 -- [469 645] (ij=77, CA3>DG)
% example 3 -- [400 411] (ij=40, VISp>VISp)
ij=40;
kidxf = find(kidx);

for j=1:length(u)
    pre_spk_id = Tlist{ijlist(ij,1)}*0;
    for k=1:length(ic)
        if ic(k)==j
            pre_spk_id = pre_spk_id+ Tlist{ijlist(ij,1)}>time_blocks(j,1) & Tlist{ijlist(ij,1)}<time_blocks(j,2);
        end
    end

    ax(j)=subplot(length(u)+2,1,j)
    fitsub = composite_model_fit(ijlist(ij,:),Tlist,params,'i>j',pre_spk_id,fitout(kidxf(ij)).feat_params);
    ylabel(strrep(u{j},'_',' '))
    box off; set(gca,'TickDir','out')

end


vel_at_spk = interp1(running_speed.start_time,running_speed.velocity,Tlist{ijlist(ij,1)},'linear');
pre_spk_id = vel_at_spk>5;
ax(length(u)+1) = subplot(length(u)+2,1,length(u)+1)
fitsub = composite_model_fit(ijlist(ij,:),Tlist,params,'i>j',pre_spk_id,fitout(kidxf(ij)).feat_params);
ylabel('Running')
box off; set(gca,'TickDir','out')
pre_spk_id = vel_at_spk<5;
ax(length(u)+2) = subplot(length(u)+2,1,length(u)+2)
fitsub = composite_model_fit(ijlist(ij,:),Tlist,params,'i>j',pre_spk_id,fitout(kidxf(ij)).feat_params);
ylabel('Resting')
box off; set(gca,'TickDir','out')
linkaxes(ax)

%%

params=[];
params.interval_centers = [-.02 .02];
params.num_bins = 101;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 6;
params.asym=false;

% fig 4d
% example 1 -- [371 365] (ij=30, CA1>CA1)
% example 2 -- [191 192] (ij=15, VISpm>VISpm)
% example 3 -- [413 417] (ij=42, VISp>VISp)

figure(7)
composite_model_fit(ijlist(42,:),Tlist,params,'i>j');
box off; set(gca,'TickDir','out')

%% Fig 4c

ijlist(ex(kidx),3)=1;
ijlist(in(kidx),3)=0;
figure(4)
subplot(1,2,1)
scatter(log10(tab_sub{ijlist(:,1),'waveform_halfwidth'}),log10(tab_sub{ijlist(:,1),'waveform_duration'}),20,ijlist(:,3),'filled')
xlabel('Waveform Half-width')
ylabel('Waveform Duration')
set(gca,'TickDir','out'); box off
xlim(log10([0.065 0.65]))
ylim(log10([0.08 1.25]))
subplot(1,2,2)
scatter(log10(tab_sub{ijlist(:,1),'firing_rate'}),log10(tab_sub{ijlist(:,1),'waveform_duration'}),20,ijlist(:,3),'filled')
xlabel('log-10 Firing Rate')
ylabel('Waveform Duration')
set(gca,'TickDir','out'); box off
xlim(log10([0.1 100]))
ylim(log10([0.08 1.25]))

% ijf=[7 77 40];
% ij = 1;
% subplot(1,2,1)
% hold on
% plot(log10(tab_sub{ijlist(ij,1),'waveform_halfwidth'}),log10(tab_sub{ijlist(ij,1),'waveform_duration'}),'r^')
% hold off
% subplot(1,2,2)
% hold on
% plot(log10(tab_sub{ijlist(ij,1),'firing_rate'}),log10(tab_sub{ijlist(ij,1),'waveform_duration'}),'r^')
% hold off

subplot(1,2,1)
hold on
text(log10(tab_sub{ijlist(:,1),'waveform_halfwidth'})+randn(size(ijlist(:,1)))/100,log10(tab_sub{ijlist(:,1),'waveform_duration'}),ijlist(:,1)*0,num2str([1:size(ijlist,1)]'))
hold off
subplot(1,2,2)
hold on
text(log10(tab_sub{ijlist(:,1),'firing_rate'})+randn(size(ijlist(:,1)))/100,log10(tab_sub{ijlist(:,1),'waveform_duration'}),ijlist(:,1)*0,num2str([1:size(ijlist,1)]'))
hold off
%%
figure(5)
subplot(3,1,1)
% edges = linspace(log10(0.065),log10(0.65),18);
% histogram(log10(tab_sub{ijlist(in(kidx),1),'waveform_halfwidth'}),edges,'EdgeColor','none')
% hold on
% histogram(log10(tab_sub{ijlist(ex(kidx),1),'waveform_halfwidth'}),edges,'EdgeColor','none')
% hold off
edges = linspace(log10(0.065),log10(0.65),256);
[fi,x]=ksdensity(log10(tab_sub{ijlist(in(kidx),1),'waveform_halfwidth'}),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(tab_sub{ijlist(ex(kidx),1),'waveform_halfwidth'}),edges,'bandwidth',range(edges)/15);
plot(x,fi,x,fe)
box off; set(gca,'TickDir','out')

subplot(3,1,2)
% edges = linspace(log10(0.1),log10(100),18);
% histogram(log10(tab_sub{ijlist(in(kidx),1),'firing_rate'}),edges,'EdgeColor','none')
% hold on
% histogram(log10(tab_sub{ijlist(ex(kidx),1),'firing_rate'}),edges,'EdgeColor','none')
% hold off
edges = linspace(log10(0.1),log10(100),256);
[fi,x]=ksdensity(log10(tab_sub{ijlist(in(kidx),1),'firing_rate'}),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(tab_sub{ijlist(ex(kidx),1),'firing_rate'}),edges,'bandwidth',range(edges)/15);
plot(x,fi,x,fe)
box off; set(gca,'TickDir','out')

subplot(3,1,3)
% edges = linspace(log10(0.08),log10(1.25),18);
% histogram(log10(tab_sub{ijlist(in(kidx),1),'waveform_duration'}),edges,'EdgeColor','none')
% hold on
% histogram(log10(tab_sub{ijlist(ex(kidx),1),'waveform_duration'}),edges,'EdgeColor','none')
% hold off
edges = linspace(log10(0.08),log10(1.25),256);
[fi,x]=ksdensity(log10(tab_sub{ijlist(in(kidx),1),'waveform_duration'}),edges,'bandwidth',range(edges)/15);
[fe,x]=ksdensity(log10(tab_sub{ijlist(ex(kidx),1),'waveform_duration'}),edges,'bandwidth',range(edges)/15);
plot(x,fi,x,fe)
box off; set(gca,'TickDir','out')

%%
same_probe = [tab_sub{ijlist(:,1),'probe_id'}==tab_sub{ijlist(:,2),'probe_id'}];
d_est = tab_sub{ijlist(:,1),'probe_vertical_position'}-tab_sub{ijlist(:,2),'probe_vertical_position'};
histogram(abs(d_est(same_probe)),100)

%%

ustructure = unique(tab_sub(reshape(ijlist(:,1:2),[],1),'ecephys_structure_acronym'));

cmat = zeros(size(ustructure,1));
for i=1:size(ustructure,1)
    for j=1:size(ustructure,1)
        prem = strcmp(tab_sub{ijlist(:,1),'ecephys_structure_acronym'},ustructure{i,'ecephys_structure_acronym'});
        posm = strcmp(tab_sub{ijlist(:,2),'ecephys_structure_acronym'},ustructure{j,'ecephys_structure_acronym'});
        cmat(i,j) = sum(prem & posm);
    end
end

figure(6)
imagesc(log10(cmat),'AlphaData',cmat>0);
% imagesc(cmat_e./cmat,'AlphaData',cmat>0);
box off; set(gca,'TickDir','out')
colorbar
set(gca,'XTick',1:size(cmat,1))
set(gca,'YTick',1:size(cmat,1))
set(gca,'XTickLabels',ustructure.ecephys_structure_acronym)
set(gca,'YTickLabels',ustructure.ecephys_structure_acronym)
