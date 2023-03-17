ijlist = [280 264 1; % excitatory examples
    22 281 1;
    124 110 1;   
    140 162 2; % inhibitory examples
    207 196 2;
    270 267 2;
    165 102 3; % excitatory pre
    165 172 3;
    165 109 3;
    165 36 3;
    165 147 3;
    165 45 3;
    176 117 4; % inhibitory pre
    176 39 4;
    176 197 4;
    176 175 4;
    176 34 4;
    176 110 4];

% correlogram parameters
params.interval_centers = [-.02 .02];
params.num_bins = 201;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

% fit parameters
params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 6;

addpath(genpath(pwd))
Tsplit = linspace(0,60*60,7);

%%

paramsall=[];
clf
for i=1:size(ijlist,1)
    for j=1:length(Tsplit)-1
        pre_spk_id = Tlist{ijlist(i,1)}>Tsplit(j) & Tlist{ijlist(i,1)}<Tsplit(j+1);
        fitout = composite_model_fit(ijlist(i,:),Tlist,params,'i>j',pre_spk_id);
        paramsall(i,j,:) = [exp(fitout.feat_params') fitout.eff fitout.effse];
        drawnow
    end
end

%% Fig 3a

figure(21)
sym = {'o','^','s'};
cmap = lines(6);
cmap = cmap([4 5 2 1],:);
subplot(2,1,1)
for i=1:6
    x = squeeze(paramsall(i,:,1))*1000;
    y = squeeze(paramsall(i,:,2))*1000;
    plot(x,y,sym{mod(i-1,3)+1},'MarkerFaceColor',cmap(ijlist(i,3),:),'MarkerEdgeColor','w')
    hold on
    error_ellipse(cov(x,y),mean([x' y'])) 
end
set(gca,'TickDir','out'); box off
hold off
ylabel('Time Constant [ms]')
xlim([0.5 3])

subplot(2,1,2)
for i=1:6
    x = squeeze(paramsall(i,:,1))*1000;
    y = squeeze(paramsall(i,:,3));
    plot(x,y,sym{mod(i-1,3)+1},'MarkerFaceColor',cmap(ijlist(i,3),:),'MarkerEdgeColor','w')
    hold on
    error_ellipse(cov(x,y),mean([x' y'])) 
end
hold off
set(gca,'TickDir','out'); box off
xlabel('Latency [ms]')
ylabel('Efficacy')
xlim([0.5 3])

%% Fig 3b

params=[];
% correlogram parameters
params.interval_centers = [-.01 .01];
params.num_bins = 101;
params.dt = range(params.interval_centers)/(params.num_bins-1);
params.interval = params.interval_centers+[-1 1]*params.dt/2;
params.t = linspace(params.interval(1),params.interval(2),params.num_bins+1);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

% fit parameters
params.mask = ones(size(params.t));
params.mask(abs(params.t)<0.001)=0;
params.mask = params.mask>0;
params.nknots = 3;

% select pair
ij=3;

figure(10)
fitout = composite_model_fit(ijlist(ij,:),Tlist,params,'i>j');
for j=1:length(Tsplit)-1
    pre_spk_id = Tlist{ijlist(ij,1)}>Tsplit(j) & Tlist{ijlist(ij,1)}<Tsplit(j+1);
    fitoutsub(j) = composite_model_fit(ijlist(ij,:),Tlist,params,'i>j',pre_spk_id,fitout.feat_params);
    drawnow
end

figure(11)
[C,bins,n] = split_xcorr(ijlist(ij,:),Tlist,params,'time',length(Tsplit)-1,1);
for j=1:length(Tsplit)-1
    pre_spk_id = Tlist{ijlist(ij,1)}>Tsplit(j) & Tlist{ijlist(ij,1)}<Tsplit(j+1);
    subplot(2,length(Tsplit)-1,j)
    hold on
    plot(params.t,fitoutsub(j).yhat)
    hold off
end


%% fit all candidate pairs

[r,c] = find(result.detected_cnx); % result from Fig1_demo.m

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

clear fitout
for i=1:length(r)
    fprintf('%03i...\n',i)
    fitout(i) = composite_model_fit([r(i) c(i)],Tlist,params,'i>j');
    drawnow
end

%% get fit features

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

x = data.x(data.ids);
y = data.y(data.ids);
ddd0 = pdist([x' y']); % distances between all possible pairs
ddd = sqrt((x(r(c))-x(c)).^2 + (y(r)-y(c)).^2); % distances between candidate pairs
ex = arrayfun(@(x) x.eff>0,fitout);
in = arrayfun(@(x) x.eff<0,fitout);

% select subset of candidate pairs
kidx = rsq>.5 & dur<2 & lat<6 & lat>.1 & eff<.5;


%% Fig 3c

[dhe,e] = histcounts(ddd(ex & kidx),linspace(0,1500,20));
[dhi,e] = histcounts(ddd(in & kidx),linspace(0,1500,20));

subplot(3,1,1)
bar(e(1:end-1)+mean(diff(e))/2,dhe,1,'EdgeColor','flat')
hold on
bar(e(1:end-1)+mean(diff(e))/2,dhi,1,'EdgeColor','flat')
hold off
box off; set(gca,'TickDir','out')

[dhe,e] = histcounts(lat(ex & kidx),linspace(0,6,20));
[dhi,e] = histcounts(lat(in & kidx),linspace(0,6,20));

subplot(3,1,2)
bar(e(1:end-1)+mean(diff(e))/2,dhe,1,'EdgeColor','flat')
hold on
bar(e(1:end-1)+mean(diff(e))/2,dhi,1,'EdgeColor','flat')
hold off
box off; set(gca,'TickDir','out')

[dhe,e] = histcounts(dur(ex & kidx),linspace(0,2,20));
[dhi,e] = histcounts(dur(in & kidx),linspace(0,2,20));

subplot(3,1,3)
bar(e(1:end-1)+mean(diff(e))/2,dhe,1,'EdgeColor','flat')
hold on
bar(e(1:end-1)+mean(diff(e))/2,dhi,1,'EdgeColor','flat')
hold off
box off; set(gca,'TickDir','out')


%% Fig 3c-inset for connection probability

[dhe,e] = histcounts(ddd(ex & kidx),linspace(0,1000,30));
[dhi,e] = histcounts(ddd(in & kidx),linspace(0,1000,30));
[dh0] = histcounts(ddd0,linspace(0,1000,30));

figure(3)
clf
plot(e(1:end-1)+mean(diff(e))/2,dhe./dh0/2,'o','MarkerFaceColor','r')
dd = e(1:end-1)+mean(diff(e))/2;
[b,dev,stats] = glmfit(log2(dd),[dhe' 2*dh0'],'binomial');
[exp(b(2)) exp(b(2)-stats.se(2)*2) exp(b(2)+stats.se(2)*2)]
x0 = log2(linspace(10,1000,512));
yhat = glmval(b,x0,'logit',stats);
hold on
plot(2.^(x0),yhat)

plot(e(1:end-1)+mean(diff(e))/2,dhi./dh0/2,'o','MarkerFaceColor','b')
[b,dev,stats] = glmfit(log2(dd),[dhi' 2*dh0'],'binomial');
[exp(b(2)) exp(b(2)-stats.se(2)*2) exp(b(2)+stats.se(2)*2)]
yhat = glmval(b,x0,'logit',stats);
plot(2.^(x0),yhat)

hold off
set(gca,'YScale','log')
box off; set(gca,'TickDir','out')
ylim([0.0003 0.025])

