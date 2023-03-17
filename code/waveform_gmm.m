
flist = dir('D:\allen_neuropixel\*_result.mat');
prefixes = arrayfun(@(x) x.name(1:9),flist,'UniformOutput',false);

dat_dir = 'D:\allen_neuropixel\';

X=[];
for i=1:length(prefixes)
    prefix = prefixes{i};
    load([dat_dir prefix '_result.mat'],'data')
    tab = readtable([dat_dir prefix '_units.csv']);
    tab_sub = tab(data.ids,:);

    X = [X; log10([tab_sub{:,'waveform_duration'} tab_sub{:,'waveform_halfwidth'} tab_sub{:,'firing_rate'}])];
end

%%

figure(3)
plot(X(:,2),X(:,1),'o')

GMModel = fitgmdist(X,3);
[idx,~,posterior] = cluster(GMModel,X);

%%
subplot(1,2,1)
scatter(X(:,2),X(:,1),20,posterior(:,2),'filled','MarkerFaceAlpha',.1)
subplot(1,2,2)
scatter(X(:,3),X(:,1),20,posterior(:,2),'filled','MarkerFaceAlpha',.1)