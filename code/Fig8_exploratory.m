flist = dir('D:\allen_neuropixel\*_result.mat');
prefixes = arrayfun(@(x) x.name(1:9),flist,'UniformOutput',false);

dat_dir = 'D:\allen_neuropixel\';
% prefix = '715093703';
prefix = prefixes{23};
load([dat_dir prefix '_result.mat'])
tab = readtable([dat_dir prefix '_units.csv']);
tab_sub = tab(data.ids,:);

X = log10([tab_sub{:,'waveform_duration'} tab_sub{:,'waveform_halfwidth'} tab_sub{:,'firing_rate'}]);
[type_est,~,posterior] = cluster(GMModel,X);
type_est = type_est-1;
type_est(type_est<0)=3;

%%
sidx = find(result.detected_cnx(:));
% 
% [r,c] = find(result.detected_cnx);
% ijlist = [r c];
% ijlist = ijlist(ijlist(:,1)==167,:);
% sidx = (ijlist(:,2)-1)*size(result.detected_cnx,1)+ijlist(:,1);

figure(1)
offset=1;
plotPairs(result,params,sidx,offset)


%% find recip examples

% sidx = find(result.detected_cnx(:));
sidx = find(result.detected_cnx'&result.detected_cnx&triu(ones(size(result.detected_cnx,1)),1));
[r,c] = find(result.detected_cnx'&result.detected_cnx&triu(ones(size(result.detected_cnx,1)),1));
ijlist = [r c];

figure(1)
offset=1;
plotPairs(result,params,sidx,offset)

%% find divergence

m = sum(result.detected_cnx');
[r,c] = find(result.detected_cnx);
ijlist = [r c m(r)'];
ijlist = ijlist(ijlist(:,3)>1,:);
[~,sidx] = sort(ijlist(:,3)*size(result.detected_cnx,1)+ijlist(:,1),'descend');
ijlist = ijlist(sidx,:);
sidx = (ijlist(:,2)-1)*size(result.detected_cnx,1)+ijlist(:,1);

figure(1)
offset=1;
plotPairs(result,params,sidx,offset)


%%
ijlist = ijlist(type_est(ijlist(:,1))==2,:);

subnets = cell(0);
upre = unique(ijlist(:,1));
for i=1:length(upre)
    subnets{i} = ijlist(ijlist(:,1)==upre(i),1:2);
    subnets{i} = [subnets{i}(1,1); subnets{i}(:,2)];
end
[~,sidx] = sort(cellfun(@length,subnets),'descend');
subnets = subnets(sidx);

%% find chain

[r,c] = find(result.detected_cnx);
ijlist_all = [r c];
chains = getChains(ijlist_all);
[~,sidx]=sort(cellfun(@length,chains),'descend');
chains = chains(sidx);
deselect = zeros(length(chains),1);
for i=1:length(chains)
    for j=1:(i-1)
        if chains{i} == chains{j}(1:length(chains{i}))
            deselect(i)=1;
        end
        if length(intersect(chains{i},chains{j}))>4
            deselect(i)=1;
        end
    end
end
chains = chains(~deselect);
chains = chains(cellfun(@(x) type_est(x(1))==2,chains));



%%
figure(2)
snum = 5;
subnet = subnets{snum};
% subnet = chains{19};
% subnet = [167 168 171 169];

plotSubnet(params,subnet,Tlist,type_est(subnet))

%% append to examples

subnet_spikes{length(subnet_spikes)+1} = Tlist(subnet);
subnet_ecephys_id{length(subnet_spikes)} = prefix;
subnet_unit_ids{length(subnet_spikes)} = tab_sub{subnet,'unit_id'};

%%


load('../data/abi_subnet_examples.mat')
load('../data/abi_waveform_gmm.mat')

%% requires ABI data

example = 7;

dat_dir = 'D:\allen_neuropixel\';
prefix = subnet_ecephys_id{example};
load([dat_dir prefix '_result.mat'])
tab = readtable([dat_dir prefix '_units.csv']);
tab_sub = tab(data.ids,:);

X = log10([tab_sub{:,'waveform_duration'} tab_sub{:,'waveform_halfwidth'} tab_sub{:,'firing_rate'}]);
[type_est,~,posterior] = cluster(GMModel,X);
type_est = type_est-1;
type_est(type_est<0)=3;

ij=[];
for i=1:length(subnet_unit_ids{example})
    ij(i) = find(tab_sub{:,'unit_id'}==subnet_unit_ids{example}(i));
end

tab_sub(ij,'ecephys_structure_acronym')
