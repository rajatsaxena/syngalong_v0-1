
% Get dataset from CRCNS SSC-3
% http://dx.doi.org/10.6080/K07D2S2F
load('../data/DataSet23.mat')
%% Pre-process spike times

Tlist=cell(0);
for i=1:length(data.spikes)
    Tlist{i}=data.spikes{i}'/1000;
	% add sub-resolution noise (0.0001s) to prevent aliasing in correlograms
    Tlist{i}=Tlist{i}+(rand(size(Tlist{i}))-.5)*0.0001;
end
Tlist=Tlist';

% exclude neurons with low firing rates
data.ids = cellfun(@numel,Tlist)>1000;
Tlist=Tlist(cellfun(@numel,Tlist)>1000);


%% get correlograms and identify candidates

params = [];

% correlogram parameters
params.trange = [-0.025 0.025];
params.tn = 102;
params.t = linspace(params.trange(1),params.trange(2),params.tn);
params.t = params.t(1:end-1)+mean(diff(params.t))/2;

% detection parameters
params.jit_min_t = 0.0011;
params.ssi_min_t = 0.001;
params.ssi_max = 0.4;
params.alph = 0.1;
params.jit_conv_window = 2;

result = detection_pipeline(Tlist,params);

%% Example putative synapses

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


%% Fig1A

data.x_subset=data.x(data.ids);
data.y_subset=data.y(data.ids);

figure(2)
clf
scatter(data.x_subset,data.y_subset,50,'filled')
axis equal
hold on
for i=1:length(Tlist)
    for j=1:length(Tlist)
        if result.detected_cnx(i,j)
            line([data.x_subset(i) data.x_subset(j)],[data.y_subset(i) data.y_subset(j)]);
        end
    end
end
cmap=lines(5); cmap=cmap(2:end,:);
for i=1:size(ijlist,1)
    line([data.x_subset(ijlist(i,1)) data.x_subset(ijlist(i,2))],[data.y_subset(ijlist(i,1)) data.y_subset(ijlist(i,2))],'Color',cmap(ijlist(i,3),:),'LineWidth',2);
end

% optionally, add labels
% for i=1:length(Tlist)
%     text(data.x_subset(i),data.y_subset(i),num2str(i))
% end
% hold off


%% Fig 1B/C

plotPairs(result,params,ijlist,1)

