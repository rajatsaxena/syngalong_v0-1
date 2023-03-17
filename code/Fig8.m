
load('../data/abi_subnet_examples.mat')

%% generate panels in Fig 8A-D using example subnetworks
% waveform classification is based on ../data/abi_waveform_gmm.mat

example = 2;

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

subnet = 1:length(subnet_unit_ids{example});
plotSubnet(params,subnet,subnet_spikes{example})