
load('/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/timelocked/01.mat')

% Plot data

cfg = [];
cfg.fontsize = 6;
cfg.xlim = [-0.1 0.4];
cfg.graphcolor = 'brk';
cfg.showlabels = 'no';
cfg.comment = 'no';
cfg.channel = 'MEG';

% Grad

%cfg.ylim = [-2e-12 9e-12];
cfg.layout = 'CTF275.lay';
figure;
ft_multiplotER(cfg, mumufes.tonal.pitch)
