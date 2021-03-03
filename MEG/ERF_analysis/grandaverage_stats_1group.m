clear all
%run /users/david/pathway_fieldtrip

%%
% This script takes timelocked data for each participant, creates a grand
% average, performs statistical analyses and makes some plots. 

% non-musicians:
selection = [1 2 4 5 6 7 8 9 10 11 12 14 15 16 17 18 19 20 21 25];  

% classical musicians:
% selection = [23 24 26 27 28 29 30 31];

% contemporary musicians:
% selection = [3 13 22];

%% set seed for reproducibility of permutations 

rng(157483,'twister')
seed = rng;
rng(seed)

%% Set up relevant directories

out_dir_figures = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/figures/';
if ~exist(out_dir_figures,'dir')
    mkdir(out_dir_figures)
end

out_dir_data = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/data/';
if ~exist(out_dir_data,'dir')
    mkdir(out_dir_data)
end

out_dir_stats = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/stats/';
if ~exist(out_dir_stats,'dir')
    mkdir(out_dir_stats)
end

%% Setup relevant variables

starting_point = 1; % where to start the storage of participants in the 
% relevant cell array.

blocks = {'tonal','atonal'}; % high entropy, low entropy
features = {'standard','pitch', 'intensity', 'timbre', 'location'}; % pitch, intensity, timbre
plot_interact = 'no'; % should we be able to interact with plots?
input_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/timelocked/data/'; 
out_prefix = ''; % prefix to save output files
%expertise = {'nonmus','classical', 'contemporary'}; % nonmusicians, classical, contemporary (short: no, cl, co)
scodes.nonmus = [1:20]; % subject codes for nonmusicians
%scodes.classical = [14 15 16]; % subject codes for musicians
%scodes.contemporary = [17 18 19];% numbers refer to the position of the selection; see above


%% load data

grand_all = []; % initialize structure for stats

for a = 1:length(selection)
    subject = selection(a);
    fprintf('loading subject %d \n',subject)
    load([input_dir,sprintf('%02.f',subject)])
    index = a + starting_point - 1;
    for hh = 1:length(features)
        feature = features{hh};
        for jj = 1:length(blocks)
            bl = blocks{jj};
            if isfield(mumufes.(bl).(feature).cfg, 'previous')
                mumufes.(bl).(feature).cfg.previous = [];
            end
            grand_all.(bl).(feature){index} = mumufes.(bl).(feature);
                
            if hh ~= 1
                grand_all.(bl).([feature '_MMN']){index}  = mumufes.(bl).([feature '_MMN']);
               
            end
        end
        if hh ~= 1
        grand_all.cond_diff.(feature){index} = grand_all.tonal.([feature '_MMN']){index};
       
        % Compute difference in MMN between tonal and atonal for stats and
       
        % plotting:
        
        grand_all.cond_diff.(feature){index}.avg = grand_all.tonal.([feature '_MMN']){index}.avg...
                                                  - grand_all.atonal.([feature '_MMN']){index}.avg;
        grand_all.cond_avg.(feature){index} = grand_all.tonal.([feature '_MMN']){index};
        grand_all.cond_avg.(feature){index}.avg = mean(cat(3,grand_all.tonal.([feature '_MMN']){index}.avg,...
                                                             grand_all.atonal.([feature '_MMN']){index}.avg),3);
        end
    end
end

%% GRAND AVERAGE

disp('grand averaging')

cfg = [];
cfg.channel = 'MEG';
cfg.method = 'across';
cfg.latency = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no'; 

blocks2 = fieldnames(grand_all);

grand_avg = [];

for jj = 1:length(blocks2)
    features2 = fieldnames(grand_all.(blocks2{jj}));
    for hh = 1:length(features2)
        grand_avg.(blocks2{jj}).(features2{hh}).all = ft_timelockgrandaverage(cfg,grand_all.(blocks2{jj}).(features2{hh}){1:end});      
    end
end


% % calculating the planar gradient
% 
% cfg                 = [];
% cfg.feedback        = 'yes';
% cfg.method          = 'template';
% cfg.neighbours = ft_prepare_neighbours(cfg,grand_all.atonal.pitch_MMN{1}); % 
% 
% cfg.planarmethod    = 'sincos';
% pitch_atonal_planar  = ft_megplanar(cfg, grand_all.atonal.pitch_MMN{1});
% 

cfg = [];
cfg.fontsize = 6;
cfg.xlim = [-0.1 0.4];
cfg.ylim = [-8e-14 4e-14]
cfg.graphcolor = 'brkm';
cfg.showlabels = 'no';
cfg.comment = 'no';
cfg.channel = 'MEG';
%cfg.ylim = [-2e-12 9e-12];
cfg.layout = 'CTF275.lay';

 % comparing conditions (=cond) in the three musician groups

figure;

ft_multiplotER(cfg, grand_avg.atonal.pitch_MMN.all, grand_avg.tonal.pitch_MMN.all)

% ft_multiplotER(cfg, grand_avg.atonal.intensity_MMN.all, grand_avg.tonal.intensity_MMN.all)
% ft_multiplotER(cfg, grand_avg.atonal.timbre_MMN.all, grand_avg.tonal.timbre_MMN.all)
% ft_multiplotER(cfg, grand_avg.atonal.location_MMN.all, grand_avg.tonal.location_MMN.all)

savefig(sprintf('%s/pitch_MMN_nonmus_only',out_dir_figures))


%% saving grand average

save([out_dir_data out_prefix 'grand_average'],'grand_avg');
%clear grand_all


%% Statistics (basically for reporting statistics; not for plotting)

% Define neighbour channels
cfg_neighb = [];
cfg_neighb.method = 'distance';
neighbours = ft_prepare_neighbours(cfg_neighb,grand_all.atonal.standard{1}); % 

disp('doing statistics')

% Configure statistical test (cluster based permutations)

cfg = [];
cfg.channel     = {'MEG'}; % I.e. only combined planar gradiometers
cfg.latency     = [0.075 0.3]; % time window
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo'; % sampling / permutation method
cfg.correctm    = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;
cfg.tail = 0;
cfg.clustertail = 0; %Two-sided T
cfg.alpha       = 0.025;
cfg.numrandomization = 10000;% Ideally 10.000

cfg.statistic   = 'depsamplesT';

cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number (only relevant for dependent samples)


%% 1. simple effects of deviant for all three groups separately - is there an MMN at all?: standard -deviant

results = [];
Nsub = 20;
comparison_name = 'MMN_check';

cfg.design = [];
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

for aa = 2:length(features)
    feature = features{aa};
    for bb = 1:length(blocks)
        block = blocks{bb};
        results.(comparison_name).([block feature]).all= ...
            ft_timelockstatistics(cfg, grand_all.(block).(feature){:},...
            grand_all.(block).standard{:});
        results.(comparison_name).([block feature]).all.cfg.previous = [];
    end
end

%% 2. t-test effect of entropy (LE vs. HE) (cond)

comparison_name = 'tonal_vs_atonal';

%%% main effect

Nsub = 20;

cfg.design = [];
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

for aa = 2:length(features)
    feature = features{aa};
    results.(comparison_name).(feature).all= ft_timelockstatistics(cfg,...
        grand_all.tonal.([feature '_MMN']){:},...
        grand_all.atonal.([feature '_MMN']){:});
    results.(comparison_name).(feature).all.cfg.previous = [];
end


save([out_dir_stats out_prefix 'cb_permutations'],'results','-v7.3')


%% topoplots

cfg = [];
cfg.parameter = 'avg';
% cfg.colorbar = 'SouthOutside';
cfg.colormap = parula;
cfg.style = 'straight';
cfg.comment = 'no';
cfg.marker = 'off';
cfg.gridscale = 268;

features = {'pitch', 'intensity', 'timbre', 'location'};
TW = {[70,300],[70,300],[70,300],[70,300]}; % in ms

sens = {'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'};

labels = grand_avg.tonal.pitch_MMN.all.label;
[ch, channs] = ismember(sens, labels);

% in one page:

lim1 = {[-10.1*10^(-14) 7.1*10^(-14)]/10^(-15),[0 15.1*10^(-13)]/10^(-13)};
positions = {[0.12,0.72,0.1,0.1];
    [0.32,0.72,0.1,0.1];
    [0.52,0.72,0.1,0.1];
    [0.72,0.72,0.1,0.1];
    [0.12,0.50,0.1,0.1];
    [0.32,0.50,0.1,0.1];
    [0.52,0.50,0.1,0.1];
    [0.72,0.50,0.1,0.1]};

close all

figure('Color','white');
gg = gcf;
set(gg, 'PaperOrientation', 'landscape')
set(gg, 'PaperUnits', 'normalized')
set(gg, 'PaperPosition', [0 0 1 1])
set(gg, 'Renderer','painters')
set(gg, 'Units', 'normalized')
set(gg, 'Position', [0 0 0.8 0.6])
n = 0;

conditions = {'MMN'};
blocks = {'tonal','atonal'};
block_names = {'tonal','atonal'};

zlim = {[-2.5*10^(-14) 2.5*10^(-14)], [-3.5*10^(-14) 3.5*10^(-14)], ...
    [-2.5*10^(-14) 2.5*10^(-14)], [-3.5*10^(-14) 3.5*10^(-14)]};

cfg.layout = 'CTF275.lay';

for bb = 1:length(blocks)
    block = blocks{bb};
    for hh = 1:length(features)
        feature = features{hh};
        
        cfg.zlim = zlim{hh}/10^(-15);
        
        tw = TW{hh};     
        n = n +1; %plot No
        val = []; idx =[];
        dataset = grand_avg.(block).([feature '_MMN']).all;
        tw_idx = find(dataset.time >= tw(1)/1000 & dataset.time <= tw(2)/1000);
        current_data = mean(dataset.avg(channs,:));
        [val, idx] = min(current_data(tw_idx));
        dataset.avg = dataset.avg/10^(-15);            
        cfg.xlim = [(dataset.time(idx-7 + tw_idx(1)-1)),(dataset.time(idx-7 + tw_idx(1)-1))];
        axes('Position',positions{n} + [0,0,0.19,0.19]) % size of the brains; the last two numbers
        cfg.style = 'straight_imsat';
        
        ft_topoplotER(cfg, dataset);
        
        title(['\fontsize{14}' num2str(round(dataset.time(idx+tw_idx(1)-1)*1000)) ' ms'])
        
        if n > 4
            c = colorbar('south');
            c.TickDirection = 'out';
            c.Ticks = [zlim{hh}(1)/10^(-15),0,zlim{hh}(2)/10^(-15)];
            c.Position = [positions{n}(1)+0.08, positions{n}(2)-0.01, 0.1, 0.02];
            annotation('textbox',[positions{n}(1)+0.18, positions{n}(2)-0.034, 0.05 .05],...
                'String','fT','EdgeColor','none','FontSize',12)
        end
        
    end
end

annotation('textbox',[.09 .82 .05 .05],'String','tonal','EdgeColor','none','FontSize',18)
annotation('textbox',[.08 .60 .05 .05],'String','atonal','EdgeColor','none','FontSize',18)
annotation('textbox',[.23 .95 .05 .05],'String','Pitch','EdgeColor','none','FontSize',18)
annotation('textbox',[.42 .95 .2 .05],'String','Intensity','EdgeColor','none','FontSize',18)
annotation('textbox',[.62 .95 .05 .05],'String','Timbre','EdgeColor','none','FontSize',18)
annotation('textbox',[.82 .95 .05 .05],'String','Location','EdgeColor','none','FontSize',18)

print([out_dir_figures 'Fig4-topo_plots'],'-dpdf','-r300');

%close all
