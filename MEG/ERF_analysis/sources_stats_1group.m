
% SOURCES - STATISTICS for one group

clear all

%%
% This script takes timelocked sources for each participant, creates a grand
% average, performs statistical analyses and makes some plots. 

%selection = [1 2 4 6 7 8] % 9 10 11 12 14 15 16 18 19 20 21 25] ; % participantID included in analyses


selection = [1 2 4 6 7 8 9 10 11 12 14 15 16 18 19 20 21 25];

%% set seed for reproducibility of permutations 

rng(157483,'twister')
seed = rng;
rng(seed)

%% Set up relevant directories

% testfigures = subjects 1 to 6
% testfigures2 = subjects 1 to 12

out_dir_fig = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/stats/figures/';
if ~exist(out_dir_fig,'dir')
    mkdir(out_dir_fig)
end

out_dir_data = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/stats/data/';
if ~exist(out_dir_data,'dir')
    mkdir(out_dir_data)
end

out_dir_stats = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/stats/statistics/';
if ~exist(out_dir_stats,'dir')
    mkdir(out_dir_stats)
end

%% Setup relevant variables

starting_point = 1; % where to start the storage of participants in the 
% relevant cell array.

blocks = {'tonal','atonal'}; % high entropy, low entropy
features = {'pitch', 'intensity', 'timbre', 'location'}; % pitch, intensity, timbre
plot_interact = 'no'; % should we be able to interact with plots?
input_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/sources/'; 
out_prefix = ''; % prefix to save output files
%expertise = {'nonmus','classical', 'contemporary'}; % nonmusicians, classical, contemporary (short: no, cl, co)
scodes.nonmus = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]; % subject codes for nonmusicians
% scodes.classical = [4 5]; % subject codes for musicians
% scodes.contemporary = [6 7];% numbers refer to the position of the selection; see above


%% load data

sources_all = []; % initialize structure for stats
sources_norm = [];

for a = 1:length(selection)
    subject = selection(a);
    fprintf('loading subject %d \n',subject)
    load([input_dir,sprintf('%02.f_sources_test3',subject)])
    index = a + starting_point - 1;
    
    for hh = 1:length(features)
        feature = features{hh};
        
        for jj = 1:length(blocks)
            bl = blocks{jj};
            
            cnames = fieldnames(mumufes_sources.(bl).(feature));
            
            for c=1:length(cnames)
                cname = cnames{c};
                
                if a == 1
                    
                    sources_norm.(feature).(cname) = mumufes_sources.tonal.(feature).(cname).pow;
                    sources_norm.(feature).(cname)(:,1,2) = mumufes_sources.atonal.(feature).(cname).pow;
                else
                    
                    sources_norm.(feature).(cname)(:,end+1,1) = mumufes_sources.tonal.(feature).(cname).pow;
                    sources_norm.(feature).(cname)(:,end+1,2) = mumufes_sources.atonal.(feature).(cname).pow;
                    
                end
                
                
                if isfield(mumufes_sources.(bl).(feature).(cname).cfg, 'previous')
                    mumufes_sources.(bl).(feature).(cname).cfg.previous = [];
                end
                
                %mumufes_sources.(bl).(feature).(cname).avg.pow = mumufes_sources.(bl).(feature).(cname).pow;
                
                sources_all.(bl).(feature).(cname){index} = mumufes_sources.(bl).(feature).(cname);
                sources_all.(bl).(feature).(cname){index}.inside(:,:,:) = 1;
                
            end
        end
        
        % Compute difference in MMN between tonal and atonal for stats and plotting:
        
        sources_all.cond_diff.(feature).MMN{index} = sources_all.tonal.(feature).MMN{index};
        
        sources_all.cond_diff.(feature).MMN{index}.pow = sources_all.tonal.(feature).MMN{index}.pow...
            - sources_all.atonal.(feature).MMN{index}.pow;
        
        sources_all.cond_avg.(feature).MMN{index} = sources_all.tonal.(feature).MMN{index};
        sources_all.cond_avg.(feature).MMN{index}.pow = mean(cat(3,sources_all.tonal.(feature).MMN{index}.pow,...
            sources_all.atonal.(feature).MMN{index}.pow),3);
        
    end
end

%% normalize data
% 
% sources_stats_norm = sources_all;
% 
% for f = 1:length(features)
%     feature = features{f};
%     mu = mean(sources_norm.(feature),'all', 'omitnan');
%     sd = std(sources_norm.(feature),0,'all','omitnan');
%     c_names = fieldnames(sources_all);
%     for c = 1:length(c_names)
%         c_name = c_names{c};
%         c_data = sources_all.(c_name).(feature);
%         
%         for i = 1:length(c_data)
%             c_data{i}.pow = (c_data{i}.pow-mu)/sd;
%             
%         end
%         sources_stats_norm.(c_name).(feature) = c_data;
%         
%     end
% end


%% Statistics

disp('hello! doing statistics');

% Configure statistical test (cluster based permutations)

cfg = [];
cfg.parameter   = 'pow';
cfg.method      = 'montecarlo'; % sampling / permutation method
cfg.correctm    = 'cluster';
cfg.clusteralpha = 0.025;
cfg.clusterstatistic = 'maxsum';
cfg.tail = 0;
cfg.clustertail = 0; % two-sided T-test
cfg.alpha       = 0.025;
cfg.numrandomization = 1000;% 1000 is fine 

cfg.statistic   = 'depsamplesT';

cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number (only relevant for dependent samples)

%% 1. MMN sources

results_sources = [];

Nsub = 18;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

% dummy = sources_all.tonal.pitch{1};
% %dummy.pwr(:,:) = 0;
% dummy.pow(:,:) = 0;
% dummy_all = repmat({dummy},1,Nsub);

for aa = 1:length(features)
    feature = features{aa};
    for jj = 1:length(blocks)
        bl = blocks{jj};
        results_sources.MMN.(bl).(feature) = ft_sourcestatistics(cfg, sources_all.(bl).(feature).std{:},...
            sources_all.(bl).(feature).dev{:});
    end
end


%% 2. conditions: tonal versus atonal: Different sources? (t-test effect of entropy) (section 2 and 3 not yet fully adapted)

% comparison_name = 'tonal_vs_atonal';
% 
% %%% main effect
% 
% Nsub = 18;
% cfg.design = [];
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% 
% 
% for aa = 1:length(features)
%     feature = features{aa};
%     results_sources.(comparison_name).(feature).pow = ft_sourcestatistics(cfg,...
%         sources_all.tonal.(feature){:},...
%         sources_all.atonal.(feature){:});
%     results_sources.(comparison_name).(feature).pow.cfg.previous = [];
% end
% 
 save([out_dir_stats out_prefix 'cb_permutations_sources_test3'],'results_sources','-v7.3')


%% 3. GRAND AVERAGE

disp('source grand averaging')

cfg = [];
%cfg.method = 'across';
%cfg.latency = 'all';
cfg.parameter = 'pow';
cfg.keepindividual = 'no';

blocks2 = fieldnames(sources_all);

grand_avg_sources = [];

for jj = 1:length(blocks2)
    features2 = fieldnames(sources_all.(blocks2{jj}));
    for hh = 1:length(features2)
        
        cnames = fieldnames(sources_all.(blocks2{jj}).(features2{hh}));
        
        for c=1:length(cnames)
            cname = cnames{c};
            
            grand_avg_sources.(blocks2{jj}).(features2{hh}).(cname) = ...
                ft_sourcegrandaverage(cfg,sources_all.(blocks2{jj}).(features2{hh}).(cname){1:end});
        end
    end
end

% saving sources grand average

save([out_dir_data out_prefix 'grand_avg_sources_test3'],'grand_avg_sources', '-v7.3');

%clear grand_all

%% PLOT RESULTS

% 1. comparing sources for atonal deviants and tonal deviants 

% 1.1 features = {'pitch', 'intensity', 'timbre', 'location'}; % pi

views = [90, -90];
%clims = {[0 2.5e-24], [0 7e-24], [0 3e-24], [0 7e-24]};

clims = {[0 4.5e-24], [0 10e-24], [0 5e-24], [0 10e-24]};

pars = {'MMN'}

for p=1:length(pars)
    par = pars{p}
  
for f=1:length(features)
    feature = features{f};
    
    for c=1:length(blocks)
        block = blocks{c};
        
        source_tonal = grand_avg_sources.(block).(feature).(par);
        source_tonal.anatomy = sources_all.(block).(feature).(par){1,1}.anatomy;
        
         %source_tonal.pow = results_sources.MMN.(block).(feature).stat;
%         %
%         source_tonal.mask = results_sources.MMN.(block).(feature).mask;
%         %
%          source_tonal.mask =  source_tonal.pow > max(source_tonal.pow(:))*.30 |...
%              source_tonal.pow < min(source_tonal.pow(:))*.30; % 30 % of maximum; positive and negative effects; larger than 80% of the max. 
        
%         s_clusters = ismember(results_sources.MMN.(block).(feature).posclusterslabelmat,...
%             find(struct2table(results_sources.MMN.(block).(feature).posclusters).prob <= 0.01));
%         %source_tonal.mask = source_tonal.mask.*s_clusters;
%         source_tonal.mask = s_clusters; 
%         
%         %  source_tonal.mask = zeros(size(results_sources.MMN.tonal.pitch.pos,1),1);
        % indVoxClust = find(results_sources.MMN.tonal.pitch.posclusterslabelmat==1);
               
        cfg = [];
        cfg.method         = 'surface';
        cfg.funparameter   = 'pow';
        %cfg.maskparameter  = cfg.funparameter;
        %cfg.maskparameter  = 'mask';
        %cfg.funcolorlim    = clims{f};
        cfg.funcolormap    = 'viridis'; % 'parula'; %'viridis'; %'parula'; %'jet';
        %cfg.opacitymap     = 'rampup';
        cfg.projmethod     = 'project'; %'project';
        cfg.surfinflated   = 'surface_inflated_both.mat'; %surface_inflated_both.mat'; %'surface_inflated_both.mat';
        %cfg.surfdownsample = 10;
        cfg.camlight       = 'no';
        
        for v=1:length(views)
            
            ft_sourceplot(cfg, source_tonal);
            
            view ([views(v) 0]);
            
           %savefig(sprintf('%s/%s_%s_%d_surface',out_dir_fig, block, feature, views(v)));
           %print(sprintf('%s/%s_%s_%d_surface',out_dir_fig, block, feature, views(v)), '-dpdf'); 
           
           %close all
           
        end
        
    end
    
end
end

%% savefig(sprintf('%s/pitch_tonal_surface3',out_dir_fig));





% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.maskparameter = cfg.funparameter;
% %cfg.funcolorlim   = [0.0 maxval];
% %cfg.opacitylim    = [0.0 maxval];
% cfg.maskparameter = 'mask';
% cfg.opacitymap    = 'rampup';
% ft_sourceplot(cfg,source_tonal);
% print(sprintf('%s/intensity_tonal_slice',out_dir_fig), '-dpdf'); 
% %savefig(sprintf('%s/pitch_tonal_slice',out_dir_fig));
% 


% MNI white-matter surface



% 
% 
% % 1.2 atonal (pitch, intensity, location, timbre)
% 
% source_atonal = grand_avg_sources.atonal.intensity;
% source_atonal.anatomy = sources_all.atonal.intensity{1,1}.anatomy;
% 
% source_atonal.mask =  source_atonal.pow > max(source_atonal.pow(:))*.40 |...
%     source_atonal.pow < min(source_atonal.pow(:))*.40 
% 
% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.maskparameter = cfg.funparameter;
% cfg.maskparameter = 'mask';
% cfg.opacitymap    = 'rampup';
% ft_sourceplot(cfg,source_atonal);
% 
% print(sprintf('%s/intensity_atonal_slice',out_dir_fig), '-dpdf');
% 
% close all
% 
% %savefig(sprintf('%s/intensity_atonal_slice',out_dir_fig));
% 
% % MNI white-matter surface
% 
% cfg = [];
% cfg.method         = 'surface';
% cfg.funparameter   = 'pow';
% cfg.maskparameter  = 'mask';
% cfg.funcolorlim    = [-1 1];
% cfg.funcolormap    = 'jet';
% cfg.opacitymap     = 'rampup';
% cfg.projmethod     = 'project';
% cfg.surfinflated   = 'surface_inflated_both.mat'; 
% cfg.camlight       = 'no'; 
% 
% ft_sourceplot(cfg, source_atonal);
% %savefig(sprintf('%s/pitch_atonal_surface1',out_dir_fig));
% ft_sourceplot(cfg, source_atonal);
% view ([90 0]);
% %savefig(sprintf('%s/pitch_atonal_surface2',out_dir_fig));
% ft_sourceplot(cfg, source_atonal);
% view ([-90 0]);
% %savefig(sprintf('%s/pitch_atonal_surface3',out_dir_fig));
% 
% close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. Interaction entropy X expertise (difference of diffenrences, here only one interaction. we need 3)
% 
% cfg.statistic   = 'indepsamplesT';
% cfg.ivar = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar = [];
% comparison_name = 'interaction';
% 
% %%% main effect
% 
% % Nonmus - Classical (non, cla)
% 
% Nsub1 = 2;
% Nsub2 = 2;
% cfg.design = [];
% cfg.design(1,1:(Nsub1+Nsub2))  = [ones(1,Nsub1) 2*ones(1,Nsub2)];
% 
% for aa = 2:length(features)
%     feature = features{aa}; 
%     
%     data_non = sources_all.cond_diff.(feature)(scodes.nonmus);
%     data_cla = sources_all.cond_diff.(feature)(scodes.classical);
%     results_sources.(comparison_name).(feature).non_cla = ft_timelockstatistics(cfg,data_non{:},data_cla{:});
%     results_sources.(comparison_name).(feature).non_cla.cfg.previous = [];
% end
% 
% % Classical - Contemporary (cla, con)
% 
% Nsub1 = 2;
% Nsub2 = 2;
% cfg.design = [];
% cfg.design(1,1:(Nsub1+Nsub2))  = [ones(1,Nsub1) 2*ones(1,Nsub2)];
% 
% for aa = 2:length(features)
%     feature = features{aa}; 
%     
%     data_cla = sources_all.cond_diff.(feature)(scodes.classical);
%     data_con = sources_all.cond_diff.(feature)(scodes.contemporary);
%     results_sources.(comparison_name).(feature).cla_con = ft_timelockstatistics(cfg,data_cla{:},data_con{:});
%     results_sources.(comparison_name).(feature).cla_con.cfg.previous = [];
% end
% 
% % Nonmus - Contemporary (non, cont)
% 
% Nsub1 = 2;
% Nsub2 = 2;
% cfg.design = [];
% cfg.design(1,1:(Nsub1+Nsub2))  = [ones(1,Nsub1) 2*ones(1,Nsub2)];
% 
% for aa = 2:length(features)
%     feature = features{aa}; 
%     
%     data_non = sources_all.cond_diff.(feature)(scodes.nonmus);
%     data_con = sources_all.cond_diff.(feature)(scodes.contemporary);
%     results_sources.(comparison_name).(feature).non_con = ft_timelockstatistics(cfg,data_non{:},data_con{:});
%     results_sources.(comparison_name).(feature).non_con.cfg.previous = [];
% end

% comparison_name = 'slide_HE_LE_two_TW';
% twindows = {[0.1 0.2],[0.25 0.35]};
% 
% for aa = 1:length(twindows)
%     cfg.latency = twindows{aa};
%     results.(comparison_name).(['TW_' num2str(aa)]) = ft_timelockstatistics(cfg,grand_all.LE.s.MMN{:},grand_all.HE.s.MMN{:});
%     results.(comparison_name).(['TW_' num2str(aa)]).cfg.previous = [];
% end

%% simple effects (HE and LE); we would need 6 comparisons here. Leave it
% %%% out for now.
% 
% for aa = 1:length(features)
%     feature = features{aa};
%     for bb = 1:length(blocks)
%         block = blocks{bb};
%         data_mus = grand_all.(block).(feature).MMN(scodes.mus);
%         data_nmus = grand_all.(block).(feature).MMN(scodes.nmus);
%         results.(comparison_name).(feature).(block) = ft_timelockstatistics(cfg,data_mus{:},data_nmus{:});
%         results.(comparison_name).(feature).(block).cfg.previous = [];
%     end
% end




% %% PLOTTING
% disp('plotting')test3
% if ~exist('grand_avg', 'var')
%     load([out_dir_data out_prefix '_grand_average']) % load grand average
% end
% 
% layouts = {'neuromag306mag.lay','neuromag306cmb.lay'};
% lay_name = {'mag','grad'};
% limits = {[-1.7e-13 1.4e-13],[-1e-12 5.5e-12]};
% cfg = [];
% cfg.fontsize = 6;
% cfg.graphcolor = 'brk';
% cfg.showlabels = 'no';
% cfg.comment = 'no';
% cfg.channel = 'MEG';
% cfg.xlim = [-0.1 0.4];
% 
% 
%     layout = lay_name{ll};
%     cfg.ylim = limits{ll};
%     cfg.layout = layouts{ll};%'neuromag306mag.lay'; %'neuromag306cmb.lay' 'vertical';
%     
%     expertise2 = {'mus','nmus','all'};
% %     %% stdDev
%     for hh = 1:length(features)
%         for jj = 1:length(blocks)
%             for ee = 1:length(expertise2)
%                 exp2 = expertise2{ee};
%                 figure('Color','white');
%                 gg = gcf;
%                 set(gg, 'PaperOrientation', 'landscape')
%                 set(gg, 'PaperUnits', 'normalized')
%                 set(gg, 'PaperPosition', [0 0 1 1])
%                 set(gg, 'Renderer','painters')
%                 ft_multiplotER(cfg, grand_avg.(blocks{jj}).(features{hh}).s.(exp2), grand_avg.(blocks{jj}).(features{hh}).d.(exp2));
%                 title([feature_names{hh} ' ' block_names{jj} ' ' exp2])
%                 if hh == 1 && jj == 1 && ee == 1
%                     print([out_dir_figures out_prefix '_all_sensors_stDev_' layout],'-dpsc');
%                 else
%                     print([out_dir_figures out_prefix '_all_sensors_stDev_' layout],'-dpsc','-append');
%                 end
%                 if strcmp(plot_interact, 'no')
%                     close all
%                 end
%             end
%         end
%     end
    
%     %% MMN
%     
%     for hh = 1:length(features)
%         for ee = 1:length(expertise2)
%             exp2 = expertise2{ee};
%             figure('Color','white');
%             gg = gcf;
%             set(gg, 'PaperOrientation', 'landscape')
%             set(gg, 'PaperUnits', 'normalized')
%             set(gg, 'PaperPosition', [0 0 1 1])
%             set(gg, 'Renderer','painters')
%             ft_multiplotER(cfg, grand_avg.LE.(features{hh}).MMN.(exp2), grand_avg.HE.(features{hh}).MMN.(exp2));
%             title([feature_names{hh} ' ' exp2])
%             if hh == 1 && ee == 1
%                 print([out_dir_figures out_prefix '_all_sensors_Diff_' layout],'-dpsc');
%             else
%                 print([out_dir_figures out_prefix '_all_sensors_Diff_' layout],'-dpsc','-append');
%             end
%             if strcmp(plot_interact, 'no')
%                 close all
%             end
%         end
%     end
%     
% %     %% MMN difference LE-HE
%     
%     for hh = 1:length(features)
%         figure;
%         gg = gcf;
%         set(gg, 'PaperOrientation', 'landscape')
%         set(gg, 'PaperUnits', 'normalized')
%         set(gg, 'PaperPosition', [0 0 1 1])
%         set(gg, 'Renderer','painters')
%         ft_multiplotER(cfg, grand_avg.entropy_diff.(features{hh}).HE_LE.mus,...
%             grand_avg.entropy_diff.(features{hh}).HE_LE.nmus);
%         title(feature_names{hh})
%         if hh == 1
%             print([out_dir_figures out_prefix '_all_sensors_HE_LE_diff_' layout],'-dpsc');
%         else
%             print([out_dir_figures out_prefix '_all_sensors_HE_LE_diff_' layout],'-dpsc','-append');
%         end
%         if strcmp(plot_interact, 'no')
%             close all
%         end
%     end
%     
% %     %% MMN average LE-HE
%     
%     for hh = 1:length(features)
%         figure;
%         gg = gcf;
%         set(gg, 'PaperOrientation', 'landscape')
%         set(gg, 'PaperUnits', 'normalized')
%         set(gg, 'PaperPosition', [0 0 1 1])
%         set(gg, 'Renderer','painters')
%         ft_multiplotER(cfg, grand_avg.entropy_avg.(features{hh}).HE_LE.mus,...
%             grand_avg.entropy_avg.(features{hh}).HE_LE.nmus);
%         title(feature_names{hh})
%         if hh == 1
%             print([out_dir_figures out_prefix '_all_sensors_HE_LE_avg_' layout],'-dpsc');
%         else
%             print([out_dir_figures out_prefix '_all_sensors_HE_LE_avg_' layout],'-dpsc','-append');
%         end
%         if strcmp(plot_interact, 'no')
%             close all
%         end
%     end

%