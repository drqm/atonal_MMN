clear all

addpath('/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip_git/fieldtrip/')
ft_defaults

%%
features = {'standard','pitch','location','intensity','timbre'};
%features = {'pitch','location','intensity','timbre'};
feature_names = {'standard','pitch','location','intensity','timbre'};
block_names = {'tonal','atonal'};
plot_interact = 'no';
out_prefix = 'entropy';

blocks = {'tonal','atonal'};
conditions = {'deviant','standard','MMN'};
input_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/timelocked/data/';

stats_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/stats/';
data_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/data/';
out_dir_figures = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/figures/';

if ~exist('grand_avg','var')
load([data_dir 'grand_average.mat'])
end

if ~exist('results','var')
load([stats_dir 'cb_permutations.mat'])
end

time = grand_avg.tonal.(features{1}).all.time*1000;
hemispheres = {'right','left'};

%new channels based on the largest P50 response: 
sens.right = {'MRP56', 'MRP57', 'MRP55', 'MRT15', 'MRT16', 'MRT14'}; %'MRT22', 'MRT33', 'MRT23', 'MRT32', 'MRT12', 'MRT42',};              
sens.left = {'MLO14','MLP55', 'MLP56', 'MLT16', 'MLT27', 'MLT15'}; %'MLT22', 'MLT33', 'MLT23', 'MLT32', 'MLT12', 'MLT42',};                 

%old channels
%sens.right = {'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'}; % channels found in the previous section.
%sens.left = {'MLO13','MLO14', 'MLO24', 'MLP53', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT14', 'MLT15', 'MLT16', 'MLT24', 'MLT25', 'MLT26', 'MLT27'};  % channels found in the previous section.

%channels.grad = {index_right,index_left};

ntests = 4; % bonferroni correcting
n_participants = 20;

%% Load individual data

starting_point = 1;

grand_all = [];

selection = 1:20;
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
        end
    end
end

%% MMN check

blank_data = grand_avg.tonal.pitch_MMN.all;
labels = blank_data.label;
features = {'pitch','location','intensity','timbre'};
feature_names = {'pitch','location','intensity','timbre'};

for ss = 1:length(hemispheres)
    plot_count = 0;
    figure('Color','white');
    gg = gcf;
    set(gg, 'PaperOrientation', 'portrait')
    set(gg, 'PaperUnits', 'normalized')
    set(gg, 'PaperPosition', [0 0 1 1])
    set(gg, 'Units', 'normalized')
    set(gg, 'Renderer','painters')
    set(gg, 'Position', [0 0 1 1])
    
    hem = hemispheres{ss};   
    [ch, chann] = ismember(sens.(hem), labels);
     
    for hh = 1:length(features)
        feature = features{hh};
        for kk = 1:length(blocks)
            block = blocks{kk};
            cluster_time = [];
            pval = [];
            if isfield(results.MMN_check.([block feature]).all,'negclusters')
                cluster_mask = results.MMN_check.([block feature]).all.negclusterslabelmat(chann,:);
                significance = results.MMN_check.([block feature]).all.negclusters(1).prob;
                [row, col] = find(cluster_mask == 1);
                if significance < 0.001
                    pval = 'p < .001';
                else
                    pval = ['p =' sprintf('%.3f',significance)];
                    pval(4) = ' ';
                end
                if significance <= 0.025
                    cluster_time = results.MMN_check.([block feature]).all.time(col)*1000;
                end
            end
            
            y_cluster = zeros(1,length(cluster_time)) - 68; % green line
            plot_count = plot_count +1;
            stand_avg = mean(grand_avg.(block).standard.all.avg(chann,:)/1e-15);
            dev_avg = mean(grand_avg.(block).(feature).all.avg(chann,:)/1e-15);
            MMN_avg = mean(grand_avg.(block).([feature '_MMN']).all.avg(chann,:)/1e-15);
            error_up = MMN_avg + 1.96*(mean(sqrt(grand_avg.(block).([feature '_MMN']).all.var(chann,:)))/1e-15)/sqrt(n_participants);
            error_down = MMN_avg - 1.96*(mean(sqrt(grand_avg.(block).([feature '_MMN']).all.var(chann,:)))/1e-15)/sqrt(n_participants);
            time_error=[time,fliplr(time)];                %#create continuous x value array for plotting
            MMN_error=[error_down,fliplr(error_up)];  %#create y values for out and then back
            
            subplot(4,2,plot_count);
            ylimit = [-70 70];
            set(gca,'fontsize',12)
           
            if ~isempty(cluster_time)
                plot(cluster_time,y_cluster,'color','g','LineWidth',4)
                hold on
            end
            
            for iss = 1:length(grand_all.(block).([feature '_MMN']))
                ind_plot_s = plot(time,mean(grand_all.(block).([feature '_MMN']){1,iss}.avg(chann,:)/1e-15),'k');
                alpha(.9)
                ind_plot_s.Color(4) = 0.1; % single subjects lines
                hold on
            end
            
            hline = refline([0,0]); % horizontal line / x-line
            hline.Color = 'k';
            line([0,0],ylimit, 'color','black','LineStyle','--'); % y-axis
            line([250,250],ylimit, 'color','black','LineStyle','--');   
            hold on;
            
            fill(time_error,MMN_error,'k','LineStyle','none'); % shaded area, std error
            alpha(.15)
            hold on;
            
            
            %             std_plt = plot(time,stand_avg,':',...
            %                 'LineWidth',3,...
            %                 'MarkerEdgeColor','b');
            %             hold on;
            %
            %             dev_plt = plot(time,dev_avg,':',...
            %                 'LineWidth',3,...
            %                 'MarkerEdgeColor','r');
            %             hold on;
            %
           
            std_plt = plot(time,stand_avg,'b','LineWidth',1.2); % standard line, blue
            hold on;
            
            dev_plt = plot(time,dev_avg,'r','LineWidth',1.2); % deviant line; red
            %dev_plt.Color(4) = 0.45;
            hold on;
            
            MMN_plt = plot(time,MMN_avg,'k','LineWidth',2.2);
            hold on;
            
            if isempty(pval)
                
                text(280,50,'\fontsize{9} No clusters found')
            else
                text(280,50,['\fontsize{9}' pval])
            end
            
            hold off
            
            ylim(ylimit)
            xlim([-100 400])
            box off
            
            if plot_count == 2
                xlabel('\fontsize{9} time (ms)')
                ylabel('\fontsize{9} amplitude (fT/cm)')
                legend([std_plt,dev_plt,MMN_plt],{'Standard','Deviant','MMN'},'Location','northwest',...
                 'Position',[0.3,0.02,0.4,0.05],'FontSize',15,'Orientation','horizontal')
                legend boxoff
            end
        end
    end
    
    annotation('textbox',[.02 .82 .05 .05],'String','Pitch','EdgeColor','none','FontSize',13)
    annotation('textbox',[.02 .6 .05 .05],'String','Location','EdgeColor','none','FontSize',13)
    annotation('textbox',[.02 .38 .05 .05],'String','Intensity','EdgeColor','none','FontSize',13)
    annotation('textbox',[.02 .16 .05 .05],'String','Timbre','EdgeColor','none','FontSize',13)
    annotation('textbox',[.27 .93 .05 .05],'String','tonal','EdgeColor','none','FontSize',20)
    annotation('textbox',[.71 .93 .05 .05],'String','atonal','EdgeColor','none','FontSize',20)
%   saveas(gcf,[out_dir_figures, out_prefix '_stDev_chansel_grad_' hem],'fig')
    print([out_dir_figures  'Fig2-stDev_chansel_grad_' hem],'-dpdf', '-r300')
  % print([out_dir_figures  'Fig3-stDev_chansel_grad_' hem],'-djpeg', '-r1200')

    %close all
end



%% Main effect Entropy / FIGURE 4

features = {'pitch','location','intensity','timbre'};

for ss = 1:length(hemispheres)
    plot_count = 0;
    figure('Color','white');
    
    gg = gcf;
    set(gg, 'PaperOrientation', 'landscape')
    set(gg, 'PaperUnits', 'normalized')
    set(gg, 'PaperPosition', [0 0 1 1])
    set(gg, 'Units', 'normalized')
    set(gg, 'Renderer','painters')
    set(gg, 'Position', [0 0 1 1])
    positions = {[0.13 0.83 0.07 0.08]; % = sizes of brains (the last two numbers); positions = first two number
         [0.57 0.83 0.07 0.08];
         [0.13 0.35 0.07 0.08];
         [0.57 0.35 0.07 0.08]};
    
     hem = hemispheres{ss}; 
    [ch, chann] = ismember(sens.(hem), labels);
    
    for hh = 1:length(features)
        feature = features{hh};
        cluster_time = [];
        pval = [];
        significance = [];
        if isfield(results.tonal_vs_atonal.(feature).all,'negclusters') && ~isempty(results.tonal_vs_atonal.(feature).all.negclusters)% find and select timepoints with the most prominent cluster
            cluster_mask = results.tonal_vs_atonal.(feature).all.negclusterslabelmat(chann,:);
            significance = results.tonal_vs_atonal.(feature).all.negclusters(1).prob; %*ntests; not multiplying the p-values by number of test; lowered the threshold 
            [row, col] = find(cluster_mask == 1);
            if significance < 0.001
                pval = 'p < .001';
            else
                pval = ['p =' sprintf('%.3f',significance)];
                pval(4) = ' ';
            end
            if significance <= 0.025
                cluster_time = results.tonal_vs_atonal.(feature).all.time(col)*1000;
            end
        end
        
        y_cluster = zeros(1,length(cluster_time)) - 78;
        plot_count = plot_count +1;
       
        tonal_avg = mean(grand_avg.tonal.([feature '_MMN']).all.avg(chann,:)/1e-15); %compute mean of selected channels
        atonal_avg = mean(grand_avg.atonal.([feature '_MMN']).all.avg(chann,:)/1e-15); 
        
        diff_avg =  mean(grand_avg.cond_diff.(feature).all.avg(chann,:)/1e-15);
        
        diff_error_up = diff_avg + 1.96*(mean(sqrt(grand_avg.cond_diff.(feature).all.var(chann,:)))/1e-15)/sqrt(n_participants);
        diff_error_down = diff_avg - 1.96*(mean(sqrt(grand_avg.cond_diff.(feature).all.var(chann,:)))/1e-15)/sqrt(n_participants);
        
        time_error=[time,fliplr(time)];                %#create continuous x value array for plotting
        diff_error=[diff_error_down,fliplr(diff_error_up)];         
        
        subplot(2,2, plot_count);
        set(gca,'fontsize',14)
        ylimit = [-80 80];
       
        for iss = 1:length(grand_all.cond_diff.(feature))
            ind_plot_s = plot(time,mean(grand_all.cond_diff.(feature){1,iss}.avg(chann,:)/1e-15),'k');
            ind_plot_s.Color(4) = 0.1; % single subjects lines
            hold on
        end
        
        hline = refline([0,0]);
        hline.Color = 'k';
        line([0,0],ylimit, 'color','black','LineStyle','--');
        line([250,250],ylimit, 'color','black','LineStyle','--');
        hold on
       
        fill(time_error,diff_error,'k','LineStyle','none'); % standard error
        alpha(.15)
        hold on
        
        atonal_plt = plot(time,atonal_avg,'b','LineWidth',1.2);
        hold on;
        
        tonal_plt = plot(time,tonal_avg,'r','LineWidth',1.2);
        hold on;
        
        diff_plot = plot(time,diff_avg,'k','LineWidth',2.2);
        hold on;
        
        if isempty(pval)
            text(300,32,'\fontsize{14} p = 1')
        elseif significance > 1
            text(300,32,'\fontsize{14} p = 1')
        else
            text(300,32,['\fontsize{14}' pval])
        end
        if ~isempty(cluster_time)
            plot(cluster_time,y_cluster,'g','LineWidth',8)
        end
        hold on
%        hold off
        
        ylim(ylimit)
        xlim([-100 400])
        
        title(['\fontsize{20} ' feature_names{hh}])
        box off
        
        if plot_count == 1
            xlabel('\fontsize{14} time (ms)')
            ylabel('\fontsize{14} amplitude (fT/cm)')
           legend([tonal_plt,atonal_plt, diff_plot],{'tonal','atonal','Difference (tonal-atonal)'},...
                  'Position',[0.2,0.02,0.65,0.05],'FontSize',20,'Orientation','horizontal')
           legend boxoff
        end
        
        tw = 52:91;
        
        cfg = [];
        cfg.parameter = 'avg';
        % cfg.xlim = 0:0.025:0.3;
        cfg.colorbar = [];%'East';
        cfg.colormap = parula;
        cfg.style = 'straight';
        cfg.comment = 'no';
        cfg.marker = 'off';
        cfg.gridscale = 268;
        cfg.layout = 'CTF275.lay';
        cfg.zlim = [-15*10^(-15) 15*10^(-15)]/10^(-15);
        
        dataset = grand_avg.cond_diff.(feature).all;
        dataset.avg = dataset.avg/10^(-15);
        
        current_data = mean(dataset.avg(chann,:));
        [val, idx] = min(current_data(tw));
        cfg.xlim = [(dataset.time(idx-7 + tw(1)-1)),(dataset.time(idx + 7 + tw(1)-1))];
        cfg.style = 'straight_imsat';
        axes('Position',positions{hh})
        
        topo = ft_topoplotER(cfg,dataset);
        
        c = colorbar;
        c.TickDirection = 'out';
        c.Position = [positions{hh}(1) + 0.08,positions{hh}(2) + 0.01,0.01,0.08];
        hold off
    end   
    
%   saveas(gcf,[out_dir_figures, out_prefix '_difference_chansel_grad_' hem],'fig')
    print([out_dir_figures 'Fig4-difference_chansel_grad_' hem],'-dpdf','-r300')
 %   print([out_dir_figures 'Fig5-difference_chansel_grad_' hem],'-djpeg','-r1200')
   %close all
end