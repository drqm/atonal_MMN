function out = epoching_MM(input)
out =[];
mumufes = [];
rawdir = input.rawdir;
ICA_dir = input.ICA_dir;

timelockon = input.timelock;
%fname = input.fname;
blocks = input.blocks; %bl
blockn = input.blockn;  %bn
blocki = input.blocki;
comps = input.comps;
ID = input.ID;
outdir = input.outdir;
data_dir = [outdir 'data/'];
%data_dir = [outdir 'testdata/'];
%data_dir = [outdir 'comparisondata/'];
fig_dir = [outdir 'figures/'];
%fig_dir = [outdir 'testfigure/'];
%fig_dir = [outdir 'comparisonfigures/'];

if ~exist(data_dir,'dir')
    mkdir(data_dir)
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir)
end

load('/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip/template/neighbours/ctf275_neighb.mat');

if timelockon == 1
    nrejected = [];
    for b = 1:length(blocks)
        bn = blockn{b};
        bl = blocks{b};
        bi = blocki{b};
        bads = split(input.bads{b},';');      
        
        load([ICA_dir, bi{1}]);
        
        %% epoching
        
        cfg = [];
        cfg.dataset = [rawdir, bn{1}];
        cfg.trialdef.eventtype = 'UPPT001';  % Specifying where markers are located.
        cfg.trialdef.eventvalue = [231:234, 241:244];  % Trigger values to include.
        cfg.trialdef.prestim = 0.100;  % baseline correction
        cfg.trialdef.poststim = 0.400;  % Note that next sound starts already after 250 ms
        cfg = ft_definetrial(cfg);
        
        cfg.trl(:,1:2) = cfg.trl(:,1:2)-60;
        
        %old_trial_list = cfg.trl; %saving trigger list for quality control
        
        % semi-automatic artifact rejection
      
        % cfg = [];
        % cfg.trl = trl;

%         cfg.datafile = [rawdir, fname{1,1} '_0', num2str(bn), '.ds'];
%         cfg.headerfile = [rawdir, fname{1,1} '_0', num2str(bn), '.ds'];

        % channel selection, cutoff and padding
        cfg.artfctdef.zvalue.channel = 'MEG';
        cfg.artfctdef.zvalue.cutoff = 30;
        cfg.artfctdef.zvalue.trlpadding = 0;
        cfg.artfctdef.zvalue.artpadding = 0;
        cfg.artfctdef.zvalue.fltpadding = 0;

        % algorithmic parameters
        cfg.artfctdef.zvalue.cumulative = 'yes';
        cfg.artfctdef.zvalue.medianfilter = 'yes';
        cfg.artfctdef.zvalue.medianfiltord = 9;
        cfg.artfctdef.zvalue.absdiff = 'yes';

        % make the process interactive
        %cfg.artfctdef.zvalue.interactive = 'yes';

        [cfg, artifact_jump] = ft_artifact_zvalue(cfg);

        % reject:
        cfg.artfctdef.reject = 'complete'; % this rejects complete trials, use 'partial' if you want to do partial artifact rejection
        [cfg] = ft_rejectartifact(cfg);
              
        nrejected.(bl) = length(cfg.trlold) - length(cfg.trl);
        trl_selection = cfg.trl;
        
        % FILTERING
        
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.6;
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 40;
        cfg.demean = 'yes'; % notch dftfilter check if it can be changed
        cfg.baselinewindow = [-0.1 0];
        data = ft_preprocessing(cfg); % data = ft_preprocessing(cfg,data);
        
        cfg=[];
        cfg.detrend = 'no';
        cfg.resamplefs = 250;
        res_data = ft_resampledata(cfg,data);
        
        if ~isempty(bads) % interpolation of bad channels. 
            
            cfg = [];
            
            cfg.layout = '/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip/template/layout/CTF275.lay';
            cfg.layout = ft_prepare_layout(cfg);
            cfg.neighbours = neighbours;
            cfg.method = 'average';
 %          cfg.missingchannel = bads;
            cfg.badchannel = bads;
            cfg.senstype = 'meg';
            res_data = ft_channelrepair(cfg, res_data);
        end
        
        % ICA component rejection
        
        cfg           = [];
        cfg.component = comps{b};
        data_clean    = ft_rejectcomponent(cfg, meg_channels_comp, res_data);
        
        
        
        %% TIMELOCKING
        
        
        triggers = {[231:234]; 241; 242; 243; 244};
        features = {'standard'; 'pitch'; 'intensity'; 'timbre'; 'location'};
        
        for c = 1:length(features)
            trig = triggers{c};
            feature = features{c};
            
            cfg = [];
            cfg.channel = 'MEG';
            cfg.covariance = 'yes';
            %cfg.keeptrials = 'yes'; % for headmovement correction; 
            cfg.trials = find(ismember(data.trialinfo(:,1),trig));
            
            t_locked = ft_timelockanalysis(cfg, data_clean);
            
            cfg = [];
            cfg.channel = 'MEG';
            cfg.baseline = [-0.1 0];
            cfg.covariance = 'yes';
            t_locked = ft_timelockbaseline(cfg, t_locked);
            
            t_locked.cfg.previous = [];
            
            %mumufes.(bl).(feature) = ft_timelockanalysis(cfg,data_clean);
            %mumufes.(bl).(feature).cfg.previous = [];
            
            % regressing out headmovements
            
            % define trials / Create trial-by-trial estimates of head movement.
%             
%             cfg = [];
%             cfg.dataset = [rawdir, fname{1,1} '_0', num2str(bn), '.ds'];
%             cfg.trialdef.eventtype = 'UPPT001';  % Specifying where markers are located.
%             cfg.trialdef.eventvalue = [231:234, 241:244];  % Trigger values to include.
%             cfg.trialdef.prestim = 0.100;  % baseline correction
%             cfg.trialdef.poststim = 0.400;  % Note that next sound starts already after 250 ms
%             cfg = ft_definetrial(cfg);
%             
%             % preprocess the headposition data
%             
%             cfg.channel = {'HLC0011','HLC0012','HLC0013',...
%                 'HLC0021','HLC0022','HLC0023',...
%                 'HLC0031','HLC0032','HLC0033'};
%             
%             %cfg.trl = trl_selection;
%             cfg.trl = trl_selection(find(ismember(trl_selection(:,4),trig)),:);
%             headpos = ft_preprocessing(cfg);
%             sinfo = headpos.sampleinfo;
%             cfg=[];
%             cfg.detrend = 'no';
%             cfg.resamplefs = 250;
%             headpos = ft_resampledata(cfg,headpos);
%             headpos.sampleinfo = sinfo;
%             % calculate the mean coil position per trial
%             
%             ntrials = length(headpos.sampleinfo)
%             coil1 = [];coil2 = []; coil3 = [];
%             
%             for t = 1:ntrials
%                 coil1(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
%                 coil2(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
%                 coil3(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
%             end
%             
%             % calculate the headposition and orientation per trial
%             cc = circumcenter(coil1, coil2, coil3);
%             
%             % demean to obtain translations and rotations from the average position and orientation
%             cc_dem = [cc - repmat(mean(cc,2),1,size(cc,2))]';
%             
%             %Fit the headmovement regressors to the data and remove variance that can be explained by these confounds.
%             % add head movements to the regressorlist. also add the constant (at the end; column 7)
%             confound = [cc_dem ones(size(cc_dem,1),1)];
%             
%             % regress out headposition confounds
%             cfg = []; 
%             cfg.confound = confound; 
%             cfg.reject = [1:6];
%             
%             % keeping the constant (nr 7)
%             t_locked_reg = ft_regressconfound(cfg, t_locked);
% 
%           
%           % headmovement correction or not:
%           if input.regconf == 1
%            % %  t_locked.avg = squeeze(t_locked_reg.beta(7,:,:)); 
%            t_locked.avg = squeeze(mean(t_locked_reg.trial,1));
%           end    
%           
             mumufes.(bl).(feature) = t_locked;
%                      
            if c~= 1
                cfg = [];
                cfg.operation = 'subtract';
                cfg.parameter = 'avg';
                mumufes.(bl).([feature '_MMN'])= ft_math(cfg,mumufes.(bl).(feature),mumufes.(bl).standard);
                mumufes.(bl).([feature '_MMN']).grad = mumufes.(bl).standard.grad;
                mumufes.(bl).([feature '_MMN']).cfg.previous = [];
            end
            
        end
        
    end
    
    mumufes.nrejected = nrejected;
    save(sprintf('%s%02.f',data_dir,ID),'mumufes','-v7.3')
    
else load(sprintf('%s%02.f',data_dir,ID))
end

% PLOT RESULTS

features = {'pitch'; 'intensity'; 'timbre'; 'location'};

for f = 1:length(features)
    feature = features{f};
    MMN = [];
    
    for b = 1:length(blocks)
        bl = blocks{b};
        
        % substract standard from deviant
        
        standard = mumufes.(bl).standard;
        deviant = mumufes.(bl).(feature);
        MMN.(bl) = mumufes.(bl).([feature '_MMN']);
        %MMN.(bl).avg = deviant.avg - standard.avg;
        
        % % Plot data
        %
        cfg = [];
        cfg.fontsize = 6;
        cfg.xlim = [-0.1 0.4];
        cfg.graphcolor = 'kr';
        cfg.showlabels = 'no';
        cfg.comment = 'no';
        cfg.channel = 'MEG';        
        %cfg.ylim = [-2e-12 9e-12];
        cfg.layout = 'CTF275.lay';
        cfg.channel = [];
        
        figure;
        ft_multiplotER(cfg, standard, deviant);
        savefig(sprintf('%s%02.f_%s_%s',fig_dir,ID,bl,feature));
        %print(sprintf('%s%02.f_%s_%s',fig_dir,ID,bl,feature),'-dpng');
        
        close all
        
    end
    
    cfg.channel = [];
    cfg.graphcolor = 'br'; 
    figure;
    ft_multiplotER(cfg, MMN.atonal, MMN.tonal);
    savefig(sprintf('%s%02.f_%s_MMN',fig_dir,ID,feature));
    
    close all
end

