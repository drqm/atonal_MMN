
% SOURCE LOCALIZATION

% uses timelocked data + headmodels (this could be turned into a function)
% for source time courses check avg.mom

ftdir = '/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip_git/fieldtrip/';
addpath(ftdir)

ft_defaults

% nonmusicians
 subjectID = [1 2 4 6 7 8 9 10 11 12 14 15 16 18 19 20 21 25]; % 17 has no complete MRI; for 5 no MRI at all;

normalWeight =1;
filedir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/scripts/analysis_MM/';
filename = fullfile(filedir,'Preprocessing_MEG_MM.csv');
subjects = readtable(filename);

%bn = {subjects.atonal(subjectID), subjects.tonal(subjectID)};
%fname = subjects.Files(subjectID);
%mri_name = subjects.MRI(subjectID);

tl_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/timelocked/data/'; % tl = timelocked
hm_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/headmodels/'; %hm = headmodel
fig_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/figures/';
out_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/sources/';

blocks = {'tonal','atonal'}; % high entropy, low entropy
features = {'pitch', 'intensity', 'timbre', 'location'}; 
latency = {[0.15, 0.25], [0.1, 0.2], [0.1, 0.2], [0.08, 0.18]};

% channels.right = {'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'}; % channels found in the previous section.
% channels.left = {'MLO13','MLO14', 'MLO24', 'MLP53', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT14', 'MLT15', 'MLT16', 'MLT24', 'MLT25', 'MLT26', 'MLT27'};  % channels found in the previous section.

% anterior channels
channels.right = {'MRT22', 'MRT33', 'MRT23', 'MRT32', 'MRT12', 'MRT42'};               %'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'}; % channels found in the previous section.
channels.left = {'MLT22', 'MLT33', 'MLT23', 'MLT32', 'MLT12', 'MLT42'}; 
directions = [1,-1]; % IMPORTANT: if you use posterior channels change it to: [-1,1]

% posterior channels
% sens.right = {'MRP56', 'MRP57', 'MRP55', 'MRT15', 'MRT16', 'MRT14'};               
% sens.left = {'MLO14','MLP55', 'MLP56', 'MLT16', 'MLT27', 'MLT15'};                 
% directions = [-1,1]; % IMPORTANT: if you use posterior channels change it to: [-1,1]

interval = [0.05, 0.3];
hems = {'right', 'left'};

%mri_dir = ['/hpc/users/iris.mencke/workspace/CCMusic/rawdata/', mri_name{1}];
%mri_file = [mri_dir '/0002_MR000004.dcm'];


for s=SLURM_ARRAY_TASK_ID+1;   % 1:length(subjectID)
    subjectXY  = subjectID(s)
    
    load(sprintf('%s%02.f_headmodel_meg',hm_dir, subjectXY),'headmodel','grid','grad','shape','seg'); % load headmodel
    load(sprintf('%s%02.f',tl_dir,subjectXY),'mumufes'); %load timelocked data
    load(sprintf('%s%02.f_leadfield_meg_N%1.1f.mat',hm_dir,subjectXY,normalWeight), 'leadfield_meg'); % all in cm; all in ctf
    
    %%covariance: loop starts
    
    mumufes_sources = []; % initialize structure for source stats
    
    for hh = 1:length(features)
        feature = features{hh};
        for jj = 1:length(blocks)
            bl = blocks{jj};
            
            standard = mumufes.(bl).standard;
            deviant = mumufes.(bl).(feature);
            
            difdat= standard;
            difdat.avg = deviant.avg - standard.avg;
            
            lat = [];
            for hm = 1:2
                hem = hems{hm}
                direction = directions(hm);
                [check1,gradn_indx] = ismember(channels.(hem),difdat.label);
                ERF_grad = mean(difdat.avg(gradn_indx,:));
                [peaks_grad, samples_grad] = findpeaks(direction*ERF_grad(difdat.time >= interval(1) &...
                    difdat.time <= interval(2)));
                
                if isempty(peaks_grad)
                    [peaks_grad, samples_grad] = max(direction*ERF_grad(difdat.time >= interval(1) &...
                        difdat.time <= interval(2)));
                end
                
                [peak_grad, indx_grad] = max(peaks_grad);
                sample_grad = samples_grad(indx_grad);
                lat(end+1) = difdat.time(sample_grad + min(find(difdat.time >= interval(1)))...
                    - 1);            
            end
            
            timeStart = min(lat) - 0.040;
            timeEnd = max(lat) + 0.040;
             
            data_sd = ft_appenddata([],standard,deviant); % conditions are put together
            
            cfg = [];
            cfg.covariance ='no'; %(previous: 'yes')
            cfg.vartrllength = 2;
            
            avg = ft_timelockanalysis(cfg,data_sd); % covariance of both cond. is calculated
            avg.cov = (standard.cov + deviant.cov)/2; % = giorgos suggestions
            
            %
            %             cfg = [];
            %             cfg.covariance='yes';
            %             avgstd = ft_timelockanalysis(cfg,standard);
            %             avgdev = ft_timelockanalysis(cfg,deviant);
            %
            % BEAMFORMING
            
            cfg=[];
            cfg.method='lcmv';
            cfg.grid= leadfield_meg; % that was calculated earlier
            %cfg.grid.leadfield = leadfield_meg.lead;
            cfg.vol = headmodel;
            cfg.lcmv.keepfilter='yes'; % keep the beamformer solution
            cfg.lcmv.lambda = '35%'; %giorgos: 35%
            cfg.lcmv.fixedori = 'yes';
            cfg.channel     = {'MEG'};
            cfg.senstype    = 'MEG';
            cfg.grad = grad;
            
            sourceavg = ft_sourceanalysis(cfg, avg);
            
            
            % projecting sensor data on source / source inversion
            % computation with the filter above
            
            Nchans = length(avg.label); % for the matrix
            Nsources = sum(sourceavg.inside); % for the matrix
            spatfilt= reshape(([sourceavg.avg.filter{sourceavg.inside}]), Nchans, Nsources)';
            
            % here MMN is computed
            
            avgdif = standard;
            avgdif.avg = deviant.avg - standard.avg;
            
            % selecting the time points (different tome windows for
            % different time points)
            
%             timeStart = latency{hh}(1);
%             timeEnd = latency{hh}(2);
%             
            indStart = nearest(avgdif.time, timeStart);
            indEnd = nearest(avgdif.time, timeEnd);
            
            srcDevMat= spatfilt*deviant.avg(:,indStart:indEnd);
            
            srcStdMat= spatfilt*standard.avg(:,indStart:indEnd);
            srcMMN_diff = spatfilt*avgdif.avg(:,indStart:indEnd);
            
            %s1 = sourceavg;
           
            sourceavg.pos = grid.pos_mni;   
            devsource = sourceavg;
            stdsource = sourceavg; 
            MMNsource = sourceavg;
            MMN_diff = sourceavg;
            
            
            %s1.avg.pow = zeros(size(s1.avg.pow));
            %s1.avg.pow(s1.inside) = mean(srcDevMat,2) - mean(srcStdMat,2);
            
            devsource.pow(devsource.inside) = mean(abs(srcDevMat),2); % source of deviant
            stdsource.pow(stdsource.inside) = mean(abs(srcStdMat),2); % source standard
            MMNsource.pow(MMNsource.inside) = mean(abs(srcDevMat - srcStdMat),2);
            MMN_diff.pow(MMN_diff.inside) = mean(abs(srcMMN_diff),2);
            
            
            %MMNsource.pow = devsource.pow - stdsource.pow;         
        
            %s1.avg.pow(s1.inside) = mean(abs(srcDevMat),2) - mean(abs(srcStdMat),2);
            
            %s1.pos = grid.pos_mni;                   
            
            mriTemplFile =[ftdir,'template/anatomy/single_subj_T1.nii'];
            mriTempl=ft_read_mri(mriTemplFile);
            
            cfg              = [];
            cfg.voxelcoord   = 'no';
            cfg.parameter    = 'avg.pow';
            cfg.interpmethod = 'nearest';
            
            dev_int  = ft_sourceinterpolate(cfg, devsource, mriTempl);
            std_int  = ft_sourceinterpolate(cfg, stdsource, mriTempl);
            MMN_int  = ft_sourceinterpolate(cfg, MMNsource, mriTempl);
            MMN_diff  = ft_sourceinterpolate(cfg, MMN_diff, mriTempl);
           
            %source_norm = s1_int;
            
            % spatially normalize the anatomy and functional data to MNI coordinates
            
            %             cfg = [];
            %             cfg.nonlinear = 'no';
            %
            %             source_norm = ft_volumenormalise(cfg, s1_int);
            
            source_norm.cfg.previous = [];
            
            dev_int.cfg.previous = [];
            std_int.cfg.previous = [];
            MMN_int.cfg.previous = [];
            MMN_diff.cfg.previous = [];
             
            mumufes_sources.(bl).(feature).dev = dev_int;
            mumufes_sources.(bl).(feature).std = std_int;
            mumufes_sources.(bl).(feature).MMN = MMN_int;
            mumufes_sources.(bl).(feature).MMN_diff = MMN_diff;
            
            %% PLOT RESULTS
            
            % fix mask for positive and negative values
            
            cfg = [];
            cfg.method        = 'ortho';
            cfg.funparameter  = 'avg.pow';
            cfg.maskparameter = cfg.funparameter;
            %cfg.funcolorlim   = [0.0 1.2];
            %cfg.opacitylim    = [0.0 1.2];
            %cfg.opacitymap    = 'rampup';
            
            ft_sourceplot(cfg, MMN_int);
            
            %savefig(sprintf('%s%02.f_%s_%s_%s',fig_dir,subjectXY,bl, feature,'slice'));
            print(sprintf('%s%02.f_%s_%s_%s_test3',fig_dir,subjectXY,bl, feature,'ortho'), '-dpdf');
            
            % Project volumetric data to an MNI white-matter surface surface
            
%             cfg = [];
%             cfg.method         = 'surface';
%             cfg.funparameter   = 'avg.pow';
%             cfg.maskparameter  = 'mask';
%             %cfg.funcolorlim    = [-1 1];
%             cfg.funcolormap    = 'jet';
%             %cfg.opacitylim     = [0.0 1.2];
%             cfg.opacitymap     = 'rampup';
%             cfg.projmethod     = 'project';
%             %cfg.surfinflated   = 'surface_pial_both.mat';
%             cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
%             %cfg.surfinflated   = 'surface_inflated_both.mat'; %string, file that contains the inflated surface (default = [])  may require specifying a point-matching (uninflated) surffile
%             cfg.camlight       = 'no'; % default is yes
%             %cfg.surfdownsample = 15;  % downsample to speed up processing; number (default = 1, i.e. no downsampling)
%             
%             ft_sourceplot(cfg, source_norm);
%             view ([90 0]);
%             print(sprintf('%s%02.f_%s_%s_%s',fig_dir,subjectXY,bl, feature,'surface_left'), '-dpdf');
%             
%             ft_sourceplot(cfg, source_norm);
%             view ([90 0]);
%             print(sprintf('%s%02.f_%s_%s_%s',fig_dir,subjectXY,bl, feature,'surface_right'), '-dpdf');
%             
            %save figures
            %savefig(sprintf('%s%02.f_%s_%s_%s',fig_dir,subjectXY,bl, feature,'surface'));
            
            close all
            
            
        end
    end
    
    % save data
    save(sprintf('%s%02.f_%s_test3',out_dir,subjectXY, 'sources'),'mumufes_sources','-v7.3');
    
    
end

