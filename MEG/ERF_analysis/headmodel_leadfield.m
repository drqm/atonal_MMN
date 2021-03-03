
%% HEADMODEL WITH COMMON GRID (Giorgos)

ftpath = '/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip_git/fieldtrip/';
addpath(ftpath);
ft_defaults;


%%
% HEADMODEL

subjectID = [1];

%nonmusicians
% subjectID = [1 2 4 5 6 7 8 9 10 11 12 14 15 16 18 19 20 21 25];
% and 17 has no complete MRI and for 5 no MRI exists.

filedir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/scripts/analysis_MM/';
filename = fullfile(filedir,'Preprocessing_MEG_MM.csv');
subjects = readtable(filename);

bn = subjects.atonal(subjectID);
ID = subjects.ID(subjectID);
fname = subjects.Files(subjectID);
mri_name = subjects.MRI(subjectID); 

mri_dir = ['/hpc/users/iris.mencke/workspace/CCMusic/rawdata/', mri_name{1}];
meg_dir = '/hpc/users/iris.mencke/workspace/CCMusic/rawdata/';
out_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/source_space/headmodels/';

% read the DICOM files
  
mri_filename = [mri_dir '/0002_MR000004.dcm']; % load([ICA_dir, fname{1,1}(1:5) '_0', num2str(bn) '_meg_channels_comp']);
mri = ft_read_mri(mri_filename);
%mri.coordsys = 'ctf';

meg_data = [meg_dir, bn{1}]; %'AZA16_MPIEA0121.CCMUSIC_20190722_05.ds'
grad    = ft_read_sens(meg_data,'senstype','meg', 'unit','cm');
shape   = ft_read_headshape(meg_data,'unit','cm');

figure;
ft_plot_headshape(shape,'fidlabel','yes');
ft_plot_sens(grad, 'style', '*b');
view([1 0 0]);

%% PREPARE COMMON MI BRAIN HEADMODEL AND GRID 
% from 
% https://www.fieldtriptoolbox.org/tutorial/sourcemodel/#performing-group-analysis-on-3-dimensional-source-reconstructed-data

% NOTE: the path to the template file is user-specific

load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d8mm'));
%load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d4mm'));

template_grid = sourcemodel;
%template_grid = sourcemodel;

clear sourcemodel;

%% Realign individual MRI to MEG space

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
mri_realigned = ft_volumerealign(cfg,mri); % marking of the fiducials: bottom ear is right ear (of the upper left picture)

% segment the anatomical MRI
%cfg        = [];
%cfg.output = 'brain';
%seg        = ft_volumesegment(cfg, mri_realigned);

cfg = [];
mri_realigned.coordsys = 'ctf';
cfg.output = {'brain','skull','scalp'};
% cfg.braintreshold = 'no';
% cfg.scalpthreshold = 'no';
seg = ft_volumesegment(cfg,mri_realigned);

% 
% cfg = [];
% cfg.funparameter = 'brain';
% ft_sourceplot(cfg, seg);
% 
% cfg.funparameter = 'skull';
% ft_sourceplot(cfg, seg);
% 
% cfg.funparameter = 'scalp';
% ft_sourceplot(cfg, seg);


% construct the volume conductor model (i.e. head model) for each subject
% this is optional, and for the purpose of this tutorial only required for
% plotting, later on

cfg        = [];
cfg.method = 'singleshell';
headmodel  = ft_prepare_headmodel(cfg, seg);

cfg           = [];
cfg.warpmni   = 'yes';
cfg.template  = template_grid;
cfg.nonlinear = 'yes';
cfg.mri       = mri_realigned;
cfg.unit      ='mm';
grid          = ft_prepare_sourcemodel(cfg);

grad = ft_convert_units(grad,'mm');


grid.pos_mni = template_grid.pos;


%% PLOT

figure; hold on;

ft_plot_sens(grad, 'style', 'ob');
ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_mesh(grid.pos(grid.inside,:));

save(sprintf('%s%02.f_headmodel_meg',out_dir,subjectID), 'headmodel','grid','grad','shape','seg','-v7.3'); % all in cm; all in ctf


%% LEADFIELDS

normalWeight = 1;

cfg             = [];
cfg.channel     = 'MEG';
cfg.grad        = grad;
cfg.sourcemodel = grid;
cfg.headmodel   = headmodel;
cfg.method      = 'singleshell';
cfg.normalize = 'yes'; % normalization is necessary for beamformer;
cfg.normalizeparam  = normalWeight ;
cfg.singleshell.batchsize = 1000;
cfg.lcmv.reducerank = 2; % beamforming (lcmv)
leadfield_meg   = ft_prepare_leadfield(cfg); % NOTE: input of the whitened data ensures the correct sensor definition to be used.

save(sprintf('%s%02.f_leadfield_meg_N%1.1f.mat',out_dir,subjectID,normalWeight), 'leadfield_meg'); % all in cm; all in ctf

close all







