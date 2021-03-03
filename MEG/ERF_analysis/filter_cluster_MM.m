
%% Import all processed FREE LISTENING CTF data to FieldTrip format and low-pass and high-pass filter the data

%% Settings

select_files  = 'yes'; % Show GUI and select which files to process
n_jobs = 1; % Number of simultanous jobs (for appropriate memory usage)
apply_cluster = 'no'; % Whether to apply cluster or run locally to avoid memory shortage ('yes' or 'no')

inputpath     = '/hpc/users/iris.mencke/workspace/CCMusic/rawdata'; % Path to input folder, which in itself or in its sub-folders contains the input fif-files
outputpath    = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/filtered'; % Path to output folder

files = dir([inputpath,'/*.ds']); % Find the rawdata (ds-files)

% files = {'MEL23_MPIEA0121.CCMUSIC_20190527_01.ds', 'MEL23_MPIEA0121.CCMUSIC_20190527_03.ds',...
% 'MEE23_MPIEA0121.CCMUSIC_20190605_07.ds', 'MEE23_MPIEA0121.CCMUSIC_20190605_08.ds',...
% 'EGA04_MPIEA0121.CCMUSIC_20190612_08.ds', 'EGA04_MPIEA0121.CCMUSIC_20190612_05.ds',...
files = {'VNK12_MPIEA0121.CCMUSIC_20190612_MPIEA0121.ds'} %, 'VNK12_MPIEA0121.CCMUSIC_20190612_06.ds',...
% 'LZN15_MPIEA0121.CCMUSIC_20190617_07.ds', 'LZN15_MPIEA0121.CCMUSIC_20190617_06.ds',...
% 'FET05_MPIEA0121.CCMUSIC_20190702_07.ds', 'FET05_MPIEA0121.CCMUSIC_20190702_08.ds',...
% 'CHS26_MPIEA0121.CCMUSIC_20190703_06.ds', 'CHS26_MPIEA0121.CCMUSIC_20190703_05.ds',...
% 'MSR13_MPIEA0121.CCMUSIC_20190708_05.ds', 'MSR13_MPIEA0121.CCMUSIC_20190708_06.ds',...
% 'MGT07_MPIEA0121.CCMUSIC_20190711_06.ds', 'MGT07_MPIEA0121.CCMUSIC_20190711_05.ds',...
% 'GBA26_MPIEA0121.CCMUSIC_20190715_05.ds', 'GBA26_MPIEA0121.CCMUSIC_20190715_06.ds',...
% 'MBA01_MPIEA0121.CCMUSIC_20190717_07.ds', 'MBA01_MPIEA0121.CCMUSIC_20190717_05.ds',...
% 'AZA16_MPIEA0121.CCMUSIC_20190722_05.ds','AZA16_MPIEA0121.CCMUSIC_20190722_07.ds',...
% 'SAA08_MPIEA0121.CCMUSIC_20190724_07.ds','SAA08_MPIEA0121.CCMUSIC_20190724_05.ds',...
% 'BGA10_MPIEA0121.CCMUSIC_20190726_05.ds', 'BGA10_MPIEA0121.CCMUSIC_20190726_06.ds',...
% 'CGA05_MPIEA0121.CCMUSIC_20190828_06.ds', 'CGA05_MPIEA0121.CCMUSIC_20190828_05.ds'; 
%  

%files = {'GBA26_MPIEA0121.CCMUSIC_20190715_07.ds', 'GBA26_MPIEA0121.CCMUSIC_20190715_06.ds'}
%files = {'VNK12_MPIEA0121.CCMUSIC_20190612_MPIEA0121.ds', 'VNK12_MPIEA0121.CCMUSIC_20190612_06.ds'}
%files = {'FAL19_MPIEA0121.CCMUSIC_20191015_MPIEA0121.ds','FAL19_MPIEA0121.CCMUSIC_20191015_05.ds'}
%exception.input = {'FAL19_MPIEA0121.CCMUSIC_20191015_MPIEA0121.ds'};
%exception.output = {'FAL19_MPIEA0121.CCMUSIC_20191015_06.ds'};

% files oct 2019

%files = {'FAL19_MPIEA0121.CCMUSIC_20191015_05.ds', 'FLA19_MPIEA0121.CCMUSIC_20191015_06.ds'}
% 'CKD23_MPIEA0121.CCMUSIC_20191016_07.ds', 'CKD23_MPI  EA0121.CCMUSIC_20191016_05.ds',...
% 'BHL10_MPIEA0121.CCMUSIC_20191017_06.ds', 'BHL10_MPIEA0121.CCMUSIC_20191017_07.ds',...
% 'GBA27_MPIEA0121.CCMUSIC_20191023_05.ds', 'GBA27_MPIEA0121.CCMUSIC_20191023_06.ds',...
% 'GCE26_MPIEA0121.CCMUSIC_20191024_06.ds', 'GCE26_MPIEA0121.CCMUSIC_20191024_05.ds';

% files nov/dec 2019

%files = {'KWO21_MPIEA0121.CCMUSIC_20191127_08.ds', 'KWO21_MPIEA0121.CCMUSIC_20191127_07.ds',...
%'XFI02_MPIEA0121.CCMUSIC_20191127_05.ds', 'XFI02_MPIEA0121.CCMUSIC_20191127_06.ds'};
%'GWN27_MPIEA0121.CCMUSIC_20191018_08.ds', 'GWN27_MPIEA0121.CCMUSIC_20191018_06.ds',...
% 'DLA24_MPIEA0121.CCMUSIC_20191126_05.ds', 'DLA24_MPIEA0121.CCMUSIC_20191126_06.ds'}

% files = {'MSA09_MPIEA0121.CCMUSIC_20191128_08.ds', 'MSA09_MPIEA0121.CCMUSIC_20191128_07.ds',...
% 'EAE03_MPIEA0121.CCMUSIC_20191128_06.ds', 'EAE03_MPIEA0121.CCMUSIC_20191128_07.ds',...
% 'RRE13_MPIEA0121.CCMUSIC_20191202_06.ds', 'RRE13_MPIEA0121.CCMUSIC_20191202_05.ds',...
% 'YBE17_MPIEA0121.CCMUSIC_20191206_07.ds', 'YBE17_MPIEA0121.CCMUSIC_20191206_05.ds'}

% files = {'expMOR08_MPIEA0121.CCMUSIC_20200207_07.ds', 'expMOR08_MPIEA0121.CCMUSIC_20200207_07.ds'};
% 
 exception.input = {};
 exception.output = {};

% files = {'IZE09_MPIEA0121.CCMUSIC_20200310_05.ds', 'IZE09_MPIEA0121.CCMUSIC_20200310_06.ds',...
%     'KGA29_MPIEA0121.CCMUSIC_20200306_05.ds', 'KGA29_MPIEA0121.CCMUSIC_20200306_06.ds',...
%     'GKA01_MPIEA0121.CCMUSIC_20200312_05.ds', 'GKA01_MPIEA0121.CCMUSIC_20200312_06.ds'};

include_dirs = true; % Set to true for CTF folders

% EXCEPTION for MOR08
%characters = [1:6,34:35]; 
characters = [1:6];
% ORIGINAL:
%characters = [1:6,34:35]; % Character indexes to include in out filename (for all 1:end-4 )
output_suffix = '05'; % Add this suffix to the output files
skip_exist = 'yes'; % Skip if filtered file already exists ('yes' or 'no')

resample = [250]; % Target sample rate in samples per second (keep empty if no resampling should be applied)
resample_first = true; % For processing data that was not already downsampled

cfg = [];
cfg.channel = {'MEG'}; % Configure the channels

cfg.lpfilter = 'yes'; % Chose low pass filter
cfg.lpfreq = 30; % Half-cutoff frequency (DEFINE)
cfg.lpfilttype = 'fir'; % Filter type
cfg.lpfiltdir = 'twopass'; % Filter direction

cfg.hpfilter = 'yes'; % Chose high pass filter 
cfg.hpfreq = 0.5; % Half-cutoff frequency (DEFINE)
cfg.hpfilttype = 'fir'; % Filter type
cfg.hpfiltdir = 'twopass'; % Filter direction


%% Prepare the jobs

if exist([outputpath,'/filter.log'],'file')==2; delete([outputpath,'/filter.log']); end % Delete old log file
diary([outputpath,'/filter.log']); % Create log file
diary on; % Record log
cfg, resample % Record configuration

% if strcmp(apply_cluster,'yes')
%     clusterconfig('scheduler','cluster');
% else
%     clusterconfig('scheduler',[]);
% end
% clusterconfig('wait',1);
% clusterconfig('long_running',1);

if ~include_dirs
    % Delete any table entry rows that are directories
    dirs = [];
    for i=1:length(files)
        if files(i).isdir==true; dirs = [dirs i]; end
    end
    files(dirs) = [];
end

% Select which files to process
if strcmp(select_files,'yes')
   inputfile = {};
%    for i=1:length(files)
%       inputfile{i} = files(i).name;
%    end
   inputfile = files;
   selection = listdlg('PromptString','Select files to process:','ListSize',[300 300],'SelectionMode','Multiple','ListString',inputfile);
   clear('inputfile');
   fprintf(['\nProcessing selected ',num2str(length(selection)),' files:\n\n'])
   for i=1:length(selection)
       % disp(files(selection(i)).name)
       disp(files{selection(i)})
   end
   disp(' ')
   files(~ismember(1:length(files),selection)) = [];
else
   fprintf('\n\nProcessing all files.\n\n')
end

%% 

%% Submit the jobs to the cluster

n = 1;
while n<=length(files)
    valid_jobs_filter = n:n+n_jobs-1<=length(files);
    next_jobs = n:n+n_jobs-1;
    next_jobs = next_jobs(valid_jobs_filter);
    filter_cfg = {};
    for i=1:length(next_jobs);
        % disp(['Submitting cluster job for file ''',files(next_jobs(i)).name,''''])
        disp(['Submitting cluster job for file ''',files{next_jobs(i)},''''])
        filter_cfg(i).files = files{next_jobs(i)}; 
        filter_cfg(i).inputpath = inputpath; 
        filter_cfg(i).outputpath = outputpath; 
        filter_cfg(i).characters = characters; 
        filter_cfg(i).output_suffix = output_suffix;
        filter_cfg(i).skip_exist = skip_exist;
        filter_cfg(i).cfg = cfg;
        filter_cfg(i).resample = resample;
        filter_cfg(i).resample_first = resample_first;
        filter_cfg(i).exception = exception;
    end
    try
%         job_id = job2cluster(@filter_cluster_job,filter_cfg);
            filter_cluster_job(filter_cfg)
    catch err
        warning(err.message)
    end
    n = n+n_jobs;
end

diary off; % Stop recording log
fprintf('Done.\n')