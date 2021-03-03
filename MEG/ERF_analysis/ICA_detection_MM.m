%% Run ICA on all the filtered MEG data

addpath('/mnt/beegfs/users/iris.mencke/workspace/CCMusic/Mumufe/filtered');

%% Settings

select_files        = 'yes'; % Show GUI and select which files to process
search_string       = '.mat'; % Process files containing this suffix and last name

method              = 'runica'; % This is the default and uses the implementation from EEGLAB
max_MEG_components  = 64; % Number of MEG components after PCA

%% Run ICA

if exist('ICA_detection.log','file')==2; delete('ICA_detection.log'); end % Delete old log file
diary('ICA_detection.log'); % Create log file
diary on; % Start recording log

MAT = dir(['*',search_string]); % Find the processed mat files

if strcmp(select_files,'yes')
   inputfile = {};
   for i=1:size(MAT,1)
      [~,inputfile{i},~] = fileparts(MAT(i).name);
   end
   selection = listdlg('PromptString','Select files to create:','ListSize',[300 300],'SelectionMode','Multiple','ListString',inputfile);
   MAT(~ismember(1:size(MAT,1),selection)) = [];
   inputfile = inputfile(selection);
   clear('selection');
   fprintf(['\nProcessing ',num2str(length(inputfile)),' files:\n\n'])
   for i=1:length(inputfile)
       disp(cell2mat(inputfile(i)))
   end
   disp(' ')
else
   fprintf('\n\nProcessing all files.\n\n')
end

n = 0; % Reset progress counter
w = waitbar(n/(size(MAT,1)),'Running ICA'); % Show progress
for i=1:size(MAT,1)
    tic; % Reset processing time
    load(MAT(i(1)).name);
    % Find the MEG channels
    cfg = [];
    cfg.method = method;
    cfg.channel = 'MEG';
    cfg.runica.pca = max_MEG_components; % Limit number of components after PCA to max.
    meg_channels_comp = ft_componentanalysis(cfg,data); % Run ICA
    save([MAT(i(1)).name(1:end-4),'_meg_channels_comp'],'meg_channels_comp'); % Save the ICA data
    clear('meg_channels_comp'); % Remove the processed data from memory
    clear('data'); % Remove the data from memory
    n = n + 1; % Add 1 to progress counter
    pt = toc; % Measure processing time
    waitbar(n/(size(MAT,1)),w,['Running ICA. ETA: ',sprintf('%0.2f',(size(MAT,1)-n)*pt/3600),' hours']); % Create progress bar
end
diary off; % Stop recording log
close(w);
fprintf('Done.\n')
