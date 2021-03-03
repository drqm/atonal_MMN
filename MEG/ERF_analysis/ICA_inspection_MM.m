%% Visual inspection of ICA

%% Settings

% Basic settings
search_string      = '*.mat'; % Find the input files by using the defined search string (e.g. '*.mat')
exclude_endings    = {'_eventlist','_ICA_rej'}; % Exclude any files with following name endings (e.g. {'_eventlist', '_ICA_rej'})
select_files       = 'yes'; % Show GUI and select which files to process

% Visualization settings
maximize           = 'yes'; % Maximize component window ('yes' or 'no')
width_factor       = 0.80;
height_factor      = 0.60;
showEOG            = 'no'; % Show data recorded with EOG electrodes? ('no' or 'yes')
showallcomp        = 'no'; % Show waveforms for all components? ('no' or 'yes')
showhighampl       = 'no'; % Show 10 components with highest amplitudes? ('no' or 'yes')
swap_x_y           = 'no';
reverse_axes       = 'no';
layout_file        = '/hpc/users/iris.mencke/workspace/toolboxes/matlab/fieldtrip/template/layout/CTF275.lay'; % Specify file with layout. Use empty [] if layout should be retrieved from data.

% Advanced settings
% (Developed 2016 for MiB by Niels Trusbak Haumann)

detect_artifacts   = 'yes'; % Automatically detect artifact components? ('yes' or 'no')

% EOG
artifact.EOG.eeg_amplitude     = 25; % Min. EEG amplitude
artifact.EOG.meg_amplitude     = 500; % Min. magnetometer amplitude
artifact.EOG.grad_amplitude    = 100; % Min. gradiometer amplitude
artifact.EOG.frequency         = [1 4]; % Frequency range in Hz
artifact.EOG.scalp             = 'frontal'; % Scalp distribution
artifact.EOG.min_duration      = 0.1; % Min. part of complete data with artifact in seconds
artifact.EOG.ioi               = []; % Inter-onset-interval range in seconds

% ECG
artifact.ECG.eeg_amplitude     = 2.5; % Min. EEG amplitude
artifact.ECG.meg_amplitude     = 50; % Min. magnetometer amplitude
artifact.ECG.grad_amplitude    = 10; % Min. gradiometer amplitude
artifact.ECG.frequency         = [18 40]; % Frequency range in Hz
artifact.ECG.scalp             = 'temporal'; % Scalp distribution
artifact.ECG.min_duration      = 0; % Min. part of complete data with artifact in seconds
artifact.ECG.ioi               = [0.5 2]; % Inter-onset-interval range in seconds
artifact.ECG.ioi_min           = 0.1; % Min. measured inter-onset-interval in seconds
artifact.ECG.ioi_hist_binsize  = 0.1; % Bin size of measured inter-onset-intervals in seconds
artifact.ECG.ioi_hist_freq_min = 2; % Minimum occurences of measured inter-onset-intervals

% Scaling and units
scalingmeg.factor              = 1e+15; % Scaling factor for the magnetometer type (1 = no rescaling. E.g. 1e15 to rescale T to fT)
scalingmeg.unit                = 'fT'; % Unit for the mangnetometer type (e.g. 'fT')
scalinggrad.factor             = 1e+13; % Scaling factor for the gradiometer type (1 = no rescaling. E.g. 1e+13 to rescale T/m to fT/cm)
scalinggrad.unit               = 'fT/cm'; % Unit for the gradiometer type (e.g. 'fT/cm')
scalingeeg.factor              = 1e+6; % Scaling factor for the EEG type (1 = no rescaling. E.g. 1e+6 to rescale V to ?V)
scalingeeg.unit                = '?V'; % Unit for the gradiometer type (e.g. '?')


%% Select the data and sensortype

files = dir(search_string);

% Delete any table entry rows that are directories or components
dirs = []; comps = []; endings = [];
for i=1:length(files)
    if files(i).isdir==true; dirs = [dirs i]; end
    if length(files(i).name)>8
        if strcmp(files(i).name(end-8:end),'_channels_comp.mat')
            comps = [comps i]; 
        end
    end
    for j=1:length(exclude_endings)
        if length(files(i).name)>length(exclude_endings{j})+4
           if strcmp(files(i).name(end-length(exclude_endings{j})-3:end),[exclude_endings{j},'.mat'])
               endings = [endings i]; 
           end
        end
    end
end
files([dirs, comps, endings]) = [];

% Select which files to process
if strcmp(select_files,'yes')
   inputfile = {};
   for i=1:length(files)
      inputfile{i} = files(i).name;
   end
   last_selection = [];
   if exist('ICA_inspection_memory','var')
       if ~isempty(ICA_inspection_memory)
           last_selection = find(ismember(inputfile,[ICA_inspection_memory,'.mat']));
       end
   end
   if ~isempty(last_selection)
      selection = listdlg('PromptString','Select the data to show:','ListSize',[300 300],'SelectionMode','Single','ListString',inputfile,'InitialValue',last_selection);
   else
       selection = listdlg('PromptString','Select the data to show:','ListSize',[300 300],'SelectionMode','Single','ListString',inputfile);
   end
   clear('inputfile');
   [~,search_string,~] = fileparts(files(selection).name);
   ICA_inspection_memory = search_string;
end

% Select sensortype
sensortype = 'meg'; % Sensortype ('eeg', 'mag', 'grad') (DEFINE)
fprintf(['\nShowing components for ''',upper(sensortype),''' channels in ''',search_string,'''\n\n'])


%% Load the data

load([search_string,'_',sensortype,'_channels_comp']); % Load the component data
load(search_string); % Load the data
cfg = []; EOG = [];
if strcmp(showEOG,'yes'); % Check whether EOG data should be loaded
    cfg.dataset = [search_string,'.fif'];
    cfg.channel = 'EOG';
    EOG = ft_preprocessing(cfg); % Load EOG channels
end


%% Show EOG channels

if ~isempty(EOG); % Check whether EOG channels were prepared
    cfg.blocksize = 10;
    ft_databrowser(cfg,EOG);
end


%% Show component topographies and waveforms

cfg = []; % Clear the configuration
if isfield(data,'elec')
    if ~strcmp(sensortype,'eeg');
        data = rmfield(data,'elec');
    end
end
if isempty(layout_file)
    try
        cfg.layout = ft_prepare_layout([],data); % Prepare the layout for topoplot
    catch
        warning('Attempting to apply quickcap64 layout for EEG')
        load('/usr/local/common/matlab_toolbox/fieldtrip/r9093/template/layout/quickcap64')
        lay_index = [];
        for i=1:length(data.label)
            lay_index(i) = find(strcmp(upper(lay.label),upper(data.label(i))));
        end
        cfg.layout.pos = lay.pos(lay_index,:);
        cfg.layout.width = lay.width(lay_index,:);
        cfg.layout.height = lay.height(lay_index,:);
        cfg.layout.label = lay.label(lay_index,:);
    end
else
    cfg.layout = layout_file;
    cfg.layout = ft_prepare_layout(cfg); % Prepare the layout for topoplot
end
if strcmp(swap_x_y,'yes')
    cfg.layout.pos = cfg.layout.pos(:,[2 1]);
end
if strcmp(reverse_axes,'yes')
    cfg.layout.pos = -cfg.layout.pos;
end

cfg.viewmode = 'component';
cfg.continuous = 'yes';
cfg.blocksize = 10;
evalc([search_string,'_',sensortype,'_channels_comp = ',sensortype,'_channels_comp']);
eval(['ft_databrowser(cfg,',search_string,'_',sensortype,'_channels_comp)']); % Show all components
clear([search_string,'_',sensortype,'_channels_comp']);
if strcmp(maximize,'yes')
    screensize = get( groot, 'Screensize' ); % Obtain screen size information
    x0=screensize(1);y0=screensize(2);width=screensize(3);height=screensize(4); % Increase the figure size to fill the whole screen
    set(gcf,'units','points','position',[x0,y0,width*width_factor,height*height_factor]);
end


%% Detect artifact components

if strcmp(detect_artifacts,'yes')
    if isempty(layout_file)
            cfg.channel = 'MEG'; 
            cfg.layout = ft_prepare_layout(cfg,data); % Prepare the layout for scalp estimates
    else
        cfg.layout = layout_file;
        cfg.layout = ft_prepare_layout(cfg); % Prepare the layout for topoplot
    end
    if strcmp(swap_x_y,'yes')
        cfg.layout.pos = cfg.layout.pos(:,[2 1]);
    end
    if strcmp(reverse_axes,'yes')
        cfg.layout.pos = -cfg.layout.pos;
    end

    components = {};
    artifact_duration_component = [];
    disp(' ')
    artifact_types = fields(artifact); % Find the artifact types
    for i=1:length(artifact_types)

        % Create filtered versions of the component data to detect segments with each defined type of artifacts
        
        % Filter the components
        fprintf(['Detecting ',artifact_types{i},' components...'])
        bpfilt = fdesign.bandpass('N,Fc1,Fc2',25,eval(cell2mat(['artifact.',artifact_types(i),'.frequency(1)'])),eval(cell2mat(['artifact.',artifact_types(i),'.frequency(2)'])),eval([sensortype,'_channels_comp.fsample'])); % Design the band pass filter
        d = design(bpfilt); % Create bandpass filter
        component_data = filtfilt(d.Numerator,1,eval([sensortype,'_channels_comp.trial{1,1}'''])); % Bandpass filter the component data
        component_data = component_data';

        % Multiply component data by the maximum backprojection
        component_data = component_data.*repmat(eval(['max(abs(pinv(',sensortype,'_channels_comp.unmixing)''),[],2)']),1,size(component_data,2));
        
        % Find segments above amplitude threshold and optionally inter-onset-intervals
        artifact_duration_component(:,i) = sum(abs(component_data) > eval(['artifact.',artifact_types{i},'.',sensortype,'_amplitude'])/eval(['scaling',sensortype,'.factor']),2)/eval([sensortype,'_channels_comp.fsample']);
        artifact_duration_component(artifact_duration_component(:,i) < eval(['artifact.',artifact_types{i},'.min_duration']),i) = 0;
        if ~isempty(eval(['artifact.',artifact_types{i},'.ioi']))
            ioi_component = [];
            ioi_match = zeros(size(component_data,1),1);
            for j=1:size(component_data,1)
                iois = diff(find(diff(abs(component_data(j,:)) > eval(['artifact.',artifact_types{i},'.',sensortype,'_amplitude'])/eval(['scaling',sensortype,'.factor']),1,2)))/eval([sensortype,'_channels_comp.fsample']);
                iois = iois(iois >= eval(['artifact.',artifact_types{i},'.ioi_min']));
                if ~isempty(iois)
                    [frequency ioi] = hist(iois,0:eval(['artifact.',artifact_types{i},'.ioi_hist_binsize']):max(iois));
                    if sum(frequency>=eval(['artifact.',artifact_types{i},'.ioi_hist_freq_min']))>0;
                        [value index] = max(frequency);
                        ioi_component(j) = ioi(index);
                    else
                        ioi_component(j) = NaN;
                    end
                else
                    ioi_component(j) = NaN;
                end
                if ioi_component(j) >= eval(['artifact.',artifact_types{i},'.ioi(1)']) & ioi_component(j) <= eval(['artifact.',artifact_types{i},'.ioi(2)']); ioi_match(j) = 1; end
            end
            clear('iois');
        else
            ioi_match = ones(size(component_data,1),1);
        end
        
        % Find the scalp distribution for each component
        scalp = {}; % Clear scalp distribution variables
        for j=1:size(component_data,1)
            scalp{j,i}.frontal = 'no'; scalp{j,i}.parietal = 'no';
            topo = eval(['pinv(',sensortype,'_channels_comp.unmixing(j,:))']);
            chan_index = ismember(cfg.layout.label, eval([sensortype,'_channels_comp.topolabel']));
            if mean(abs(topo(find(cfg.layout.pos(chan_index,2)>0.45/3))))>mean(abs(topo(find(cfg.layout.pos(chan_index,2)<0.45/3)))); scalp{j,i}.frontal = 'yes'; scalp{j,i}.parietal = 'no'; end; % Check whether frontal amplitude is stronger than parietal
            if mean(abs(topo(find(cfg.layout.pos(chan_index,2)<-0.45/3))))>mean(abs(topo(find(cfg.layout.pos(chan_index,2)>-0.45/3)))); scalp{j,i}.parietal = 'yes'; scalp{j,i}.frontal = 'no'; end; % Check whether parietal amplitude is stronger than frontal
            if mean(abs(topo(find(abs(cfg.layout.pos(chan_index,1))>0.45/3))))>mean(abs(topo(find(abs(cfg.layout.pos(chan_index,1))<0.45/3)))); scalp{j,i}.temporal = 'yes'; scalp{j,i}.vertex = 'no'; end; % Check whether temporal amplitude is stronger than vertex
            if mean(abs(topo(find(abs(cfg.layout.pos(chan_index,1))<0.45/3))))>mean(abs(topo(find(abs(cfg.layout.pos(chan_index,1))>0.45/3)))); scalp{j,i}.vertex = 'yes'; scalp{j,i}.temporal = 'no'; end; % Check whether vertex amplitude is stronger than temporal
            clear('topo');
        end
        
        % Detect artifact components and show results
        
        scalp_match = zeros(size(component_data,1),1);
        for j=1:size(component_data,1)
            if strcmp(eval(['scalp{j,i}.',eval(['artifact.',artifact_types{i},'.scalp'])]),'yes'); scalp_match(j) = 1; end
        end
        evalc(['components.',artifact_types{i},' = find(artifact_duration_component(:,i)>0 & scalp_match>0 & ioi_match>0)']);
        [value order] = sort(eval(['artifact_duration_component(components.',artifact_types{i},',i)']),'descend');
        evalc(['components.',artifact_types{i},' = components.',artifact_types{i},'(order)']);
        
        message_to_user{i} = {};
        if length(eval(['components.',artifact_types{i}]))==0
            fprintf(' Did not detect any components')
            message_to_user{i} = ['Did not detect any ',artifact_types{i},' components'];
        elseif length(eval(['components.',artifact_types{i}]))==1
            fprintf(' Detected component ')
            message_to_user{i} = ['Detected ',artifact_types{i},' component '];
        else
            fprintf(' Detected components ')
            message_to_user{i} = ['Detected ',artifact_types{i},' components '];
        end 
        for j=1:length(eval(['components.',artifact_types{i}]))
            fprintf(num2str(eval(['components.',artifact_types{i},'(j)'])))
            message_to_user{i} = [message_to_user{i} num2str(eval(['components.',artifact_types{i},'(j)']))];
            if j~=length(eval(['components.',artifact_types{i}]))
               fprintf(', ')
               message_to_user{i} = [message_to_user{i} ', '];
            end
        end
        fprintf('.\n')
        message_to_user{i} = [message_to_user{i} '.'];
    end
    msgbox(message_to_user,[search_string,' [',upper(sensortype),']']);
end
clear('component_data');
disp(' ')


%% Calculate the amplitude peaks for each component and order the peaks

[peak,sample] = max(abs(eval([sensortype,'_channels_comp.trial{1,1}'])),[],2); % Find the peak amplitudes for each component
[peak_ordered,order_peak] = sort(peak,'descend'); sample_ordered = sample(order_peak); % Sort the peaks and samples according to amplitude
sampleblocksize = round(eval([sensortype,'_channels_comp.fsample'])*10); % Find the size of a sample block of 10 seconds


%% Show the amplitude peaks of all 64 components

if strcmp(showallcomp,'yes'); % Check whether all components should be shown
    figure;
    for i=1:64;
        samples = sample(i(1))-round(sampleblocksize/2):sample(i(1))+round(sampleblocksize/2); % Find the samples to show
        samples = samples(samples >= 1 & samples <= eval(['length(',sensortype,'_channels_comp.time{1,1})'])); % Filter out samples which exceed recording time
        component = i(1); % Find the component to show
        subplot(8,8,i(1))
        plot(eval([sensortype,'_channels_comp.time{1,1}(samples)']),eval([sensortype,'_channels_comp.trial{1,1}(component,samples)']))
        title(['Component ',num2str(component)])
    end
    screensize = get( groot, 'Screensize' ); % Obtain screen size information
    x0=screensize(1);y0=screensize(2);width=screensize(3);height=screensize(4); % Increase the figure size to fill the whole screen
    set(gcf,'units','points','position',[x0,y0,width,height]);
    figure_all = gcf;
    set(figure_all,'Name','Amplitude peaks of all 64 components');
end


%% Show 10 components with highest peak amplitudes

if strcmp(showhighampl,'yes'); % Check whether high amplitudes should be shown
    topographies = 1:2:19;
    waveforms = 2:2:20;
    cfg.comment = 'no';
    figure;
    for i=1:10;
        samples = sample_ordered(i(1))-round(sampleblocksize/2):sample_ordered(i(1))+round(sampleblocksize/2); % Find the samples to show
        samples = samples(samples >= 1 & samples <= eval(['length(',sensortype,'_channels_comp.time{1,1})'])); % Filter out samples which exceed recording time
        component = order_peak(i(1)); % Find the component to show
        cfg.component = component;
        subplot(5,4,topographies(i(1)))
        ft_topoplotIC(cfg,eval([sensortype,'_channels_comp']));
        subplot(5,4,waveforms(i(1)))
        plot(eval([sensortype,'_channels_comp.time{1,1}(samples)']),eval([sensortype,'_channels_comp.trial{1,1}(component,samples)']))
        title(['Component ',num2str(component)])
    end
    screensize = get( groot, 'Screensize' ); % Obtain screen size information
    x0=screensize(1);y0=screensize(2);width=screensize(3);height=screensize(4); % Increase the figure size to fill the whole screen
    set(gcf,'units','points','position',[x0,y0,width,height]);
    figure_sorted = gcf;
    set(figure_sorted,'Name','10 components with highest peak amplitudes');
end


%% Show in the command window the amplitude peak times of all components

for i=1:length(sample);
   peak_time = eval([sensortype,'_channels_comp.time{1,1}(sample(i(1)))']);
   peak_segment = floor(peak_time/10)+1;
   disp(['Component ',num2str(i(1)),' peaks at time = ',num2str(peak_time),' seconds (in segment ',num2str(peak_segment),') with amplitude = ',num2str(round(peak(i(1))*eval(['scaling',sensortype,'.factor']))),' ',eval(['scaling',sensortype,'.unit']),'.']) 
end

clear('data');
clear([sensortype,'_channels_comp']);
clear('EOG');
