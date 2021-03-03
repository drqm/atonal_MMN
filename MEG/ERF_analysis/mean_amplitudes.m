clear all

%%
% This script extracts the mean gradient amplitudes from the ERFs and time-
%windows of relevance.

%% Find the channels with the largest M50 response. Run only once.

load('/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/results/data/grand_average.mat')
channs = cellfun(@(x) ~strcmp('1',x(end)),grand_avg.tonal.pitch.all.label);
grad_labs = grand_avg.tonal.pitch.all.label(channs);
grand_grand = grand_avg.tonal.standard.all.avg(channs,:);
timewindow = find(grand_avg.tonal.standard.all.time >= 0.04 & grand_avg.tonal.standard.all.time <= 0.09);
[ordered, order] = sort(mean(grand_grand(:,timewindow),2));
ordered = num2cell(ordered);
M50 = [grad_labs(order),ordered]; % an ordered list is produced.
                                  % check visually the 4 channels with
                                  % the largest peak in each hemisphere
                                                   
%%
format bank % avoid scientific notation
subjects = [1 2 4 5 6 7 9 10 11 12 14 15 16 17 18 19 20 21 25]; % select the subjects
%subjects = 1:50;
sampling_rate = 250;

input_dir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/timelocked/data/';
filedir = '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/scripts/analysis_MM/';
sub_data_file = fullfile(filedir,'Preprocessing_MEG_MM.csv');

sub_data = readtable(sub_data_file,'Delimiter',';');
sub_data = sub_data(ismember(sub_data.ID,subjects),:);

sens = [];

sens.right_ant = {'MRT22', 'MRT33', 'MRT23', 'MRT32', 'MRT12', 'MRT42'};               %'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'}; % channels found in the previous section.
sens.left_ant = {'MLT22', 'MLT33', 'MLT23', 'MLT32', 'MLT12', 'MLT42'};                  %'MLO13','MLO14', 'MLO24', 'MLP53', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT14', 'MLT15', 'MLT16', 'MLT24', 'MLT25', 'MLT26', 'MLT27'};  % channels found in the previous section.

sens.right_post = {'MRP56', 'MRP57', 'MRP55', 'MRT15', 'MRT16', 'MRT14'};                       % 'MRO13','MRO14', 'MRO24', 'MRP53', 'MRP54', 'MRP55', 'MRP56', 'MRP57', 'MRT14', 'MRT15', 'MRT16', 'MRT24', 'MRT25', 'MRT26', 'MRT27'}; % channels found in the previous section.
sens.left_post = {'MLO14','MLP55', 'MLP56', 'MLT16', 'MLT27', 'MLT15' };                         %'MLO13','MLO14', 'MLO24', 'MLP53', 'MLP54', 'MLP55', 'MLP56', 'MLP57', 'MLT14', 'MLT15', 'MLT16', 'MLT24', 'MLT25', 'MLT26', 'MLT27'};  % channels found in the previous section.

flip = {[1,-1],[-1,1]};

hems = {'right','left'};
time_windows = {[70,300]};
time_around_peak = 25; %e.g. 25ms before and 25ms after

for s = 1:length(subjects)
    subject = subjects(s);
    
    fprintf('loading subject %d \n',subject)
    load(sprintf([input_dir, '%02.f.mat'],subject))
    conds = {'tonal', 'atonal'}; %select MMN conditions
   
    for l = 1:length(time_windows)
        interval = time_windows{l};
        for h = 1:length(hems)
            hem = hems{h};
            
            for c = 1:length(conds)
                cond = conds{c};
                feats = fieldnames(mumufes.(cond));
                feats = feats(~cellfun(@isempty,strfind(feats, '_MMN'))); %select MMN conditions
                
                for ff = 1:length(feats)
                    feat = feats{ff};
                    col_name = [hems{h},'_',feat,'_',cond];
                    if s == 1
                        sub_data.(col_name) = zeros(length(subjects),1);
                    end
                    
                    dataset = mumufes.(cond).(feat);
                    [check1,gradn_indx_ant] = ismember(sens.([hem,'_ant']),dataset.label);
                    [check1,gradn_indx_post] = ismember(sens.([hem,'_post']),dataset.label);
                    
                    ERF_grad_ant = mean(dataset.avg(gradn_indx_ant,:))*flip{h}(1,1);
                    ERF_grad_post = mean(dataset.avg(gradn_indx_post,:))*flip{h}(1,2);
                    
                    ERF_grad = mean([ERF_grad_ant; ERF_grad_post]);
                    
                    [peaks_grad, samples_grad] = findpeaks(ERF_grad(dataset.time >= (interval(1)/1000) &...
                        dataset.time <= (interval(2)/1000)));
                    
                    if isempty(peaks_grad)
                        [peaks_grad, samples_grad] = max(ERF_grad(dataset.time >= (interval(1)/1000) &...
                            dataset.time <= (interval(2)/1000)));
                    end
                    
                    [peak_grad, indx_grad] = max(peaks_grad);
                    sample_grad = samples_grad(indx_grad);
                    lat = dataset.time(sample_grad + min(find(dataset.time >= (interval(1)/1000)))...
                            - 1);
                    ERF_samples = (sample_grad + min(find(dataset.time >= (interval(1)/1000)))...
                            - 1 - ceil(time_around_peak*sampling_rate/1000)):...
                            (sample_grad + min(find(dataset.time >= (interval(1)/1000))) -...
                            1 + floor(time_around_peak*sampling_rate/1000));
                    
                     MA =  round((mean(ERF_grad(ERF_samples))/1e-15),2)*-1;
                     sub_data.(col_name)(sub_data.ID == subject) = MA;
                     sub_data.([col_name '_lat'])(sub_data.ID == subject) = lat*1000;

                end
            end
        end
    end
end

sub_data.bad_channels_atonal = [];
sub_data.bad_channels_tonal = [];
sub_data.MRI = [];
sub_data(:,1:5) = [];

writetable(sub_data, '/hpc/users/iris.mencke/workspace/CCMusic/Mumufe/scripts/analysis_MM/Peak_MGAs_four_channels.csv')