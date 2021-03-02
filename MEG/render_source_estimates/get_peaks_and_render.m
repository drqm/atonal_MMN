ft_defaults;

vol_path = 'volumes';
surface = 'cortex_20484.surf.gii'; %Higher res, slower plotting
M = gifti(surface);
grand_avg_file = 'grand_avg_sources_test3.mat';
source_stats_file = 'cb_permutations_sources_test3.mat';
load(grand_avg_file)
load(source_stats_file)
load('anatomy_labels')

%% First, get mask for ROIs (if not already done)
% atlas = niftiread('AAL3V1.nii');
% info = readtable('ROI_MNI_V7_vol.txt');

% cpos = grand_avg_sources.tonal.pitch.std.pos;
% origin = [45 63 36];
% cpos(:,1) = origin(1)-cpos(:,1)/2;
% cpos(:,2) = cpos(:,2)/2 + origin(2);
% cpos(:,3) = cpos(:,3)/2 + origin(3);
% cpos = cpos + 1;
% labels = cell(length(cpos),1);
% codes = zeros(length(cpos),1);
% origin = [45 63 36];
% for r = 1:length(cpos)
%     codes(r) = atlas(cpos(r,1),cpos(r,2),cpos(r,3));
%     if codes(r) == 0
%         labels(r) = {'NaN'};
%     else
%         labels(r) = info.nom_c(info.color == codes(r));
%     end  
% end
% save('anatomy_labels','codes','labels')
%%
%%
ROIs = {'HESCHLL','HESCHLR','T1L','T1R','T1AL','T1AR','T2L','T2AR',...
    'T3L','T3R','GSML','GSMR','F3OPL','F3OPR','F3TL','F3TR',...
    'ORL','ORR'};
ROIidx = ismember(labels,ROIs);
hem_idx = [1,-1];
cnames = {'tonal','atonal'}; %ieldnames(grand_avg_sources);
hems = {'right','left'};
cond = {};
feat = {};
pow ={};
hem ={};
x = {};
y = {};
z = {};
par = 'MMN_diff';
for c = 1:length(cnames)
    cname = cnames{c}; % current condition name
    ccond = grand_avg_sources.(cname); % current condition data
    fnames = fieldnames(ccond); % feature names
    for f = 1:length(fnames)
        fname = fnames{f}; % current feature name
        cfeat = ccond.(fname).(par); % current feature data
        for h = 1:length(hems)
            hidx = cfeat.pos(:,1)*hem_idx(h) > 0 & ROIidx == 1;
            cpos = cfeat.pos(hidx,:);
            cpow = cfeat.pow(hidx);
            [maxp, max_idx] = max(cpow);
            cond{end+1} = cname;
            feat{end+1} = fname;
            hem{end+1} = hems{h};
            pow{end+1} = maxp;
            x{end+1} = cpos(max_idx,1);
            y{end+1} = cpos(max_idx,2);
            z{end+1} = cpos(max_idx,3);
        end
    end
end
peaks = table(cond',feat',hem',pow',x',y',z');
peaks.Properties.VariableNames = {'condition','feature','hemisphere','power','x','y','z'};
writetable(peaks, 'peak_coords.csv')
%% Convert sources to nifti

cfg = [];
cfg.filetype = 'nifti';
cfg.parameter = 'pow';
par = 'MMN_diff';
for c = 1:length(cnames)
    cname = cnames{c}; % current condition name
    ccond = grand_avg_sources.(cname); % current condition data
    fnames = fieldnames(ccond); % feature names
    for f = 1:length(fnames)
        fname = fnames{f}; % current feature name
        cfeat = ccond.(fname).(par); % current feature data
        clustermask = int8(results_sources.MMN.(cname).(fname).negclusterslabelmat == 1);
        min(results_sources.MMN.(cname).(fname).stat)
        clustermask = clustermask.*int8(results_sources.MMN.(cname).(fname).stat <= -3);
        clusterstat = results_sources.MMN.(cname).(fname).stat*-1;
        cmask = cfeat;
        cstat = cfeat;
        cmask.pow = clustermask;
        cstat.pow = clusterstat;
        cfg.filename = [vol_path, '/', cname,'_',fname];
        ft_sourcewrite(cfg,cfeat);
        cfg.filename = [vol_path, '/', cname,'_',fname,'_mask'];
        ft_sourcewrite(cfg,cmask);
        cfg.filename = [vol_path, '/', cname,'_',fname,'_stat'];
        ft_sourcewrite(cfg,cstat);
    end
end

%%
addpath C:\Users\au571303\Documents\MATLAB\spm12\toolbox\cat12
spm('defaults', 'fMRI','cat_12');
cm = colormap('parula');
%cm = repmat([1,0,0],length(cm),1);
figure;
set(gcf,'color','white')
set(gcf,'renderer','painters')
set(gcf,'PaperOrientation','portrait')
set(gcf,'PaperUnits','centimeters')
%set(gcf, 'PaperSize',[21,29.7])
set(gcf, 'PaperPosition',[0,0,18,16])
set(gcf,'Units','normalized')
set(gcf,'Position',[0.1,0.05,0.5,0.8])

views = [-90,90];
lights = {[30,0],[-30,0]};
clims = {[0,4e-24],[0,8e-24],[0,5e-24],[0,8e-24]};
tlabels = {[0,2e-24,4e-24],[0,4e-24,8e-24],[0,2.5e-24,5e-24],[0,4e-24,8e-24]};
for c = 1:length(cnames)
    cname = cnames{c}; % current condition name
    ccond = grand_avg_sources.(cname); % current condition data
    fnames = fieldnames(ccond); % feature names
    for f = 1:length(fnames)
        fname = fnames{f}; % current feature name
        %results = [vol_path, '/', cname,'_',fname,'.nii'];
        results = [vol_path, '/', cname,'_',fname,'.nii'];
        mask =  [vol_path, '/', cname,'_',fname,'_mask.nii'];
        P = spm_mesh_project(M, results);
        %         mk = spm_mesh_project(M, mask);
        %         P = P.*mk;
        for v = 1:length(views)
            Y = subplot(4,4, 2*(c-1)+ v + 4*(f-1));
            H = cat_surf_render(surface, 'Parent',Y);
            cat_surf_render('Overlay',H,P);
            cat_surf_render('ColourMap',H,cm);
            cat_surf_render('Clim',H,clims{f});
            if v == 1 & c == 1
                %cat_surf_render('ColourBar',H,'on')%{['on'],'off'});
                cbar = colorbar();
                cbar.Position = [0.92,Y.Position(2),0.02,Y.Position(4)];
                cbar.Ticks = tlabels{f};
                cbar.FontSize = 10;
                aa = annotation('textarrow',[0,0],[0,0],'String',fname,...
                    'HeadStyle','none','LineStyle', 'none',...
                    'TextRotation',90,'FontSize',20);
                set(aa,'Position',[0.08,Y.Position(2)+0.08,0.02,0.2]);
                set(aa,'VerticalAlignment','middle')
                set(aa,'HorizontalAlignment','center')
            end
            if f == 1 & v == 1
                aa = annotation('textarrow',[0,0],[0,0],'String',cname,...
                    'HeadStyle','none','LineStyle', 'none',...
                    'TextRotation',0,'FontSize',24);
                set(aa,'Position',[Y.Position(1)+0.18,Y.Position(2)+0.2,0.02,0.2]);
                set(aa,'VerticalAlignment','middle')
                set(aa,'HorizontalAlignment','center')
            end
            spm_mesh_inflate(H.patch,Inf,1);
            spm_mesh_inflate(H.patch,Inf,1);
            view(views(v),0)
            camlight(lights{v}(1),lights{v}(2))
        end
    end
end

print('MMN_sources.pdf','-dpdf','-r300')
print('MMN_sources.png','-dpng','-r300')

close all