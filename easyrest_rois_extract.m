function easyrest_rois_extract(ER_)
%EASYREST_ROIS_EXTRACT extracts measure values in specified ROIs and store 
%data in a matlab file to be open with easyrest_rois_plot. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDITABLE FIELDS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ER.project = '';            % the complete path of the project's mat file (ER_project_name.mat)
ER.output_name = '';        % output name for the extracted data

% eclude specific measures to save time: 
ER.psd_extract = 1;
ER.psdw_extract = 1;
ER.psdm_extract = 1;
ER.gFC_HIST_extract = 1;
ER.coh_ma_extract = 1;


% ROI class (experimental)
% ER.roi_classifier{1} = {'PCC','MPFC','Broca'};      %only nii are allowed. DO NOT PUT the extension. %comment for no classification
% ER.roi_classifier{2} = {'ACC','PCC'};
% ER.roi_classifier{3} = {'V1'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW, UNLESS YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -global opt
global opt

version = version_er;      
%%---------------------------------
% Welcome screen
welcome_screen(version)
%%---------------------------------

if nargin
    fprintf('\nReading options from input variable...');
    clear ER;
    ER = ER_;
    fprintf('done.\n');
end
load (ER.project)

%%---------------------------------
% Check for update
% if ER.check_update
%     if exist('check_updates')~=2
%         error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''/home/user/MATLAB/ER''))');
%     end
%     exit = check_updates(version);
%     if exit;return;end
% end
%%---------------------------------
measure = measure_definition;
indx = arrayfun( @(s) isempty(s.do),measure);       %mi assicuro che non ci siano campi empty
measure(indx) = [];
measure(~[measure.do]) = [];    %elimino measures non eseguite

if  ~isfield(ER,'roi_classifier')
    roi_classes = 1;
    list_roi = dir([opt.folders.seed,'/*.nii']);        %TODO ADD ALSO HDR
    % add brain mask as first roi
    roi_class(1).roi(1).name = 'BrainMask';
    if opt.bm.use_ss
        warning off backtrace
        warning(sprintf('\nThis project uses subject-specific brainmask, be aware:\n1) Voxel-wise analysis is not possible in this modality.\n2) ROI-based analysis might be baised.\n3) Comparison of standardized measures might be biased.'));
        warning on backtrace
        roi_class(1).roi(1).img = 0.*uint8(opt.AUX.bmask_img{1,1}) + 1; % creating a fake brain mask composed of only ones. This will be intersected with subject specific brainmask.
        roi_class(1).roi(1).path = [];
    else
        roi_class(1).roi(1).img = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path{1,1})));
        roi_class(1).roi(1).path = opt.AUX.bmask_path{1,1};
    end
    %------------------------------
    for l = 1:length(list_roi)
        roi_class(1).roi(l+1).name = list_roi(l).name(1:end-4);
        roi_class(1).roi(l+1).img = spm_read_vols(spm_vol([opt.folders.seed,'/',list_roi(l).name]));
        roi_class(1).roi(l+1).path = [opt.folders.seed,'/',list_roi(l).name];
    end
    roi_max = length(roi_class(1).roi);
else    %this case seems that doesnt have the brain mask...It could be the case to remove this modality
    roi_classes = length(ER.roi_classifier);
    roi_max = 0;
    for c = 1:roi_classes
        for l = 1:length(ER.roi_classifier{c})
            roi_class(c).roi(l).name = ER.roi_classifier{c}{l};
            roi_class(c).roi(l).img = spm_read_vols(spm_vol([opt.folders.seed,'/',ER.roi_classifier{c}{l},'.nii']));
            roi_class(c).roi(l).path = [opt.folders.seed,'/',ER.roi_classifier{c}{l},'.nii'];
            roi_max_temp = length(roi_class(c).roi);
            if roi_max_temp > roi_max
                roi_max = roi_max_temp;
            end
        end
    end
end

%BM.img = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));

fprintf('\nExtracting data variable\n');
reverseStr  = '';
count = 0;
total_roi = 0;
for z = 1:roi_classes
    total_roi =  length(roi_class(z).roi) + total_roi;
end
total_subj = 0;
for z = 1:opt.group_number 
    total_subj =  opt.subject_number(z)+total_subj;
end
total = total_subj*opt.session_number*length(measure)*total_roi;

data = nan(opt.group_number,max(opt.subject_number),opt.session_number,length(measure),roi_classes,roi_max);

if ~opt.bm.use_ss;
    try
        bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path{1,1})));
    catch  %for retrocompatibility
        bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    end
end

for g = 1:opt.group_number 
    for v = 1:opt.subject_number(g)
        if opt.bm.use_ss
            bmask = uint8(opt.AUX.bmask_img{g,v});
        end
        for s = 1: opt.session_number
            for m = 1: length(measure)
                vol_str = [opt.folders.firstlevel,'/',measure(m).dir_name,'/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_',measure(m).name,'.nii'];
                volume = spm_read_vols(spm_vol(vol_str));
                if ndims(volume) == 4
                    volume = squeeze(volume(:,:,:,measure(m).subbrick));
                end
                
                for r = 1: roi_classes
                    for l = 1: length(roi_class(r).roi)
                        count = count + 1;
                        
                        tmp = volume(roi_class(r).roi(l).img == 1 & bmask == 1);
                        tmp(tmp==0) = [];
                        data(g,v,s,m,r,l) = er_nanmean(tmp);
                        
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                        
                    end
                end
                
            end
        end
    end
end

subject_number = opt.subject_number;

roi_backup = roi_class;         % for restoring roi_class after cleaning for saving
for z = 1:roi_classes
    for l = 1:length(roi_class(z).roi)
        roi_class(z).roi(l).img = [];
    end
end

info_session_group.session_n = opt.session_number;
info_session_group.session_name = opt.AUX.session_names_prefix;
info_session_group.group_n = opt.group_number;
info_session_group.group_name = opt.AUX.group_names_prefix;

try 
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data','roi_class','measure','roi_class','subject_number','info_session_group','-append');
catch
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data','roi_class','measure','roi_class','subject_number','info_session_group');
end
fprintf('\n');

roi_class = roi_backup; % reload roi_class 


% DA data
if opt.MEASURES.DA % I don't use a variable to exclude it because it will be fast.
    fprintf('\nPointing to Denoising Analysis (DA) data\n');
    reverseStr  = '';
    total = 1;
    DA.path = opt.DA.extracted_data_path;
    DA.proj_name = opt.folders.prject_name;
    %DA.reg_names = opt.DA_extracted.name_contract;
    DA.folder_name = opt.folders.denoising_analysis;
        percentDone = 100 * 1 / total;
        msg = sprintf('Percent done: %3.1f', percentDone);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'DA','-append');
end


%PSD data
if opt.MEASURES.PSD && ER.psd_extract
    fprintf('\nExtracting data_psd variable\n');
    tr = opt.tr;
    filter_band = opt.filter_band;
    nfft = opt.AUX.psd_nfft;
    reverseStr  = '';
    count = 0;
    total = total_subj*opt.session_number*total_roi;
    data_ps = nan(opt.group_number,max(opt.subject_number),opt.session_number,roi_classes,roi_max,nfft/2);
    for g = 1:opt.group_number 
        for v = 1:opt.subject_number(g)
            if opt.bm.use_ss
                bmask = uint8(opt.AUX.bmask_img{g,v});
            end
            for s = 1: opt.session_number
                vol_str = [opt.folders.firstlevel,'/PSD/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_PSD.nii'];
                volume = spm_read_vols(spm_vol(vol_str));
                dime = size(volume);
                for r = 1: roi_classes
                    for l = 1: length(roi_class(r).roi)
                        media = [];
                        for k=1:dime(4)
                            tmp = volume(:,:,:,k);
                            media(k) = mean(tmp(roi_class(r).roi(l).img == 1 & bmask == 1));
                        end
                        %[~,results] = system(['3dROIstats -quiet -mask ',roi_class(r).roi(l).path,' ',vol_str]);
                        %data_ps(g,v,s,r,l,:) = str2num(results);
                        data_ps(g,v,s,r,l,:) = media;
                        count = count + 1;
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data_ps','tr','filter_band','nfft','-append');
end

% GFC HISTOGRAM
if opt.MEASURES.gFC && ER.gFC_HIST_extract
    fprintf('\nExtracting data_gFC_HIST variable\n');
    reverseStr  = '';
    count = 0;
    total = total_subj*opt.session_number*total_roi;
    data_gFC_HIST = nan(opt.group_number,max(opt.subject_number),opt.session_number,roi_classes,roi_max,200);
    for g = 1:opt.group_number 
        for v = 1:opt.subject_number(g)
            if opt.bm.use_ss
                bmask = uint8(opt.AUX.bmask_img{g,v});
            end
            for s = 1: opt.session_number
                vol_str = [opt.folders.firstlevel,'/gFC/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_gFC_HIST.nii'];
                volume = spm_read_vols(spm_vol(vol_str));
                dime = size(volume);
                for r = 1: roi_classes
                    for l = 1: length(roi_class(r).roi)
                        media = [];
                        for k=1:dime(4)
                            tmp = volume(:,:,:,k);
                            media(k) = mean(tmp(roi_class(r).roi(l).img == 1 & bmask == 1));
                        end
                        %[~,results] = system(['3dROIstats -quiet -mask ',roi_class(r).roi(l).path,' ',vol_str]);
                        %data_gFC_HIST(g,v,s,r,l,:) = str2num(results);
                        data_gFC_HIST(g,v,s,r,l,:) = media;
                        count = count + 1;
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data_gFC_HIST','-append');
end


%COH MAGNITUDE data
if opt.MEASURES.COH && ER.coh_ma_extract
    if ~exist('filter_band','var')
        filter_band = opt.filter_band;
    end
    if isfield(opt.coherency,'freq')
        freq_coh = opt.coherency.freq;
    else   %for Backward compatibility
        if exist('tr','var')
            tr = opt.tr;
        end
    end
    seed_number = size(opt.AUX.COH.roi_names,1);
    name_coh_ma = cell(seed_number,1);
    fprintf('\nExtracting data_coh_ma (coherency''s magnitude) variable\n');
    reverseStr  = '';
    count = 0;
    total = total_subj*opt.session_number*total_roi*seed_number;
    data_coh_ma = nan(opt.group_number,max(opt.subject_number),opt.session_number,seed_number,roi_classes,roi_max,length(freq_coh));
    for g = 1:opt.group_number 
        for v = 1:opt.subject_number(g)
            if opt.bm.use_ss
                bmask = uint8(opt.AUX.bmask_img{g,v});
            end
            for s = 1: opt.session_number
                for m = 1:seed_number
                    seed_name = c_b(opt.AUX.COH.roi_names(m,:));
                    seed_name = seed_name(1:end-4);
                    name_coh_ma{m} = seed_name;
                    vol_str = [opt.folders.firstlevel,'/COH/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_COH_',seed_name,'_MA.nii'];
                    volume = spm_read_vols(spm_vol(vol_str));
                    dime = size(volume);
                    for r = 1: roi_classes
                        for l = 1: length(roi_class(r).roi)
                            %{
                            [~,results] = system(['3dROIstats -quiet -mask ',roi_class(r).roi(l).path,' ',vol_str]);
                            str_pattern = 'errors';
                            indx = strfind(results,'errors');
                            if ~isempty(indx)
                                indx = indx + length(str_pattern) -1;
                                results(1:indx) = [];
                            end
                            data_coh_ma(g,v,s,m,r,l,:) = str2num(results);
                            %}
                            media = [];
                            for k=1:dime(4)
                                tmp = volume(:,:,:,k);
                                media(k) = mean(tmp(roi_class(r).roi(l).img == 1 & bmask == 1));
                            end
                            data_coh_ma(g,v,s,m,r,l,:) = media;
                            count = count + 1;
                            percentDone = 100 * count / total;
                            msg = sprintf('Percent done: %3.1f', percentDone);
                            fprintf([reverseStr, msg]);
                            reverseStr = repmat(sprintf('\b'), 1, length(msg));
                        end
                    end
                end
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data_coh_ma','name_coh_ma','freq_coh','filter_band','-append');
end


%PSDw data
if opt.MEASURES.PSDW && ER.psdw_extract
    fprintf('\nExtracting data_psd (welch''s method) variable\n');
    freqw = opt.psdw.freq ;
    if ~exist('filter_band','var')
        filter_band = opt.filter_band;
    end
    reverseStr  = '';
    count = 0;
    total = total_subj*opt.session_number*total_roi;
    data_psw = nan(opt.group_number,max(opt.subject_number),opt.session_number,roi_classes,roi_max,length(freqw));
    for g = 1:opt.group_number 
        for v = 1:opt.subject_number(g)
            if opt.bm.use_ss
                bmask = uint8(opt.AUX.bmask_img{g,v});
            end
            for s = 1: opt.session_number
                vol_str = [opt.folders.firstlevel,'/PSDW/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_PSDW.nii'];
                volume = spm_read_vols(spm_vol(vol_str));
                dime = size(volume);
                for r = 1: roi_classes
                    for l = 1: length(roi_class(r).roi)
                        media = [];
                        for k=1:dime(4)
                            tmp = volume(:,:,:,k);
                            media(k) = mean(tmp(roi_class(r).roi(l).img == 1 & bmask == 1));
                        end
                        data_psw(g,v,s,r,l,:) = media;
                        %[~,results] = system(['3dROIstats -quiet -mask ',roi_class(r).roi(l).path,' ',vol_str]);
                        %data_psw(g,v,s,r,l,:) = str2num(results);
                        count = count + 1;
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data_psw','freqw','filter_band','-append');
end

%PSDm data

if opt.MEASURES.PSDM && ER.psdm_extract
    freqm = opt.psdm.freq ;
    if ~exist('filter_band','var')
        filter_band = opt.filter_band;
    end
    fprintf('\nExtracting data_psd (multitaper method) variable\n');
    reverseStr  = '';
    count = 0;
    total = total_subj*opt.session_number*total_roi;
    data_psm = nan(opt.group_number,max(opt.subject_number),opt.session_number,roi_classes,roi_max,length(freqm));
    for g = 1:opt.group_number 
        for v = 1:opt.subject_number(g)
            if opt.bm.use_ss
                bmask = uint8(opt.AUX.bmask_img{g,v});
            end
            for s = 1: opt.session_number
                vol_str = [opt.folders.firstlevel,'/PSDM/',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,v},opt.AUX.session_names_prefix{s},'_PSDM.nii'];
                volume = spm_read_vols(spm_vol(vol_str));
                dime = size(volume);
                for r = 1: roi_classes
                    for l = 1: length(roi_class(r).roi)
                        media = [];
                        for k=1:dime(4)
                            tmp = volume(:,:,:,k);
                            media(k) = mean(tmp(roi_class(r).roi(l).img == 1 & bmask == 1));
                        end
                        data_psm(g,v,s,r,l,:) = media;
                        %[~,results] = system(['3dROIstats -quiet -mask ',roi_class(r).roi(l).path,' ',vol_str]);
                        %data_psm(g,v,s,r,l,:) = str2num(results);
                        count = count + 1;
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_ROIs_',ER.output_name,'.mat'],'data_psm','freqm','filter_band','-append');
end
        
fprintf('\n\nRoi extraction complete. Bye!\n\n');

return
end

function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end