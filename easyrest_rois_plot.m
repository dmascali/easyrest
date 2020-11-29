function easyrest_rois_plot(ER_)
%EASYREST_ROIS_PLOT plots measure values previosuly extracted with easyrest_rois_extract 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDITABLE FIELDS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ER.roi_file = '/data/danielem/visanet/ER_analysis/ER_an_ROIs.mat';  %specify the full path
% of the extracted rois values (from easyrest_rois_extract)

%%%%%%%% MEASURE SELECTOR %%%%%%%%%%%%%
% comment unwanted measure

%-------AMPLITUDE----------
% raw
ER.MEASURE{1} = 'ALFF';
ER.MEASURE{2} = 'fALFF';
ER.MEASURE{3} = 'RSFA';
ER.MEASURE{4} = 'fRSFA';
% z-score
ER.MEASURE{5} = 'zALFF';
ER.MEASURE{6} = 'zfALFF';
ER.MEASURE{7} = 'zRSFA';
ER.MEASURE{8} = 'zfRSFA';
% m-normalized
ER.MEASURE{9} = 'mALFF';
ER.MEASURE{10} = 'mRSFA';
% m-1 normalized
ER.MEASURE{11} = 'mALFF-1';
ER.MEASURE{12} = 'mRSFA-1';
%-------REHO----------------
ER.MEASURE{12} = 'REHO7';
ER.MEASURE{13} = 'REHO19';
ER.MEASURE{14} = 'REHO27';
%-------SPEN----------------
ER.MEASURE{15} = 'SpEn';
%-------GFC-----------------
ER.MEASURE{16} = 'gFC_Zr';
ER.MEASURE{17} = 'gFC_Zz';
ER.MEASURE{18} = 'gFC_P';
ER.MEASURE{19} = 'gFC_VT';
%-------DeC-----------------
ER.MEASURE{20} = 'DeC';
ER.MEASURE{21} = 'zDeC';
%-------COEN----------------
ER.MEASURE{22} = 'zCoEn';
%-------3dDELAY-------------
ER.MEASURE{23} = 'DELAY';
%-------COHERENCY-----------
ER.MEASURE{24} = 'COH_M';
ER.MEASURE{25} = 'COH_D';
%-------STV-----------------
ER.MEASURE{26} = 'StV';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ER.indx_roi = [1 9];                 % For selecting rois
ER.indx_measure_seed = [1];          % For selecting seeds
ER.indx_measure_thr = [];            % For selecting thresholds in DeC, gFC_VT

ER.subj_selection = 0;               % 1/0; 1 allow for subject selection. It requires definition of ER.sample(s).good_subjs

ER.stats.paired = 0;                 % 1/0; 1: perform a paired t-test (when it is possible). ATTENTION: Affect both BAR plot and ROI-to-ROI stats

ER.save_data = 0;                    % save current plotted data in current folder

%%%%%%%  MANAGE OUTLIERS  %%%%%%%%%%%%%%
ER.remove_outliers = 0;                 %this modality works only for data in the bar plot (no correlation, ROI to ROI, DA ecc..)
ER.remove_outliers_plot_hist = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    BAR PLOT     %%%%%%%%%%%%%%%
ER.bar_plot = 1;                               % 1/0; 1: do bar plots
ER.bar_plot_options.plot_points = 0;           % 1/0; 1: plot single subject means over istograms 
ER.bar_plot_options.plot_type = 1;             % [1,2,3]; 1: bar plot; 2: box plot; 3: line plot
ER.bar_plot_options.disable_stat = 0;          % 1/0; 1: turn off statisticall analysis
ER.bar_plot_options.regress_out_covariate = 0; % if 1 regress out covariate (Y-X*B). Use ER.bar_plot_options.cova
% ER.bar_plot_options.stim.vect = [0 0 1];
% ER.bar_plot_options.stim.col = [0.3 0.3 0.3];
% ER.bar_plot_options.stim.name = {'stim'};
ER.bar_plot_options.savedata_suffixname = '';        % suffixname for save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    DA PLOT      %%%%%%%%%%%%%%%
ER.DA_plot = 1;
ER.DA_plot_movements = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    PSD PLOT     %%%%%%%%%%%%%%%
ER.psd_plot_options.plot = 1;                % 1/0; 1: plot Power spectral density (periodogram)
ER.psdw_plot_options.plot = 1;               % 1/0; 1: plot Power spectral density (welch's method)
ER.psdm_plot_options.plot = 1;               % 1/0; 1: plot Power spectral density (multitaper method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    CORRELATION PLOT     %%%%%%%
ER.correlation_plot = 1;                     % 1/0; 1: do correlation plots
ER.correlation_partial = 1;                  % 1/0; 1: compute partial correlation (you must add a covariate matrix (.cova) in the correlation structure)
ER.correlation_only_below_thr = 0.001;       % If not empty, only correlations below this threshold will be plotted  
ER.correlation_type = 'Pearson';             % Choose among: 'Pearson', 'Kendall', 'Spearman'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  COHERENCE MAGNITUDE PLOT  %%%%%%
ER.coh_ma_plot = 1;                          % 1/0; 1: plot magitude of coherency (as a function of frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  gFC HISTOGRAM PLOT  %%%%%%%%%
ER.gFC_HIST_plot = 1;                        % 1/0; 1: plot histogram of global functional connectivity (gFC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  ROI TO ROI matrix  %%%%%%%%%%
ER.RtoR_file = '';      % Specify the full path of the correlation matrix .mat (ie., RtR_*_correlation.mat');
ER.RtoR.plot = 1;       % 1/0; Do ROI to ROI plots
ER.RtoR.ostt.do = 1;    % 1/0; Do One-sample t-tests
ER.RtoR.tstt.do = 1;    % 1/0; Do Two-sample t-tests. You need to specify the samples to use (tstt.samples)
ER.RtoR.tstt.samples = [1,2];       % Select 2 sample to be compared: the test will be first sample > second sample
ER.RtoR.corr.do = 1;    % 1/0; Do correlation plots. Use the same field in ER.sample (i.e,. corr(1).vect ...) used by the above correlation PLOT 
ER.RtoR.savedata_suffixname = '';        % name for save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ER.bar_plot_options.cova = [];      % A covariate matrix for bar plot only (rows: subjs for each sample; columns: nuisance variables)
                                    % The regression is made across samples
                                    % NaN are NOT allowed in this matrix
                                    
% DEFINE sample for comparisons:
                                    
ER.sample(1).selector = [1,2,1;1,1,2];     %i: group, j = session, k = roi_class. Add raw for merging different sample. Different roi_class must have the same number of elements (ie., rois) to be compared
ER.sample(1).name = 'C';             % Assign a name to the sample
ER.sample(1).session_average = 1;   % 1/0; if 1 different sessions will be averaged. if 0 simple merging.
ER.sample(1).corr(1).name = '';     % corr(m) computes m correlations with between measures and corr(m).vect. Give a name for each corr(m) in corr(m).name
ER.sample(1).corr(1).vect = [];
ER.sample(1).corr(1).cova = [];     % only for PARTIAL CORRELATION (also for ROItoROI). If you want to remove the effect of confounds (e.g., age, sex, score) put in the matrix corr(m).cova. cova is a mxn matrix= m runs on subjects while n on covariates. 
ER.sample(1).corr(2).name = '';
ER.sample(1).corr(2).vect = [];
ER.sample(1).corr(2).cova = [];
ER.sample(1).good_subjs = [];        % A vector of ones and zeros for selecting good subjects only. Put it empty or comment for no selection.
ER.sample(1).RtoR_cova = [];        % covariates only for ROI to ROI plots (one-sample and two-sample, for correlation use corr(l).cova)

ER.sample(2).selector = [1,2,2;1,1,1];
ER.sample(2).name = 'NC';



% Other examples:

% ER.sample(1).selector = [1,1,1;1,1,2];     %i: group, j = session, k = roi_class. Add raw for merging different sample. Different roi_class must have the same number of elements (ie., rois) to be compared
% ER.sample(1).name = 'ATT LEFT';
% 
% ER.sample(2).selector = [1,2,1;1,2,2];
% ER.sample(2).name = 'ATT RIGTH';

% ER.sample(1).selector = [1,1,1;1,2,1];     %i: group, j = session, k = roi_class. Add raw for merging different sample. Different roi_class must have the same number of elements (ie., rois) to be compared
% ER.sample(1).name = 'HEMI LEFT';
% 
% ER.sample(2).selector = [1,1,2;1,2,2];
% ER.sample(2).name = 'HEMI RIGTH';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW, UNLESS YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global sample_number sample_names temp_rois filter_band cut_subjs
global flag_type_plot flag_paired flag_plot_poits sample_number sample_names temp_rois subject_number filter_band cut_subjs RtoR corr_type do_stat plot_stim

flag_check_update = 1;      % 1/0; automatic check for update. DO NOT DISABLE (usefull only in the case the repository is down)         

version = version_er;      
%%---------------------------------
% Welcome screen
welcome_screen(version)
%%---------------------------------

%%---------------------------------
% Check for update
% if flag_check_update
%     if exist('check_updates')~=2
%         error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''/home/user/MATLAB/ER''))');
%     end
%     exit = check_updates(version);
%     if exit;return;end
% end
%%---------------------------------

if nargin
    fprintf('\nReading options from input variable...');
    clear ER;
    ER = ER_;
    fprintf('done.\n');
end
load (ER.roi_file)

%%---------------------------------
% Check for missing information due to update of ER

if ~exist('info_session_group','var')
    fprintf('\n-> The ROI file comes from an older version of ER, fixing missing information...');
    t = 1;
    while t
    [FILENAME, PATHNAME] = uigetfile('ER*.mat', 'Select the .mat of the project');
    load ([PATHNAME,FILENAME]);
    if exist('opt','var')
        t = 0;
        info_session_group.session_n = opt.session_number;
        info_session_group.session_name = opt.AUX.session_names_prefix;
        info_session_group.group_n = opt.group_number;
        info_session_group.group_name = opt.AUX.group_names_prefix;
        clear opt;
        save(ER.roi_file,'data','info_session_group','-append');
    else
        fprintf('\nWrong mat file. Please insert the .mat of the project, the one with the structure ''opt''');
    end
    end
fprintf('done!\n');
end
%%---------------------------------


flag_type_plot = ER.bar_plot_options.plot_type ;
flag_paired = ER.stats.paired ;
flag_plot_poits = ER.bar_plot_options.plot_points;


% measure selector------------------
indx = zeros(length(measure),1);
for j=1:length(measure)
    if any(strcmp(measure(j).name_selector,ER.MEASURE))
        indx(j) = 1;
    end
end
if ~isempty(ER.indx_measure_seed) 
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'StV')
            count = count + 1;
            if isempty(ER.indx_measure_seed(ER.indx_measure_seed == count))
            indx(j) = 0;
            end
        end
    end
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'ZStV')
            count = count + 1;
            if isempty(ER.indx_measure_seed(ER.indx_measure_seed == count))
            indx(j) = 0;
            end
        end
    end
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'COH_M')
            count = count + 1;
            if isempty(ER.indx_measure_seed(ER.indx_measure_seed == count))
            indx(j) = 0;
            end
        end
    end
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'COH_D')
            count = count + 1;
            if isempty(ER.indx_measure_seed(ER.indx_measure_seed == count))
            indx(j) = 0;
            end
        end
    end
end

if ~isempty(ER.indx_measure_thr) 
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'gFC_VT')
            count = count + 1;
            if isempty(ER.indx_measure_thr(ER.indx_measure_thr == count))
            indx(j) = 0;
            end
        end
    end
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'DeC')
            count = count + 1;
            if isempty(ER.indx_measure_thr(ER.indx_measure_thr == count))
            indx(j) = 0;
            end
        end
    end
    count = 0;
    for j=1:length(measure)
        if strcmp(measure(j).name_selector,'zDeC')
            count = count + 1;
            if isempty(ER.indx_measure_thr(ER.indx_measure_thr == count))
            indx(j) = 0;
            end
        end
    end
end


measure(~indx) = [];
data(:,:,:,~indx,:,:) = [];
%---------------------------------- 


% roi stringer fixer---------------   
for z=1:length(roi_class)
    roi_class(z).roi = str_fixer(roi_class(z).roi);
end
%----------------------------------   

% sample, subject number + sample names 
sample_number = length(ER.sample);
for z = 1:sample_number
    sample_names{z} = ER.sample(z).name;
end
%----------------------------------  

cut_subjs = zeros(sample_number,1);
pairedbility_check = [];

fprintf('\nBuilding samples for statistical comparisons...');

for s = 1:sample_number
    if ~isfield(ER.sample(s),'session_average') || isempty(ER.sample(s).session_average)
        ER.sample(s).session_average = 0;
    end
    ER.sample(s).war_average = 0;
    if length(unique(ER.sample(s).selector(:,1))) ~=1   % Se i soggetti vengono dallo stesso gruppo i colori tra sessioni saranno uguali. Altrimenti saranno diversi
        ER.sample(s).same_group_subject = 0;
        if ER.sample(s).session_average
            ER.sample(s).war_average = 1;
            ER.sample(s).session_average = 0;
        end
    else
        ER.sample(s).same_group_subject = 1;
    end
    temp = [];
    check_nrois = [];
    sub_merged_temp = [];
    sub_tot_temp = 0;
    n_cutted_subjs = 0;
    if ~ER.sample(s).session_average
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(2,temp,data(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),:,ER.sample(s).selector(ss,3),:)); 
            check_nrois(ss) = length(roi_class(ER.sample(s).selector(ss,3)).roi);
            sub_tot_temp = sub_tot_temp + subject_number(ER.sample(s).selector(ss,1));
            sub_merged_temp(ss) = subject_number(ER.sample(s).selector(ss,1));
            if ss == 1
                subj_merged_indx{ss} = 1: subject_number(ER.sample(s).selector(ss,1));
            else
                subj_merged_indx{ss} = subj_merged_indx{ss-1}(end) +1 : subj_merged_indx{ss-1}(end) + subject_number(ER.sample(s).selector(ss,1));
            end
        end
    else
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(7,temp,data(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),:,ER.sample(s).selector(ss,3),:)); 
            check_nrois(ss) = length(roi_class(ER.sample(s).selector(ss,3)).roi);
        end
        temp = mean(temp,7);
        sub_merged_temp(1) = subject_number(ER.sample(s).selector(ss,1));
        sub_tot_temp = sub_merged_temp(1);
        subj_merged_indx{1} = 1:sub_tot_temp;
        
    end
    if isfield(ER,'subj_selection') && ER.subj_selection && isfield(ER.sample(s),'good_subjs') && ~isempty(ER.sample(s).good_subjs)  % for cutting bad subjects
        if length(ER.sample(s).good_subjs) ~= sub_tot_temp;
            error('good_subjs selector does not match the length of data. The vector must have the same length of data after concatenation of sessions or groups.');
        end
        cut_subjs(s) = 1;
        n_cutted_subjs = sum(~ER.sample(s).good_subjs);
        temp(:,~ER.sample(s).good_subjs,:,:,:,:) = [];
        if ~ER.sample(s).session_average
            for ss = 1: size(ER.sample(s).selector,1)
                merg_tmp = ER.sample(s).good_subjs(subj_merged_indx{ss});
                n_rem = sum(~merg_tmp);
                sub_merged_temp(ss) = sub_merged_temp(ss) -n_rem;
                if ss == 1
                    subj_merged_indx{ss} = 1: (subject_number(ER.sample(s).selector(ss,1))-n_rem);
                else
                    subj_merged_indx{ss} = subj_merged_indx{ss-1}(end) +1 : (subj_merged_indx{ss-1}(end) + subject_number(ER.sample(s).selector(ss,1)) -n_rem);
                end
            end
        else
            merg_tmp = ER.sample(s).good_subjs(subj_merged_indx{1});
            n_rem = sum(~merg_tmp);
            sub_merged_temp(1) = sub_merged_temp(1) -n_rem;
            subj_merged_indx{1} = 1:sub_merged_temp(1);
        end
    else
        cut_subjs(s) = 0;
    end
    if length(unique(check_nrois)) ~= 1
        error('Roi classes cannot be merged due to different size');
    end
    ER.sample(s).data = temp;
    ER.sample(s).tot_subj = sub_tot_temp - n_cutted_subjs;
    ER.sample(s).merged_subj = sub_merged_temp;
    ER.sample(s).merged_subj_index = subj_merged_indx;
    if ~ER.sample(s).session_average
        ER.sample(s).merging_number = ss;
    else
        ER.sample(s).merging_number = 1;
    end
    pairedbility_check = [pairedbility_check,ER.sample(s).tot_subj];
end

% printing sample information
fprintf('done!\n\n')
for s = 1:sample_number
    fprintf('Sample name: %s\n',ER.sample(s).name);
    for ss = 1: size(ER.sample(s).selector,1)
	       fprintf('-> Group: %s\tSession: %s\t RoiClass: %d\n',info_session_group.group_name{ER.sample(s).selector(ss,1)}, info_session_group.session_name{ER.sample(s).selector(ss,2)},ER.sample(s).selector(ss,3));
    end
    if ER.sample(s).session_average
        ses_ave_str = 'ON';
    elseif ~ER.sample(s).session_average && ER.sample(s).same_group_subject
        ses_ave_str = 'OFF';
    elseif ER.sample(s).same_group_subject == 0
        ses_ave_str = 'NOT ALLOWED';
    end
    fprintf('   Session Average: %s\n',ses_ave_str);
    if ER.sample(s).war_average
        warning('Averaging sessions is only possible among the same group of subjects. Switching to merge mode.');
    end
    fprintf('\n');
end

warning off backtrace
if length(unique(pairedbility_check)) ~=1 % i sample hanno un numero diverso di soggetti
    if ER.stats.paired
        warning('Data are NOT compatible with PAIRED t-test. Two-sample t-test will be performed.');
    end
    flag_paired = 0;
end

if flag_paired 
    warning('Paired Two-Sample T-Test will be performed. Pay attention not mixing apples with oranges!');
end

if isfield(ER.bar_plot_options,'disable_stat') && ER.bar_plot_options.disable_stat == 1
    do_stat = 0;
    warning('T-tests have been disabled.');
else
    do_stat = 1;
end
warning on backtrace    

% stim plotting

if isfield(ER.bar_plot_options,'stim')
    plot_stim.do = 1;
    stim = ER.bar_plot_options.stim;
    u_stim = unique(stim.vect);
    if u_stim(1) == 0
        u_stim(1) = [];
    end
    n_stim = length(u_stim);
    for st = 1:n_stim
        indx_st = find(stim.vect == u_stim(st));
        continuity = [indx_st(2:end) inf] - indx_st;
        pos = [];
        wei = [];
        tmp_pos = [];
        tmp_wei = 1;
        for in = 1:length(indx_st)
            if isempty(tmp_pos)
                tmp_pos = indx_st(in);
                pos = [pos, tmp_pos];
            end
            if continuity(in) == 1
                tmp_wei  = tmp_wei  + 1;
            else
                wei = [wei, tmp_wei];
                tmp_wei = 1;
                tmp_pos = [];
            end
        end
        pos = pos-0.5;
        plot_stim.POS{st} = pos;
        plot_stim.WEI{st} = wei;
        plot_stim.COL{st} = stim.col(st,:); 
        plot_stim.NAME{st} = stim.name{st};
    end
    plot_stim.n_stim = n_stim;
else
    plot_stim.do = 0;
end


for s = 1: sample_number
    siz{s} = size(ER.sample(s).data); % 1 = group; 2 = subjs; 3 = session; 4 = measure; 5 = roi_class: 6 = rois
end

temp_rois = roi_class(ER.sample(1).selector(1,3)).roi;
if ~isempty(ER.indx_roi)
    roi_number = length(ER.indx_roi);
else
    roi_number = length(temp_rois);
    ER.indx_roi = 1:roi_number;
end

measure_number = siz{1}(4);

s_p_row = roi_number;
s_p_col = measure_number;

for z = 1:sample_number
     ER.sample(z).data = reshape(ER.sample(z).data,[siz{z}(2), siz{z}(4), siz{z}(6)]);
end

if ER.remove_outliers   %put NaN at the place of outliers
    fprintf('\nRemoving ouliers using peirce''s criteria...');
    %tot_bad_sub = cell(sample_number);
    for s=1:sample_number
        fprintf('\n---Sample: %s',sample_names{s}');
        fprintf('\nN\tMeasure\tROI\tSubjs');
        for i = 1:measure_number
            count = 0;
            for l = ER.indx_roi
                count = count +1;
                ER.sample(s).data(:,i,l) = apply_peirce(ER.sample(s).data(:,i,l));
                n_out = sum(isnan(ER.sample(s).data(:,i,l)));
                if  n_out > 0
                    indx_sub = find(isnan(ER.sample(s).data(:,i,l)));
                    fprintf('\n%d %s\t%s\t%s',n_out,[str_fixer_measure(measure(i).output_dir)],temp_rois(l).name,num2str(indx_sub'));
                    %tot_bad_sub{s} = [tot_bad_sub{s};indx_sub];
                end
            end
        end
    end
    fprintf('\n');
    if ER.remove_outliers_plot_hist % every counting assumes that original data has no nan
        figure;
        set(gcf,'name','Outliers histogram')
        subplot(3,sample_number,s)
        for s=1:sample_number
            subplot(3,sample_number,s)
            tmp = isnan(ER.sample(s).data);
            tmp = sum(sum(tmp,2),3);
            norm_count = 100*tmp/(measure_number*length(ER.indx_roi));
            bar([1:ER.sample(s).tot_subj],norm_count);
            title(['Sample: ',sample_names{s},' SUBJECT effect'])
            xlim([1 ER.sample(s).tot_subj]);
            %ylim([0 100]);
            xlabel('Subject number');
            ylabel('Normalized count of outliers (%)');
            
            subplot(3,sample_number,sample_number + s)
            tmp = isnan(ER.sample(s).data);
            tmp = sum(sum(tmp,1),3);
            norm_count = 100*tmp/(ER.sample(s).tot_subj*length(ER.indx_roi));
            bar([1:measure_number],norm_count);
            title(['Sample: ',sample_names{s},' MEASURE effect'])
            %xlim([1 ER.sample(s).tot_subj]);
            %ylim([0 100]);
            xlabel('Measure number');
            ylabel('Normalized count of outliers (%)');
            
            subplot(3,sample_number,2*sample_number + s)
            tmp = isnan(ER.sample(s).data);
            tmp = squeeze(sum(sum(tmp,1),2));
            tmp = tmp(ER.indx_roi);
            norm_count = 100*tmp/(ER.sample(s).tot_subj*measure_number);
            bar([ER.indx_roi],norm_count);
            title(['Sample: ',sample_names{s},' ROI effect'])
            %xlim([1 ER.sample(s).tot_subj]);
            %ylim([0 100]);
            xlabel('ROI number');
            ylabel('Normalized count of outliers (%)');
        end
    end
        
end

%------------------------------------------------------
%regress out covariate across samples for bar plot only
%------------------------------------------------------
if ~isfield(ER.bar_plot_options,'regress_out_covariate'); ER.bar_plot_options.regress_out_covariate = 0; end;
if ~isfield(ER.bar_plot_options,'cova'); ER.bar_plot_options.cova = []; end;

if ER.bar_plot && ER.bar_plot_options.regress_out_covariate  
   if isempty(ER.bar_plot_options.cova)
       warning off backtrace
       warning('No covariate matrix has been provided for bar plot regression. Disabling regression out of covariate across samples.');
       warning on backtrace
       ER.bar_plot_options.regress_out_covariate = 0;
   else %initialize cova
       s_cova = size(ER.bar_plot_options.cova);
       if s_cova(2) > s_cova(1)
           ER.bar_plot_options.cova = ER.bar_plot_options.cova';
       end
       if sum(cut_subjs) > 0
           %todo: to create a vector merging samples
%            ER.bar_plot_options.cova(~ER.sample(s).good_subjs,:) = [];
%            s_cova = size(ER.sample(s).corr(c).cova);
       end       
       ER.bar_plot_options.cova = detrend(ER.bar_plot_options.cova,'constant');
       %todo: enabling NaN presence?
       ER.bar_plot_options.cova_N = max(size(ER.bar_plot_options.cova));
       ER.bar_plot_options.cova = [ER.bar_plot_options.cova ones(ER.bar_plot_options.cova_N,1)]; % adding costant term
%        if flag_paired
%         % subject specific regressors
%             subjs = ER.bar_plot_options.cova_N/sample_number;
%             tmp = [];
%             for s = 1:sample_number 
%                 tmp = [tmp; 1; zeros(subjs-1,1)];
%             end
%             for ss = 1:subjs-1
%                tmp = [tmp,circshift(tmp(:,1),ss)];
%                 
%             end
%             ER.bar_plot_options.cova = [ER.bar_plot_options.cova(:,end-1) ,tmp];  %removal of the costant term
%        end       
       %create data matrix (concatenating samples)
       data = [];
       for z = 1:sample_number
           data = cat(1,data,ER.sample(z).data);
       end
       %reshape data to have a 2D matrix
       s_data = size(data); %subj,measure,roi
       data = reshape(data,[s_data(1),s_data(2)*s_data(3)]);
       beta = ER.bar_plot_options.cova\data;
       data = data-ER.bar_plot_options.cova(:,1:end-1)*beta(1:end-1,:);
       data = reshape(data,[s_data(1),s_data(2),s_data(3)]);
       %rebuilding samples
       last_index = 0;
       for z = 1:sample_number
           a = last_index+1;
           b = a + ER.sample(z).tot_subj -1;
           ER.sample(z).data_reg_out_cova = data(a:b,:,:);
           last_index = b;
       end
       clear data;
   end
end
%------------------------------------------------------
   
%------------------------------------------------------
%save data check (common to all plots)
if ~isfield(ER,'save_data'); ER.save_data  = 0; end;
%------------------------------------------------------

%------------------------------------------------------
%save data output name
if ER.save_data
    if ~isfield(ER.bar_plot_options,'savedata_suffixname'); ER.bar_plot_options.savedata_suffixname  = ''; end;% datestr(now,'dd-mmm-yyyy_HH:MM:SS'); end
%     if isempty(ER.bar_plot_options.savedata_suffixname)
%         ER.bar_plot_options.savedata_suffixname = datestr(now,'dd-mmm-yyyy_HH:MM:SS');
%     end
    [~,outp,~] = fileparts(ER.roi_file);
    output_name_barplot = ['_ER_plotted_',outp(9:end),'_',str_fixer_measure([measure.output_dir]),'_',ER.bar_plot_options.savedata_suffixname];
end
%------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Bar plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ER.bar_plot
    figure;
    if ER.bar_plot_options.regress_out_covariate
        set(gcf,'name',['Bar plot with COVA. ',datestr(now,'HH:MM:SS')]);
        out_name_cova = 'COVA';
    else
        set(gcf,'name',['Bar plot. ',datestr(now,'HH:MM:SS')]);
        out_name_cova = '';
    end
    if ER.save_data
        plotted_data = [];
        plotted_data_labels_measure = [];
        plotted_data_labels_rois = [];
        plotted_data_stats = [];
    end
    count = 0;
    for l = ER.indx_roi
        for i = 1:measure_number
            count = count +1;
            subplot(s_p_row, s_p_col, count)
            isto_cell = cell(1,sample_number);
            for z = 1:sample_number
                if ER.bar_plot_options.regress_out_covariate
                    isto_cell{z} = ER.sample(z).data_reg_out_cova(:,i,l);
                else
                    isto_cell{z} = ER.sample(z).data(:,i,l);
                end
            end
            stats = isto(isto_cell, ER.sample, sample_number, sample_names);
            if ER.save_data
                plotted_data{count} = isto_cell;
                plotted_data_labels_measure{count} = str_fixer_measure(measure(i).output_dir);
                plotted_data_labels_rois{count} =temp_rois(l).name;
                plotted_data_stats = [plotted_data_stats;stats];
            end
            if mod(count-1,s_p_col) == 0
                ylabel([temp_rois(l).name],'fontsize',9,'FontWeight','bold');
            end
            if count <= s_p_col
                title([str_fixer_measure(measure(count).output_dir)],'fontsize',9,'FontWeight','bold');
            end
        end
    end
    if ER.save_data
        save([output_name_barplot,out_name_cova,'.mat'],'plotted_data','plotted_data_labels_measure','plotted_data_labels_rois','sample_names','plotted_data_stats');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Correlation plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ER.correlation_plot || (ER.RtoR.plot && ER.RtoR.corr.do)
    
    %prepare data for correlation. This part is in common with ROI to ROI
    for s =1:sample_number
        if ~isfield(ER.sample(s),'corr')
            error('Define correlation structure in ER.sample().');
        end
        indx = arrayfun( @(z) isempty(z.vect),ER.sample(s).corr);       %mi assicuro che non ci siano campi empty
        ER.sample(s).corr(indx) = [];
        bad_indx = [];
        for c = 1:length(ER.sample(s).corr)
            if cut_subjs(s)
                ER.sample(s).corr(c).vect(~ER.sample(s).good_subjs) = [];
            end
            if ER.sample(s).tot_subj ~= length(ER.sample(s).corr(c).vect)
                warning off backtrace
                warning(['Skipping correlation plot number ',num2str(c),':not matching length.']);
                warning on backtrace  
                bad_indx = [bad_indx,c];
                continue
            end
            if isempty(ER.sample(s).corr(c).name)
                ER.sample(s).corr(c).name = 'N/D';
                warning off backtrace
                warning('Undefined name of correlation vector');
                warning on backtrace
            end
            if ~iscolumn(ER.sample(s).corr(c).vect)
                ER.sample(s).corr(c).vect = ER.sample(s).corr(c).vect';
            end
            % check for NaN in the vector
            ER.sample(s).corr(c).nan_index = isnan(ER.sample(s).corr(c).vect);
            %----------------------------
            ER.sample(s).corr(c).do_partial = 0;
            if ER.correlation_partial
                if isfield(ER.sample(s).corr(c),'cova') && ~isempty(ER.sample(s).corr(c).cova)
                   ER.sample(s).corr(c).do_partial = 1; 
                   s_cova = size(ER.sample(s).corr(c).cova);
                   if s_cova(2) > s_cova(1)
                       ER.sample(s).corr(c).cova = ER.sample(s).corr(c).cova';
                   end
                   if cut_subjs(s)
                       ER.sample(s).corr(c).cova(~ER.sample(s).good_subjs,:) = [];
                   end
                   % check for NaN in cova
                   nan_index_cova = sum(isnan(ER.sample(s).corr(c).cova),2);
                   ER.sample(s).corr(c).nan_index = logical(sum([ER.sample(s).corr(c).nan_index,nan_index_cova],2));   %updating NaN index
                   % remove NaN from cova
                   ER.sample(s).corr(c).cova(ER.sample(s).corr(c).nan_index,:) = []; %removing NaNs (bad subjs cut already done)
                   %----------------------------
                   ER.sample(s).corr(c).cova_N = max(size(ER.sample(s).corr(c).cova));
                end
            end
            % now I can remove NaN from vector (taking into account cova)
            ER.sample(s).corr(c).vect(ER.sample(s).corr(c).nan_index) = []; %removing NaNs (bad subjs cut already done)
            ER.sample(s).corr(c).N = length(ER.sample(s).corr(c).vect); %real N after bad subj and NaN removal
        end
        ER.sample(s).corr(bad_indx) = [];
        ER.sample(s).corr_number = length(ER.sample(s).corr);
    end
    
    %start plot
    if ER.correlation_plot
        corr_type = ER.correlation_type;
        for s =1:sample_number
            if ER.sample(s).corr_number == 0
                continue
            end
            for c = 1:ER.sample(s).corr_number
                count = 0;
                data_x = ER.sample(s).corr(c).vect;
                deltax = max(data_x) - min(data_x); %adding extra space at boundary
                xlimit = [(min(data_x) -5*deltax/100),(max(data_x) + 5*deltax/100)]; 
                if ER.sample(s).corr(c).do_partial && ~ER.remove_outliers
                    %compute RESx
                    cova = zscore(ER.sample(s).corr(c).cova); % non influenza le correlazioni ma permette di plottare i punti sulla stessa scala di quelli senza cova
                    co.cova = [cova ones(ER.sample(s).corr(c).cova_N,1)]; % adding costant term
                    co.pre = inv(co.cova'*co.cova)*co.cova';
                    Bx = co.pre*data_x;
                    co.RESx = data_x-co.cova(:,1:end-1)*Bx(1:end-1);
                    %co.RESx = data_x;  % <---- remove
                    co.do =1;
                else
                    co.do = 0;
                end
                if isempty(ER.correlation_only_below_thr)       %caso senza threshold
                    figure;
                    set(gcf,'name',[ER.correlation_type,' correlation Plot of ',sample_names{s},' with: ',ER.sample(s).corr(c).name,'. N=',num2str(ER.sample(s).corr(c).N)]);
                    for l = ER.indx_roi
                        for i = 1:measure_number
                            count = count +1;
                            subplot(s_p_row, s_p_col, count)
                            
                            data_y = ER.sample(s).data(:,i,l);
                            data_y(ER.sample(s).corr(c).nan_index) = [];
                            if ER.sample(s).corr(c).do_partial && ER.remove_outliers
                                correlation_outliers(data_y,data_x,ER.sample(s).corr(c).name,xlimit,ER.sample(s).corr(c).cova);
                            else
                                correlation(data_y,data_x,ER.sample(s).corr(c).name,xlimit,co);
                                pbaspect([1 1 1]);
                            end
                            if mod(count-1,s_p_col) == 0
                                ylabel([temp_rois(l).name],'fontsize',9,'FontWeight','bold');
                            end
                            if count <= s_p_col
                                title([str_fixer_measure(measure(count).output_dir)],'fontsize',9,'FontWeight','bold');
                            end
                        end
                    end
                    % to remove
%                     set(gcf,'units','centimeters','position',[0 0 20 20]);
%                     set(findall(gcf, '-property', 'FontSize'), 'FontSize', 4, 'fontWeight', 'normal')
%                     set(findobj(gca,'ylabel','text'),'FontSize',7)
%                     set(findall(gca, '-property', 'Linewidth'), 'Linewidth', 1)
% 
%                     print('-dpng','-r400',[sample_names{s},'_',ER.sample(s).corr(c).name,'.png'])
                    %
                else                                            %caso CON threshold
                    index_roi = [];
                    index_measure = [];
                    total = measure_number*length(ER.indx_roi);
                    h = waitbar(0,'Looking for significant correlations...','Name',['Correlation of ',sample_names{s},' with ',ER.sample(s).corr(c).name,'. N=',num2str(ER.sample(s).corr(c).N)]);
                    count = 0;
                    for l = ER.indx_roi
                        for i = 1:measure_number  
                            count = count +1;
                            data_y = ER.sample(s).data(:,i,l);
                            data_y(ER.sample(s).corr(c).nan_index) = [];
                            if ER.sample(s).corr(c).do_partial && ER.remove_outliers
                                outcome = below_thr_correlation_outliers(data_y,data_x,ER.sample(s).corr(c).cova,ER.correlation_only_below_thr);
                            else
                                outcome = below_thr_correlation(data_y,data_x,co,ER.correlation_only_below_thr);
                            end
                            if outcome
                                index_roi = [l,index_roi];
                                index_measure = [i,index_measure];
                            end
                            waitbar(count/ total);
                        end
                    end
                    close(h)
                    corr_n = length(index_roi);
                    if corr_n == 0; fprintf(['\nNo significant correaltion found for ',sample_names{s},' with ',ER.sample(s).corr(c).name,'\n']);continue; end
                        figure;
                        set(gcf,'name',[ER.correlation_type,' correlation Plot of ',sample_names{s},' with: ',ER.sample(s).corr(c).name,'. Threshold p <= ',num2str(ER.correlation_only_below_thr)])
                    s_p_row_c = ceil(corr_n/8);   %<- to change
                    s_p_col_c = 8;      % <- to change
                    count = 0;
                    for ll = 1:corr_n
                        count = count +1;
                        subplot(s_p_row_c, s_p_col_c, count)
                        data_y = ER.sample(s).data(:,index_measure(ll),index_roi(ll));
                        data_y(ER.sample(s).corr(c).nan_index) = [];
                        if ER.sample(s).corr(c).do_partial && ER.remove_outliers
                            correlation_outliers(data_y,data_x,ER.sample(s).corr(c).name,xlimit,ER.sample(s).corr(c).cova);
                        else
                            correlation(data_y,data_x,ER.sample(s).corr(c).name,xlimit,co);
                        end
                        ylabel([temp_rois(index_roi(ll)).name],'fontsize',9,'FontWeight','bold');
                        title([str_fixer_measure(measure(index_measure(ll)).output_dir)],'fontsize',9,'FontWeight','bold');
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DA PLOTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('DA','var') && ~isempty(DA) && ER.DA_plot
    % create samples
    for s = 1:sample_number
        selectors{s} = ER.sample(s).selector;
        if cut_subjs(s)
            good_subjs{s} = ER.sample(s).good_subjs;
        else
            good_subjs{s} = [];
        end
    end
    data = da_create_sample(DA,selectors,good_subjs);
    %---------------
    %plot:
    wd = pwd;
    cd (DA.folder_name)
    %da_2ndl_plot(data,'USER_DEFINED',sample_names,'USER_DEFINED',flag_paired,1,DA.reg_names,DA.proj_name)
    % TO FIX: reg_names is missing, it also affects da_2nd_plot first lines
    da_2ndl_plot(data,'USER_DEFINED',sample_names,'USER_DEFINED',flag_paired,1,'',DA.proj_name)
    if ER.DA_plot_movements
        da_2ndl_plot_only_movements(data,'USER_DEFINED',sample_names,'USER_DEFINED_ONLY_MOVEMENTS',flag_paired,1,'',DA.proj_name)
    end
    cd (wd);
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% gFC HISTOGRAM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if exist('data_gFC_HIST','var') && ~isempty(data_gFC_HIST) && ER.gFC_HIST_plot  
    s_p_row = ceil(length(ER.indx_roi)/2);
    for s = 1:sample_number
        temp = [];
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(2,temp,data_gFC_HIST(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),ER.sample(s).selector(ss,3),:,:)); 
        end
        if cut_subjs(s)
            temp(:,~ER.sample(s).good_subjs,:,:,:,:) = [];
        end
        ER.sample(s).data_gFC_HIST = temp;
    end
    figure;
    set(gcf,'name','gFC HISTOGRAM')
    count = 0;
     x = -1:2/(200-1):1;
    for i=ER.indx_roi %ciclo sulle roi
        count = count +1;
        subplot(s_p_row,2,count)
        temp = cell(1,sample_number);
        for s = 1:sample_number
            for j=1:size(ER.sample(s).data_gFC_HIST,2) %subjs
                temp{s}(j,:) = squeeze(ER.sample(s).data_gFC_HIST(:,j,:,:,i,:));
            end
        end
        hold on
        stats = cell(1,sample_number);
        colors = hsv(sample_number);
        colors_m = abs(colors - (60/255));
        hs = cell(sample_number,1);
        hp = hs;
        leg_indx = [];
        for s = 1:sample_number
            stats{s}(1,:) = mean(temp{s},1);
            stats{s}(2,:) = std(temp{s},1)./sqrt(ER.sample(s).tot_subj-1);
            hs{s} = shadedplot(x', (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
            hp{s} = plot(x',stats{s}(1,:),'linewidth',2);
            set(hp{s},'Color',colors_m(s,:))
            leg_indx = [leg_indx, hp{s}];
        end
        legend(leg_indx,sample_names)
%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
        title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
        xlabel('Pearson Correlation');
        ylabel('Count');
        xlim([-1 1]);
        box on
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PSD estimates (periodogram)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if exist('data_ps','var') && ~isempty(data_ps) && ER.psd_plot_options.plot
    sampling = 1/tr{1,1};
    freq = sampling/2*linspace(0,1,nfft/2); 
    s_p_row = ceil(length(ER.indx_roi)/2);
    for s = 1:sample_number
        temp = [];
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(2,temp,data_ps(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),ER.sample(s).selector(ss,3),:,:)); 
        end
        if cut_subjs(s)
            temp(:,~ER.sample(s).good_subjs,:,:,:,:) = [];
        end
        ER.sample(s).data_ps = temp;
    end
    figure;
    set(gcf,'name','PSD (periodogram)')
    count = 0;
    for i=ER.indx_roi %ciclo sulle roi
        count = count +1;
        subplot(s_p_row,2,count)   
        temp = cell(1,sample_number);
        for s = 1:sample_number
            for j=1:size(ER.sample(s).data_ps,2) %subjs
                temp{s}(j,:) = squeeze(ER.sample(s).data_ps(:,j,:,:,i,:));
            end
        end
        hold on
        stats = cell(1,sample_number);
        colors = hsv(sample_number);
        colors_m = abs(colors - (60/255));
        hs = cell(sample_number,1);
        hp = hs;
        leg_indx = [];
        for s = 1:sample_number
            stats{s}(1,:) = mean(temp{s},1);
            stats{s}(2,:) = std(temp{s},1)./sqrt(ER.sample(s).tot_subj-1);
            hs{s} = shadedplot(freq, (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
            hp{s} = plot(freq,stats{s}(1,:),'linewidth',2);
            set(hp{s},'Color',colors_m(s,:))
            leg_indx = [leg_indx, hp{s}];
        end
        legend(leg_indx,sample_names)
%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
        title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
        xlabel('Frequency (Hz)');
        ylabel('PSD');
        
        if filter_band(2) == 99999;
            x_lim_top = 1/(2*tr{1,1});
        else
            x_lim_top = filter_band(2) + (10*filter_band(2)/100);
        end
        xlim([0 x_lim_top]);
        box on
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PSD estimates (Welch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('data_psw','var') && ~isempty(data_psw) && ER.psdw_plot_options.plot
    freq = freqw;
    s_p_row = ceil(length(ER.indx_roi)/2);
    for s = 1:sample_number
        temp = [];
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(2,temp,data_psw(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),ER.sample(s).selector(ss,3),:,:)); 
        end
        if cut_subjs(s)
            temp(:,~ER.sample(s).good_subjs,:,:,:,:) = [];
        end
        ER.sample(s).data_psw = temp;
    end
    figure;
    set(gcf,'name','PSD with Welch''s method')
    count = 0;
    for i=ER.indx_roi %ciclo sulle roi
        count = count +1;
        subplot(s_p_row,2,count)
        temp = cell(1,sample_number);
        for s = 1:sample_number
            for j=1:size(ER.sample(s).data_psw,2) %subjs
                temp{s}(j,:) = squeeze(ER.sample(s).data_psw(:,j,:,:,i,:));
            end
        end
        hold on
        stats = cell(1,sample_number);
        colors = hsv(sample_number);
        colors_m = abs(colors - (60/255));
        hs = cell(sample_number,1);
        hp = hs;
        leg_indx = [];
        for s = 1:sample_number
            stats{s}(1,:) = mean(temp{s},1);
            stats{s}(2,:) = std(temp{s},1)./sqrt(ER.sample(s).tot_subj-1);
            hs{s} = shadedplot(freq, (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
            hp{s} = plot(freq,stats{s}(1,:),'linewidth',2);
            set(hp{s},'Color',colors_m(s,:))
            leg_indx = [leg_indx, hp{s}];
        end
        legend(leg_indx,sample_names)
%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
        title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
        xlabel('Frequency (Hz)');
        ylabel('PSD');
        xlim([0 (filter_band(2) + (10*filter_band(2)/100))]);
        box on
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PSD estimates (Multitaper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('data_psm','var') && ~isempty(data_psm) && ER.psdm_plot_options.plot
    freq = freqm;
    s_p_row = ceil(length(ER.indx_roi)/2);
    for s = 1:sample_number
        temp = [];
        for ss = 1: size(ER.sample(s).selector,1)
            temp = cat(2,temp,data_psm(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),ER.sample(s).selector(ss,3),:,:)); 
        end
        if cut_subjs(s)
            temp(:,~ER.sample(s).good_subjs,:,:,:,:) = [];
        end
        ER.sample(s).data_psm = temp;
    end
    figure;
    set(gcf,'name','PSD with multitaper')
    count = 0;
    for i=ER.indx_roi %ciclo sulle roi
        count = count +1;
        subplot(s_p_row,2,count)
        temp = cell(1,sample_number);
        for s = 1:sample_number
            for j=1:size(ER.sample(s).data_psm,2) %subjs
                temp{s}(j,:) = squeeze(ER.sample(s).data_psm(:,j,:,:,i,:));
            end
        end
        hold on
        stats = cell(1,sample_number);
        colors = hsv(sample_number);
        colors_m = abs(colors - (60/255));
        hs = cell(sample_number,1);
        hp = hs;
        leg_indx = [];
        for s = 1:sample_number
            stats{s}(1,:) = mean(temp{s},1);
            stats{s}(2,:) = std(temp{s},1)./sqrt(ER.sample(s).tot_subj-1);
            hs{s} = shadedplot(freq, (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
            hp{s} = plot(freq,stats{s}(1,:),'linewidth',2);
            set(hp{s},'Color',colors_m(s,:))
            leg_indx = [leg_indx, hp{s}];
        end
        legend(leg_indx,sample_names)
%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
        title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
        xlabel('Frequency (Hz)');
        ylabel('PSD');
        xlim([0 (filter_band(2) + (10*filter_band(2)/100))]);
        box on
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COH Magnitude plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('data_coh_ma','var') && ~isempty(data_coh_ma) && ER.coh_ma_plot
    if isempty(freq_coh) %for Backward compatibility
       freq_coh = linspace(0,1/(2*tr),size(data_coh_ma,7));
    end
    for v = 1:size(data_coh_ma,4)
        coh_ma = data_coh_ma(:,:,:,v,:,:,:);
        plot_coh_ma(freq_coh,coh_ma,name_coh_ma{v},ER);         %se servisse far tornare dalla funzione ER modificato
    end
end

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ROI to ROI connectivity plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(ER,'RtoR_file') || isempty(ER.RtoR_file)
    ER.RtoR.plot = 0;
end

if ER.RtoR.plot
    %------------------------------------------------------
    %save data output name
    if ER.save_data
        if ~isfield(ER.RtoR,'savedata_suffixname'); ER.RtoR.savedata_suffixname  = ''; end;% datestr(now,'dd-mmm-yyyy_HH:MM:SS'); end
        %     if isempty(ER.bar_plot_options.savedata_suffixname)
        %         ER.bar_plot_options.savedata_suffixname = datestr(now,'dd-mmm-yyyy_HH:MM:SS');
        %     end
        [~,outp,~] = fileparts(ER.roi_file);
        output_name_RtoR = ['_ER_plotted_',outp(9:end),'_RtoR_',ER.RtoR.savedata_suffixname];
    end
    %------
    load (ER.RtoR_file)
    atlas_name =RtoR_info.RtoR.name;
    RtoR.n_rois = size(RtoR_info.RtoR.roi_names,1);
    RtoR.roi_names = cell(RtoR.n_rois,1);
    for r = 1:RtoR.n_rois
        str = c_b(RtoR_info.RtoR.roi_names(r,:));
        RtoR.roi_names{r} = str_fixer_measure(str(1:end));
    end
    for s = 1:sample_number
        temp = [];
        % -- merging sessions
       if ~ER.sample(s).session_average
            for ss = 1: size(ER.sample(s).selector,1)
                a = RtoR_corr{ER.sample(s).selector(ss,1),ER.sample(s).selector(ss,2)};
                temp = cat(1,temp,a); 
            end
        else
            for ss = 1: size(ER.sample(s).selector,1)
                a = RtoR_corr{ER.sample(s).selector(ss,1),ER.sample(s).selector(ss,2)};
                temp = cat(5,temp,a);
            end
            temp = mean(temp,5);
       end
        % -- -----------------
        if isfield (ER.sample(s),'RtoR_cova') && ~isempty(ER.sample(s).RtoR_cova)
           ER.sample(s).RtoR.cova_do = 1; 
           s_cova = size(ER.sample(s).RtoR_cova);
           if s_cova(2) > s_cova(1)
               ER.sample(s).RtoR_cova = ER.sample(s).RtoR_cova';
           end
           ER.sample(s).RtoR.cova = ER.sample(s).RtoR_cova;
        else
            ER.sample(s).RtoR.cova_do = 0;
        end
        if cut_subjs(s) == 1;
            temp(~ER.sample(s).good_subjs,:,:,:) = [];
            if ER.sample(s).RtoR.cova_do
                ER.sample(s).RtoR.cova(~ER.sample(s).good_subjs,:) = [];
                if size(temp,1) ~= size(ER.sample(s).RtoR.cova,1)
                    error('No matching covariates (check RtoR_cova).');
                end
            end
        end
        ER.sample(s).RtoR.data = temp;  % 4d i = subj; j = r,z,p; lxl roi
    end

    if ER.RtoR.ostt.do
        for s =1:sample_number
            z_corr = squeeze(ER.sample(s).RtoR.data(:,2,:,:));
            ostt(s).r_corr_mean = tanh(squeeze(mean(z_corr,1)));
            X = ones(size(z_corr,1),1);
            if ER.sample(s).RtoR.cova_do
                % z-score covariate. Because mean centering is mandatory in
                % this type of model
                zcova = zscore(ER.sample(s).RtoR.cova);
                %zcova = bsxfun(@minus,ER.sample(s).RtoR.cova,mean(ER.sample(s).RtoR.cova));    % mean-centering and z-score are equivalent
                X = [X zcova];
                cova_str = ' with cova';
            else
                cova_str = '';
            end
            ostt(s).row_name = [sample_names{s},cova_str];
            [ostt(s).T,ostt(s).P] = RtoR_glm(z_corr,X,1);
        end
        RtoR_plot(ostt,['ATLAS: ',atlas_name,' One-Sample t-test'],1);
        if ER.save_data
            try 
                save(output_name_RtoR,'ostt','ER','-append');
            catch
                save(output_name_RtoR,'ostt','ER');
            end
        end
        
    end
    
    if ER.RtoR.tstt.do
        selector = ER.RtoR.tstt.samples;
        z_corr1 = squeeze(ER.sample(selector(1)).RtoR.data(:,2,:,:));
        z_corr2 = squeeze(ER.sample(selector(2)).RtoR.data(:,2,:,:));
        z_corr = cat(1,z_corr1,z_corr2);
        x0_1 = [ones(size(z_corr1,1),1);zeros(size(z_corr2,1),1)];
        x0_2 = ~x0_1;
        X = [x0_1 x0_2];
        C = [1 -1];
        if flag_paired 
            % check if it is possible
            size1 = size(z_corr1,1);
            size2 = size(z_corr2,1);
            if size1~=size2
                warning backtrace off
                warning('Sorry, something went wrong! Paired t-test cannot be performed, switching to unpaired modality.');
                warning backtrace on
                flag_paired = 0;
            else
                % OLD WAY (as FLS does). PROBLEM HERE IS COVARIATE. IF YOU
                % ADD A COVARIATE X BECOMES SINGULAR.  NEW:
%                 non era invertible perchè subject-specifc. Le covariate
%                 devono essere scan-specific.        
                X_pairs = [eye(size1);eye(size1)];
                x0_1 = [ones(size(z_corr1,1),1);-1*ones(size(z_corr1,1),1)];  %   con due regressori la matrice non era invertibile. 
                X = [x0_1, X_pairs];
                C = 1;
                % WORKAROUND FRO COVARIATE ISSUE. MAKE A ONE SAMPLE T-TEST
                % ON THE DIFFERENCE. IN THIS WAY IN CASE OF COVARIATE
                % MODALITY NO PROBLEM AT ALL
%                 z_corr = z_corr1 -z_corr2;
%                 X = ones(size(z_corr,1),1);
%                 C = 1;
                % NOTE: the two modalities are matematically equivalent, no
                % differnce neither in t-values nor in p-values
            end
        end 
        if (ER.sample(selector(1)).RtoR.cova_do + ER.sample(selector(2)).RtoR.cova_do) == 2
            cova = cat(1,ER.sample(selector(1)).RtoR.cova,ER.sample(selector(2)).RtoR.cova);
            if flag_paired 
                %cova = ER.sample(selector(1)).RtoR.cova;   % OLD, was done
                %to allow fro subject-specfic covariate (like age). No
                %sense.
                 %cova = ER.sample(selector(1)).RtoR.cova -ER.sample(selector(2)).RtoR.cova;
                  
                  %cova = cova -mean(cova);    %it seems mean centering is
                  %never appropriate in this case

            end
            % no need for converting in z-score for this type of model
            X = [X cova];
            cova_str = ' with cova';
        else
            cova_str = '';
        end    
        [tstt(1).T,tstt(1).P] = RtoR_glm(z_corr,X,C);

        RtoR_plot(tstt,['ATLAS: ',atlas_name,' Two-Sample t-test',cova_str,': ',sample_names{selector(1)},'>',sample_names{selector(2)}],1);
        tstt.name = ['t-test',cova_str,': ',sample_names{selector(1)},'>',sample_names{selector(2)}];
        if ER.save_data
            try 
                save(output_name_RtoR,'tstt','ER','-append');
            catch
                save(output_name_RtoR,'tstt','ER');
            end
        end
    end
    
    if ER.RtoR.corr.do
        for s =1:sample_number
            if ER.sample(s).corr_number == 0
                continue
            end
            for c = 1:ER.sample(s).corr_number
                data_x = ER.sample(s).corr(c).vect;
                X = ones(size(data_x,1),1);
                X = [X, data_x];
                z_corr = squeeze(ER.sample(s).RtoR.data(:,2,:,:));
                z_corr(ER.sample(s).corr(c).nan_index,:,:) = []; %removing NaNs (bad subjs cut already done)
                if ER.sample(s).corr(c).do_partial
                    cova = ER.sample(s).corr(c).cova;
                    % no need for converting in z-score for this type of model
                    X = [X, cova];
                    [corr_res(1).T,corr_res(1).P] = RtoR_glm(z_corr,X,[0 1]);
                    cova_str = '(with COVA) ';
                else
                    [corr_res(1).T,corr_res(1).P] = RtoR_glm(z_corr,X,[0 1]);
                    cova_str = '';
                end
                figure_title = ['ATLAS: ',atlas_name,' Correlation ',cova_str,'of ',sample_names{s},' with: ',ER.sample(s).corr(c).name];
                RtoR_plot(corr_res,figure_title,1);
                if ER.save_data
                    try
                        save(output_name_RtoR,'corr_res','ER','-append');
                    catch
                        save(output_name_RtoR,'corr_res','ER');
                    end
                end
                corr_res = [];
            end
        end
    end
   
end

return
end

function [T_2d,P_2d] = RtoR_glm(Y,X,C)
siz = size(Y);
n = siz(1);
r = siz(2);
Y = reshape(Y,[n, r*r]);
rang = rank(X);
warning off; C = [C, zeros(1,(size(X,2)-size(C,2)))];warning on;  %add missing zeros
C_number = size(C,1);

B = (X\Y);
RES = Y-X*B;
sigma2 = sum(RES.^2,1)/(n-rang);

if C_number == 1
    T = C*B./sqrt(sigma2*(C*inv(X'*X)*C'));
else
    T = C*B./sqrt( repmat(sigma2,[C_number,1]) .* diag((C*inv(X'*X)*C')) );
end

P =  2 * tcdf(-abs(T), n-rang);

%reshape
T_2d = nan(size(C,1),r,r);
P_2d = nan(size(C,1),r,r);
for j=1:size(C,1)
    T_2d(j,:,:) = reshape(T(j,:),r,r);
    P_2d(j,:,:) = reshape(P(j,:),r,r);
end
return
end

function RtoR_plot(data,fig_title,contrast)
global RtoR
n_row = length(data);

if isfield(data,'r_corr_mean');
    n_col = 3;
else
    n_col = 2;
end

figure;
set(gcf,'name',fig_title)

%let's find the max and min of tvalues
maxx = -100;
minn = 100;
for s = 1:n_row
    tmp = squeeze(data(s).T(contrast,:,:));
    maxxx = max(tmp(tmp ~= inf));
    minnn = min(tmp(tmp ~= inf));
    if maxxx > maxx
        maxx = maxxx;
    end
    if minnn < minn
        minn = minnn;
    end
end

% if abs(maxx-minn) < 0.001
%     if maxx > 0
%         minn = 0;
%     else
%         minn = maxx;
%         maxx = 0;
%     end
% end

% for simmetric colorbar
if abs(maxx) > abs(minn)
    minn = -abs(maxx);
    maxx = abs(maxx);
else
    minn = -abs(minn);
    maxx = abs(minn);
end
    

for s = 1:n_row  %run on samples
    if n_col == 3
        s2 = subplot(n_row,n_col, (n_col*(s-1)+1));
        imagesc(data(s).r_corr_mean,[-1,1]);
        pbaspect([1 1 1])
        %
        title('Pearson Correlation (r)')
        set(gca,'xtick',[])
        
        if s == n_row
            pos = get(s2,'position');
            hb = colorbar('location','southoutside');
            set(s2,'position',pos);

        end
        add_indx = 1;
        if mod(((n_col*(s-1))+ (add_indx))-1,n_col) == 0
            set(gca,'ytick',[1:RtoR.n_rois])
            set(gca,'yTickLabel',RtoR.roi_names);
            set(gca,'FontSize',9)
            if isfield(data(s),'row_name');
                ht = text(-0.3,1.15,data(s).row_name,'Units','normalized','HorizontalAlignment','Center');
                set(ht,'fontsize',11,'FontWeight','bold');
            end
        end
    else
        add_indx = 0;
    end
    s2 = subplot(n_row,n_col, (n_col*(s-1))+ (1+add_indx));
    imagesc(squeeze(data(s).T(contrast,:,:)),[minn,maxx]);
    colormap jet;
    pbaspect([1 1 1])
    title('T-values (t)')
    set(gca,'xtick',[])
    if mod(((n_col*(s-1))+ (1+add_indx))-1,n_col) == 0
        set(gca,'ytick',[1:RtoR.n_rois])
        set(gca,'yTickLabel',RtoR.roi_names);
         set(gca,'FontSize',9)
    else
        set(gca,'ytick',[])
    end
    if s == n_row
        pos = get(s2,'position');
        hb = colorbar('location','southoutside');
        set(s2,'position',pos);
    end
    
    s2 = subplot(n_row,n_col, (n_col*(s-1))+ (2+add_indx));
    imagesc(squeeze(data(s).P(contrast,:,:)),[0,0.05]);
    pbaspect([1 1 1])
    title('P-values (p)')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    if s == n_row
        pos = get(s2,'position');
        hb = colorbar('location','southoutside');
        set(s2,'position',pos);
    end
end


return
end

function test_stat = isto(dati,sample,sample_number,sample_names)
  
global flag_type_plot flag_paired flag_plot_poits do_stat plot_stim

    markers = {'o','s','d','+','*','>','<','^'};        %add more marker here if u have merged more groups
    
    media = zeros(sample_number,1);
    sem = media;
    for s = 1: sample_number
        media(s) = er_nanmean(dati{s});
        sem(s) = er_nanstd(dati{s})/sqrt(sum(~isnan(dati{s}))-1);
    end
     
    hold on
    switch flag_type_plot
        case 2  %box plot
            temp_data = [];
            temp_g = [];
            for s = 1:sample_number
                temp_data = [temp_data;dati{s}];
                temp_g = [temp_g; s*ones(length(dati{s}),1)];
            end
            boxplot(temp_data,temp_g)
            if flag_plot_poits
                if flag_paired; line_x = cell(2,1); line_y = line_x; count_s = 0;end
                for s = 1:sample_number
                    if sample(s).same_group_subject
                        colors = hsv(sample(s).merged_subj(1));
                        colors = repmat(colors,[sample(s).merging_number,1]);
                    else
                        colors = hsv(sample(s).tot_subj);
                    end
                    a = 1:sample(s).merging_number; % for separing group
                    a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
                    if flag_paired; line_x_ = []; line_y_ = [];end
                    cl_indx = 0;
                    for l = 1:sample(s).merging_number
                            tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
                            x = repmat(s,1,length(tmp)); %the x axis location
                            x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
                            for j=1:sample(s).merged_subj(l)
                                cl_indx = cl_indx +1;
                                plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',5);
                            end
                            if flag_paired 
                                line_x_ = [line_x_,x];
                                line_y_ = [line_y_,tmp'];
                            end
                    end
                    if flag_paired;
                        line_x{s} = line_x_;
                        line_y{s} = line_y_;
                        count_s = count_s +1;
                    end
                end
                if flag_paired && count_s == 2;
                    x = [line_x{1};line_x{2}];
                    y = [line_y{1};line_y{2}];
                    %plot(x,y,'--k','MarkerFaceColor','none')
                    plot(x,y,':','Color',[0.3,0.3,0.3],'MarkerFaceColor','none')
                end
            end
        case 1  % bar plot
            %bar_col = bone(sample_number);
            for z =1: sample_number
                %bar([z],media(z),'FaceColor',bar_col(z,:)); %'EdgeColor',[0.3,0.3,0.3]); 'EdgeColor','none');
                bar([z],media(z),'FaceColor',[0.3,0.3,0.3],'EdgeColor','none'); %'EdgeColor',[0.3,0.3,0.3]); 'EdgeColor','none');
            end
            if flag_plot_poits
                for s = 1:sample_number
                        if sample(s).same_group_subject
                            colors = hsv(sample(s).merged_subj(1));
                            colors = repmat(colors,[sample(s).merging_number,1]);
                        else
                            colors = hsv(sample(s).tot_subj);
                        end
                        a = 1:sample(s).merging_number; % for separing group
                        a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
                        cl_indx = 0;
                        for l = 1:sample(s).merging_number
                                tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
                                x = repmat(s,1,length(tmp)); %the x axis location
                                x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
                                for j=1:sample(s).merged_subj(l)
                                    cl_indx = cl_indx +1;
                                    plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',4)
                                end
                        end
    %                 else
    %                     subjs_number = sample(s).tot_subj;
    %                     colors = hsv(subjs_number);
    %                     start_index = 1;
    %                     final_index = sample(s).merged_subj(1);
    %                     a = 1:sample(s).merging_number; % for separing group
    %                     a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
    %                     for l = 1:sample(s).merging_number
    %                             tmp = dati{s}(start_index:final_index); %temporarily store data in variable "tmp"
    %                             x = repmat(s,1,length(tmp)); %the x axis location
    %                             x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
    %                             c = 0;
    %                             for j=start_index:final_index
    %                                 c = c + 1;
    %                                 plot(x(c),tmp(c),markers{l},'Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',4)
    %                             end
    %                             start_index = final_index + 1;
    %                             if start_index > subjs_number
    %                                 break
    %                             end
    %                             final_index = start_index + sample(s).merged_subj(l+1) -1;
    %                     end
    %                 end
                end
            end
            h=errorbar(media, sem,'k');
            set(h,'LineStyle','none');
            set(h,'LineWidth',4);
        case 3 % line plot
            if flag_plot_poits
                for s = 1:sample_number
                        if sample(s).same_group_subject
                            colors = hsv(sample(s).merged_subj(1));
                            colors = repmat(colors,[sample(s).merging_number,1]);
                        else
                            colors = hsv(sample(s).tot_subj);
                        end
                        a = 1:sample(s).merging_number; % for separing group
                        a = flipdim((zscore(a)./(sample(s).merging_number+3)),2);
                        cl_indx = 0;
                        for l = 1:sample(s).merging_number
                                tmp = dati{s}([sample(s).merged_subj_index{l}]); %temporarily store data in variable "tmp"
                                x = repmat(s,1,length(tmp)); %the x axis location
                                x = (x+(rand(size(x))-0.5)*0.2) -a(l); %add a little random "jitter" to aid visibility
                                for j=1:sample(s).merged_subj(l)
                                    cl_indx = cl_indx +1;
                                    plot(x(j),tmp(j),markers{l},'Color',colors(cl_indx,:),'MarkerFaceColor',colors(cl_indx,:),'MarkerSize',4)
                                end
                        end
                end
            end
            h=errorbar(media, sem,'k');
            set(h,'LineWidth',2);
    end
    
    %-------------testing for significance-----------------------
    if do_stat
    
        if sample_number == 1
            MM = 0;
            [~,p,~,stat] = ttest(dati{1},MM,0.05,'both');    
            combs = [1];
        else
            combs = combs_no_rep(sample_number,2);
            if flag_paired == 1
                for l = 1:size(combs,1)   
                    [~,p(l),~,stat(l)] = ttest(dati{combs(l,1)},dati{combs(l,2)},0.05,'both');  
                end
            else
                for l = 1:size(combs,1)   
                    [~,p(l),~,stat(l)] = ttest2(dati{combs(l,1)},dati{combs(l,2)},0.05,'both');  
                end

            end
        end

        for l = 1:size(combs,1)
            sigstar_cell{l} = combs(l,:);
        end
        hsig = sigstar(sigstar_cell,p);
        set(hsig(:,1), 'LineWidth',0.8);
        set(hsig(:,2), 'FontSize',11);
        %set(hsig(:,2), 'FontWeight','bold');
        test_stat.sigstar = sigstar_cell;
        test_stat.p = p;
    else
        test_stat = 'Not performed';
    end
    %------------------------------------------------------------
        
    %set(gca,'YTick',[0:0.1:0.4])
    %set(gca,'YTickLabel',['0';' ';'0.1';' ';'0.2';' ';'0.3';' ';'0.4'])
    
    %a = ylim;
    %text(-1.1,a(2)+0.01/a(2),['\fontsize{12px}',letter]);
    %text(xpos,1.070,['\fontsize{12px}',letter],'Units', 'Normalized', 'VerticalAlignment', 'Top')
    
   
    set(gca,'xtick',[1:sample_number]);
    set(gca,'XTickLabel',sample_names);
%     set(gca,'FontSize',9);     
    %set(gca,'FontWeight','bold'); 
%     set(gca, 'LineWidth'   , 0.7     );
    
          
    if plot_stim.do
        y_lim = ylim;
        wei_y = y_lim(2) - y_lim(1);
        for n = 1:plot_stim.n_stim
           for nn = 1:length(plot_stim.POS{n})
               pos = [plot_stim.POS{n}(nn) y_lim(1) plot_stim.WEI{n}(nn) (wei_y)];
               r = rectangle('Position',pos,'FaceColor', plot_stim.COL{n},'LineStyle','none');
               uistack(r,'bottom')
           end
       end
        xlim([0.5 (sample_number+0.5)]);
        ylim(y_lim);
    else
        xlim([0 (sample_number+1)]);
    end

    box on
    hold off
    
    
    
end

function correlation(y,x,x_name,xlimit,co)
global corr_type
hs = scatter(x,y,20,'filled','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);  
hold on
a=~isnan(x);b=~isnan(y);a=a.*b;
if sum(a)>3	    
    [r,p] = corr(x(a==1),y(a==1),'type',corr_type);
    if p <= 0.05
        [cc2,~]=fit(x(a==1),y(a==1),'poly1');
        c2= plot(cc2);
        set(c2,'Color',[180/255,0,0],'Linewidth',2);
        ht = text(0.05,0.90,sprintf(['r: ',num2str(r,2),'\np: ',num2str(p,2)]),'Units','normalized');
        set(ht, 'FontSize',9,'FontWeight','bold','Color',[180/255,0,0]);
        legend off
        ylabel([]);
    end
end
if co.do
    %compute RESy
    By = co.pre*y;
    RESy = y-co.cova(:,1:end-1)*By(1:end-1);
    hs2 = scatter(co.RESx,RESy,20,'filled','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0,0,180/255]);
    set(hs,'MarkerEdgeColor',[180/255,0,0]);
    a=~isnan(co.RESx);b=~isnan(RESy);a=a.*b;
    if sum(a)>3	    
        [r,p] = corr(co.RESx(a==1),RESy(a==1),'type',corr_type);
        if p <= 0.05
            [cc2,~]=fit(co.RESx(a==1),RESy(a==1),'poly1');
            c2= plot(cc2);
            set(c2,'Color',[0,0,180/255],'Linewidth',2);
            ht = text(0.6,0.90,sprintf(['r: ',num2str(r,2),'\np: ',num2str(p,2)]),'Units','normalized');
            set(ht, 'FontSize',9,'FontWeight','bold','Color',[0,0,180/255]);
            legend off
            ylabel([]);
        end
    end
end
xlabel(x_name);    
% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim(xlimit);
hold off
grid on
box on
return
end

function ou = rm_first_char(in)
in(1:2) = '';
ou = in;
return
end

function correlation_outliers(y,x,x_name,xlimit,cova)
%serve the case of the partial correlation plus outliers
global corr_type
hs = scatter(x,y,20,'filled','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0.4 0.4 0.4]);
hold on
a=~isnan(x);b=~isnan(y);a=a.*b;
if sum(a)>3	    
    [r,p] = corr(x(a==1),y(a==1),'type',corr_type);
    if p <= 0.05
        [cc2,~]=fit(x(a==1),y(a==1),'poly1');
        c2= plot(cc2);
        set(c2,'Color',[180/255,0,0],'Linewidth',2);
        ht = text(0.05,0.90,sprintf(['r: ',num2str(r,2),'\np: ',num2str(p,2)]),'Units','normalized');
        set(ht, 'FontSize',9,'FontWeight','bold','Color',[180/255,0,0]);
        legend off
        ylabel([]);
    end
end
x(~a) = [];
y(~a) = [];
cova(~a,:) = [];
%compute RESy
cova = zscore(cova);
cova = [cova ones(size(cova,1),1)]; % adding costant term
pre = inv(cova'*cova)*cova';
Bx = pre*x;
By = pre*y;
RESx = x-cova(:,1:end-1)*Bx(1:end-1);
RESy = y-cova(:,1:end-1)*By(1:end-1);
hs2 = scatter(RESx,RESy,20,'filled','o','MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor',[0,0,180/255]);
set(hs,'MarkerEdgeColor',[180/255,0,0]);
[r,p] = corr(RESx,RESy,'type',corr_type);
if p <= 0.05
    [cc2,~]=fit(RESx,RESy,'poly1');
    c2= plot(cc2);
    set(c2,'Color',[0,0,180/255],'Linewidth',2);
    ht = text(0.6,0.90,sprintf(['r: ',num2str(r,2),'\np: ',num2str(p,2)]),'Units','normalized');
    set(ht, 'FontSize',9,'FontWeight','bold','Color',[0,0,180/255]);
    legend off
    ylabel([]);
end


xlabel(x_name);
xlim(xlimit);
hold off
grid on
box on
return
end

function [outcome] = below_thr_correlation(y,x,co,thr)
% this serve only to find the indx of rois e measure significant
% correlation
global corr_type
p = 1;
p_cova = 1;
a=~isnan(x);b=~isnan(y);a=a.*b;
if sum(a)>3	    
    [~,p] = corr(x(a==1),y(a==1),'type',corr_type);
end
if co.do
    %compute RESy
    By = co.pre*y;
    RESy = y-co.cova(:,1:end-1)*By(1:end-1);
    a=~isnan(co.RESx);b=~isnan(RESy);a=a.*b;
    if sum(a)>3	    
        [~,p_cova] = corr(co.RESx(a==1),RESy(a==1),'type',corr_type);
    end
end
if p <= thr || p_cova <= thr
    outcome = 1;
else
    outcome = 0;
end
return
end

function [outcome] = below_thr_correlation_outliers(y,x,cova,thr)
% this serve only to find the indx of rois e measure significant
% correlation. Case plus outliers
global corr_type
p = 1;
p_cova = 1;
a=~isnan(x);b=~isnan(y);a=a.*b;
if sum(a)>3	    
    [~,p] = corr(x(a==1),y(a==1),'type',corr_type);
end
x(~a) = [];
y(~a) = [];
cova(~a,:) = [];
%compute RESy
cova = [cova ones(size(cova,1),1)]; % adding costant term
pre = inv(cova'*cova)*cova';
Bx = pre*x;
By = pre*y;
RESx = x-cova(:,1:end-1)*Bx(1:end-1);
RESy = y-cova(:,1:end-1)*By(1:end-1);
[~,p_cova] = corr(RESx,RESy,'type',corr_type);
if p <= thr || p_cova <= thr
    outcome = 1;
else
    outcome = 0;
end

return
end

function plot_coh_ma(freq,coh_ma,seed_name,ER)
global sample_number sample_names temp_rois subject_number filter_band cut_subjs
 
s_p_row = ceil(length(ER.indx_roi)/2);
for s = 1:sample_number
    temp = [];
    for ss = 1: size(ER.sample(s).selector,1)
        temp = cat(2,temp,coh_ma(ER.sample(s).selector(ss,1),1:subject_number(ER.sample(s).selector(ss,1)),ER.sample(s).selector(ss,2),1,ER.sample(s).selector(ss,3),:,:)); 
    end
    if cut_subjs(s)
        temp(:,~ER.sample(s).good_subjs,:,:,:,:,:) = [];
    end
    ER.sample(s).coh_ma = temp;
end
figure;
set(gcf,'name',['Coherency Magnitude SEED = ',seed_name])
count = 0;
for i=ER.indx_roi %ciclo sulle roi
    count = count +1;
    subplot(s_p_row,2,count)
    temp = cell(1,sample_number);
    for s = 1:sample_number
        for j=1:size(ER.sample(s).coh_ma,2) %subjs
            temp{s}(j,:) = squeeze(ER.sample(s).coh_ma(:,j,:,:,:,i,:));
        end
    end
    hold on
    stats = cell(1,sample_number);
    colors = hsv(sample_number);
    colors_m = abs(colors - (60/255));
    hs = cell(sample_number,1);
    hp = hs;
    leg_indx = [];
    for s = 1:sample_number
        stats{s}(1,:) = mean(temp{s},1);
        stats{s}(2,:) = std(temp{s},1)./sqrt(ER.sample(s).tot_subj-1);
        hs{s} = shadedplot(freq, (stats{s}(1,:)-stats{s}(2,:)), (stats{s}(1,:)+stats{s}(2,:)), colors(s,:));
        hp{s} = plot(freq,stats{s}(1,:),'linewidth',2);
        set(hp{s},'Color',colors_m(s,:))
        leg_indx = [leg_indx, hp{s}];
    end
    legend(leg_indx,sample_names)
%         %s1 = shadedplot(freq, (media_REST-sd_REST), (media_REST+sd_REST), [153/255 255/255 162/255]);
%         s2 = shadedplot(freq, (media_L-sd_L), (media_L+sd_L), [255/255 102/255 51/255]);
%         s3 = shadedplot(freq, (media_R-sd_R), (media_R+sd_R), [153/255 153/255 255/255]);
    title([temp_rois(i).name],'fontsize',10,'FontWeight','bold');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    xlim([filter_band(1) filter_band(2)])
    box on
    hold off
end

return
end

function [data] = da_create_sample(DA,selectors,good_subjs)

load (DA.path);
DA = DA_extracted;
data.fd_TRnorm = da_create_sample_aux(selectors,DA.fd_TRnorm,good_subjs);
data.rms = da_create_sample_aux(selectors,DA.rms,good_subjs);
data.rms_dvars = da_create_sample_aux(selectors,DA.rms_dvars,good_subjs);
data.rms_sd = da_create_sample_aux(selectors,DA.rms_sd,good_subjs);
for r = 1:size(DA.vTv,4)
    data.vTv{r} = da_create_sample_aux(selectors,DA.vTv,good_subjs,r);
    data.vTv_mean{r} = da_create_sample_aux(selectors,DA.vTv_mean,good_subjs,r);
    data.var_exp{r}= da_create_sample_aux(selectors,DA.var_exp,good_subjs,r);
    for s = 1:size(DA.specifcity,5)
        data.specificity{r,s}= da_create_sample_aux(selectors,DA.specifcity,good_subjs,r,s);
    end
end
return
end

function Y = da_create_sample_aux(selectors,y,good_subjs,reg_n,spe)
global sample_number subject_number

for s = 1:sample_number
    temp = [];
    for ss = 1: size(selectors{s},1)
        group_sel = selectors{s}(ss,1);
        sess_sel = selectors{s}(ss,2);
        for k = 1:subject_number(group_sel)
            if nargin == 4
                temp = [temp;y{group_sel,sess_sel,k,reg_n}];
            elseif nargin == 5
                temp = [temp;y{group_sel,sess_sel,k,reg_n,spe}];
            else
                temp = [temp;y{group_sel,sess_sel,k}];
            end
        end
    end
    if ~isempty(good_subjs{s})
        temp(~good_subjs{s},:) = [];
    end
    Y{s} = temp;
end

return
end

function [roi] = str_fixer(roi)

for i=1:length(roi)
    roi(i).name = strrep(roi(i).name(1:end), '_', ' ');
    if ~isempty(strfind(roi(i).name,'bmi'))
        roi(i).name = roi(i).name(5:end);
    end
        
end

return
end

function [str] = str_fixer_measure(str)
str = strrep(str, '_', ' ');
str = strrep(str, '.', '');
a = strfind(str,'bmi');
if ~isempty(a)
    str(a:a+2) = [];
end
        
return
end






