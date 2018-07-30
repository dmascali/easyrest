function easyrest_analysis(ER_)
% cov is now implemented only in OSTT TSTT and TSTTp
% subj selector is implemented in OSTT, TSTT and TSTTp
% extra contrasts implemented only in OSTT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDITABLE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ER.project = '/data/danielem/visanet/ER_analysis/ER_analysis.mat';

ER.an.cov(1).name = 'age_all';
ER.an.cov(1).c = [69;72;59;83;66;71;70;65;58;62;75;68;74;71;70;70;70;74;73;74;76;55;78;80;68;71;77;70;73;72;62;63;76;81;70;79;88;76;72;76;85;68;76;85;82;57;72;50;50;82;61;73;66;76;54;56;59;59;68;58;67;61;79;73;65;74;61;70;76;69;59;74;67];
ER.an.cov(2).name = 'sex_all';
ER.an.cov(2).c = [-1;-1;1;-1;-1;-1;-1;-1;1;-1;-1;-1;-1;-1;-1;-1;-1;1;-1;1;-1;-1;1;1;-1;-1;-1;-1;1;-1;-1;1;-1;-1;-1;-1;-1;1;1;-1;-1;1;-1;1;-1;1;1;1;1;1;1;1;1;1;-1;-1;-1;-1;1;-1;1;-1;1;1;-1;1;-1;-1;1;1;-1;1;1];

ER.an.model.one_sample_ttest.do = 1;
ER.an.model.one_sample_ttest.model(1).name = 'ADL';
ER.an.model.one_sample_ttest.model(1).selector = [1,1];                                 % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session         
ER.an.model.one_sample_ttest.model(1).covariate = [];                                   % Indexs of the covariates you want to include in the model
ER.an.model.one_sample_ttest.model(1).extra_contrast(1).tcon.name = 'cov1-effect';      % If you need extra contrasts (besides the usual [1 0 ...])
ER.an.model.one_sample_ttest.model(1).extra_contrast(1).tcon.convec = [0 1];
ER.an.model.one_sample_ttest.model(1).no_zero_contrast = 0;                             % if 1, the usual contrast [1 0 ...] (i.e., the group/session contrast) will NOT be calculated. This situation is useful if you are interested only in the covairate effect and not in the beta zero regressor 
ER.an.model.one_sample_ttest.model(1).subj_selector = [];                               % A vector of ones and zeros for selecting good subjects only

ER.an.model.one_sample_ttest.model(2).name = 'ADNL';
ER.an.model.one_sample_ttest.model(2).selector = [2,1];
ER.an.model.one_sample_ttest.model(2).covariate = [];
ER.an.model.one_sample_ttest.model(3).name = 'HS';
ER.an.model.one_sample_ttest.model(3).selector = [3,1];
ER.an.model.one_sample_ttest.model(3).covariate = [];
ER.an.model.one_sample_ttest.model(4).name = 'AD';
ER.an.model.one_sample_ttest.model(4).selector = [1,1;2,1];
ER.an.model.one_sample_ttest.model(4).covariate = [];


ER.an.model.two_sample_ttest.do = 1;
ER.an.model.two_sample_ttest.model(1).name = 'HS_vs_ADL';
ER.an.model.two_sample_ttest.model(1).sample1 = [3,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(1).sample2 = [1,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session   
ER.an.model.two_sample_ttest.model(1).session_average = 1;  % 1/0; if 1 different sessions will be averaged. Sessions must belong to the same group. If 0 simple merging.
ER.an.model.two_sample_ttest.model(1).paired = 0;
ER.an.model.two_sample_ttest.model(1).covariate = [1,2];    % Indexs of the covariates you want to include in the model (NB In case of paired use the same covariate ordering of NOT paired. ER takes care of reordering the vector)
ER.an.model.two_sample_ttest.model(1).subj_selector = [];   % A vector of ones and zeros for selecting good subjects only (NB 1D vector that runs on [sample1,sample2]. If paired it must run only on [sample1])

ER.an.model.two_sample_ttest.model(2).name = 'HS_vs_ADNL';
ER.an.model.two_sample_ttest.model(2).sample1 = [3,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(2).sample2 = [2,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(2).paired = 0;
ER.an.model.two_sample_ttest.model(2).covariate = [];

ER.an.model.two_sample_ttest.model(3).name = 'ADNL_vs_ADL';
ER.an.model.two_sample_ttest.model(3).sample1 = [2,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(3).sample2 = [1,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(3).paired = 0;
ER.an.model.two_sample_ttest.model(3).covariate = [];

ER.an.model.two_sample_ttest.model(4).name = 'HS_vs_AD';
ER.an.model.two_sample_ttest.model(4).sample1 = [3,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(4).sample2 = [1,1;2,1];      % 2D vector: first colum, group selecotr; second colum, session selector. Add extra row to combine group/session        
ER.an.model.two_sample_ttest.model(4).paired = 0;
ER.an.model.two_sample_ttest.model(4).covariate = [];

ER.an.model.one_way_anova.do = 1;
ER.an.model.one_way_anova.model(1).name = '';
ER.an.model.one_way_anova.model(1).selector = [1,1;2,1;3,1];  % nx2 matrix. Each row is a group. Columns are group and condition indexes respectively.
ER.an.model.one_way_anova.model(1).covariate = [];

ER.an.measure.all = 0;  % 1/0. 1: do analysis on all available measure. 0: use only measures specified in the selector below;
ER.an.measure.selector = {'StV'};   % Select measure to analyize (only if ER.an.measure.all = 0) 
                                      % {'ALFF','fALFF','mALFF','RSFA','mRSFA','fRSFA','gFC','SpEn','DELAY',StV','REHO19'};
ER.an.seed.all = 1;
ER.an.seed.selector = {'PCC'};

%contrast thresholds (only two different thresholds can be used). they are
%applied to all defined contrasts
ER.an.model.contrast(1).titlestr = 'All contrasts, p<0.001 k=10';
ER.an.model.contrast(1).threshdesc = 'none';    %multiple comparison correction
ER.an.model.contrast(1).thresh = 0.001; % p-values 
ER.an.model.contrast(1).extent = 10;    % minum volume
ER.an.model.contrast(2).titlestr = 'All contrasts, p<0.05 FWE, k=0';
ER.an.model.contrast(2).threshdesc = 'FWE';    %multiple comparison correction
ER.an.model.contrast(2).thresh = 0.05; % p-values 
ER.an.model.contrast(2).extent = 0;    % minum volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW, UNLESS YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global opt
global opt

flag_check_update = 1;      % 1/0; automatic check for update. DO NOT DISABLE (usefull only in the case the repository is down)         

version = version_er;      
%%---------------------------------
% Welcome screen
welcome_screen(version)
%%---------------------------------

%%---------------------------------
% Check for update
if flag_check_update
    if exist('check_updates')~=2
        error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''/home/user/MATLAB/ER''))');
    end
    exit = check_updates(version);
    if exit;return;end
end
%%---------------------------------

if nargin
    fprintf('\nReading options from input variable...');
    clear ER;
    ER = ER_;
end
load (ER.project)
opt.an = ER.an;
fprintf('done');


fprintf('\nInitializing analysis.');
warning_str = prepare_and_check_consistency;
opt.an.letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
fprintf('done\n');
warning off backtrace
if ~isempty(warning_str)
    for l = 1:size(warning_str,1)
        warning(warning_str(l,:));
    end
end
warning on backtrace

spm('defaults','fMRI')
spm_jobman('initcfg');

if ~isfield(opt.an.model,'contrast')
    opt.an.model.contrast(1).titlestr = 'All contrasts, p<0.001 k=10';
    opt.an.model.contrast(1).threshdesc = 'none';    %multiple comparison correction
    opt.an.model.contrast(1).thresh = 0.001; % p-values 
    opt.an.model.contrast(1).extent = 10;    % minum volume
    opt.an.model.contrast(2).titlestr = 'All contrasts, p<0.05 FWE, k=0';
    opt.an.model.contrast(2).threshdesc = 'FWE';    %multiple comparison correction
    opt.an.model.contrast(2).thresh = 0.05; % p-values 
    opt.an.model.contrast(2).extent = 0;    % minum volume
end

for l = 1:length(opt.an.measure)    %on measure

    if ~opt.an.measure(l).do
        continue
    end
     
    %[~,~] = system(['rm -r ',opt.folders.secondlevel,'/',opt.an.measure(l).name]);
    warning off
    mkdir(opt.folders.secondlevel,opt.an.measure(l).output_dir);
    warning on
    dest_measure_dir = [opt.folders.secondlevel,'/',opt.an.measure(l).output_dir];
    
    %------------------- One sample t-test-------------------------
    if opt.an.model.one_sample_ttest.do
        
        model = opt.an.model.one_sample_ttest.model;
        
        for m = 1:length(model)   %on model type
            
            [~,~] = system(['rm -r ',dest_measure_dir,'/OSTT_',model(m).name]);
            mkdir(dest_measure_dir,['OSTT_',model(m).name]);
            opt.an.current_model_dir = [dest_measure_dir,'/OSTT_',model(m).name];
            [ostt] = prepare_data_ostt(m,l,opt.an.measure(l).subbrick);
            one_sample_ttest(ostt)
            
        end
    end
    %--------------------------------------------------------------
    
    %------------------- Two sample t-test-------------------------
    if opt.an.model.two_sample_ttest.do       
        model = opt.an.model.two_sample_ttest.model;
        for m = 1:length(model)   %on model type
            if model(m).paired
                [~,~] = system(['rm -r ',dest_measure_dir,'/TSTTp_',model(m).name]);
                mkdir(dest_measure_dir,['TSTTp_',model(m).name]);
                opt.an.current_model_dir = [dest_measure_dir,'/TSTTp_',model(m).name];
                tsttp = prepare_data_tstt_paired(m,l,opt.an.measure(l).subbrick);
                two_sample_ttest_paired(tsttp)
            else
                [~,~] = system(['rm -r ',dest_measure_dir,'/TSTT_',model(m).name]);
                mkdir(dest_measure_dir,['TSTT_',model(m).name]);
                opt.an.current_model_dir = [dest_measure_dir,'/TSTT_',model(m).name];
                tstt = prepare_data_tstt(m,l,opt.an.measure(l).subbrick);
                two_sample_ttest(tstt)
            end
        end
    end
    %--------------------------------------------------------------
    
    %--------------------one_way_anova-----------------------------
    if opt.an.model.one_way_anova.do
        
        model = opt.an.model.one_way_anova.model;
        
        for m = 1:length(model)   %on model type
            [~,~] = system(['rm -r ',dest_measure_dir,'/OWANOVA_',model(m).name]);
            mkdir(dest_measure_dir,['OWANOVA_',model(m).name]);
            opt.an.current_model_dir = [dest_measure_dir,'/OWANOVA_',model(m).name];
            
            anv = prepare_data_one_way_anova(m,l,opt.an.measure(l).subbrick);
            
            one_way_anova(anv)
            
        end
    end
    %--------------------------------------------------------------
    

end
fprintf('\nEasyRest Analysis done. Bye!\n');
return
end

function warning_str = prepare_and_check_consistency

global opt

warning_str = [];

if opt.an.model.one_sample_ttest.do
    for l = 1:length(opt.an.model.one_sample_ttest.model)
%         if any(size(opt.an.model.one_sample_ttest.model(l).selector) ~= [1,2])
%             error('Badly defined group/session selector.');
%         end        
%         if size(opt.an.model.two_sample_ttest.model(l).selector,2) ~= 2 
%             error(['TSTT #',num2str(l),'. Selector variable must have two columns (group,session).']);
%         end
        if ~isfield(opt.an.model.one_sample_ttest.model(l),'session_average') || isempty(opt.an.model.one_sample_ttest.model(l).session_average)
            opt.an.model.one_sample_ttest.model(l).session_average = 0;
        end
        if opt.an.model.one_sample_ttest.model(l).selector(1,1) > opt.group_number
            error(['OSTT #',num2str(l),'. Group selector is greater than the available group number.']);
        elseif opt.an.model.one_sample_ttest.model(l).selector(1,2) > opt.session_number
            error(['OSTT #',num2str(l),'. Condition/session selector is greater than the available condition/session number.']);
        end
    end
end

if opt.an.model.two_sample_ttest.do
    for l = 1:length(opt.an.model.two_sample_ttest.model)
        if size(opt.an.model.two_sample_ttest.model(l).sample1,2) ~= 2  || size(opt.an.model.two_sample_ttest.model(l).sample2,2) ~= 2 
            error(['TSTT #',num2str(l),'. Sample variables must have two columns (group,session).']);
        end
        %now is checking only the first row, change to check multiple rows
        if opt.an.model.two_sample_ttest.model(l).sample1(1,1) > opt.group_number
            error(['TSTT #',num2str(l),'. Group selector in sample1 is greater than the available group number.']);
        elseif opt.an.model.two_sample_ttest.model(l).sample1(1,2) > opt.session_number
            error(['TSTT #',num2str(l),'. Condition/session selector in sample1 is greater than the available condition/session number.']);
        elseif opt.an.model.two_sample_ttest.model(l).sample2(1,1) > opt.group_number
            error(['TSTT #',num2str(l),'. Group selector in sample2 is greater than the available group number.']);
        elseif opt.an.model.two_sample_ttest.model(l).sample2(1,2) > opt.session_number
            error(['TSTT #',num2str(l),'. Condition/session selector in sample2 is greater than the available condition/session number.']);
        end
        if ~isfield(opt.an.model.two_sample_ttest.model(l),'session_average')  || isempty(opt.an.model.two_sample_ttest.model(l).session_average)
            opt.an.model.two_sample_ttest.model(l).session_average = 0;
        end
        %there is no check on "pairedbility". Add it?        
        if opt.an.model.two_sample_ttest.model(l).session_average 
            if length(unique(opt.an.model.two_sample_ttest.model(l).sample1(:,1))) ~=1 || length(unique(opt.an.model.two_sample_ttest.model(l).sample2(:,1))) ~=1
                warning_str = strvcat(warning_str,['TSTT #',num2str(l),'. Session average is only possible among the same group of subjects. Switching to merge mode.']);
                opt.an.model.two_sample_ttest.model(l).session_average = 0;
            end
        end
        
        
    end
end

fprintf('.');
measure = measure_definition;
fprintf('.');

measure_count = 0;
if ~opt.an.measure.all
    for l = 1:length(measure)
        if ~any(strcmpi(opt.an.measure.selector,measure(l).name_selector))
            measure(l).do = 0;
        end
        measure_count = measure_count + measure(l).do;
    end
else
    for l = 1:length(measure)
        measure_count = measure_count + measure(l).do;
    end
end

if ~opt.an.seed.all
    for l = 1:length(measure)
        if isempty(measure(l).seed_name)
            continue
        end
        if ~any(strcmpi(opt.an.seed.selector,measure(l).seed_name))
           measure(l).do = 0;
        end
    end
end

   
opt.an.measure = measure;

if measure_count == 0
    error('No measure to analyze.')
end

return
end


function [ostt] = prepare_data_ostt(m,l,subbrick) %m model, l measure
global opt
[data,N] = build_data_cell(opt.an.model.one_sample_ttest.model(m),l,subbrick,'selector');

cov_selctor = opt.an.model.one_sample_ttest.model(m).covariate;
cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
if ~isempty(cov_selctor)
    for ll = 1:length(cov_selctor)
        if length(opt.an.cov(cov_selctor(ll)).c) ~= total_rows
            error(['Covariate ',num2str(cov_selctor(ll)),' in OSTT model ',num2str(m),' has a wrong number of elements']);
        end
        cov(ll).c = opt.an.cov(cov_selctor(ll)).c;
        cov(ll).cname = opt.an.cov(cov_selctor(ll)).name;
        cov(ll).iCFI = 1;
        cov(ll).iCC = 1;
    end
end

contrasts = [];
if isfield(opt.an.model.one_sample_ttest.model(m),'extra_contrast')
    contrasts = opt.an.model.one_sample_ttest.model(m).extra_contrast;
    indx = arrayfun( @(s) isempty(s.tcon),contrasts);       %mi assicuro che non ci siano campi empty
    contrasts(indx) = [];
    bad_indx = [];
    for ll = 1:length(contrasts)
        if ~isempty(contrasts(ll).tcon.convec)
            if isempty(contrasts(ll).tcon.name) 
                error(['Contrast ',num2str(ll),' in OSTT model ',num2str(m),' has no name definition.']);
            end
            contrasts(ll).tcon.sessrep = 'none';
        else
            bad_indx = [bad_indx,ll];
        end
    end
    contrasts(bad_indx) = [];
    if isempty(contrasts)
        contrasts =[];
    end
end
   
% Select subj
if isfield(opt.an.model.one_sample_ttest.model(m),'subj_selector') && ~isempty(opt.an.model.one_sample_ttest.model(m).subj_selector)
    selector = opt.an.model.one_sample_ttest.model(m).subj_selector;
    if N ~=  length(selector) 
        error('Subj_selector is bad defined.');
    end
    data(~selector) = [];
    
    for l = 1:length(cov)
        cov(l).c(~selector) = [];
    end
end

ostt.data = data;
ostt.cov = cov;
ostt.contrasts = contrasts;
ostt.no_zero_contrast = opt.an.model.one_sample_ttest.model(m).no_zero_contrast;
return
end

function [tstt] = prepare_data_tstt(m,l,subbrick) %m model, l measure

global opt

[data1,N1] = build_data_cell(opt.an.model.two_sample_ttest.model(m),l,subbrick,'sample1');
[data2,N2] = build_data_cell(opt.an.model.two_sample_ttest.model(m),l,subbrick,'sample2');

N = N1 + N2;

cov_selctor = opt.an.model.two_sample_ttest.model(m).covariate;
cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
if ~isempty(cov_selctor)
    for ll = 1:length(cov_selctor)
        if length(opt.an.cov(cov_selctor(ll)).c) ~= N
            error(['Covariate ',num2str(cov_selctor(ll)),' in TSTT model ',num2str(m),' has a wrong number of elements']);
        end
        cov(ll).c = opt.an.cov(cov_selctor(ll)).c;
        cov(ll).cname = opt.an.cov(cov_selctor(ll)).name;
        cov(ll).iCFI = 1;
        cov(ll).iCC = 1;
    end
end

% Select subj
if isfield(opt.an.model.two_sample_ttest.model(m),'subj_selector') && ~isempty(opt.an.model.two_sample_ttest.model(m).subj_selector)
    selector = opt.an.model.two_sample_ttest.model(m).subj_selector;
    if N ~=  length(selector) 
        error('Subj_selector in tstt is bad defined.');
    end
    sel1 = selector(1:N1);
    sel2 = selector(N1+1:end);
    data1(~sel1) = [];
    data2(~sel2) = [];
    for l = 1:length(cov)
        cov(l).c(~selector) = [];
    end
end

tstt.data1 = data1;
tstt.data2 = data2;
tstt.cov = cov;

return
end

function tsttp = prepare_data_tstt_paired(m,l,subbrick) %m model, l measure

global opt

[data1,N1] = build_data_cell(opt.an.model.two_sample_ttest.model(m),l,subbrick,'sample1');
[data2,N2] = build_data_cell(opt.an.model.two_sample_ttest.model(m),l,subbrick,'sample2');


if N1 ~= N2
    error (['TSTTp #',num2str(m),'. You can''t perform a paired t-test on this data']);
end

for s = 1:N1
    pair(s).scans = {data1{s,1};data2{s,1}};
end

cov_selctor = opt.an.model.two_sample_ttest.model(m).covariate;
cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
if ~isempty(cov_selctor)
    for ll = 1:length(cov_selctor)
        if length(opt.an.cov(cov_selctor(ll)).c) ~= (N1 + N2)
            error(['Covariate ',num2str(cov_selctor(ll)),' in TSTTp model ',num2str(m),' has a wrong number of elements']);
        end
        %Important, I assume covariate ordering is [N1;N2]. So I reorder it
        % to be [S01n1, S01n2, S02n1, S02n2...]
        c =opt.an.cov(cov_selctor(ll)).c; if ~isrow(c); c = c';end
        tmp1 = c(1:N1); 
        tmp2 = c(N1+1:end);
        c = [tmp1;tmp2]; c= c(:)';
        cov(ll).c = c;
        cov(ll).cname = opt.an.cov(cov_selctor(ll)).name;
        cov(ll).iCFI = 1;
        cov(ll).iCC = 1;
    end
end

% Select subj
if isfield(opt.an.model.two_sample_ttest.model(m),'subj_selector') && ~isempty(opt.an.model.two_sample_ttest.model(m).subj_selector)
    selector = opt.an.model.two_sample_ttest.model(m).subj_selector;
    if N1 ~=  length(selector)
        error('Subj_selector in tstt is bad defined.');
    end
    pair(~selector) = [];
    %for covariate I need a selector vector interleaved
    if ~isrow(selector); selector = selector';end
    selector = [selector;selector]; selector = selector(:)';
    for l = 1:length(cov)
        cov(l).c(~selector) = [];
    end
end


tsttp.pair = pair;
tsttp.cov = cov;

return
end


function one_sample_ttest(ostt)
global opt
matlabbatch{1}.spm.stats.factorial_design.dir = {opt.an.current_model_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = ostt.data;
matlabbatch{1}.spm.stats.factorial_design.cov = ostt.cov;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
if ~ ostt.no_zero_contrast
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'X > 0';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'X < 0';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    if ~isempty(ostt.contrasts)
        for ll = 1:length(ostt.contrasts)
            matlabbatch{3}.spm.stats.con.consess{2+ll} = ostt.contrasts(ll);
        end
    end
else
    if ~isempty(ostt.contrasts)
        for ll = 1:length(ostt.contrasts)
            matlabbatch{3}.spm.stats.con.consess{ll} = ostt.contrasts(ll);
        end
    end
end
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = opt.an.model.contrast(1).titlestr;
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = opt.an.model.contrast(1).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(1).thresh = opt.an.model.contrast(1).thresh;
matlabbatch{4}.spm.stats.results.conspec(1).extent = opt.an.model.contrast(1).extent;
matlabbatch{4}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.conspec(2).titlestr = opt.an.model.contrast(2).titlestr;
matlabbatch{4}.spm.stats.results.conspec(2).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = opt.an.model.contrast(2).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(2).thresh = opt.an.model.contrast(2).thresh;
matlabbatch{4}.spm.stats.results.conspec(2).extent = opt.an.model.contrast(2).extent;
matlabbatch{4}.spm.stats.results.conspec(2).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.print = true;

save ([opt.an.current_model_dir,'/batch_one_sample_ttest.mat'],'matlabbatch');

spm_jobman('run_nogui',matlabbatch);


return
end

function two_sample_ttest(tstt)
global opt
matlabbatch{1}.spm.stats.factorial_design.dir = {opt.an.current_model_dir};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = tstt.data1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = tstt.data2;
matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = tstt.cov;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat =  {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat =  {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '1 > 2';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '1 < 2';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = opt.an.model.contrast(1).titlestr;
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = opt.an.model.contrast(1).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(1).thresh = opt.an.model.contrast(1).thresh;
matlabbatch{4}.spm.stats.results.conspec(1).extent = opt.an.model.contrast(1).extent;
matlabbatch{4}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.conspec(2).titlestr = opt.an.model.contrast(2).titlestr;
matlabbatch{4}.spm.stats.results.conspec(2).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = opt.an.model.contrast(2).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(2).thresh = opt.an.model.contrast(2).thresh;
matlabbatch{4}.spm.stats.results.conspec(2).extent = opt.an.model.contrast(2).extent;
matlabbatch{4}.spm.stats.results.conspec(2).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.print = true;

save ([opt.an.current_model_dir,'/batch_two_sample_ttest.mat'],'matlabbatch');

spm_jobman('run_nogui',matlabbatch);


return 
end


function two_sample_ttest_paired(tsttp)
global opt
matlabbatch{1}.spm.stats.factorial_design.dir = {opt.an.current_model_dir};
matlabbatch{1}.spm.stats.factorial_design.des.pt.pair = tsttp.pair;
matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = tsttp.cov;
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat =  {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat =  {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '1 > 2';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '1 < 2';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = opt.an.model.contrast(1).titlestr;
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = opt.an.model.contrast(1).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(1).thresh = opt.an.model.contrast(1).thresh;
matlabbatch{4}.spm.stats.results.conspec(1).extent = opt.an.model.contrast(1).extent;
matlabbatch{4}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.conspec(2).titlestr = opt.an.model.contrast(2).titlestr;
matlabbatch{4}.spm.stats.results.conspec(2).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = opt.an.model.contrast(2).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(2).thresh = opt.an.model.contrast(2).thresh;
matlabbatch{4}.spm.stats.results.conspec(2).extent = opt.an.model.contrast(2).extent;
matlabbatch{4}.spm.stats.results.conspec(2).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.print = true;

save ([opt.an.current_model_dir,'/batch_two_sample_ttest_paired.mat'],'matlabbatch');

spm_jobman('run_nogui',matlabbatch);


return 
end

function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end

function anv = prepare_data_one_way_anova(m,l,subbrick)

global opt


levels = size(opt.an.model.one_way_anova.model(m).selector,1);

for ll = 1:levels
    group(ll) = opt.an.model.one_way_anova.model(m).selector(ll,1);
    session(ll) = opt.an.model.one_way_anova.model(m).selector(ll,2);
    subjects{ll} = opt.AUX.subject_names{group(ll)};
end

if opt.group_number == 1
    for ll = 1: levels 
        group_prefix{ll} = '';
    end
else
    for ll = 1: levels
        group_prefix{ll} = opt.AUX.group_names_prefix{group(ll)};
    end
end

if opt.session_number == 1
    for ll = 1: levels
        session_prefix{ll} = '';
    end
else
    for ll = 1: levels 
        session_prefix{ll} = opt.AUX.session_names_prefix{session(ll)};
    end
end

total_rows = 0;
for ll = 1:levels
    scans = cell(opt.subject_number(group(ll)),1);
    for s = 1:opt.subject_number(group(ll))
        scans{s,1} =  [opt.folders.firstlevel,'/',opt.an.measure(l).dir_name,'/',group_prefix{ll},c_b(subjects{ll}(s,:)),session_prefix{ll},'_',opt.an.measure(l).name,'.nii,',num2str(subbrick)];
    end
    icell(ll).scans = scans;
    total_rows = total_rows + length(icell(ll).scans);
end

cov_selctor = opt.an.model.one_way_anova.model(m).covariate;
cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
if ~isempty(cov_selctor)
    for ll = 1:length(cov_selctor)
        if length(opt.an.cov(cov_selctor(ll)).c) ~= total_rows
            error(['Covariate ',num2str(cov_selctor(ll)),' in TSTT model ',num2str(m),' has a wrong number of elements']);
        end
        cov(ll).c = opt.an.cov(cov_selctor(ll)).c;
        cov(ll).cname = opt.an.cov(cov_selctor(ll)).name;
        cov(ll).iCFI = 1;
        cov(ll).iCC = 1;
    end
end
anv.icell = icell;
anv.cov = cov;
        
return
end

function one_way_anova(anv)
global opt
matlabbatch{1}.spm.stats.factorial_design.dir = {opt.an.current_model_dir};
%
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell = anv.icell;
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
%
matlabbatch{1}.spm.stats.factorial_design.cov = anv.cov;%struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};


% TODO use comb to generalize this
if length(anv.icell) == 3
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Global Pattern';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.convec = {
                                                       [1 -1 0
                                                       0 1 -1]
                                                       }';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'g1<g2';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'g2<g3';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
end

if length(anv.icell) == 4
    
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Global Pattern';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.convec = {
                                                       [1 -1 0 0
                                                       0 1 -1 0
                                                       0 0 1 -1]
                                                       }';
    matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'g2<g1';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [1 -1 0 0];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'g3<g1';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.convec = [1 0 -1 0];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'g4<g1';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.convec = [1 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'g3<g2';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.convec = [0 1 -1 0];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'g4<g2';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.convec = [0 1 0 -1 ];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'g4<g3';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.convec = [0 0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

end
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat = {[opt.an.current_model_dir,'/SPM.mat']};
matlabbatch{4}.spm.stats.results.conspec(1).titlestr = opt.an.model.contrast(1).titlestr;
matlabbatch{4}.spm.stats.results.conspec(1).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(1).threshdesc = opt.an.model.contrast(1).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(1).thresh = opt.an.model.contrast(1).thresh;
matlabbatch{4}.spm.stats.results.conspec(1).extent = opt.an.model.contrast(1).extent;
matlabbatch{4}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.conspec(2).titlestr = opt.an.model.contrast(2).titlestr;
matlabbatch{4}.spm.stats.results.conspec(2).contrasts = Inf;
matlabbatch{4}.spm.stats.results.conspec(2).threshdesc = opt.an.model.contrast(2).threshdesc;
matlabbatch{4}.spm.stats.results.conspec(2).thresh = opt.an.model.contrast(2).thresh;
matlabbatch{4}.spm.stats.results.conspec(2).extent = opt.an.model.contrast(2).extent;
matlabbatch{4}.spm.stats.results.conspec(2).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.print = true;

save ([opt.an.current_model_dir,'/batch_anova.mat'],'matlabbatch');

spm_jobman('run_nogui',matlabbatch);

return
end


function [data,N] = build_data_cell(model,l,subbrick,string_eval) %m model, l measure
global opt
N = 0;
comb = eval(['size(model.',string_eval,',1)']);
for ll = 1: comb
    group(ll) = eval(['model.',string_eval,'(ll,1)']);
    session(ll) = eval(['model.',string_eval,'(ll,2)']);
    tmp = {opt.AUX.subject_names{group(ll),:}};
    subjects{ll} = tmp(~cellfun(@isempty,tmp));
    N = size(subjects{ll},2) +N;
end

if opt.group_number == 1
    for ll = 1:comb
        group_prefix{ll} = '';
    end
else
    for ll = 1:comb
        group_prefix{ll} = opt.AUX.group_names_prefix{group(ll)};
    end
end

if opt.session_number == 1
    for ll = 1: comb
        session_prefix{ll} = '';
    end
else
    for ll = 1: comb
        session_prefix{ll} = opt.AUX.session_names_prefix{session(ll)};
    end
end

if ~model.session_average
    count = 0;
    data = cell(N,1);
    for ll = 1:comb
        for s = 1:opt.subject_number(group(ll))
            count = count +1;
            data{count,1} = [opt.folders.firstlevel,'/',opt.an.measure(l).dir_name,'/',group_prefix{ll},subjects{ll}{s},session_prefix{ll},'_',opt.an.measure(l).name,'.nii,',num2str(subbrick)];
        end
    end
    N = count;
elseif   model.session_average && comb > 1 %session average
    data = [];
    output_name = [group_prefix{1},opt.an.measure(l).name];
    for ll = 1:comb
        for s = 1:opt.subject_number(group(ll))
            data{ll,s} = [opt.folders.firstlevel,'/',opt.an.measure(l).dir_name,'/',group_prefix{ll},subjects{ll}{s},session_prefix{ll},'_',opt.an.measure(l).name,'.nii']; %for AFNI
        end
        output_name = [output_name,session_prefix{ll}];
    end
    [data,N] = average_maps(data,output_name,subjects,subbrick);
end
return
end


function [data_new,N] = average_maps (data,output_name,subjects,subbrick) 
global opt

s = size(data);
data_new = [];
for ll = 1:s(2) % subjs
    str1 = [];
    str2 = [];
    for l = 1:s(1) %session number
        str1 = [str1,'-',opt.an.letters{l},' ',data{l,ll},' ']; 
        str2 = [str2,opt.an.letters{l},'+'];
    end
    str2(end) = []; str2 = ['(',str2,')/',num2str(s(1))];
    [~,~]= system(['3dcalc ',str1,'-prefix ',opt.an.current_model_dir,'/',output_name,'_',subjects{l}{ll},'.nii -expr ''(',str2,')''']);
    data_new{ll,1} = [opt.an.current_model_dir,'/',output_name,'_',subjects{l}{ll},'.nii,',num2str(subbrick)];  %for SPM
end

N = s(2);
return
end