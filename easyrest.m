function easyrest(ER_)
%EASYREST   Easy-to-use toolbox for analyzing resting-state data.
%_______________________________________________________________________
% ____   __   ____  _  _ 
%(  __) / _\ / ___)( \/ )
% ) _) /    \\___ \ )  / 
%(____)\_/\_/(____/(__/  
% ____  ____  ____  ____ 
%(  _ \(  __)/ ___)(_  _)  AFNI/SPM based toolbox for resting/steady
% )   / ) _) \___ \  )(    state data processing and analyzing.
%(__\_)(____)(____/ (__)     
%_______________________________________________________________________
%
%   All options must be put in the editable fields of the script.
%   Alternatively you can pass the ER variable as input, i.e.,
%   EASYREST(ER).
%  
%   EASYREST requires AFNI and SPM.
%
%   To install easyrest on your machine add the ER folders in your matlab
%   path (i.e.: addpath(genpath('/home/user/MATLAB/ER')) )
%   
%   An implementation of single subject processing is now available: easyrest_ss.
%   However, you should use easyrest since all its outputs are fully comparable 
%   between subjects. Output from eastrest_ss are not optimized for a second 
%   level analysis.

% daniele.mascali@gmail.com  Developed in 2015

%TODO: 1) To reduce hard disk space requirement save files as nii.gz (afni
%has this functionality) Problem: SPM can't read nii.gz 
% 2) Allow the case of non matching sessions (subjects with different
% number of sessions; see prepare_data function for some extra details).
% 3) Plotting of RtoR must be runnable even without ER_ROIs file.
% 4) removes first timepoints from regressors to have them to mach the
% LFF.nii 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDITABLE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------
% Data definitions
%---------------------------------------
    % FOLDER defintions
ER.folders.data{1} = '/home/mallow/Scrivania/easyrest/DATA/group1';    % index run on groups
%ER.folders.data{2} = '/home/mallow/Scrivania/easyrest/DATA/group2';   % comment extra groups
ER.folders.rois = '/home/mallow/Scrivania/easyrest/rois';              % roi folder for SEED based measures (eg, seed to voxel; only .nii or hdr/img rois are allowed). Mandatory only if seed-based measures are selceted.                                                    
ER.folders.results = '/data/danielem/visanet';                         % result directory in which the new project will be created
ER.folders.project = 'analysis';                                       % name of the project 

    % Wildcard FILE selector
%ER.file_str.wc1_file = 'wc1v_*.nii';              % GM tissue probability or binary mask
%ER.file_str.wc2_file = 'wc2v_*.nii';              % WM tissue probability or binary mask 
%ER.file_str.wc3_file = 'wc3v_*.nii';              % CSF tissue probability or binary mask 
ER.file_str.functional_file{1} = 'bw_avg_RS_1.nii'; % index run on sessions/conditions
ER.file_str.functional_file{2} = 'bw_avg_RS_2.nii'; % comment extra sessions/conditions

    % RP regression
ER.file_str.rp_file{1} = 'rp_mean_rs1.txt';        % Wildcard file selector for rp files. Put a cell for each session/condtion rp files (i.e., index run on sessions/conditions).
ER.file_str.rp_file{2} = 'rp_mean_rs2.txt';        % comment extra session/condition.

    % subjects-,condition-, group-name definition
ER.file_str.subject_names = 'VA*';                 % Non ?? fondamentale, se non specificato i soggetti vengono nominati con subj_###.
ER.file_str.condition_names = {'LEFT','RIGHT'};    % Condition name of each scan (as many as the number of sessions). Mandatory if # scans > 1
%ER.file_str.group_names = {'AD','HS'};             % group name of data folder (as many as the number of data folders). Mandatory if # data folders > 1. Not applicable for # data folders == 1
ER.file_str.condition_assignment{1} = [1 2; 1 2; 2 1; 1 2; 2 1; 1 2; 2 1; 1 2; 2 1; 2 1];         % If scan names are in sequential order but conditions are randomized use this entry to assign to each subject the right conditon. (rows are subjects, columns are sessions; use integer to indicate the condition_names; e.g.: [1,2;2,1;...]
%ER.file_str.condition_assignment{2} = []; 

%---------------------------------------
% preprocessing options
%---------------------------------------

ER.others.volume_selector = [];           % Optional, put the starting and ending volumes to be considered for the analysis. 
                                          % If volume_selector differs among sessions or groups use volume_selector{group,session} = [N1 N2].
                                          % Index starts from 0. Eg, first 20 volumes; [0 19];
                                          % ATTENTION: Issues using data with different lengths: 1)COH: the # of averaged segments will be different 2) PSDM: this will not work. 3) PSDW: the # of averaged segments will be different
                                          % ATTENTION2: Do not use the volume selector if there are different volume numbers in the same group and session (solution: precut your data)
ER.others.tr = 3.1;                       % Repetion time of EPI. If TR differs among sessions or groups use tr{group,session} = [tr].
ER.others.filter_band = [0.008, 0.09];
ER.others.psc = 0;                             % 1/0 [0]; 1 psc; 0 raw ; DON'T USE IF YOUR DATA IS ALREADY MEAN CENTERED (AT ZERO)
ER.others.psc_save = 0;                        % 1/0 [0]; save psc epi (Usually, you don't need them..so save disk space) 
ER.others.tredrsfc_extraoption = '';           % See AFNI 3dRSFC (or 3dTproject if censoring is on) for extra options (-blur is not permitted any more [replaced by new smoothing function])


ER.others.rp_regression.do = 1;                 % 1/0; 1: for rp regression. NB: No need to have zero-centerd rp, indeed, the mean is removed by (3dRSFC (w0 freq) or 3dTproject). It MIGHT make sense to have no-centered no-detrended rp for denoising analyisis
ER.others.rp_regression.derivate = 1;           % compute first derivates, if not already present
ER.others.rp_type = 'FSL';                      % {'FSL','SPM','AFNI','HCP'}. Specifiy how rp are calculated (it changes the unit of rotations and the ordering). 

ER.others.aCompCor.do = 1;                      % aCompCor method for noise mitigation.It requieres tissue probability maps
                                                % NOTE: aCompCor will perform better on UNsmoothed data. So pass data not smoothed and
                                                % smooth them afterwards
 %Optimized aCompCor parameters, you should not change them
ER.others.aCompCor.rpOrtogonalize = 1;          % [def = 1] extract PCA/mean over data already ortogonalized to rp.
                                                %           In this way the model is maximally predictive.                                                
ER.others.aCompCor.asCONN = 0;                  % [def = 0] CONN extracts PCA over mean regressed data, and the first component is equal to straight average. 
ER.others.aCompCor.TvarianceNormalise = 1;      % [def = 1] whether or not to varaince normalize the data before running PCA ( Bezhadi = 1, CONN = 0);
    
ER.others.aCompCor.ROI(1).tissue = 2;           % 1=GM, 2=WM, 3=CSF, 4=GSR (Global Signal Regression)
ER.others.aCompCor.ROI(1).dime = 5;             % Number of PCA components. 0 for straight average.
ER.others.aCompCor.ROI(1).deri = 0;             % 0/n>1. Regress also the first n derivates. O: no derivatives 
ER.others.aCompCor.ROI(1).erode = 1;            % Erode by 1 voxel to avoid tissue contamination
ER.others.aCompCor.ROI(1).thr = 0.5;            % Threshold of probability maps

ER.others.aCompCor.ROI(2).tissue = 3;           % 1=GM, 2=WM, 3=CSF, 4=GSR (Global Signal Regression)
ER.others.aCompCor.ROI(2).dime = 5;             % Number of PCA components. 0 for straight average.
ER.others.aCompCor.ROI(2).deri = 0;             % 0/n>1. Regress also the first n derivates. O: no derivatives
ER.others.aCompCor.ROI(2).erode = 1;            % Erode by 1 voxel to avoid tissue contamination
ER.others.aCompCor.ROI(2).thr = 0.5;            % Threshold of probability maps

ER.others.aCompCor.ROI(3).tissue = 4;           % 1=GM, 2=WM, 3=CSF, 4=GSR (Global Signal Regression)
ER.others.aCompCor.ROI(3).dime = 0;             % Number of PCA components. 0 for straight average.
ER.others.aCompCor.ROI(3).deri = 1;             % 0/n>1. Regress also the first n derivates. O: no derivatives
ER.others.aCompCor.ROI(3).erode = 0;            % Erode by 1 voxel to avoid tissue contamination
ER.others.aCompCor.ROI(3).thr = 0;              % ATT: in this case this paramenter is override (0 can't be changed)

ER.others.censoring.do = 1;                     % Censoring is not compatible with amplitude/spectral quantities (no ALFF will be produced)                       
ER.others.censoring.mode = 'ZERO';              % ZERO, KILL, NTRP, see afni 3dTproject
ER.others.censoring.value = 0.3;                % FD (Power) above which to censor 
ER.others.censoring.pre_TR = [];                % cens n previous TR 
ER.others.censoring.post_TR = [];               % cens n next TR 

%smoothing
ER.others.FWHM = 4;                                 %FWHM of smoothing. Smoothing is applied inside the choosen brainmask only. Put 0 or empty for no smoothing
ER.others.smoothing_mode = '3dBlurInMask';          %two possible strings: '3dBlurInMask' or '3dBlurToFWHM' (see AFNI help for details)
ER.others.BlurToFWHM_extraoption = '-nbhd NULL';    % Add the option '-nbhd NULL' to avoid estimation of local smoothness (it should preserve it)

%---------------------------------------
% Group mask defintion
%---------------------------------------
ER.others.bm_method = 5;                       % 1: all brain (GM + WM + CSF) from structural. 
                                               % 2: GM only.
                                               % 3: a priori group brainmask (deprecated). 
                                               % 4: from single-subject mask files. No automask procedure will be applied
                                               % 5: 3dAutomask 
                                               % 6: Define a groupbm based only on the EPI time series of all subjs.
                                                  % It removes voxels which timeseries have zero std (usefull when data already masked)
                                               % 7: 3dSkullStrip, useful when the skull is still present in the data   
% case 1 and 2
ER.others.bm.tissue_thr = [0.15 0.15 0.50];    % mandatory for bm_method = 1  default values: [0.15 0.15 0.50]. If you don't want a tissue simply put a threshold > 1 (i.e, [0.90,0.50,2])
ER.others.bm.gm_thr = 0.75;                    % mandatory for bm_method = 2  GM threshold
ER.others.bm.group = [];                       % select groups of subject (for example , 1 or [1,2]) in which perform the bmask computation. Let it empty for selecting all groups (i.e., all subjects) [this is the default behaviour]
ER.others.bm.modality = 'soft';                % you can choose: 'hard' or 'soft'. HARD: threshold and then intersection. SOFT: mean and then threshold. Defautl is 'hard'.
ER.others.bm.doautomask = 0;                   % If you dont perform 3dautomask, simply zero-signal voxels will be removed.  
% case 3:
ER.others.bm.apriori = '';                     % mandatory for bm_method = 3  complete path of a priori brain mask
% case 4:
ER.others.bm.ss_masks = '';                    % mandatory for bm_method = 4  wildecard for single subject mask files (NB they are supposed to be inside data folder) 
ER.others.bm.ss_masks_thr = [];                % Number 0 < x <=1. if not empty, the final brain mask will be computed as the sum of all masks, divided by the number of masks and the thresholded at the specificed values. If empty or commented the final mask will be simply the intersection of all single masks.
% case 5:
ER.others.bm.m5_dilate = 0;                    % for bm_method = 5. 1/0: 1 dilate the mask by 1 voxel. 0 no dilatation 
% other options:
ER.others.bm.automask_session_selector = [];   % Most brainmask method use automask on the EPI series to avoid inclusion of non brain voxels. Use this selector to apply automask only on certain sessions. 
ER.others.bm.save_ss = 1;                      % 1/0; 1: save single subject masks when possible (method 1,2 (hard mode) and 5)
ER.others.bm.use_ss = 0;                       % 1/0; 1: instead of computing a group brainmask a different brainmask is used for each subject. Work for all modalities but bm=3,5 and 7.
                                               %      Be aware that normalized measures (z-score, "m"-score) should be considered only with a group brainmask. 
                                               %      Also, voxel wise comparison is discourage with single subject brainmasks. 
                                               %      When comparing ROIs be sure that the number of voxels are similar.
%---------------------------------------
% Measure selector 
%---------------------------------------                  

%------DENOSING ANALYSIS (DA)-------
ER.others.measure.DA = 1;                  % Performs an extensive DENOISING analysis (DA). Reallingnment parameters must be present.
                                           % It requires GM maps
ER.others.DA.force_overwrite = 0;          % if 0, it doesn't recompute already calculated subjects/sessions.
%-----------------------------------

% 3dRSFC (i.e., ALFF) is mandatory. In addition you can compute
ER.others.measure.RSFA = 0;           % 0) RSFA: Resting-state Fluctuation Amplitude (equivalent to ALFF). Normally, you don't need this quantity.
ER.others.measure.REHO = 0;           % 1) REHO: Regional Homogeneity [ref: Regional homogeneity approach to fMRI data analysis by Zang, Jiang, Lu, He, and Tiana (2004, NeuroImage)] 
ER.others.measure.gFC  = 0;           % 2) gFC: Global Functional Connectivity (this can take a VERY LONG TIME)
ER.others.measure.CoE  = 0;           % 3) CoE: Correlation Entropy. Computes entropy of the global correlation histogram (requires gFC)   
ER.others.measure.DeC  = 0;           % 4) DeC: Degree Connectivity (graph theory). Using each voxel as a node an adjacency matrix (undirected and unweighted) is difened. The deegre is defined as the sum of the thresholded adjacency matrix. One or more thresholds can be defined in ER.others.DeC.Thrs
ER.others.measure.SpE  = 0;           % 5) SpE: Spectral Entropy  
        % seed-based measure 
ER.others.measure.StoV = 0;          % 6) SEED-to-VOXEL connectivity (fastconn). YOU MUST PUT THE ROIS IN folder.rois (only .nii or .hdr/img are allowed)
ER.others.measure.DEL  = 0;           % 7) DEL: 3dDelay (a seed-based measure). YOU MUST PUT THE ROIS IN folder.rois (only .nii or .hdr/img are allowed)  
ER.others.measure.COH  = 0;           % 8) COH: Coherency (a seed-based measure). YOU MUST PUT THE ROIS IN folder.rois (only .nii or .hdr/img are allowed)  
        % ROI to ROI
ER.others.measure.RtoR = 0;          % 9) All scans must have the same number of volumes
        % Spectral measure:
ER.others.measure.PSD  = 0;           % 10) PSD: Power Spectral Density (periodogram) (output name _PS.nii) 
ER.others.measure.PSDW = 0;          % 11) PSDW: Power Spectral Density estimated with Welch's method (it reduces noise in the estimated power spectra in exchange for reducing the frequency resolution) (output name _PSW.nii)     
ER.others.measure.PSDM = 0;          % 12) PSDM: Power Spectral Density estimated with multitaper method (it reduces noise in the estimated power spectra) (output name _PSM.nii)     
        % Special alff estimation:
ER.others.measure.ALFFW = 0;         % 13) ALFFW: ALFF computed with Welch's method for spectral estimation
ER.others.measure.ALFFM = 0;         % 14) ALFFM: ALFF computed with multitaper for spectral estimation

%---------------------------------------
% Measure options
%---------------------------------------  
ER.others.new_rois = 0;              % 1/0; 1: In the case you add new Seed in the roi folders (only if you re-run easyrest on the same project) 
ER.others.gFC.Vthr =[0.25 0.5 0.05]; % [t0,t1,dt] Save the COUNT of how many voxels survive thresholding at several levels abs(correlation) >= thr, for thr = t0, t0+dt, ..., t1. [0.25 0.5 0.05].  
ER.others.DeC.Thrs =[0.25 0.5 0.05]; % [t0,t1,dt] For DeC. Save the COUNT of how many voxels survive thresholding at several levels correlation >= thr, for thr = t0, t0+dt, ..., t1. [0.25 0.5 0.05].  
ER.others.coh.band_division = [];    % Specify frequency bands for band-averaged coherence and delay (boh from coherency). Add multiple raws for multi-band analysis. If empty, only the filter band of the preprocessing step will be used.
ER.others.coh.window_length = [];    % Specify the number of time points (in volumes) for welch's method segmentation (for COHERENCE). You should let it empty for default segmentation (i.e., npoints = 1/(tr*filter_band(1)) ). Change the default only if you don't have enough time points (or if you want a better spectral density estimation to the detriment of low frequency discrimination)  
ER.others.psdw.window_length = [];   % Specify the number of time points (in volumes) for welch's method segmentation (for PSD estimation). You should let it empty for default segmentation (i.e., npoints = 1/(tr*filter_band(1)) ). Change the default only if you don't have enough time points (or if you want a better spectral density estimation to the detriment of low frequency discrimination)  
ER.others.psdm.dpss_nw = 3;          % time-halfbandwidth product (nw) for the discrete prolate spheroidal (Slepian) sequences (dpss). This number define the final spectral risolution. Typical values are: 2,3 or 4.   
ER.others.psdm.dpss_k = [];          % number of dpss sequence to use (ie. number of tapers) (K must be a positive integer in the range 1:N). You should let it empty for default: k = 2nw-1. greater k are usually less energy concetred.
ER.others.RtoR.roi_folder_or_file = ''; % for ROI to ROI connectivity. Two modalities: 
                                        %    1) Put the path of a directory containing only roi files (both .nii or .hdr); a mask dataset will be constructed from these files. 
                                        %    2) Put the path of a mask data file (i.e., a 3d dataset with values from 1 to N, describing N different ROIs). If a txt file with the same name of the mask data file is present the labels will be taken form this file.
ER.others.RtoR.name = '';               % The name of the ROI to ROI connectivity matrix. Change the name to avoid overwriting
ER.others.RtoR.COH = 1;             % Compute also ROI to ROI connectivity using magnitude of coherency (at the momement works only for one band (i.e., the filter band))

%smoothing order (ONLY for REHO, gFC, DeC)
ER.others.special_smoothing_order = 1;  %[-1,0,1,2]:-1: Use smoothed timeseries (if selected in preprocessing)
                                                   % 0: Use UNsmoothed timeseries
                                                   % 1: Smooth after measure computation (FWHM set below option)
                                                   % 2: (not implemented yet) As 1, but preserve unsmoothed maps
ER.others.special_smoothing_FWHM = [];  % FWHM for map smoothing in case 1. Empty value assume the same smoothing            
%additional note on smoothing: 
% if smoothing is selected during preprocessing:
%   StV is computed on the smoothed data but the seed average timeseries is
%   extracted from UNsmoothed data.
%   RtoR is computed from UNsmoothed data
%   if not specified other measures are computed on smoothed data (if
%   available)

%---------------------------------------
% Advanced options
%--------------------------------------- 
ER.others.The_IKnowWhatIamDoing_Condition = 0; % 1/0; if 1, there is no pause for checking data arrangement (generally, deprecated ...
                                               % unless you want to use nohup (e.g.: nohup matlab -nojvm -r name_of_your_script_without_extension -logfile name_of_logfile.out </dev/null & )
ER.others.update = 1;                          % 1/0; automatic check for update. DO NOT DISABLE                                             
ER.others.sendstatus = '';                     % Specify an e-mail address where to send the exit status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT BELOW, UNLESS YOU KNOW WHAT YOU ARE DOING!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TODO: add check on data dimensionality (wc1, functional)

if nargin
    %fprintf('\nReading options from input variable...');
    ER = ER_;
    %fprintf('done.\n');
end

% sendstatus if required
if isfield(ER.others,'sendstatus') && ~isempty(ER.others.sendstatus)
    sendstatus(ER.others.sendstatus)
end

clearvars -global opt
global opt
global DA_extracted
DA_extracted = [];

wd = pwd;

warning('off', 'MATLAB:MKDIR:DirectoryExists');

version = version_er;      
%%---------------------------------
% Welcome screen
welcome_screen(version)
%%---------------------------------

%%---------------------------------
% Check for update
if ER.others.update && ~ER.others.The_IKnowWhatIamDoing_Condition 
    if exist('check_updates')~=2
        error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''/home/user/MATLAB/ER''))');
    end
    exit = check_updates(version);
    if exit;return;end
end
%%---------------------------------

%%---------------------------------
%Set AFNI env variables
setenv('AFNI_AUTOGZIP','NO');
setenv('AFNI_ENVIRON_WARNINGS','NO');
%%---------------------------------


%%---------------------------------
%   Prepare data
if exist([ER.folders.results,'/ER_',ER.folders.project,'/ER_',ER.folders.project,'.mat'],'file')
    load ([ER.folders.results,'/ER_',ER.folders.project,'/ER_',ER.folders.project,'.mat'])
    if opt.AUX.preprocessing_done
        redo_mode_prepare(ER.folders,ER.others)
    else
        opt = [];
        opt.ver = version;
        exit = prepare_data(ER.folders,ER.file_str,ER.others);
        if exit;cd(wd);return;end
    end
else
    exit = prepare_data(ER.folders,ER.file_str,ER.others);
    if exit;cd(wd);return;end
end
%%---------------------------------

%%---------------------------------
project_summary;
%%---------------------------------

cd (opt.folders.preprocessing)

%%---------------------------------
%   RP initialization
    % this step might exit with error; so it's been put as first step, also
    % considering that it doesn't require brainmask
if ~opt.AUX.preprocessing_done && (opt.rp_regression.do || opt.censoring.do)
    if opt.rp_regression.do == 0
        % no need for derivative
        opt.rp_regression.derivate = 0;
    end
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Realignment parameters initialization\n']);
    total = sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count + 1;
                rp_initialization(g,j,k)
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,opt.AUX.output_name{g,j,k});
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    str = [];mean_ = []; std_ = [];
    str = [str,sprintf('\n-> Summary of Framewise displacement:')];
    str = [str,sprintf('\nGroup\tSession\t\tValues\n')];
    for g =1:opt.group_number
        %str = [str,sprintf('Group%d\t',g)];
        for j=1:opt.session_number
            str = [str,sprintf('Group%d\tSess%d:\t\t',g,j)];
            for k=1:opt.subject_number(g)
                str = [str,sprintf('%.3f\t',opt.QC.FD_mean{g,j}(k,1))];
            end
            mean_(g,j) = mean(opt.QC.FD_mean{g,j});
            std_(g,j) = std(opt.QC.FD_mean{g,j});
            str = [str,sprintf('\n')];
        end
    end
    str = [str,sprintf('-> Stats:\n')];
    for g =1:opt.group_number
        str = [str,sprintf('Group%d',g)];
        for j=1:opt.session_number
            str = [str,sprintf('\tSess%d: %.3f+/-%.3f\t',j,mean_(g,j),std_(g,j))];
        end
        str = [str,sprintf('\n')];
    end
    % actual printing
    fprintf(str);
    opt.AUX.log = strvcat(opt.AUX.log,str);
    fprintf(opt.AUX.fid,'%s\n',str');
    %
    opt.QC.FD_mean_mean = mean_;
    opt.QC.FD_mean_std  =  std_;         
    clear mean_ std_;
    fprintf('\n');
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
%   Censoring initialization
    %does not require brainmask
if ~opt.AUX.preprocessing_done && opt.censoring.do
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Censoring initialization\n']);
    total = sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
%     if isfield(opt.DATA,'RP_1D')
%         cen.data = opt.DATA.RP_1D;
%     else
%         cen.data = opt.DATA.RP; %this is for the case rp regression is set to zero, but RP are provided
%     end
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count + 1;
                censoring_initialization(g,j,k)
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,'Temporal mask definition...');
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    % PRINT STATS 
    str = [];mean_ = []; std_ = [];
    str = [str,sprintf('\n-> Summary of censored time points:')];
    str = [str,sprintf('\nGroup\tSession\t\tValues\n')];
    for g =1:opt.group_number
        str = [str,sprintf('Group%d',g)];
        for j=1:opt.session_number
            str = [str,sprintf('\tSess%d:\t\t',j)];
            for k=1:opt.subject_number(g)
                str = [str,sprintf('%d\t',opt.censoring.censor_number{g,j}(k,1))];
            end
            mean_(g,j) = mean(opt.censoring.censor_number{g,j});
            std_(g,j) = std(opt.censoring.censor_number{g,j});
            str = [str,sprintf('\n')];
        end
    end
    str = [str,sprintf('-> Stats:\n')];
    for g =1:opt.group_number
        str = [str,sprintf('Group%d\t',g)];
        for j=1:opt.session_number
            str = [str,sprintf('sess%d: %.2f+/-%.2f\t',j,mean_(g,j),std_(g,j))];
        end
        str = [str,sprintf('\n')];
    end
    % actual printing
    fprintf(str);
    opt.AUX.log = strvcat(opt.AUX.log,str);
    fprintf(opt.AUX.fid,'%s\n',str');
    %
    opt.censoring.censor_number_mean = mean_;
    opt.censoring.censor_number_std  =  std_;         
    clear cen mean_ std_;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
%   Group BRAINMASK definition
if ~opt.AUX.preprocessing_done 
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Definition of group brain mask\n']);
    brainmask(ER.others.bm_method);
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   Group BRAINMASK interection with ROIs
if opt.AUX.seedmeasure && ~opt.AUX.bmi_done
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Intersection of group brain mask with ROIs\n']);
    switch opt.bm.use_ss
        case 0 
            brainmask_ROI_intersection;
        case 1
            brainmask_ROI_intersection_SS; %in this modality no bmi file is produeced. This step only serves to know how many voxels wil survay the intersection for each subject.
    end
    opt.AUX.bmi_done = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
%   Group BRAINMASK interection with ROIs for ROI TO ROI 
if opt.AUX.TODO_RtoR 
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Intersection of group brain mask with ROIs for ROI to ROI connectivity\n']);
    switch opt.RtoR.modalitiy 
        case 'roi_list'
            brainmask_ROI_intersection_RtoR;
        case 'mask_file'
            brainmask_ROI_intersection_RtoR_mf;
    end
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

cd (opt.folders.preprocessing)

%%---------------------------------
%   CompCor Calculation
if opt.aCompCor.do && ~opt.AUX.preprocessing_done
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': aCompCor calculation\n']);
    total = sum(opt.subject_number) + sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
    [~,~] = system('rm _aCompCor_*.1D _temp.txt');   % in case of previous failure
    for g =1:opt.group_number
        for k=1:opt.subject_number(g)
            count = count + 1;
            acompcor_mask_definition(g,k)
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,'Mask definition...');
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count + 1;
                acompcor(g,k,j);
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,['PCA decomposition, ',opt.AUX.output_name{g,j,k}]);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    fprintf('\n');
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
% Create Xort regressors
if ~opt.AUX.preprocessing_done
    switcher = num2str([opt.rp_regression.do, opt.aCompCor.do]);
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
            str_out = ['_Xort_',opt.AUX.output_name{g,j,k},'.1D'];
            switch switcher
                case {'0  0'}
                % no Xort
                case {'0  1'}
                    system(['1dcat ',c_b(opt.aCompCor.aCC1D{g,j,k}),' > ',str_out]);
                case {'1  0'}
                    system(['1dcat ',c_b(opt.DATA.RP_1D{g,j,k}),' > ',str_out]);
                case {'1  1'}
                    system(['1dcat ',c_b(opt.DATA.RP_1D{g,j,k}),' ',c_b(opt.aCompCor.aCC1D{g,j,k}),' > ',str_out]);
            end
            opt.X1D{g,j,k} = [opt.folders.preprocessing,'/',str_out];
            end
        end
    end
end
%%---------------------------------

%%---------------------------------
%   3dRSFC
if ~opt.AUX.preprocessing_done && strcmp(opt.prepro_mode,'3dRSFC')
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': 3dRSFC (preprocessing + ALFF computation)\n']);

    mkdir(opt.folders.firstlevel,'ALFF');
    dest_folder = [opt.folders.firstlevel,'/ALFF'];

    total = sum(opt.subject_number*opt.session_number)+ 2 + opt.AUX.do_smoothing;
    reverseStr  = '';
    count = 0;
    if opt.psc
        opt.AUX.bmask_tmp_hdr = spm_vol(opt.AUX.bmask_path);
        opt.AUX.bmask_tmp = uint8(spm_read_vols(opt.AUX.bmask_tmp_hdr));
        opt.AUX.bmask_vxl_err = 0;
    end
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count + 1;
                do_3dRSFC(g,k,j);
                list_head = dir('*.HEAD');
                for z=1:length(list_head)
                    [status,string] = system(['3dAFNItoNIFTI ',list_head(z).name]);control(status,string);
                    [status,string] = system(['rm ',list_head(z).name,' ',list_head(z).name(1:(end-4)),'BRIK']);control(status,string);
                end
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,opt.AUX.output_name{g,j,k});
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    %-------Smooting maps--------------
    if opt.AUX.do_smoothing
        count = count + 1;
        percentDone = 100 * (count) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['Smoothing amplitude-based maps...']);  
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        smooth_map_ALFF('_ALFF.nii',g,k); 
        smooth_map_ALFF('_mALFF.nii',g,k); 
        smooth_map_ALFF('_fALFF.nii',g,k); 
        if opt.AUX.TODO_RSFA
            smooth_map_ALFF('_RSFA.nii',g,k);
            smooth_map_ALFF('_mRSFA.nii',g,k); 
            smooth_map_ALFF('_fRSFA.nii',g,k); 
        end
    end
    %----------------------------------
        count = count + 1;
        percentDone = 100 * (count) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['z-score conversion...']);  
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
%    bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    zscore_conversion('_ALFF','_zALFF','_mALFF-1');    %additional mALFF-1 output
    zscore_conversion('_fALFF','_zfALFF'); 
    [status,string] = system(['mv *_ALFF.nii *_fALFF.nii *_mALFF.nii *-1.nii *zALFF.nii *zfALFF.nii ',dest_folder]);control(status,string);
    if opt.AUX.TODO_RSFA
        zscore_conversion('_RSFA','_zRSFA','_mRSFA-1');    %additional mRSFA-1 output
        zscore_conversion('_fRSFA','_zfRSFA'); 
        [status,string] = system(['mv *_fRSFA.nii *_mRSFA.nii *_RSFA.nii *-1.nii *zRSFA.nii *zfRSFA.nii ',dest_folder]);control(status,string);
        opt.AUX.TODO_RSFA = 0;
        opt.MEASURES.RSFA = 1;
    end
        count = count + 1;
        percentDone = 100 * (count) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['All done!']);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf('\n');
    if opt.psc
        opt.AUX.bmask_vxl_err(1) = [];
        if ~isempty(opt.AUX.bmask_vxl_err)
            occurrence = length(opt.AUX.bmask_vxl_err);
            str = ['Group-brain mask has been reshaped ',num2str(occurrence),' time(s) to allow for correct percent signal change transformation. Voxel lost each time: '];
            for zz = 1:occurrence 
                str = [str,num2str(opt.AUX.bmask_vxl_err(zz)),' '];
            end
            warning off backtrace
            warning(str);
            warning on backtrace
        end
        fields = {'bmask_tmp_hdr','bmask_tmp'};
        rmfield(opt.AUX,fields);
    end
    %opt.AUX.preprocessing_done = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
%   3dTproject
if ~opt.AUX.preprocessing_done && strcmp(opt.prepro_mode,'3dTproject') 
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': 3dTproject (preprocessing with censoring)\n']);

    total = sum(opt.subject_number*opt.session_number)+1;
    reverseStr  = '';
    count = 0;
    if opt.psc
        opt.AUX.bmask_tmp_hdr = spm_vol(opt.AUX.bmask_path);
        opt.AUX.bmask_tmp = uint8(spm_read_vols(opt.AUX.bmask_tmp_hdr));
        opt.AUX.bmask_vxl_err = 0;
    end
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count + 1;
                do_3dTproject(g,k,j)
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,opt.AUX.output_name{g,j,k});
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
%     %------convert to nifti from brick
%         percentDone = 100 * (count +1) / total;
%         msg = sprintf('Percent done: %3.1f (%s)', percentDone,'AFNItoNIFTI conversion...');
%         fprintf([reverseStr, msg]);
%         reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     list_head = dir('*.HEAD');    
%     for j=1:length(list_head);
%         [status,string] = system(['3dAFNItoNIFTI ',list_head(j).name]);control(status,string);
%         [status,string] = system(['rm ',list_head(j).name,' ',list_head(j).name(1:(end-4)),'BRIK']);control(status,string);
%     end
%     %----------------------------------
        percentDone = 100 * (count +1) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,'All done!');
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        fprintf('\n');
    if opt.psc
        opt.AUX.bmask_vxl_err(1) = [];
        if ~isempty(opt.AUX.bmask_vxl_err)
            occurrence = length(opt.AUX.bmask_vxl_err);
            str = ['Group-brain mask has been reshaped ',num2str(occurrence),' time(s) to allow for correct percent signal change transformation. Voxel lost each time: '];
            for zz = 1:occurrence 
                str = [str,num2str(opt.AUX.bmask_vxl_err(zz)),' '];
            end
            warning off backtrace
            warning(str);
            warning on backtrace
        end
        fields = {'bmask_tmp_hdr','bmask_tmp'};
        rmfield(opt.AUX,fields);
    end
    % PRINT DOF STATS 
    str = [];mean_ = []; std_ = [];
    str = [str,sprintf('\n-> Degrees of freedom:')];
    str = [str,sprintf('\nGroup\tSession\t\tValues\n')];
    for g =1:opt.group_number
        str = [str,sprintf('Group%d',g)];
        for j=1:opt.session_number
            str = [str,sprintf('\tSess%d:\t\t',j)];
            for k=1:opt.subject_number(g)
                str = [str,sprintf('%d\t',opt.censoring.DOF{g,j}(k,1))];
            end
            mean_(g,j) = mean(opt.censoring.DOF{g,j}(:,1));
            std_(g,j) = std(opt.censoring.DOF{g,j}(:,1));
            str = [str,sprintf('\n')];
        end
    end
    str = [str,sprintf('-> Stats:\n')];
    for g =1:opt.group_number
        str = [str,sprintf('Group%d\t',g)];
        for j=1:opt.session_number
            str = [str,sprintf('sess%d: %.2f+/-%.2f\t',j,mean_(g,j),std_(g,j))];
        end
        str = [str,sprintf('\n')];
    end
    % actual printing
    fprintf(str);
    opt.AUX.log = strvcat(opt.AUX.log,str);
    fprintf(opt.AUX.fid,'%s\n',str');
    %
    opt.censoring.DOF_mean = mean_;
    opt.censoring.DOF_std =  std_;         
    clear mean_ std_;
    %opt.AUX.preprocessing_done = 1; 
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
end
%%---------------------------------

%%---------------------------------
%   SMOOTHING
if ~opt.AUX.preprocessing_done && opt.AUX.do_smoothing
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Smoothing ',num2str(opt.FWHM),' mm (',opt.smoothing_mode,')\n']);
 
    total = sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count +1;
                base_name = opt.AUX.output_name{g,j,k};
                switch opt.smoothing_mode
                    case {'3dBlurToFWHM'}
                        %it doesn't need the preserve option
                        [status,string] = system([opt.smoothing_mode,' -input ',base_name,'_LFF.nii -FWHM ',num2str(opt.FWHM),' -prefix ',base_name,'_sLFF.nii ',brainmask_selector(g,k,2),opt.AUX.BlurToFWHM_extraoption,' -nodetrend']);control(status,string);
                    case {'3dBlurInMask'}
                        [status,string] = system([opt.smoothing_mode,' -input ',base_name,'_LFF.nii -float -FWHM ',num2str(opt.FWHM),' -prefix ',base_name,'_sLFF.nii ',brainmask_selector(g,k,2),' -preserve']);control(status,string);
                end
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    fprintf('\n');
end
opt.AUX.preprocessing_done = 1; 
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');

%%---------------------------------
%   DENOISING ANALYSIS
if opt.AUX.TODO_DA
    close all; % In SS mode a figure is open, it interferes with plotting rusults.
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Denoising Analysis\n']);

    mkdir(opt.folders.preprocessing,'denoising_analysis');
    dest_folder = [opt.folders.preprocessing,'/denoising_analysis'];  %main folder with .mat
    opt.folders.denoising_analysis = dest_folder;
    mkdir(opt.folders.denoising_analysis,'results') %folder with plots of the results
    opt.folders.denoising_results = [opt.folders.denoising_analysis,'/results'];
    mkdir(opt.folders.denoising_results,'specificity');  %folder with plots of the results of specificity
    opt.folders.denoising_results_spe = [opt.folders.denoising_results,'/specificity'];
    mkdir(opt.folders.denoising_results,'roi_to_roi');  %folder with plots of the results of specificity
    opt.folders.denoising_results_rtr = [opt.folders.denoising_results,'/roi_to_roi'];
    cd (dest_folder);
    total = sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
    %----------------------------------------------------------------------
    % define colors for various regressions and others varius stuff
    opt.DA.colors = [0 0 0; ...
          204 0 0; ...
          76 153 0; ...
          0 153 153; ...
          76 0 153];
    opt.DA.colors = opt.DA.colors./255;
    opt.DA.fd_ylim = [0 2]; %mm
    %----------------------------------------------------------------------
    % For specificity computation lets define a couple (or three ROIs to
    % compare)
    da_specificity_roi_definition;
    %----------------------------------------------------------------------
    % For ROI-to-ROI distance define atlas (likely the one of Power)
    prepare_atlas;
    %----------------------------------------------------------------------
    opt.DA_extracted = [];
    opt.AUX.DA_skip_SS = 0;      % this is a debbugin variable. Must be 0. Use 1 if you have already extracted SS analysis.
    if ~opt.AUX.DA_skip_SS
        opt.DA.fid = fopen('ER_DA_voxel_count.txt','w');
        fprintf(opt.DA.fid,'%s\t\t','GM');
        for s = 1:length(opt.DA.roi_specificity)
            for l = 1:length(opt.DA.roi_specificity(s).roi)
                fprintf(opt.DA.fid,'%s\t',opt.DA.roi_specificity(s).roi(l).name_contract);
            end
            fprintf(opt.DA.fid,'\t');
        end
    end
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count +1;
                base_name = opt.AUX.output_name{g,j,k};
                denoising_analysis(g,k,j)
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    if ~opt.AUX.DA_skip_SS
        fclose(opt.DA.fid);
    end
    % load extracted data  (fast...no need for %).
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                denoising_analysis_extractSS(g,k,j);
            end
        end
    end
    %now lets do the 2nd level analysis
    %denoising_analysis_plot_2ndl;
    opt.AUX.TODO_DA = 0;
    opt.MEASURES.DA = 1;
    opt.DA.extracted_data_path = [opt.folders.results,'/ER_',opt.folders.prject_name,'_DA_extracted.mat'];
    save(opt.DA.extracted_data_path ,'DA_extracted','-v7.3');
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
    cd (opt.folders.preprocessing);
end

%%---------------------------------

%%---------------------------------
%   3dREHO
if opt.AUX.TODO_REHO
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Regional Homogeneity (ReHo) computation\n']);

    mkdir(opt.folders.firstlevel,'ReHo');
    dest_folder = [opt.folders.firstlevel,'/ReHo'];
    [~,~] = system(['rm *_REHO19.nii *_zREHO19.nii *_REHO7.nii *_zREHO7.nii *_REHO27.nii *_zREHO27.nii']);   % in case of previous failure
    if opt.AUX.special_smoothing_order < 0
        LFF_string = '_sLFF.nii';
    else
        LFF_string = '_LFF.nii';
    end
    total = sum(opt.subject_number*opt.session_number) + 1 + sum(opt.AUX.special_smoothing_order(opt.AUX.special_smoothing_order > 0));
    reverseStr  = '';
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count +1;
                base_name = [opt.AUX.output_name{g,j,k}];
                input_name = [base_name,LFF_string];
                [status,string] = system(['3dReHo -prefix ',base_name,'_REHO19.nii -inset ',input_name,' -nneigh 19 -chi_sq ',brainmask_selector(g,k,2)]);control(status,string);
                [status,string] = system(['3dReHo -prefix ',base_name,'_REHO7.nii -inset ',input_name,' -nneigh 7 -chi_sq ',brainmask_selector(g,k,2)]);control(status,string);
                [status,string] = system(['3dReHo -prefix ',base_name,'_REHO27.nii -inset ',input_name,' -nneigh 27 -chi_sq ',brainmask_selector(g,k,2)]);control(status,string);
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    if opt.AUX.special_smoothing_order == 1
            count = count +1;
            percentDone = 100 * (count) / total;
            msg = sprintf('Percent done: %3.1f (%s)', percentDone,['Smoothing maps...']);  
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            smooth_map('_REHO19.nii',g,k); 
            smooth_map('_REHO7.nii',g,k); 
            smooth_map('_REHO27.nii',g,k); 
    end
        count = count +1;
        percentDone = 100 * (count) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['z-score conversion...']);  
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    zscore_conversion('_REHO19','_zREHO19'); 
    zscore_conversion('_REHO7','_zREHO7'); 
    zscore_conversion('_REHO27','_zREHO27'); 
    [status,string] = system(['mv *_REHO19.nii *_zREHO19.nii *_REHO7.nii *_zREHO7.nii *_REHO27.nii *_zREHO27.nii ',dest_folder]);control(status,string); 
    opt.AUX.TODO_REHO = 0;
    opt.MEASURES.REHO = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   Power Spectral Densitiy (PS) (a crude estimate)
if opt.AUX.TODO_PSD
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Power Spectral Density (PSD) computation\n']);
    mkdir(opt.folders.firstlevel,'PSD');
    dest_folder = [opt.folders.firstlevel,'/PSD'];
    [~,~] = system('rm *_PSD.nii');   % in case of previous failure
    if opt.AUX.do_smoothing   %use smoothed data if possible
        LFF_string = '_sLFF.nii';
    else
        LFF_string = '_LFF.nii';
    end
    total = sum(opt.subject_number*opt.session_number);
    reverseStr  = '';
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count +1;
                base_name = [opt.AUX.output_name{g,j,k}];
                input_name = [base_name,LFF_string];
                [status,string] = system (['3dPeriodogram -nfft ',num2str(opt.AUX.psd_nfft),' -taper 0 -prefix ',base_name,'_PSD.nii ',input_name]);control(status,string);
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end  
        end
    end
    [status,string] = system(['mv *_PSD.nii ',dest_folder]);control(status,string);
    opt.AUX.TODO_PSD = 0;
    opt.MEASURES.PSD = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   Spectral Entropy (SpE)
if opt.AUX.TODO_SpE
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Spectra Entropy (SpE) computation\n']);
    mkdir(opt.folders.firstlevel,'SpEn');
    dest_folder = [opt.folders.firstlevel,'/SpEn'];
    cd ([opt.folders.firstlevel,'/PSD'])
    PS_list = dir(['*_PSD.nii']);
    total = length(PS_list);
    reverseStr  = '';
    for i=1:length(PS_list)
        spectral_entropy(PS_list(i).name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,PS_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
    [status,string] = system(['mv *_SpEn.nii ',dest_folder]);control(status,string); 
    cd (opt.folders.preprocessing)
    opt.AUX.TODO_SpE = 0;
    opt.MEASURES.SpE = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   Power Spectral Densitiy (PSW) (with Welch's method)
if opt.AUX.TODO_PSDW
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Power Spectral Density (Welch''s Method - PSDW) computation\n']);
    mkdir(opt.folders.firstlevel,'PSDW');
    dest_folder = [opt.folders.firstlevel,'/PSDW'];
    LFF_list = dir(['*_LFF.nii']);
    total = length(LFF_list);
    reverseStr  = '';
    for i=1:length(LFF_list)
        output_name = [LFF_list(i).name(1:end-7),'PSDW.nii'];
        psd_welch(LFF_list(i).name,output_name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,LFF_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
    [status,string] = system(['mv *_PSDW.nii ',dest_folder]);control(status,string);
    [status,string] = system('rm *_PSDW.mat');control(status,string);
    opt.AUX.TODO_PSDW = 0;
    opt.MEASURES.PSDW = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   ALFFW (ALFF with Welch's method)
if opt.AUX.TODO_ALFFW
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': ALFF (Welch''s Method - ALFFW) computation\n']);
    mkdir(opt.folders.firstlevel,'ALFFW');
    dest_folder = [opt.folders.firstlevel,'/ALFFW'];
    cd ([opt.folders.firstlevel,'/PSDW'])
    PSDW_list = dir('*_PSDW.nii');
    total = length(PSDW_list);
    reverseStr  = '';
    bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    for i=1:length(PSDW_list)
        output_name = [PSDW_list(i).name(1:end-8)];
        ALFF_welch(PSDW_list(i).name,bmask,output_name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,PSDW_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
    [status,string] = system(['mv *ALFFw.nii *mALFFw-1.nii ',dest_folder]);control(status,string);
    cd (opt.folders.preprocessing)
    opt.AUX.TODO_ALFFW = 0;
    opt.MEASURES.ALFFW = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   Power Spectral Densitiy (PSM) (with multitaper method)
if opt.AUX.TODO_PSDM
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Power Spectral Density (Multitaper Method - PSDM) computation\n']);
    mkdir(opt.folders.firstlevel,'PSDM');
    dest_folder = [opt.folders.firstlevel,'/PSDM'];
    LFF_list = dir('*_LFF.nii');
    total = length(LFF_list);
    reverseStr  = '';
    for i=1:length(LFF_list)
        output_name = [LFF_list(i).name(1:end-7),'PSDM.nii'];
        psd_multitaper(LFF_list(i).name,output_name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,LFF_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
    [status,string] = system(['mv *_PSDM.nii ',dest_folder]);control(status,string);
    [status,string] = system('rm *_PSDM.mat'); % no control here, beacause it's not essential, if it doesnt find mat file job can proceed
    opt.AUX.TODO_PSDM = 0;
    opt.MEASURES.PSDM = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   ALFFM (ALFF with multitaper method)
if opt.AUX.TODO_ALFFM
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': ALFF (Multitaper Method - ALFFM) computation\n']);
    mkdir(opt.folders.firstlevel,'ALFFM');
    dest_folder = [opt.folders.firstlevel,'/ALFFM'];
    cd ([opt.folders.firstlevel,'/PSDM'])
    PSDM_list = dir('*_PSDM.nii');
    total = length(PSDM_list);
    reverseStr  = '';
    bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    for i=1:length(PSDM_list)
        output_name = PSDM_list(i).name(1:end-8);
        ALFF_multitaper(PSDM_list(i).name,bmask,output_name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,PSDM_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
    [status,string] = system(['mv *ALFFm.nii *mALFFm-1.nii ',dest_folder]);control(status,string);
    cd (opt.folders.preprocessing)
    opt.AUX.TODO_ALFFM = 0;
    opt.MEASURES.ALFFM = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end

%%---------------------------------

%%---------------------------------
%   3dDELAY 
if opt.AUX.TODO_DEL && opt.AUX.seedmeasure
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Delay (3ddelay) computation\n']);
    
    mkdir(opt.folders.firstlevel,'DELAY');
    dest_folder = [opt.folders.firstlevel,'/DELAY'];    
   
    LFF_list = dir('*_LFF.nii');
    
    total = length(LFF_list)*size(opt.DATA.ROIS,1);
    reverseStr  = '';
    count = 0;

    for i=1:length(LFF_list)   %ciclo subjs
        for l=1:size(opt.DATA.ROIS,1)   %ciclo rois
           count = count + 1;
           roi_str =  c_b(opt.AUX.roi_names(l,:));
           output_name = [LFF_list(i).name(1:end-7),'DEL_',roi_str(1:(end-4))];
           do_3ddelay(LFF_list(i).name,c_b(opt.DATA.ROIS(l,:)),output_name);
           
           percentDone = 100 * count / total;
           msg = sprintf('Percent done: %3.1f (%s)', percentDone,output_name);
           fprintf([reverseStr, msg]);
           reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    [~,~] = system(['rm ',dest_folder,'/*']);       % in case you are performing a second time a seed-based measure
    system(['mv *DEL_* ',dest_folder]);  
    opt.AUX.TODO_DEL = 0;
    opt.MEASURES.DEL = 1;
    opt.AUX.DEL.roi_names = opt.AUX.roi_names;
    opt.AUX.DEL.roi_path = opt.DATA.ROIS;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   SEED TO VOXEL Correlation (FASTCONN)
if opt.AUX.TODO_StoV  && opt.AUX.seedmeasure
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Seed to Voxel Connectivity (StV) computation\n']);
    
    mkdir(opt.folders.firstlevel,'StV');
    dest_folder = [opt.folders.firstlevel,'/StV'];    
    [~,~] = system('rm *_StV_*');   % in case of previous failure
    
    ndetrend = 3; % if 3drsfc was run without -nodetrend
%     if opt.aCompCor.do
%         tmp = load(c_b(opt.aCompCor.X1D{1,1,1}));   %assuming same matrix for each subject
%         nreg = size(tmp,2) + ndetrend;
%     elseif opt.rp_regression.do == 1 && ~opt.aCompCor.do
%         tmp = load(c_b(opt.DATA.RP_1D{1,1}(1,:)));  %assuming same matrix for each subject
%         nreg = size(tmp,2) + ndetrend;
%     else
%         nreg = ndetrend;
%     end
    
    if opt.aCompCor.do || opt.rp_regression.do
        tmp = load(c_b(opt.X1D{1,1,1}));   %assuming same matrix for each subject
        nreg = size(tmp,2) + ndetrend;
    else
        nreg = ndetrend;
    end
    
    total = sum(opt.subject_number*opt.session_number)*size(opt.DATA.ROIS,1);
    reverseStr  = '';
    count = 0;

    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                if opt.AUX.do_smoothing   %use smoothed for correlation data if possible (but not for seed extraction)
                    LFF_string = '_sLFF.nii';
                else
                    LFF_string = '_LFF.nii';
                end
                base_name = opt.AUX.output_name{g,j,k};
                file_name= [base_name,LFF_string];
                hdr_o = spm_vol(file_name);
                epi = single(spm_read_vols(hdr_o));
                s = size(epi);
                %work for the seed now (always from not smoothed data)
                if opt.AUX.do_smoothing
                    LFF_string = '_LFF.nii';
                    file_name= [base_name,LFF_string];
                    epi_unsm = single(spm_read_vols(spm_vol(file_name)));
                end
                %----------------------------------------------------
                bm = double(uint8(brainmask_selector(g,k,3)));
                for l=1:size(opt.DATA.ROIS,1)   %ciclo rois
                    seed = spm_read_vols(spm_vol(c_b(opt.DATA.ROIS(l,:))));
                    seed = seed.*bm;        %redundant when use_ss == 0
                    ss = nan(s(4),1);
                    for i=1:s(4)
                        if opt.AUX.do_smoothing
                            tmp = epi_unsm(:,:,:,i);
                        else
                            tmp = epi(:,:,:,i);
                        end
                        ss(i,1) = mean(tmp(seed == 1));
                    end
                    roi_str = c_b(opt.AUX.roi_names(l,:));
                    output_name = [base_name,'_StV_',roi_str(1:(end-4))];
                    fastconn(epi,hdr_o(1),ss,opt.tr{g,j},opt.filter_band,nreg,output_name,g,k,j);
                    count = count +1;
                        percentDone = 100 * count / total;
                        msg = sprintf('Percent done: %3.1f (%s)', percentDone,output_name);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
            end
        end
    end
    clear epi epi_unsm hdr_0;
    [~,~] = system(['rm ',dest_folder,'/*']);       % in case you are performing a second time a seed-based measure
    system(['mv *StV_* ',dest_folder]);  
    opt.AUX.TODO_StoV = 0;
    opt.MEASURES.StoV = 1;
    opt.AUX.StoV.roi_names = opt.AUX.roi_names;
    opt.AUX.StoV.roi_path = opt.DATA.ROIS;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   Coherency (COH)
if opt.AUX.TODO_COH  && opt.AUX.seedmeasure
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Coherency (COH) computation\n']);
    
    mkdir(opt.folders.firstlevel,'COH');
    dest_folder = [opt.folders.firstlevel,'/COH'];    
   
    LFF_list = dir('*_LFF.nii');
    
    total = length(LFF_list)*size(opt.DATA.ROIS,1);
    reverseStr  = '';
    count = 0;

    for i=1:length(LFF_list)   %ciclo subjs
        hdr_o = spm_vol(LFF_list(i).name);
        epi = spm_read_vols(hdr_o);
        for l=1:size(opt.DATA.ROIS,1)   %ciclo rois
           count = count +1;
           roi_str =  c_b(opt.AUX.roi_names(l,:));
           output_name = [LFF_list(i).name(1:end-7),'COH_',roi_str(1:(end-4))];
           coherency(epi,hdr_o,c_b(opt.DATA.ROIS(l,:)),output_name);
           
           percentDone = 100 * count / total;
           msg = sprintf('Percent done: %3.1f (%s)', percentDone,output_name);
           fprintf([reverseStr, msg]);
           reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
    [~,~] = system(['rm ',dest_folder,'/*']);       % in case you are performing a second time a seed-based measure
    system(['mv *COH_* ',dest_folder]);  
    opt.AUX.TODO_COH = 0;
    opt.MEASURES.COH = 1;
    opt.AUX.COH.roi_names = opt.AUX.roi_names;
    opt.AUX.COH.roi_path = opt.DATA.ROIS;    
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   ROI to ROI Correlation 
if opt.AUX.TODO_RtoR  
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': ROI to ROI Connectivity (RtoR) computation\n']);
    LFF_list = dir('*_LFF.nii');    %use always unsmoothed data (even if sLFF is present)
    total = length(LFF_list);
    reverseStr  = '';
    count = 0;
    RtoR_data = cell(opt.group_number,opt.session_number);
    RtoR_corr = RtoR_data;
    if opt.RtoR.do_coh
        RtoR_corr_coh = RtoR_corr;
    end
    for g =1:opt.group_number
        for j=1:opt.session_number
            datamatrix = nan(opt.subject_number(g),opt.N{g,j},opt.RtoR.n_rois); % il problema ?? che N potrebbe non essere lo stesso per tutti gli scan % TO FIX now N is a cell
            corrmatrix = single(nan(opt.subject_number(g),3,opt.RtoR.n_rois,opt.RtoR.n_rois));  % 1=r, 2= z, 3 = p;
            if opt.RtoR.do_coh
                corrmatrix_coh = single(corrmatrix);
            end
            if strcmp(opt.prepro_mode,'3dTproject') %add here also the case N is different among subjs
                corrmatrix_zscore = single(corrmatrix);
            end
            for k=1:opt.subject_number(g)
                count = count +1;
                [~,string] = system(['3dROIstats -quiet -mask ',opt.AUX.RtoR.masks_path,' ',opt.AUX.output_name{g,j,k},'_LFF.nii']);
                subj = str2num(string);
                [corrmatrix(k,1,:,:),corrmatrix(k,3,:,:)] = corrcoef(subj);
                corrmatrix(k,2,:,:) = atanh(corrmatrix(k,1,:,:));
                if strcmp(opt.prepro_mode,'3dTproject') %add here also the case N is different among subjs
                    corrmatrix_zscore = corrmatrix;
                    corrmatrix_zscore(k,2,:,:) = atanh(corrmatrix(k,1,:,:))*sqrt(opt.censoring.DOF{g,j}(k,1) -3);
                end
                try
                    datamatrix(k,:,:) = subj;
                catch 
                    datamatrix = 'Error: the number of volumes is different among sessions no datamatrix can be extracted';
                end
                if opt.RtoR.do_coh
                    [corrmatrix_coh(k,1,:,:),corrmatrix_coh(k,2,:,:),corrmatrix_coh(k,3,:,:)] = RtoR_coherency_matrix(subj);
                end
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,opt.AUX.output_name{g,j,k});
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            RtoR_data{g,j} = datamatrix;
            RtoR_corr{g,j} = corrmatrix;
            if opt.RtoR.do_coh
                RtoR_corr_coh{g,j} = corrmatrix_coh;
            end
            if strcmp(opt.prepro_mode,'3dTproject') %add here also the case N is different among subjs
                RtoR_corr_zscore{g,j} = corrmatrix_zscore;
            end
        end
    end
    RtoR_info.subject_number = opt.subject_number;
    RtoR_info.session_number = opt.session_number;
    RtoR_info.group_number = opt.group_number;
    RtoR_info.group_names_prefix = opt.AUX.group_names_prefix;
    RtoR_info.session_names_prefix = opt.AUX.session_names_prefix;
    RtoR_info.RtoR = opt.RtoR;
    save([opt.RtoR.destination_dir,'/RtR_',opt.RtoR.name,'_extraction.mat'],'RtoR_data','RtoR_info','-v7.3');
    save([opt.RtoR.destination_dir,'/RtR_',opt.RtoR.name,'_correlation.mat'],'RtoR_corr','RtoR_info','-v7.3');
    if opt.RtoR.do_coh
        RtoR_corr = RtoR_corr_coh;
        save([opt.RtoR.destination_dir,'/RtR_',opt.RtoR.name,'_correlation_COH.mat'],'RtoR_corr','RtoR_info');
    end
    if strcmp(opt.prepro_mode,'3dTproject') %add here also the case N is different among subjs
        RtoR_corr = RtoR_corr_zscore;
        save([opt.RtoR.destination_dir,'/RtR_',opt.RtoR.name,'_correlation_Zscore.mat'],'RtoR_corr','RtoR_info');
    end
    opt.AUX.TODO_RtoR = 0;
    opt.MEASURES.RtoR = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   gFC
if opt.AUX.TODO_gFC
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Global Functional Connectivity (gFC) computation\n']);
    mkdir(opt.folders.firstlevel,'gFC');
    dest_folder = [opt.folders.firstlevel,'/gFC'];    
    [~,~] = system(['rm *_gFC_Zr.nii *_gFC_Zz.nii *_gFC_P.nii *_gFC_VT.nii *_gFC_HIST.nii']);   % in case of previous failure
    if opt.AUX.special_smoothing_order < 0
        LFF_string = '_sLFF.nii';
    else
        LFF_string = '_LFF.nii';
    end
    total = sum(opt.subject_number*opt.session_number) + sum(opt.AUX.special_smoothing_order(opt.AUX.special_smoothing_order > 0));
    reverseStr  = '';
    %t1 = num2str(opt.AUX.gFC.thr);
    t2 = [num2str(opt.AUX.gFC.Vthr(1)),' ',num2str(opt.AUX.gFC.Vthr(2)),' ',num2str(opt.AUX.gFC.Vthr(3))];
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                count = count +1;
                base_name = opt.AUX.output_name{g,j,k};
                input_name = [base_name,LFF_string];
                out_1_str = ['-Zmean ',base_name,'_gFC_Zr.nii'];
                out_2_str = ['-Pmean ',base_name,'_gFC_P.nii'];
        %         out_3_str = ['-Thresh ',t1,' ',LFF_list(i).name(1:end-7),'gFC_T.nii'];
                out_4_str = ['-VarThresh ',t2,' ',base_name,'_gFC_VT.nii'];
                if opt.AUX.gFC.hist; hist_str = [' -Hist 200 ',base_name,'_gFC_HIST.nii'];else hist_str = '';end
                [status,string] = system(['3dTcorrMap -input ',input_name,' ',brainmask_selector(g,k,2),' -polort -1 ',out_1_str,' ',out_2_str,' ',out_4_str,hist_str]);control(status,string);
                % conversion of gFC_Z to z-fisher
                [status,string] = system(['3dcalc -a ',base_name,'_gFC_Zr.nii -expr ''atanh(a)'' -prefix ',base_name,'_gFC_Zz.nii']);control(status,string);
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
    [~,string] = system(['3dinfo -nv ',base_name,'_gFC_VT.nii']);
    opt.AUX.gFC.num_thr = str2num(string);        
    if opt.AUX.gFC.num_thr > 1
        list_VT = dir('*_gFC_VT.nii');
        for l = 1:length(list_VT)
            [status,string] = system(['3drefit -TR 1 ',list_VT(l).name]);control(status,string);
        end
    end
    if opt.AUX.gFC.hist % this would fix the problem of the extraction but would also remove the bin scale
        list_H = dir('*_gFC_HIST.nii');
        for l = 1:length(list_H)
            [status,string] = system(['3drefit -TR 1 ',list_H(l).name]);control(status,string);
        end
    end
    if opt.AUX.special_smoothing_order == 1
            count = count +1;
            percentDone = 100 * (count) / total;
            msg = sprintf('Percent done: %3.1f (%s)', percentDone,['Smoothing maps...']);  
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            smooth_map('_gFC_Zr.nii',g,k); 
            smooth_map('_gFC_Zz.nii',g,k); 
            smooth_map('_gFC_P.nii',g,k); 
            smooth_map('_gFC_VT.nii',g,k); 
            % I am not smoothing the histogram!
    end
    [status,string] = system(['mv *_gFC_Zr.nii *_gFC_Zz.nii *_gFC_P.nii *_gFC_VT.nii ',dest_folder]); control(status,string);
    if opt.AUX.gFC.hist;  [status,string] = system(['mv *_gFC_HIST.nii ',dest_folder]); control(status,string); end;
    opt.AUX.TODO_gFC = 0;
    opt.MEASURES.gFC= 1;
    opt.AUX.lastcall_gFC_smooth_order = opt.AUX.special_smoothing_order;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   Correlation Entropy (CoE)
if opt.AUX.TODO_CoE
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Correlation Entropy (CoE) computation\n']);
    mkdir(opt.folders.firstlevel,'CoEn');
    dest_folder = [opt.folders.firstlevel,'/CoEn'];
    cd ([opt.folders.firstlevel,'/gFC'])
    hist_list = dir('*_gFC_HIST.nii');
    total = length(hist_list)+2;
    reverseStr  = '';
    [status,string] = system(['3dTstat -sum -prefix _sum.nii ',hist_list(1).name]);control(status,string);
    for i=1:length(hist_list)
        correlation_entropy(hist_list(i).name);
        percentDone = 100 * i / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,hist_list(i).name);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end  
        percentDone = 100 * (i +1) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['z-score conversion...']);  
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    zscore_conversion('_CoEn','_zCoEn');   
        percentDone = 100 * (i +2) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['done!']);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    [status,string] = system(['mv *_CoEn.nii *_zCoEn.nii ',dest_folder]);control(status,string); 
    [status,string] = system('rm _sum.nii');control(status,string); 
    cd (opt.folders.preprocessing)
    opt.AUX.TODO_CoE = 0;
    opt.MEASURES.CoE = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

%%---------------------------------
%   Degree Correlation (DeC)
if opt.AUX.TODO_DeC
    opt.AUX.steps_current = opt.AUX.steps_current +1;
    fprintf(['\nStep ',num2str(opt.AUX.steps_current),'/',num2str(opt.AUX.steps_total),': Degree Correlation (DeC) computation\n']);
    mkdir(opt.folders.firstlevel,'DeC');
    dest_folder = [opt.folders.firstlevel,'/DeC'];
    cd ([opt.folders.firstlevel,'/gFC'])
    [~,~] = system(['rm *_sum.nii']);   % in case of previous failure
    if opt.AUX.lastcall_gFC_smooth_order == -1 && opt.AUX.special_smoothing_order == 1
        warning off backtrace;
        warning('To avoid a double smoothing of DeC, special_smoothing_order will be temporarily set to 0.');
        warning on backtrace;
        sso = 0;
    else
        sso = opt.AUX.special_smoothing_order;
    end
    total = sum(opt.subject_number*opt.session_number)*opt.AUX.DeC.num_thr + sum(opt.AUX.special_smoothing_order(sso > 0));
    reverseStr  = '';
    %bmask = uint8(spm_read_vols(spm_vol(opt.AUX.bmask_path)));
    count = 0;
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)  
                base_name = opt.AUX.output_name{g,j,k};
                for l = 1:opt.AUX.DeC.num_thr
                    [status,string] = system(['3dTstat -sum -prefix ',num2str(l,'%.3d'),'_sum.nii ',base_name,'_gFC_HIST.nii''[',num2str(opt.AUX.DeC.Thrs_afni_indx(l)),'..199]''']);control(status,string); 
                    count = count +1;
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,base_name);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
            zscore_conversion_list('_sum','_zsum',uint8(brainmask_selector(g,k,3)));  %uint8 added but not tested!
            [status,string] = system(['3dTcat -prefix ',base_name,'_DeC.nii *_sum.nii']);control(status,string); 
            [status,string] = system(['3dTcat -prefix ',base_name,'_zDeC.nii *_zsum.nii']);control(status,string); 
            [status,string] = system('rm *_sum.nii *_zsum.nii');control(status,string); 
            end
        end
    end
    if sso == 1
        count = count +1;
        percentDone = 100 * (count) / total;
        msg = sprintf('Percent done: %3.1f (%s)', percentDone,['Smoothing maps...']);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        smooth_map('_DeC.nii',g,k);
        smooth_map('_zDeC.nii',g,k);
        % I am not smoothing the histogram!
    end
    [status,string] = system(['mv *_DeC.nii *_zDeC.nii ',dest_folder]);control(status,string); 
    cd (opt.folders.preprocessing)
    opt.AUX.TODO_DeC = 0;
    opt.MEASURES.DeC = 1;
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    fprintf('\n');
end
%%---------------------------------

war_count = 0;
for k = 1:size(opt.AUX.log,1)
    a = strfind(opt.AUX.log(k,:),'WARNING');
    war_count = war_count + length(a);
end
fprintf('\nO error found in AFNI/bash execution.');
if war_count ~=0
    fprintf('\n%d warning found in AFNI execution.',war_count);
    fprintf('\nCheck AFNI/bash log file in the project folder.\n');
else
    fprintf('\n0 warning found in AFNI esecution');
end

% close log file
fclose(opt.AUX.fid);
%
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');

cd (wd)

fprintf('\nProject Complete. Bye!\n');
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
%%%%%%%%%%%%%%%%%%%%%% END OF MAIN FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exit] = prepare_data(fo,fs,ot)

global opt
exit = 0;

if exist('conn_dir')~=2
    error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''/home/user/MATLAB/ER''))');
end

group_number = length(fo.data);
if group_number == 1
    fs.group_names{1} = '';
elseif (isfield(fs,'group_names') && length(fs.group_names) ~= group_number ) || ~isfield(fs,'group_names') 
    error('Please, define group_names for each data folder.');
end

session_number = length(fs.functional_file);
if session_number > 1 
    if (isfield(fs,'condition_names') && length(fs.condition_names) ~= session_number) || ~isfield(fs,'condition_names')
        error('Please, define condition_names for each scan session.');
    end
end

opt.rp_regression = ot.rp_regression;
if opt.rp_regression.do
    if length(fs.rp_file) ~= session_number
        error('Please, define rp_file for each scan session.');
    end
end

if ~isfield(ot.measure,'DA');ot.measure.DA =0;end  
if ~isfield(ot.measure,'RSFA');ot.measure.RSFA =0;end  
if ~isfield(ot.measure,'REHO');ot.measure.REHO =0;end  
if ~isfield(ot.measure,'gFC');ot.measure.gFC =0;end 
if ~isfield(ot.measure,'CoE');ot.measure.CoE =0;end 
if ~isfield(ot.measure,'DeC');ot.measure.DeC =0;end 
if ~isfield(ot.measure,'SpE');ot.measure.SpE =0;end 
if ~isfield(ot.measure,'StoV');ot.measure.StoV =0;end 
if ~isfield(ot.measure,'DEL');ot.measure.DEL =0;end 
if ~isfield(ot.measure,'COH');ot.measure.COH =0;end 
if ~isfield(ot.measure,'PSD');ot.measure.PSD =0;end 
if ~isfield(ot.measure,'PSDW');ot.measure.PSDW =0;end 
if ~isfield(ot.measure,'PSDM');ot.measure.PSDM =0;end 
if ~isfield(ot.measure,'ALFFW');ot.measure.ALFFW =0;end 
if ~isfield(ot.measure,'ALFFM');ot.measure.ALFFM =0;end 
if ~isfield(ot.measure,'RtoR');ot.measure.RtoR =0;end 

if (ot.measure.StoV || ot.measure.DEL || ot.measure.COH) && ( ~isfield(fo,'rois') || isempty(fo.rois) )
    error('Please, define folders.roi (i.e., the folder containing ROIs for seed-to-voxel analysis).');
end


if ot.measure.RtoR && ot.bm.use_ss
    warning('ROI-to-ROI is not implemented in single subject brain mask.');
    ot.measure.RtoR = 0;
elseif ot.measure.RtoR 
    if ~isfield(ot,'RtoR')
        error('Please, define ER.others.RtoR in order to perform ROI to ROI analysis.');
    elseif ~isfield(ot.RtoR,'roi_folder_or_file') || isempty(ot.RtoR.roi_folder_or_file)
        error('Please, define ER.others.RtoR.roi_folder_or_file. It must contain the path of the folder with the parcellation of the brain.');
    elseif ~isfield(ot.RtoR,'name') || isempty(ot.RtoR.name) 
        error('Please, define ER.others.RtoR.name. It must contain a name to indentify the parcellation chosen.');
    else
        opt.RtoR.name = r_b(ot.RtoR.name);
        opt.RtoR.roi_folder_or_file = ot.RtoR.roi_folder_or_file;
    end
    if ~isfield(ot.RtoR,'COH')
        ot.RtoR.COH = 0;
    end
    opt.RtoR.do_coh = ot.RtoR.COH;
end

wc1 = 0; 
if isfield(fs,'wc1_file') && ~isempty(fs.wc1_file)
    wc1 = 1;
elseif ot.bm_method == 2
    error('The brain mask method 2 requires GM maps');
end
wc23 = 0; 
if (isfield(fs,'wc2_file') && ~isempty(fs.wc2_file)) && (isfield(fs,'wc3_file') && ~isempty(fs.wc3_file)) && wc1 == 1 
    wc23 = 1;
elseif ot.bm_method == 1
    error('The brain mask method 1 requires GM, WM and CSF maps');
end

if opt.rp_regression.do % DA requires rp regression|| ot.measure.DA
    if ~isfield(ot,'rp_type')
        error('The field ER.others.rp_type must be specified ({''FSL'',''SPM'',''AFNI'',''HCP''}).');
    end 
    opt.AUX.rp_tye = ot.rp_type;
    rp_define_ordering_and_unit;
end
% check for denoising analysis (requires rp and wc1)
if ot.measure.DA && ~opt.rp_regression.do
    warning off backtrace
    warning('Denoising analysis requires rp parameters to be performed. Turning off this step.');
    ot.measure.DA = 0;
    warning on backtrace
elseif ot.measure.DA && opt.rp_regression.do
    opt.DA.force_overwrite = ot.DA.force_overwrite;
    if ~wc1
        error('Denoising analysis requires GM probability maps.');
    end
end

% ----- set some default values: ----------
if ~isfield(ot,'tredrsfc_extraoption');  ot.tredrsfc_extraoption = '';end
if ~isfield(ot,'psc_save');         ot.psc_save = 0;end
if ~isfield(ot,'psc');              ot.psc = 0;end
if ~isfield(ot,'volume_selector');  ot.volume_selector = [];end
%if ~isfield(ot,'N');                ot.N = [];end
if ~isfield(ot.bm,'group');         ot.bm.group = [];end
if ~isfield(ot.bm,'modality');      ot.bm.modality = 'hard';end
if ~isfield(ot.bm,'ss_masks_thr');  ot.bm.modality = [];end
if ~isfield(ot.bm,'use_ss');        ot.bm.use_ss = 0;end
if ~isfield(ot.bm,'doautomask');    ot.bm.doautomask  = 0;end
if ~isfield(ot.bm,'automask_session_selector');    ot.bm.automask_session_selector = [];end
if ~isfield(ot,'aCompCor');         ot.aCompCor.do = 0; end
if ot.aCompCor.do == 1 && ~isfield(ot.aCompCor,'ROI')
    ER.others.aCompCor.ROI(1).tissue = 2;           % 1=GM, 2=WM, 3=CSF
    ER.others.aCompCor.ROI(1).dime = 5;             % Number of PCA components. 0 for straight average.
    ER.others.aCompCor.ROI(1).deri = 0;             % 0/n>1. Regress also the first n derivates. O: no derivatives 
    ER.others.aCompCor.ROI(1).erode = 2;            % Erode by 1 voxel to avoid tissue contamination
    ER.others.aCompCor.ROI(1).thr = 0.5;            % Threshold of probability maps
    ER.others.aCompCor.ROI(2).tissue = 3;           % 1=GM, 2=WM, 3=CSF
    ER.others.aCompCor.ROI(2).dime = 5;             % Number of PCA components. 0 for straight average.
    ER.others.aCompCor.ROI(2).deri = 0;             % 0/n>1. Regress also the first n derivates. O: no derivatives
    ER.others.aCompCor.ROI(2).erode = 1;            % Erode by 1 voxel to avoid tissue contamination
    ER.others.aCompCor.ROI(2).thr = 0.5;            % Threshold of probability maps
end
if ot.aCompCor.do == 1 && ~isfield(ot.aCompCor,'rpOrtogonalize')
    ER.others.aCompCor.rpOrtogonalize = 1;
end
if ot.aCompCor.do == 1 && ~isfield(ot.aCompCor,'asCONN')
    ER.others.aCompCor.asCONN = 0;
end
if ot.aCompCor.do == 1 && ~isfield(ot.aCompCor,'TvarianceNormalise')
    ER.others.aCompCor.TvarianceNormalise = 0;
end
if ~isfield(ot,'censoring');                ot.censoring.do = 0;end
if ~isfield(ot,'FWHM');                     ot.FWHM = 0; end
if ~isfield(ot, 'special_smoothing_order'); ot.special_smoothing_order = 1; end
if ~isfield(ot, 'special_smoothing_FWHM');  ot.special_smoothing_FWHM  = []; end
if ~isfield(ot,'BlurToFWHM_extraoption');   ot.BlurToFWHM_extraoption = ''; end;
opt.censoring = ot.censoring;

% CENSORING modality
if opt.censoring.do
    opt.prepro_mode = '3dTproject';  
    %TODO,put here warning about alff ....
    
    %
    % checks and defaults
    if ~isfield(opt.censoring,'mode'); opt.censoring.mode = 'KILL';end
    allowed_str = {'KILL','ZERO','NTRP'};
    if sum(strcmp(allowed_str,opt.censoring.mode)) == 0
        error('Provide a correct censor modality, ER.others.censoring.mode = ''KILL'' or ''ZERO'' or ''NTRP.''');
    end
    if ~isfield(opt.censoring,'extraoption'); opt.censoring.extraoption = '';end
    if ~isempty(opt.censoring.extraoption);
        opt.censoring.extraoption = [' ',opt.censoring.extraoption];%added here a blank space
    end
    if ~isfield(opt.censoring,'value'); opt.censoring.value = 0.3;end
else
    opt.prepro_mode = '3dRSFC';   
end

% COMPCOR checks
if ot.aCompCor.do == 1
    for l = 1:length(ot.aCompCor.ROI)
       if  (ot.aCompCor.ROI(l).tissue ~= 1 && ot.aCompCor.ROI(l).tissue ~= 2 && ot.aCompCor.ROI(l).tissue ~= 3 && ot.aCompCor.ROI(l).tissue ~= 4)
           error('ER.others.aCompCor.ROI(1).tissue must be 1, 2, 3 or 4 (GM, WM, CSF or GSR).');
       end
       switch ot.aCompCor.ROI(l).tissue
           case 1
               if ~isfield(fs,'wc1_file') || isempty(fs.wc1_file)
                   error('aCompCor error. Grey matter probability maps required.');
               end
           case 2
               if ~isfield(fs,'wc2_file') || isempty(fs.wc2_file)
                   error('aCompCor error. White matter probability maps required.');
               end     
           case 3
               if ~isfield(fs,'wc3_file') || isempty(fs.wc3_file)
                   error('aCompCor error. CSF probability maps required.');
               end   
           case 4 % makes the others useless...it should be imporoved
               if (~isfield(fs,'wc1_file') || isempty(fs.wc1_file)) || (~isfield(fs,'wc2_file') || isempty(fs.wc2_file)) || (~isfield(fs,'wc3_file') || isempty(fs.wc3_file))
                   error('GSR requirs GM, WM and CSF probability maps.');
               end
       end
    end
end

%smoothing
opt.FWHM = ot.FWHM;
opt.AUX.special_smoothing_order = ot.special_smoothing_order;
if ~isempty(opt.FWHM) && opt.FWHM > 0
    if ~isfield(ot,'smoothing_mode');  ot.smoothing_mode = '3dBlurInMask';  end;%set default
    switch ot.smoothing_mode
        case {'3dBlurInMask'}
        case {'3dBlurToFWHM'}
            opt.AUX.BlurToFWHM_extraoption = ot.BlurToFWHM_extraoption;
            if ~isempty(opt.AUX.BlurToFWHM_extraoption)
                opt.AUX.BlurToFWHM_extraoption = [' ',opt.AUX.BlurToFWHM_extraoption];%added here a blank space
            end
        otherwise
            error('Unrecognized smoothing mode! Choose between: ''3dBlurInMask'' or ''3dBlurToFWHM''');
    end
    opt.smoothing_mode = ot.smoothing_mode;
    opt.AUX.do_smoothing = 1;
else
    opt.FWHM = 0;
    opt.smoothing_mode = 'none';
    opt.AUX.do_smoothing = 0;
    if opt.AUX.special_smoothing_order == -1
        warning off backtrace
        warning('Special_smoothing_order = -1 is not possible when preprocessing does not include smoothing. Changin value to 0.');
        warning on backtrace
        opt.AUX.special_smoothing_order = 0;
    end
end
% special smoothing FWHM (it refers to reho dc...)
if opt.AUX.special_smoothing_order == 1
    if isempty(ot.special_smoothing_FWHM)
        error ('No specified smoothing in ER.others.special_smoothing_FWHM.');
    end
    if ot.special_smoothing_FWHM == 0
        opt.AUX.special_smoothing_order = 0;
    end
end
opt.AUX.special_smoothing_FWHM = ot.special_smoothing_FWHM;
% make sure -blur is not among options of 3drsfc or 3dTproject (same
% variable)
if ~isempty(strfind(ot.tredrsfc_extraoption,'-blur'))
    error ('Smoothing via ER.others.tredrsfc_extraoption is not supported anymore. Use new ER.others.smoothing function.');
end
    

for g = 1:group_number
    
    cd(fo.data{g})
    if wc1
        WC1{g} = conn_dir(fs.wc1_file); 
        NWC1(g) = size(WC1{g},1);
    end
    if wc23
        WC2{g} = conn_dir(fs.wc2_file); 
        WC3{g} = conn_dir(fs.wc3_file);
        NWC2(g) = size(WC2{g},1);
        NWC3(g) = size(WC3{g},1);
    end
    for j=1:session_number
       FUNCTIONALS{g,j} = conn_dir(fs.functional_file{j});
       if opt.rp_regression.do || opt.censoring.do
           RP{g,j} = conn_dir(fs.rp_file{j});
           NRP(g,j) = size(RP{g,j},1);
       end
    end
    MASKS{g} = conn_dir(ot.bm.ss_masks);

    for j=1:session_number
        NFILES(g,j)=size(FUNCTIONALS{g,j},1);   % Al momento ogni soggetto 
        % ha lo stesso numero di sessioni. Non ?? contemplata la possibilit??
        % che ne manchi qualcuna. In futuro questa variabile NFILES, sar??
        % necessaria a questo scopo. Problema: esempio, nel caso di 10
        % soggetti nel gruppo 1, con 2 sessioni, la prima eseguita da tutti
        % mentre la seconda solo da 4 soggetti, bisogna assegnare la
        % seconda ai giusti soggetti.
    end
    NSUBJECTS(g) = NFILES(g,1);
    NMASKS(g) = size(MASKS{g},1);

    % Subject names
    if isfield(fs,'subject_names') && ~isempty(fs.subject_names)
        list_subj = dir(fs.subject_names);
        if strcmp(fs.subject_names,'*')
            list_subj(1:2) = [];
        end
        for i=1:length(list_subj)
             if list_subj(i).isdir
                 subject_names{g,i} = list_subj(i).name;
             end
        end
        checker = ~cellfun(@isempty,subject_names);
        if sum(checker(g,:)) ~= NSUBJECTS(g) 
            error('Subject list name and functional images number mismatch in group %s.',fs.group_names{g});
        end   
    else 
        for i = 1:NSUBJECTS(g)
            if i < 10
                subject_names{g,i} = ['subj_00',num2str(i)];
            elseif i < 100 && i >= 10
                subject_names{g,i} = ['subj_0',num2str(i)];
            elseif i < 1000 && i >= 100
                subject_names{g,i} = ['subj_',num2str(i)];
            end
        end
    end
    %----------------

    % Output names
    if group_number == 1      %    group prefix
        group_prefix = '';
        tmp_group_prefix{g} = group_prefix;
    else
        group_prefix = [fs.group_names{g},'_']; %add _
        tmp_group_prefix{g} = group_prefix;
    end
    if ~isfield(fs,'condition_assignment') || (isfield(fs,'condition_assignment') && isempty(fs.condition_assignment{g}))
        for j = 1:session_number
            for i =1:NSUBJECTS(g)
                if session_number == 1
                    output_name{g,j,i} = [group_prefix,subject_names{g,i}];
                else
                    output_name{g,j,i} = [group_prefix,subject_names{g,i},'_',fs.condition_names{j}];
                end
            end
        end
    else 
        cas = fs.condition_assignment{g};
        size_cond = size(cas);
        if size_cond(1) ~= NSUBJECTS(g)
            error('Badly defined condition_assignment: not enough subjects.');
        end
        if size_cond(2) ~= session_number
            error('Badly defined condition_assignment: not enough conditions.');
        end
        for j = 1:session_number
            tmp = [];
            for i =1:NSUBJECTS(g) 
                 output_name{g,j,i} = [group_prefix,subject_names{g,i},'_',fs.condition_names{cas(i,j)}];
            end
        end
    end
    %----------------

    if opt.rp_regression.do || opt.censoring.do
        if  any(NRP(g,:) ~= NSUBJECTS(g)) && ~isempty(fs.rp_file) 
            error('Anatomical and rp (motion covariate) number mismatch in group %s.',fs.group_names{g});
        end
    end

    if (wc1 && wc23) && ( (NWC1(g) ~= NSUBJECTS(g)) || (NWC1(g) ~= NWC2(g) || NWC1(g) ~= NWC3(g)) )
        error('Segmented images number mismatch in group %s.',fs.group_names{g});
    end
    
    if (wc1) && (NWC1(g) ~= NSUBJECTS(g))
        error('Segmented images number mismatch in group %s.',fs.group_names{g});
    end

    if (NMASKS(g) ~= NSUBJECTS(g)) && ~isempty(ot.bm.ss_masks)
        error('Mask files number mismatch in group %s.',fs.group_names{g});
    end

    fprintf('\nGROUP %d: %s\n',g,fs.group_names{g});
        
    for i=1:NSUBJECTS(g)
         fprintf('\nSubject: %s',subject_names{g,i});
         if wc1
             %Segmented GM
             aux = WC1{g}(i,1:end);
             fprintf('\n%s',aux);
         end
         if wc23
             %Segmented WM
             aux = WC2{g}(i,1:end);
             fprintf('\n%s',aux);       
             %Segmented CSF
             aux = WC3{g}(i,1:end);
             fprintf('\n%s',aux);    
         end
         %RS
         for j = 1:session_number
            aux = FUNCTIONALS{g,j}(i,1:end);
            aux2 = output_name{g,j,i};
            fprintf('\n%s --> %s',aux,aux2);    
         end
         %RP covariates
         if opt.rp_regression.do || opt.censoring.do   
             for j = 1:session_number
                aux = RP{g,j}(i,1:end);
                fprintf('\n%s',aux);
             end
         end
    end
    fprintf('\n');
    
    % -----------------------------------------
warning off backtrace
if (ot.bm_method == 3 || ot.bm_method == 5 || ot.bm_method == 7) && ot.bm.use_ss == 1
    warning('The use of subject-specific brainmask (ER.others.bm.use_ss) is not available in bm_method 3, 5 and 7. Group-brainmask is turned on.');
    ot.bm.use_ss = 0;
end

if ~isempty(ot.bm.ss_masks_thr) && (ot.bm.ss_masks_thr > 1 || ot.bm.ss_masks_thr < 0)
    error('The threshold for brain mask method 4 (ER.others.bm.ss_masks_thr) must be a floating value comprised between 0 and 1');
end

if ot.measure.SpE && ~ot.measure.PSD
    warning('Spectral entropy calculation requires the Power Spectral Density. Changing measure.PSD from 0 to 1.');
    ot.measure.PSD = 1;
end
if ot.measure.CoE && ~ot.measure.gFC
    warning('Correlation entropy calculation requires global Functionl Connectivity (gFC). Changing measure.gFC from 0 to 1.');
    ot.measure.gFC = 1;
end
if ot.measure.DeC && ~ot.measure.gFC
    warning('Degree Connectivity requires global Functionl Connectivity (gFC). Changing measure.gFC from 0 to 1.');
    ot.measure.gFC = 1;
end
if ot.measure.ALFFW && ~ot.measure.PSDW
    warning('ALFFW calculation requires the Power Spectral Density with Welch''s method. Changing measure.PSDW from 0 to 1.');
    ot.measure.PSDW = 1;
end
if ot.measure.ALFFM && ~ot.measure.PSDM
    warning('ALFFM calculation requires the Power Spectral Density with multitaper method. Changing measure.PSDM from 0 to 1.');
    ot.measure.PSDM = 1;
end
warning on backtrace
    
    % rp regression
    if opt.rp_regression.do || opt.censoring.do   
        opt.DATA.RP = RP;
    end
    % OLD
%     if (opt.rp_regression.do || opt.censoring.do) && ~ot.aCompCor.do
%          for j = 1:session_number 
%             tmp2 = [];
%             for i = 1:NSUBJECTS(g)
%                 tmp = c_b(RP{g,j}(i,:));
%                 if 0 %test the effect of detrending
%                     rp = load(tmp);
%                     rp = detrend(rp);
%                     save('_temp.txt','rp','-ascii');
%                     [~,~] = system(['1dcat _temp.txt > ',tmp(1:end-4),'.1D']);
%                     system('rm _temp.txt');
%                 else
%                     [~,~] = system(['1dcat ',tmp,' > ',tmp(1:end-4),'.1D']);
%                 end
%                 tmp2 = strvcat(tmp2,[tmp(1:end-4),'.1D']);
%             end
%             RP_1D{g,j} = tmp2;
%          end
%          opt.DATA.RP_1D = RP_1D;
%     end
%     if (opt.rp_regression.do || opt.censoring.do)
%         opt.DATA.RP = RP;
%     end
    %----------------

end

warning off backtrace
if  ot.bm.use_ss == 1
    warning(sprintf('\nSubject-specific brainmask has been chosen. Be aware:\n1) Voxel-wise analysis is not possible in this modality.\n2) ROI-based analysis might be baised.\n3) Comparison of standardized measures might be biased.'));
end
warning on backtrace
    
if ot.The_IKnowWhatIamDoing_Condition 
    fprintf('\nIt seems that you are sure that data are correctly ordered.');
    fprintf('\nWell, it''s your responsability. Now, let''s start the analysis...\n');
    pause(0.5);
else
    fprintf('\nPlease, carefully check if data are correctly ordered!');
    [opzione] = scelta_yn ('Ready to start (type ''n'' to exit)','y');
    if opzione == 'n'
        exit = 1;
        return
    end
end


% ----Control and AUXiliary variables----
% make sure there are no spaces in project name
fo.project = r_b(fo.project);

opt.folders.data = fo.data;
opt.folders.results = [fo.results,'/ER_',fo.project];
opt.folders.preprocessing = [fo.results,'/ER_',fo.project,'/preprocessing'];
opt.folders.firstlevel = [fo.results,'/ER_',fo.project,'/firstlevel'];
opt.folders.secondlevel = [fo.results,'/ER_',fo.project,'/secondlevel'];
opt.folders.seed = [opt.folders.results,'/seed'];
opt.folders.prject_name = fo.project;


if exist([fo.results,'/ER_',fo.project],'dir') 
    dirtoremove = [fo.results,'/ER_',fo.project];
    system(['rm -r ',dirtoremove]);
end

mkdir(fo.results,['/ER_',fo.project]);
mkdir(opt.folders.results,'preprocessing');
mkdir(opt.folders.results,'firstlevel');
mkdir(opt.folders.results,'secondlevel');
mkdir(opt.folders.results,'seed');

% start log
opt.AUX.fid = fopen([opt.folders.results,'/ER_',opt.folders.prject_name,'.log'],'w');
%

if wc1; opt.DATA.WC1 = WC1; end
if wc23; opt.DATA.WC2 = WC2; opt.DATA.WC3 = WC3; end

opt.DATA.FUNCTIONALS = FUNCTIONALS;

opt.aCompCor = ot.aCompCor;
if opt.aCompCor.do 
    opt.aCompCor.masknumb = length(opt.aCompCor.ROI);
end

opt.DATA.MASKS = MASKS;

opt.MEASURES.DA = 0;
if opt.censoring.do 
    opt.MEASURES.ALFF = 0;
else
    opt.MEASURES.ALFF = 1;  %automatically produced by 3drsfc
end
opt.MEASURES.RSFA = 0;
opt.MEASURES.REHO = 0;
opt.MEASURES.PSD = 0;
opt.MEASURES.gFC = 0;
opt.MEASURES.CoE = 0;
opt.MEASURES.DeC = 0;
opt.MEASURES.StoV = 0;
opt.MEASURES.SpE = 0;
opt.MEASURES.DEL = 0;
opt.MEASURES.COH = 0;
opt.MEASURES.PSDW = 0;
opt.MEASURES.ALFFW = 0;
opt.MEASURES.PSDM = 0;
opt.MEASURES.ALFFM = 0;
opt.MEASURES.RtoR = 0;

opt.AUX.TODO_DA = ot.measure.DA;
opt.AUX.TODO_RSFA = ot.measure.RSFA;
opt.AUX.TODO_REHO = ot.measure.REHO;
opt.AUX.TODO_PSD = ot.measure.PSD;
opt.AUX.TODO_gFC = ot.measure.gFC;
opt.AUX.TODO_CoE = ot.measure.CoE;
opt.AUX.TODO_DeC = ot.measure.DeC;
opt.AUX.TODO_StoV = ot.measure.StoV;
opt.AUX.TODO_SpE = ot.measure.SpE;
opt.AUX.TODO_DEL = ot.measure.DEL;
opt.AUX.TODO_COH = ot.measure.COH;
opt.AUX.TODO_PSDW = ot.measure.PSDW;
opt.AUX.TODO_ALFFW = ot.measure.ALFFW;
opt.AUX.TODO_PSDM = ot.measure.PSDM;
opt.AUX.TODO_ALFFM = ot.measure.ALFFM;
opt.AUX.TODO_RtoR = ot.measure.RtoR;

opt.AUX.seedmeasure = logical (opt.AUX.TODO_StoV || opt.AUX.TODO_DEL || opt.AUX.TODO_COH);  % add here seed based measure
if opt.AUX.seedmeasure 
    if isfield(fo,'rois') && ~isempty(fo.rois)
        opt.folders.rois = fo.rois;
    else
        error('You need to specify a ROIs folder (ER.folders.rois).');
    end
end
opt.AUX.steps_total = (2 + opt.censoring.do + opt.rp_regression.do + opt.aCompCor.do + opt.AUX.do_smoothing + opt.AUX.TODO_DA + opt.AUX.seedmeasure + opt.AUX.TODO_REHO + opt.AUX.TODO_PSD +opt.AUX.TODO_gFC + opt.AUX.TODO_CoE + opt.AUX.TODO_StoV +opt.AUX.TODO_SpE +opt.AUX.TODO_DEL + opt.AUX.TODO_COH + opt.AUX.TODO_PSDW + opt.AUX.TODO_ALFFW + opt.AUX.TODO_PSDM + opt.AUX.TODO_ALFFM +opt.AUX.TODO_DeC+ 2*opt.AUX.TODO_RtoR );
opt.AUX.steps_current = 0;

opt.subject_number = NSUBJECTS;
opt.session_number = session_number;
opt.group_number = group_number;

% NEW- TR variability across samples

if ~iscell(ot.tr)
    for g =1:opt.group_number
        for j=1:opt.session_number
            opt.tr{g,j} = ot.tr;
        end
    end
else
    % check consistency
    if size(ot.tr,1) ~= opt.group_number; error('You must specify a TR for each group of subjects (ER.others.tr{g,s}). If TR is the same for all datasets do not use a cell, just a number.');end
    if size(ot.tr,2) ~= opt.session_number; error('You must specify a TR for each session (ER.others.tr{g,s}). If TR is the same for all datasets do not use a cell, just a number.');end
    opt.tr = ot.tr;
end
    
%opt.tr = ot.tr;
opt.AUX.fs = 1/opt.tr{1,1};   % MUST BE FIXED SINCE NOW TR is a cell. For the moment I just use the first value.
opt.filter_band = ot.filter_band;
opt.psc = ot.psc;
opt.AUX.psc_save = ot.psc_save;
opt.volume_selector = ot.volume_selector;

% Volume selector
if isempty(opt.volume_selector)
    for g =1:opt.group_number
        for j=1:opt.session_number
            opt.AUX.v_s_rp_str{g,j} = '';
            opt.AUX.v_s_str{g,j} = '';
            opt.AUX.v_s_matlab{g,j} = [];
        end
    end
    opt.AUX.v_s_do = 0;
elseif ~iscell(opt.volume_selector)   % ok it's not empty but it's not a cell (old case)
        if numel(opt.volume_selector) == 2 && (opt.volume_selector(1) < opt.volume_selector(2))
            for g =1:opt.group_number
                for j=1:opt.session_number
                    opt.AUX.v_s_rp_str{g,j} = ['''{',num2str(opt.volume_selector(1)),'..',num2str(opt.volume_selector(2)),'}'''];
                    opt.AUX.v_s_str{g,j} = ['''[',num2str(opt.volume_selector(1)),'..',num2str(opt.volume_selector(2)),']'''];
                    opt.AUX.v_s_matlab{g,j} = [(opt.volume_selector(1) +1),(opt.volume_selector(2) +1)];
                    opt.N{g,j} = opt.volume_selector(2) - opt.volume_selector(1) +1;    
                end
            end
            opt.AUX.v_s_do = 1;
        else
            error('Wrongly defined volume selector. It must be integer (starting and ending volume index).');
        end
else %new case. It's a cell
    if size(opt.volume_selector,1) ~= opt.group_number; error('You must specify a volume_selector for each group of subjects (ER.others.volume_selector{g,s}). If volume_selector is the same for all datasets do not use a cell, just two element vector.');end
    if size(opt.volume_selector,2) ~= opt.session_number; ('You must specify a volume_selector for each session (ER.others.volume_selector{g,s}). If volume_selector is the same for all datasets do not use a cell, just two element vector.');end
    for g =1:opt.group_number 
        for j=1:opt.session_number 
            if numel(opt.volume_selector{g,j}) == 2 && (opt.volume_selector{g,j}(1) < opt.volume_selector{g,j}(2))
                opt.AUX.v_s_rp_str{g,j} = ['''{',num2str(opt.volume_selector{g,j}(1)),'..',num2str(opt.volume_selector{g,j}(2)),'}'''];
                opt.AUX.v_s_str{g,j} = ['''[',num2str(opt.volume_selector{g,j}(1)),'..',num2str(opt.volume_selector{g,j}(2)),']'''];
                opt.AUX.v_s_matlab{g,j} = [(opt.volume_selector{g,j}(1) +1),(opt.volume_selector{g,j}(2) +1)];
                opt.N{g,j} = opt.volume_selector{g,j}(2) - opt.volume_selector{g,j}(1) +1;   
            else
                error('Wrongly defined volume selector. It must be integer (starting and ending volume index).');
            end
        end;
    end;
    opt.AUX.v_s_do = 1;
end
    
%----------------

% Calculating the number of volume (only in the case volume selector is
% off)
if ~isfield(opt,'N') 
        for g =1:opt.group_number 
            for j=1:opt.session_number 
                [~,N] = system(['3dinfo -nv ',c_b(opt.DATA.FUNCTIONALS{g,j}(1,:))]);        %attention: assuming that subjs form the same group and session have the same N
                opt.N{g,j} = str2num(N);   
            end;
        end;
end

% 3drsfc extra option
if isempty(ot.tredrsfc_extraoption)
    opt.tredrsfc_extraoption = ot.tredrsfc_extraoption;
else
    opt.tredrsfc_extraoption = [' ',ot.tredrsfc_extraoption];%added here a blank space
end
%----------------

% Brain mask 
opt.bm_method = ot.bm_method;
opt.bm = ot.bm;
%----------------

% Rois check
if opt.AUX.seedmeasure
    if ~exist(opt.folders.rois,'dir') 
        error('Folders.rois doesn''t exist.')
    end
    [R_NII,list_R_NII] = conn_dir([opt.folders.rois,'/*.nii']);
    [R_HDR,list_R_HDR] = conn_dir([opt.folders.rois,'/*.hdr']);
    opt.DATA.ROIS = strvcat(R_NII,R_HDR);
    list_ROIS = cat(2,list_R_NII,list_R_HDR);
    if isempty(opt.DATA.ROIS)
        warning('Roi folder is empty. Seed-based measures will NOT run.');
        opt.AUX.seedmeasure = 0;
    else
        %--------------check roi size and binarity-------------------------
        [~,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
        data_info = str2num(data_info);
        for z = 1:length(list_ROIS)
            [~,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.DATA.ROIS(z,:)]);
            roi_info = fix_float_errors(roi_info);
            roi_info = str2num(roi_info);
            if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' has not the same size of functional data.']);
            elseif roi_info(4)~=1
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' has more than one sub-brick.']);
            elseif roi_info(5)~=0 || roi_info(6)~=1
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' is not binary.']);
            end
        end
        %------------------------------------------------------------------
        roi_names = [];
        for k = 1:length(list_ROIS)
            roi_names = strvcat(roi_names,list_ROIS(k).name);
        end
        opt.AUX.roi_names = roi_names;
    end
end

% Rois check for RtoR
if opt.AUX.TODO_RtoR
    switch exist(opt.RtoR.roi_folder_or_file)
        case 0  % no valid input
            error('Roi folder/file for ROI to ROI connectivity doesn''t exist.')
        case 7  % it's a directory
            opt.RtoR.modalitiy = 'roi_list';
            [R_NII,list_R_NII] = conn_dir([opt.RtoR.roi_folder_or_file,'/*.nii']);
            [R_HDR,list_R_HDR] = conn_dir([opt.RtoR.roi_folder_or_file,'/*.hdr']);
            opt.RtoR.ROIS = strvcat(R_NII,R_HDR);
            list_ROIS = cat(2,list_R_NII,list_R_HDR);
            if isempty(opt.RtoR.ROIS)
                warning off backtrace
                warning('RtoR: Roi folder for ROI to ROI connectivity is empty. ROI to ROI connectivity will NOT run.');
                opt.AUX.TODO_RtoR = 0;
                warning on backtrace
            else
                %--------------check roi size and binarity-------------------------
                [~,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
                data_info = str2num(data_info);
                for z = 1:length(list_ROIS)
                    [~,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.RtoR.ROIS(z,:)]);
                    roi_info = fix_float_errors(roi_info);
                    roi_info = str2num(roi_info);
                    if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' has not the same size of functional data.']);
                    elseif roi_info(4)~=1
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' has more than one sub-brick.']);
                    elseif roi_info(5)~=0 || roi_info(6)~=1
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' is not binary.']);
                    end
                end
                %------------------------------------------------------------------
                rois_string = [];
                for l = 1:size(opt.RtoR.ROIS,1)
                    tmp_str = c_b(opt.RtoR.ROIS(l,:));
                    rois_string = [rois_string,' ',tmp_str];
                end
                [~,~] = system(['3dMean -prefix _temp_rtor.nii -sum',rois_string]);
                rtor_sum = spm_read_vols(spm_vol('_temp_rtor.nii'));
                n_overlapping = length(rtor_sum(rtor_sum > 1));
                system('rm _temp_rtor.nii');
                opt.RtoR.n_ovelapping = 0;
                if n_overlapping > 0
                    warning off backtrace
                    warning('RtoR: Overllaping ROIs found (%.0d voxels). Overllaping voxels will be removed unless they don''t fall in the group brain mask.',n_overlapping);
                    warning on backtrace
                    opt.RtoR.n_ovelapping = n_overlapping;
                    rtor_sum(rtor_sum == 1) = 0;
                    rtor_sum(rtor_sum > 1) = 1;
                    opt.RtoR.ovelapping = rtor_sum;
                else
                    clear rtor_sum;
                end
                roi_names = [];
                for k = 1:length(list_ROIS)
                    roi_names = strvcat(roi_names,list_ROIS(k).name(1:end-4));
                end
                opt.RtoR.roi_names = roi_names;
            end
        case 2  % it's a file mask
            opt.RtoR.modalitiy = 'mask_file';
            %--------------check roi size ------------------------------------
            [~,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
            data_info = str2num(data_info);
            [~,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.RtoR.roi_folder_or_file]);
            roi_info = str2num(roi_info);
            if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                error('RtoR: The mask dataset (or parecellation atlas) has not the same size of functional data.');
            end
            %------------------------------------------------------------------
            mask_dataset = spm_read_vols(spm_vol(opt.RtoR.roi_folder_or_file));
            max_roi = max(mask_dataset(:));
            mask_unique = unique(mask_dataset(:));
            checker = max_roi - length(mask_unique) +1;
            if checker ~= 0
                error('RtoR: The mask dataset has missing roi indexes. Create a mask dataset (or parcellation atlas) without "holes".');
            end
            opt.RtoR.mf_max_roi = max_roi;
            if exist([opt.RtoR.roi_folder_or_file(1:end-3),'txt'],'file') == 2
                fid = fopen([opt.RtoR.roi_folder_or_file(1:end-3),'txt']);
                line = fgets(fid);
                roi_names = [];
                while ischar(line)    
                    roi_names = strvcat(roi_names,line); 
                    line = fgets(fid);
                end
                opt.RtoR.roi_names = roi_names;
                if size(roi_names,1) ~= max_roi
                    error('RtoR: The .txt file with labels for the mask dataset has a wrong number of lines.');
                end
                opt.RtoR.mf_txt_exist = 1;
                opt.RtoR.mf_txt_file = [opt.RtoR.roi_folder_or_file(1:end-3),'txt'];
            else
                warning off backtrace;
                warning('RtoR: No labels file (txt) found for ROI to ROI connectivity. ROIs will be named as sequential numbers');
                warning on backtrace;
                indexs = 1:1:opt.RtoR.mf_max_roi;
                indexs = indexs';
                opt.RtoR.roi_names = num2str(indexs);
                opt.RtoR.mf_txt_exist = 0;
            end
    end
end

if opt.AUX.TODO_COH || (opt.AUX.TODO_RtoR && opt.RtoR.do_coh)
    coherency_initialization(ot);
end
if opt.AUX.TODO_PSDW; psdw_initialization(ot); end
if opt.AUX.TODO_ALFFW; ALFF_welch_initialization; end
if opt.AUX.TODO_PSDM; psdm_initialization(ot); end
if opt.AUX.TODO_ALFFM; ALFF_multitaper_initialization; end
if opt.AUX.TODO_PSD; psd_initialization; end
if opt.AUX.TODO_gFC; gFC_initialization(ot); end
if opt.AUX.TODO_DeC; DeC_initialization(ot); end

opt.AUX.output_name = output_name;
opt.AUX.subject_names = subject_names;
opt.AUX.group_names_prefix = tmp_group_prefix; 
if opt.session_number == 1
    opt.AUX.session_names_prefix{1} = '';
else
    for l = 1:opt.session_number 
        opt.AUX.session_names_prefix{l} = ['_',fs.condition_names{l}];
    end
end

opt.AUX.log = '';

opt.AUX.preprocessing_done = 0;
opt.AUX.bmi_done = 0;

opt.date_creation = date;

save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');

%--------------------------

return
end

function redo_mode_prepare(fo,ot)
global opt
opt.AUX.fid = fopen([opt.folders.results,'/ER_',opt.folders.prject_name,'.log'],'a');


% check for denoising analysis (requires rp and wc1)
if ot.measure.DA && ~opt.rp_regression.do
    warning off backtrace
    warning('Denoising analysis requires rp parameters to be performed. Turning off this step.');
    ot.measure.DA = 0;
    warning on backtrace
elseif ot.measure.DA && opt.rp_regression.do
    opt.DA.force_overwrite = ot.DA.force_overwrite;
    if ~isfield(opt.DATA,'WC1')  
    %missing GM maps. Let's try to load them
    if isfield(fs,'wc1_file') && ~isempty(fs.wc1_file)
        tmpdir = pwd;  
        for g = 1:opt.group_number
            cd(fo.data{g})
            WC1{g} = conn_dir(fs.wc1_file); 
            NWC1(g) = size(WC1{g},1);
            if NWC1(g) ~= opt.subject_number(g)
                error('Segmented images number mismatch in group %s.',fs.group_names{g});
            end
            %print
            for i=1:NSUBJECTS(g)
                 fprintf('\nSubject: %s',subject_names{g}(i,:));
                 %Segmented GM
                 aux = WC1{g}(i,1:end);
                 fprintf('\n%s',aux);
            end
            fprintf('\n');
        end
        opt.DATA.WC1 = WC1;
        cd(tmpdir);
    else
        error('Denoising analysis requires GM probability maps. You can add them using the field ER.file_str.wc1_file.');
    end
    end
end

opt.AUX.TODO_DA = ot.measure.DA;
opt.AUX.TODO_REHO = ot.measure.REHO;
opt.AUX.TODO_PSD = ot.measure.PSD;
opt.AUX.TODO_gFC = ot.measure.gFC;
opt.AUX.TODO_CoE = ot.measure.CoE;
opt.AUX.TODO_DeC = ot.measure.DeC;
opt.AUX.TODO_StoV = ot.measure.StoV;
opt.AUX.TODO_SpE = ot.measure.SpE;
opt.AUX.TODO_DEL = ot.measure.DEL;
opt.AUX.TODO_COH = ot.measure.COH;
opt.AUX.TODO_PSDW = ot.measure.PSDW;
opt.AUX.TODO_ALFFW = ot.measure.ALFFW;
opt.AUX.TODO_PSDM = ot.measure.PSDM;
opt.AUX.TODO_ALFFM = ot.measure.ALFFM;
opt.AUX.TODO_RtoR = ot.measure.RtoR;

if (ot.measure.StoV || ot.measure.DEL || ot.measure.COH) && ( ~isfield(fo,'rois') || isempty(fo.rois) )
    error('Please, define folders.roi (i.e., the folder containing ROIs for seed-to-voxel analysis).');
end

opt.AUX.new_rois = ot.new_rois;
%smoothing (special)
opt.AUX.special_smoothing_order = ot.special_smoothing_order;
if opt.AUX.special_smoothing_order == 1
    if isempty(ot.special_smoothing_FWHM)
        error ('No specified smoothing in ER.others.special_smoothing_FWHM.');
    end
    if ot.special_smoothing_FWHM == 0
        opt.AUX.special_smoothing_order = 0;
    end
end
opt.AUX.special_smoothing_FWHM = ot.special_smoothing_FWHM;

if ot.measure.RtoR 
    if ~isfield(ot,'RtoR')
        error('Please, define ER.others.RtoR in order to perform ROI to ROI analysis.');
    elseif ~isfield(ot.RtoR,'roi_folder_or_file') || isempty(ot.RtoR.roi_folder_or_file)
        error('Please, define ER.others.RtoR.roi_folder_or_file. It must contain the path of the folder with the parcellation of the brain.');
    elseif ~isfield(ot.RtoR,'name') || isempty(ot.RtoR.name) 
        error('Please, define ER.others.RtoR.name. It must contain a name to indentify the parcellation chosen.');
    else
        opt.RtoR.name = r_b(ot.RtoR.name);
        opt.RtoR.roi_folder_or_file = ot.RtoR.roi_folder_or_file;
    end
    if ~isfield(ot.RtoR,'COH')
        ot.RtoR.COH = 0;
    end
    opt.RtoR.do_coh = ot.RtoR.COH;
end

warning off backtrace
if opt.AUX.TODO_ALFFW && ( ~opt.AUX.TODO_PSDW && ~opt.MEASURES.PSDW) 
    warning('ALFFW calculation requires the Power Spectral Density with Welch''s method. Changing measure.PSDW from 0 to 1.');
    opt.AUX.TODO_PSDW = 1;
end
if opt.AUX.TODO_ALFFM && ( ~opt.AUX.TODO_PSDM && ~opt.MEASURES.PSDM) 
    warning('ALFFM calculation requires the Power Spectral Density with multitaper method. Changing measure.PSDM from 0 to 1.');
    opt.AUX.TODO_PSDM = 1;
end
if opt.AUX.TODO_SpE && ( ~opt.AUX.TODO_PSD && ~opt.MEASURES.PSD)
    warning('Spectral entropy calculation requires the Power Spectral Density. Changing measure.PSD from 0 to 1.');
    opt.AUX.TODO_PSD = 1;
end
if opt.AUX.TODO_CoE && ( ~opt.AUX.TODO_gFC && ~opt.MEASURES.gFC)
    warning('Correlation entropy calculation requires global Functionl Connectivity (gFC). Changing measure.gFC from 0 to 1.');
    opt.AUX.TODO_gFC = 1;
end
if opt.AUX.TODO_DeC && ( ~opt.AUX.TODO_gFC && ~opt.MEASURES.gFC)
    warning('Degree Connectivity requires global Functionl Connectivity (gFC). Changing measure.gFC from 0 to 1.');
    ot.measure.gFC = 1;
end
warning on backtrace

opt.AUX.seedmeasure = logical (opt.AUX.TODO_StoV || opt.AUX.TODO_DEL || opt.AUX.TODO_COH);  % add here seed based measure
bmi_step = opt.AUX.seedmeasure && (~opt.AUX.bmi_done || opt.AUX.new_rois );

opt.AUX.steps_total = (bmi_step + opt.AUX.TODO_DA + opt.AUX.TODO_REHO + opt.AUX.TODO_PSD +opt.AUX.TODO_gFC + opt.AUX.TODO_CoE + opt.AUX.TODO_StoV +opt.AUX.TODO_SpE +opt.AUX.TODO_DEL + opt.AUX.TODO_COH + opt.AUX.TODO_PSDW + opt.AUX.TODO_ALFFW + opt.AUX.TODO_PSDM + opt.AUX.TODO_ALFFM + opt.AUX.TODO_DeC + 2*opt.AUX.TODO_RtoR);
opt.AUX.steps_current = 0;

if opt.AUX.seedmeasure 
    if isfield(fo,'rois') && ~isempty(fo.rois)
        opt.folders.rois = fo.rois;
    else
        error('You need to specify a ROIs folder (ER.folders.rois).');
    end
end

% Rois check
if opt.AUX.seedmeasure && (~opt.AUX.bmi_done || opt.AUX.new_rois )
    if ~exist(opt.folders.rois,'dir') 
        error('Folders.rois doesn''t exist.')
    end
    [R_NII,list_R_NII] = conn_dir([opt.folders.rois,'/*.nii']);
    [R_HDR,list_R_HDR] = conn_dir([opt.folders.rois,'/*.hdr']);
    opt.DATA.ROIS = strvcat(R_NII,R_HDR);
    list_ROIS = cat(2,list_R_NII,list_R_HDR);
    if isempty(opt.DATA.ROIS)
        warning('Roi folder is empty. Seed-based measures will NOT run.');
    else
        %--------------check roi size and binarity-------------------------
        %[status,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]); %assuming data have all the same size
        [status,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]); %assuming data have all the same size
        data_info = str2num(data_info);
        for z = 1:length(list_ROIS)
            %[status,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.DATA.ROIS(z,:)]);
            [status,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.DATA.ROIS(z,:)]);
            roi_info = fix_float_errors(roi_info);
            roi_info = str2num(roi_info);
            if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' has not the same size of functional data.']);
            elseif roi_info(4)~=1
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' has more than one sub-brick.']);
            elseif roi_info(5)~=0 || roi_info(6)~=1
                error(['ROI: ', c_b(opt.DATA.ROIS(z,:)),' is not binary.']);
            end
        end
        %------------------------------------------------------------------
        roi_names = [];
        for k = 1:length(list_ROIS)
            roi_names = strvcat(roi_names,list_ROIS(k).name);
        end
        opt.AUX.roi_names = roi_names;
    end
    if opt.AUX.bmi_done 
        system(['rm ',opt.folders.results,'/seed/*']);
    end
    opt.AUX.bmi_done = 0;
end

% Rois check for RtoR
if opt.AUX.TODO_RtoR
    switch exist(opt.RtoR.roi_folder_or_file)
        case 0  % no valid input
            error('Roi folder/file for ROI to ROI connectivity doesn''t exist.')
        case 7  % it's a directory
            opt.RtoR.modalitiy = 'roi_list';
            [R_NII,list_R_NII] = conn_dir([opt.RtoR.roi_folder_or_file,'/*.nii']);
            [R_HDR,list_R_HDR] = conn_dir([opt.RtoR.roi_folder_or_file,'/*.hdr']);
            opt.RtoR.ROIS = strvcat(R_NII,R_HDR);
            list_ROIS = cat(2,list_R_NII,list_R_HDR);
            if isempty(opt.RtoR.ROIS)
                warning off backtrace
                warning('RtoR: Roi folder for ROI to ROI connectivity is empty. ROI to ROI connectivity will NOT run.');
                opt.AUX.TODO_RtoR = 0;
                warning on backtrace
            else
                %--------------check roi size and binarity-------------------------
                [~,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
                data_info = str2num(data_info);
                for z = 1:length(list_ROIS)
                    [~,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.RtoR.ROIS(z,:)]);
                    roi_info = fix_float_errors(roi_info);
                    roi_info = str2num(roi_info);
                    if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' has not the same size of functional data.']);
                    elseif roi_info(4)~=1
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' has more than one sub-brick.']);
                    elseif roi_info(5)~=0 || roi_info(6)~=1
                        error(['RtoR. ROI: ', c_b(opt.RtoR.ROIS(z,:)),' is not binary.']);
                    end
                end
                %------------------------------------------------------------------
                rois_string = [];
                for l = 1:size(opt.RtoR.ROIS,1)
                    tmp_str = c_b(opt.RtoR.ROIS(l,:));
                    rois_string = [rois_string,' ',tmp_str];
                end
                [~,~] = system(['3dMean -prefix _temp_rtor.nii -sum',rois_string]);
                rtor_sum = spm_read_vols(spm_vol('_temp_rtor.nii'));
                n_overlapping = length(rtor_sum(rtor_sum > 1));
                system('rm _temp_rtor.nii');
                opt.RtoR.n_ovelapping = 0;
                if n_overlapping > 0
                    warning off backtrace
                    warning('RtoR: Overllaping ROIs found (%.0d voxels). Overllaping voxels will be removed unless they don''t fall in the group brain mask.',n_overlapping);
                    warning on backtrace
                    opt.RtoR.n_ovelapping = n_overlapping;
                    rtor_sum(rtor_sum == 1) = 0;
                    rtor_sum(rtor_sum > 1) = 1;
                    opt.RtoR.ovelapping = rtor_sum;
                else
                    clear rtor_sum;
                end
                roi_names = [];
                for k = 1:length(list_ROIS)
                    roi_names = strvcat(roi_names,list_ROIS(k).name(1:end-4));
                end
                opt.RtoR.roi_names = roi_names;
            end
        case 2  % it's a file mask
            opt.RtoR.modalitiy = 'mask_file';
            %--------------check roi size ------------------------------------
            [~,data_info] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
            data_info = str2num(data_info);
            [~,roi_info] = system(['3dinfo -ni -nj -nk -nv -dminus -dmaxus ',opt.RtoR.roi_folder_or_file]);
            roi_info = str2num(roi_info);
            if roi_info(1)~=data_info(1) || roi_info(2)~=data_info(2) || roi_info(3)~=data_info(3) 
                error('RtoR: The mask dataset (or parecellation atlas) has not the same size of functional data.');
            end
            %------------------------------------------------------------------
            mask_dataset = spm_read_vols(spm_vol(opt.RtoR.roi_folder_or_file));
            max_roi = max(mask_dataset(:));
            mask_unique = unique(mask_dataset(:));
            checker = max_roi - length(mask_unique) +1;
            if checker ~= 0
                error('RtoR: The mask dataset has missing roi indexes. Create a mask dataset (or parcellation atlas) without "holes".');
            end
            opt.RtoR.mf_max_roi = max_roi;
            if exist([opt.RtoR.roi_folder_or_file(1:end-3),'txt'],'file') == 2
                fid = fopen([opt.RtoR.roi_folder_or_file(1:end-3),'txt']);
                line = fgets(fid);
                roi_names = [];
                while ischar(line)    
                    roi_names = strvcat(roi_names,line); 
                    line = fgets(fid);
                end
                opt.RtoR.roi_names = roi_names;
                if size(roi_names,1) ~= max_roi
                    error('RtoR: The .txt file with labels for the mask dataset has a wrong number of lines.');
                end
                opt.RtoR.mf_txt_exist = 1;
                opt.RtoR.mf_txt_file = [opt.RtoR.roi_folder_or_file(1:end-3),'txt'];
            else
                warning off backtrace;
                warning('RtoR: No labels file (txt) found for ROI to ROI connectivity. ROIs will be named as sequential numbers');
                warning on backtrace;
                indexs = 1:1:opt.RtoR.mf_max_roi;
                indexs = indexs';
                opt.RtoR.roi_names = num2str(indexs);
                opt.RtoR.mf_txt_exist = 0;
            end
    end
end

if opt.AUX.TODO_COH || (opt.AUX.TODO_RtoR && opt.RtoR.do_coh)
    coherency_initialization(ot);
end
if opt.AUX.TODO_PSDW; psdw_initialization(ot); end
if opt.AUX.TODO_ALFFW; ALFF_welch_initialization; end
if opt.AUX.TODO_PSDM; psdm_initialization(ot); end
if opt.AUX.TODO_ALFFM; ALFF_multitaper_initialization; end
if opt.AUX.TODO_PSD; psd_initialization; end
if opt.AUX.TODO_gFC; gFC_initialization(ot); end
if opt.AUX.TODO_DeC; DeC_initialization(ot); end

warning off backtrace
if  opt.bm.use_ss == 1
    warning(sprintf('\nThis prject uses subject-specific brainmask. Be aware:\n1) Voxel-wise analysis is not possible in this modality.\n2) ROI-based analysis might be baised.\n3) Comparison of standardized measures might be biased.'));
end
warning on backtrace

opt.date_modification = date;

save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');

return
end

function [bmask_string] = brainmask(flag_bmask)
global opt
%-------first hendle the simplest cases: method = 3,4 and 5----------------
if flag_bmask == 3  % apriori brainmask 
    warning off backtrace;warning('BM modality 3 has no zero-voxel checking'); warning on backtrace
    mkdir(opt.folders.results,'brainmask_m3');
    dest_dir = [opt.folders.results,'/brainmask_m3'];
    system(['cp ',opt.bm.apriori,' ',dest_dir,'/group_bm.nii']);
    bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
    bmask_path{1,1} = [dest_dir,'/group_bm.nii'];
    bmask_img{1,1} = spm_read_vols(spm_vol(bmask_path{1,1}));
    fprintf ('Percent done: 100.0');
elseif flag_bmask == 4 %brain mask from mask_file 
    warning off backtrace;warning('BM modality 3 has no zero-voxel checking'); warning on backtrace
    mkdir(opt.folders.results,'brainmask_m4');
    dest_dir = [opt.folders.results,'/brainmask_m4'];
    cd (dest_dir)
    total = sum(opt.subject_number);
    reverseStr = '';
    count = 0;
    if ~opt.bm.use_ss   %definition of GROUP bm
        for g = 1:opt.group_number
            for i=1:opt.subject_number(g)
                count = count +1;
                temp = spm_read_vols(spm_vol(opt.DATA.MASKS{g}(i,:)));
                if count == 1
                    hdr_output = spm_vol(opt.DATA.MASKS{g}(i,:));
                    s = size(temp);
                    group_bm = ones(s(1),s(2),s(3));
                end
                if isempty(opt.bm.ss_masks_thr)
                    group_bm = group_bm.*temp;
                else
                    group_bm = group_bm + temp;
                end

               percentDone = 100 * count / total;
               msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
               fprintf([reverseStr, msg]);
               reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
        if ~isempty(opt.bm.ss_masks_thr)
            group_bm = group_bm./count;
            group_bm(group_bm <= opt.bm.ss_masks_thr) = 0;
            group_bm(group_bm > opt.bm.ss_masks_thr) = 1;
        end
        hdr_output.fname = 'group_bm.nii';
        hdr_output.private.dat.fname = 'group_bm.nii';
        spm_write_vol(hdr_output,group_bm);

        bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
        bmask_path{1.1} = [dest_dir,'/group_bm.nii'];
        bmask_img{1,1} = group_bm;
    else  % NEW: USE SS mask as brainmask, NO GROUP MASK. NOTE: I could have use directly DATA.MASKS. But in this way masks remain even if data delated.
        for g = 1:opt.group_number
            for i=1:opt.subject_number(g)
                count = count +1;
                temp = spm_read_vols(spm_vol(opt.DATA.MASKS{g}(i,:)));
                if i == 1 &&  g== 1
                    hdr_output = spm_vol(opt.DATA.MASKS{g}(i,:));
                    s = size(temp);
                    %group_bm = ones(s(1),s(2),s(3));
                end
                hdr_output.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
                hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
                spm_write_vol(hdr_output,temp);
                
                bmask_string{g,i} = ['-mask ',dest_dir,'/',hdr_output.fname];
                bmask_path{g,i} = [dest_dir,'/',hdr_output.fname];
                bmask_img{g,i} = temp;
                
               percentDone = 100 * count / total;
               msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
               fprintf([reverseStr, msg]);
               reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
elseif flag_bmask == 5 %automask (3dautomask)
    mkdir(opt.folders.results,'brainmask_m5');
    dest_dir = [opt.folders.results,'/brainmask_m5'];
    cd (dest_dir)
    
    if opt.bm.m5_dilate 
        dilate = '-dilate 1 ';
    else
        dilate = '';
    end
    if ~isempty(opt.bm.automask_session_selector)
        session = opt.bm.automask_session_selector;
    else
        session = 1:opt.session_number;
    end
    total = sum(opt.subject_number)*length(session);
    reverseStr = '';
    count = 0;
    
    for g = 1:opt.group_number
        for i=1:opt.subject_number(g)
            for k=session
                count = count +1;
                [status, string] = system(['3dAutomask -prefix ',opt.AUX.output_name{g,k,i},'_AUTOMASK.nii -clfrac 0.1 ',dilate,c_b(opt.DATA.FUNCTIONALS{g,k}(i,:))]);control(status,string);
                temp = spm_read_vols(spm_vol([opt.AUX.output_name{g,k,i},'_AUTOMASK.nii']));
                if count == 1
                    hdr_output = spm_vol([opt.AUX.output_name{g,k,i},'_AUTOMASK.nii']);
                    s = size(temp);
                    group_bm = ones(s(1),s(2),s(3));
                end
                group_bm = group_bm.*temp;
                
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
       
    hdr_output.fname = 'group_bm.nii';
    hdr_output.private.dat.fname = 'group_bm.nii';
    spm_write_vol(hdr_output,group_bm);
    
    bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
    bmask_path{1,1} = [dest_dir,'/group_bm.nii']; 
    bmask_img{1,1} = group_bm;
%-------now hendle the dirty cases: method = 1 and 2 ----------------------
elseif flag_bmask == 1 || flag_bmask == 2 %%brain mask creation from structural tissues
    mkdir(opt.folders.results,['brainmask_m',num2str(flag_bmask)]);
    dest_dir = [opt.folders.results,'/brainmask_m',num2str(flag_bmask)];
    cd (dest_dir)
    
    if flag_bmask == 1
        wc1_thres = opt.bm.tissue_thr(1);
        wc2_thres = opt.bm.tissue_thr(2);
        wc3_thres = opt.bm.tissue_thr(3);
    else
        wc1_thres = opt.bm.tissue_thr(1);
    end
        
    if ~isempty(opt.bm.automask_session_selector) && opt.bm.doautomask 
        session = opt.bm.automask_session_selector;
    else
        session = 1:opt.session_number;
    end
    if isempty(opt.bm.group) || opt.bm.use_ss
        total = sum(opt.subject_number) + sum(opt.subject_number)*length(session);  % 2 step: gm thresholding + interesection with 3dautomask
        groups = 1:opt.group_number;
    else
        total = sum(opt.subject_number(opt.bm.group)) + sum(opt.subject_number)*length(session);
        groups = opt.bm.group;
    end
    %-------------3dautomask-----------------------------------------------
    if opt.bm.doautomask 
        [count,reverseStr,groupAUTOMASK] = bm_method_1_2_automask(total);  % TOO SLOW, maybe I can use just the meanEPI, CHECK
    else
        if ~opt.bm.use_ss
            [count,reverseStr,groupAUTOMASK,~] = remove_EPI_zerovoxels(total);
        else
            count = 0;
            reverseStr = '';
        end
    end
    %----------------------------------------------------------------------
    %-------------initialize-----------------------------------------------
    if flag_bmask == 1
        if ~strcmpi(opt.bm.modality,'hard');WC1 = [];WC2 = [];WC3 = [];end
    else
        if ~strcmpi(opt.bm.modality,'hard');WC1 = [];end
    end
    hdr_output =spm_vol(opt.DATA.WC1{1}(1,:));
    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{1}(1,:)));
    s = size(wc1);
    group_bm = ones(s(1),s(2),s(3));
    %----------------------------------------------------------------------
    if ~opt.bm.use_ss  
        for g = groups
            for i=1:opt.subject_number(g)
                count = count +1;
                if flag_bmask == 1
                    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
                    wc2 = spm_read_vols(spm_vol(opt.DATA.WC2{g}(i,:)));  
                    wc3 = spm_read_vols(spm_vol(opt.DATA.WC3{g}(i,:)));   
                else
                    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
                end
                if strcmpi(opt.bm.modality,'hard')
                    if i == 1;single_group_mask = ones(s(1),s(2),s(3));end
                    if flag_bmask == 1
                        wc1(wc1<=wc1_thres) = 0;
                        wc2(wc2<=wc2_thres) = 0;
                        wc3(wc3<=wc3_thres) = 0;
                        wc1(wc1>0) = 1;
                        wc2(wc2>0) = 1;
                        wc3(wc3>0) = 1;
                        bm = wc1 +wc2 +wc3;
                        bm(bm >0) = 1;
                    else
                        wc1(wc1<=wc1_thres) = 0;
                        wc1(wc1>0) = 1;
                        bm = wc1;
                    end
                    bm = bm.*groupAUTOMASK;
                    if opt.bm.save_ss
                        hdr_output.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_HARDbm.nii'];
                        hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_HARDbm.nii'];
                        spm_write_vol(hdr_output,bm);
                    end
                    group_bm = group_bm.*bm;
                    single_group_mask = single_group_mask.*bm;
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                else  %SOFT
                    if flag_bmask == 1
                        if i == 1; WC1_sg = [];WC2_sg = [];WC3_sg = []; end
                        WC1 = cat(4,WC1,wc1);
                        WC2 = cat(4,WC2,wc2);
                        WC3 = cat(4,WC3,wc3);
                        WC1_sg = cat(4,WC1_sg,wc1);
                        WC2_sg = cat(4,WC2_sg,wc2);
                        WC3_sg = cat(4,WC3_sg,wc3);
                    else
                        if i == 1; WC1_sg = []; end;
                        WC1 = cat(4,WC1,wc1);
                        WC1_sg = cat(4,WC1_sg,wc1);
                        
                    end
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
            end
            if opt.group_number > 1
                if strcmpi(opt.bm.modality,'hard')
                %------------------save single group mask -------------------------
                    hdr_output.fname = [opt.AUX.group_names_prefix{g},'HARDbm.nii'];
                    hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},'HARDbm.nii'];
                    spm_write_vol(hdr_output,single_group_mask);
                else
                    if flag_bmask == 1
                        WC1_sg = mean(WC1_sg,4);
                        WC2_sg = mean(WC2_sg,4);
                        WC3_sg = mean(WC3_sg,4);
                        WC1_sg(WC1_sg<=wc1_thres) = 0;
                        WC2_sg(WC2_sg<=wc2_thres) = 0;
                        WC3_sg(WC3_sg<=wc3_thres) = 0;
                        WC1_sg(WC1_sg>0) = 1;
                        WC2_sg(WC2_sg>0) = 1;
                        WC3_sg(WC3_sg>0) = 1;
                        sgm = WC1_sg +WC2_sg +WC3_sg;
                        sgm(sgm >0) = 1;
                    else
                        WC1_sg = mean(WC1_sg,4);
                        WC1_sg(WC1_sg<=wc1_thres) = 0;
                        WC1_sg(WC1_sg>0) = 1;
                        sgm = WC1_sg;  
                    end
                    sgm = sgm.*groupAUTOMASK;
                    hdr_output.fname = [opt.AUX.group_names_prefix{g},'SOFTbm.nii'];
                    hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},'SOFTbm.nii'];
                    spm_write_vol(hdr_output,sgm);
                end
                %------------------------------------------------------------------ 
            end
        end
        if ~strcmpi(opt.bm.modality,'hard')
            if flag_bmask == 1
                WC1 = mean(WC1,4);
                WC2 = mean(WC2,4);
                WC3 = mean(WC3,4);
                WC1(WC1<=wc1_thres) = 0;
                WC2(WC2<=wc2_thres) = 0;
                WC3(WC3<=wc3_thres) = 0;
                WC1(WC1>0) = 1;
                WC2(WC2>0) = 1;
                WC3(WC3>0) = 1;
                bm = WC1 +WC2 +WC3;
                bm(bm >0) = 1;
            else
                WC1 = mean(WC1,4);
                WC1(WC1<=wc1_thres) = 0;
                WC1(WC1>0) = 1;
                bm = WC1;
            end
            group_bm = bm.*groupAUTOMASK;
        end
        hdr_output.fname = 'group_bm.nii';
        hdr_output.private.dat.fname = 'group_bm.nii';
        spm_write_vol(hdr_output,group_bm);
        bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
        bmask_path{1,1} = [dest_dir,'/group_bm.nii']; 
        bmask_img{1,1} = group_bm;
   else  % USE SS mask as brainmask, NO GROUP MASK.
        for g = groups
            for i=1:opt.subject_number(g)
                count = count +1;
                if flag_bmask == 1
                    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
                    wc2 = spm_read_vols(spm_vol(opt.DATA.WC2{g}(i,:)));
                    wc3 = spm_read_vols(spm_vol(opt.DATA.WC3{g}(i,:)));
                    wc1(wc1<=wc1_thres) = 0;
                    wc2(wc2<=wc2_thres) = 0;
                    wc3(wc3<=wc3_thres) = 0;
                    wc1(wc1>0) = 1;
                    wc2(wc2>0) = 1;
                    wc3(wc3>0) = 1;
                    bm = wc1 +wc2 +wc3;
                    bm(bm >0) = 1;
                else
                    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
                    wc1(wc1<=wc1_thres) = 0;
                    wc1(wc1>0) = 1;
                    bm = wc1;
                end
                bmZ = remove_EPI_zerovoxels_SS(g,i);
                bm = bm.*bmZ; 
                %bm = bm.*groupAUTOMASK;  %this must be done subject by subject...TOFIX
                hdr_output.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
                hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
                spm_write_vol(hdr_output,bm);
                
                    percentDone = 100 * count / total;
                    msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));

                bmask_string{g,i} = ['-mask ',dest_dir,'/',hdr_output.fname];
                bmask_path{g,i} = [dest_dir,'/',hdr_output.fname];
                bmask_img{g,i} = bm;
            end
        end
    end
%--------------------------------------------------------------------------
% elseif flag_bmask == 2 %%brain mask creation from structural ONLY GM
%     mkdir(opt.folders.results,'brainmask_m2');
%     dest_dir = [opt.folders.results,'/brainmask_m2'];
%     cd (dest_dir)
%     wc1_thres = opt.bm.gm_thr;
%     if ~isempty(opt.bm.automask_session_selector) && opt.bm.doautomask 
%         session = opt.bm.automask_session_selector;
%     else
%         session = 1:opt.session_number;
%     end
%     if isempty(opt.bm.group) || opt.bm.use_ss
%         total = sum(opt.subject_number) + sum(opt.subject_number)*length(session);  % 2 step: gm thresholding + interesection with 3dautomask
%         groups = 1:opt.group_number;
%     else
%         total = sum(opt.subject_number(opt.bm.group)) + sum(opt.subject_number)*length(session);
%         groups = opt.bm.group;
%     end
%     %-------------3dautomask-----------------------------------------------
%     if opt.bm.doautomask 
%         [count,reverseStr,groupAUTOMASK] = bm_method_1_2_automask(total);  % TOO SLOW, maybe I can use just the meanEPI, CHECK
%     else
%         if ~opt.bm.use_ss
%             [count,reverseStr,groupAUTOMASK,~] = remove_EPI_zerovoxels(total);
%         else
%             count = 0;
%             reverseStr = '';
%         end
%     end
%     %----------------------------------------------------------------------
%     %-------------initialize-----------------------------------------------
%     if ~strcmpi(opt.bm.modality,'hard');WC1 = [];end
%     hdr_output =spm_vol(opt.DATA.WC1{1}(1,:));
%     wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{1}(1,:)));
%     s = size(wc1);
%     group_bm = ones(s(1),s(2),s(3));
%     %----------------------------------------------------------------------
%     if ~opt.bm.use_ss   %old case, definition of GROUP bm
%         for g = groups
%             for i=1:opt.subject_number(g)
%                 count = count +1;
%                 wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
%                 if strcmpi(opt.bm.modality,'hard')
%                     if i == 1;single_group_mask = ones(s(1),s(2),s(3));end
%                     wc1(wc1<=wc1_thres) = 0;
%                     wc1(wc1>0) = 1;
%                     bm = wc1;
%                     bm = bm.*groupAUTOMASK;
%                     if opt.bm.save_ss
%                         hdr_output.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_HARDbm.nii'];
%                         hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_HARDbm.nii'];
%                         spm_write_vol(hdr_output,bm);
%                     end
%                     group_bm = group_bm.*bm;
%                     single_group_mask = single_group_mask.*bm;
%                     percentDone = 100 * count / total;
%                     msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
%                     fprintf([reverseStr, msg]);
%                     reverseStr = repmat(sprintf('\b'), 1, length(msg));
%                 else
%                     if i == 1; WC1_sg = [];end
%                     WC1 = cat(4,WC1,wc1);
%                     WC1_sg = cat(4,WC1_sg,wc1);
%                     percentDone = 100 * count / total;
%                     msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
%                     fprintf([reverseStr, msg]);
%                     reverseStr = repmat(sprintf('\b'), 1, length(msg));
%                 end
%             end
%             if opt.group_number > 1
%                 if strcmpi(opt.bm.modality,'hard')
%                 %------------------save single group mask -------------------------
%                     hdr_output.fname = [opt.AUX.group_names_prefix{g},'HARDbm.nii'];
%                     hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},'HARDbm.nii'];
%                     spm_write_vol(hdr_output,single_group_mask);
%                 else
%                     WC1_sg = mean(WC1_sg,4);
%                     WC1_sg(WC1_sg<=wc1_thres) = 0;
%                     WC1_sg(WC1_sg>0) = 1;
%                     sgm = WC1_sg;
%                     sgm = sgm.*groupAUTOMASK;
%                     hdr_output.fname = [opt.AUX.group_names_prefix{g},'SOFTbm.nii'];
%                     hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},'SOFTbm.nii'];
%                     spm_write_vol(hdr_output,sgm);
%                 end
%                 %------------------------------------------------------------------    
%             end
%         end
%         if ~strcmpi(opt.bm.modality,'hard')
%             WC1 = mean(WC1,4);
%             WC1(WC1<=wc1_thres) = 0;
%             WC1(WC1>0) = 1;
%             bm = WC1;
%             group_bm = bm.*groupAUTOMASK;
%         end
%         hdr_output.fname = 'group_bm.nii';
%         hdr_output.private.dat.fname = 'group_bm.nii';
%         spm_write_vol(hdr_output,group_bm);
%         bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
%         bmask_path{1,1} = [dest_dir,'/group_bm.nii']; 
%         bmask_img{1,1} = group_bm;
%     else  % USE SS mask as brainmask, NO GROUP MASK.
%         for g = groups
%             for i=1:opt.subject_number(g)
%                 count = count +1;
%                 wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(i,:)));
%                 wc1(wc1<=wc1_thres) = 0;
%                 wc1(wc1>0) = 1;
%                 bm = wc1;
%                 bmZ = remove_EPI_zerovoxels_SS(g,i);
%                 bm = bm.*bmZ; 
%                 %bm = bm.*groupAUTOMASK;  %this must be done subject by subject...TOFIX
%                 hdr_output.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
%                 hdr_output.private.dat.fname = [opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i},'_bm.nii'];
%                 spm_write_vol(hdr_output,bm);
%                 
%                     percentDone = 100 * count / total;
%                     msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
%                     fprintf([reverseStr, msg]);
%                     reverseStr = repmat(sprintf('\b'), 1, length(msg));
% 
%                 bmask_string{g,i} = ['-mask ',dest_dir,'/',hdr_output.fname];
%                 bmask_path{g,i} = [dest_dir,'/',hdr_output.fname];
%                 bmask_img{g,i} = bm;
%             end
%         end
%     end
    
elseif flag_bmask == 6 %%brain mask creation eliminating only zero voxels
    mkdir(opt.folders.results,'brainmask_m6');
    dest_dir = [opt.folders.results,'/brainmask_m6'];
    cd (dest_dir)

%     if ~isempty(opt.bm.automask_session_selector)
%         session = opt.bm.automask_session_selector;
%     else
%         session = 1:opt.session_number;
%     end
    session = 1:opt.session_number;
    total = sum(opt.subject_number)*length(session);
    reverseStr = '';
    count = 0;
    
    [~,~,group_bm,hdr] = remove_EPI_zerovoxels(total);
    
    hdr_output = hdr(1);
    
    hdr_output.fname = 'group_bm.nii';
    hdr_output.private.dat.fname = 'group_bm.nii';
    spm_write_vol(hdr_output,group_bm);
    
    bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
    bmask_path{1,1} = [dest_dir,'/group_bm.nii']; 
    bmask_img{1,1} = group_bm;
elseif flag_bmask == 7 %% 3dskullstrip
    mkdir(opt.folders.results,'brainmask_m7');
    dest_dir = [opt.folders.results,'/brainmask_m7'];
    cd (dest_dir)
    
    session = 1:opt.session_number;
    total = sum(opt.subject_number)*length(session);
    reverseStr = '';
    count = 0;
    
    for g = 1:opt.group_number
        for i=1:opt.subject_number(g)
            for k=session
                count = count +1;
                [status, string] = system(['3dSkullStrip -input ',c_b(opt.DATA.FUNCTIONALS{g,k}(i,:)),' -mask_vol -prefix _temp.nii']);control(status,string);
                [status, string] = system(['3dcalc -a _temp.nii -prefix _temp2.nii -expr ''step(a-3)''']);control(status,string);
                [status, string] = system(['3dmask_tool -input _temp2.nii -prefix ',opt.AUX.output_name{g,k,i},'_SKULLSTRIPPED.nii -dilate_input 1 -fill_holes']);control(status,string);
                [~, ~] = system('rm _tem*.nii');
                temp = spm_read_vols(spm_vol([opt.AUX.output_name{g,k,i},'_SKULLSTRIPPED.nii']));
                if count == 1
                    hdr_output = spm_vol([opt.AUX.output_name{g,k,i},'_SKULLSTRIPPED.nii']);
                    s = size(temp);
                    group_bm = ones(s(1),s(2),s(3));
                end
                group_bm = group_bm.*temp;
                
                percentDone = 100 * count / total;
                msg = sprintf('Percent done: %3.1f (%s)', percentDone,[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,i}]);
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
        end
    end
       
    hdr_output.fname = 'group_bm.nii';
    hdr_output.private.dat.fname = 'group_bm.nii';
    spm_write_vol(hdr_output,group_bm);
    
    bmask_string{1,1} = ['-mask ',dest_dir,'/group_bm.nii'];
    bmask_path{1,1} = [dest_dir,'/group_bm.nii']; 
    bmask_img{1,1} = group_bm;
    
end

opt.AUX.bmask_string = bmask_string;
opt.AUX.bmask_path = bmask_path;
opt.AUX.bmask_img = bmask_img;
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');

return
end

function [count,reverseStr,group_bm,hdr] = remove_EPI_zerovoxels(total)
% Define a groupbm based only on the EPI time series of all subjs. It
% removes voxels which timeseries have zero std and NaN voxels
global opt
reverseStr = '';
count = 0;
% if ~isempty(opt.bm.automask_session_selector)
%     session = opt.bm.automask_session_selector;
% else
%     session = 1:opt.session_number;
% end
session = 1:opt.session_number;
for g = 1:opt.group_number
    for k=1:opt.subject_number(g)
        for j=session
            count = count +1;
            %[status, string] = system(['3dAutomask -prefix _temp_AUTOMASK.nii -clfrac 0.1 ',c_b(opt.DATA.FUNCTIONALS{g,k}(i,:))]);control(status,string);
            temp = spm_read_vols(spm_vol(c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))));
            if count == 1
                s = size(temp);
                group_bm = ones(s(1),s(2),s(3));
                hdr = spm_vol(c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)));
            end
            stdv = std(temp,[],4);
            a = zeros(size(stdv));
            a(stdv == 0) = 1;
            a(isnan(stdv)) = 1;   %also remove NaNs
            if opt.bm.save_ss
                bm = ones(size(a));
                bm(logical(a)) = 0;
                hdr_output = hdr(1);
                hdr_output.fname = [opt.AUX.output_name{g,j,k},'_noSignalVoxelRemoved.nii'];
                hdr_output.private.dat.fname = [opt.AUX.output_name{g,j,k},'_noSignalVoxelRemoved.nii'];
                spm_write_vol(hdr_output,bm);
            end
            group_bm(logical(a)) = 0;
            percentDone = 100 * count / total;
            msg = sprintf('Percent done: %3.1f (%s)', percentDone,['Removing zero signal voxels']);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
end
return
end

function [bmZ] = remove_EPI_zerovoxels_SS(g,k)
% It removes voxels which timeseries have zero std %also remove NaNs. Works for single
% subject (intersection across sessions)
global opt

session = 1:opt.session_number;
for j=session
    temp = spm_read_vols(spm_vol(c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))));
    if j == 1
        s = size(temp);
        bmZ = ones(s(1),s(2),s(3));
    end
    stdv = std(temp,[],4);
    a = zeros(size(stdv));
    a(stdv == 0) = 1;
    a(isnan(stdv)) = 1;   %also remove NaNs
    bmZ(logical(a)) = 0;
end
return
end

function [count,reverseStr,group_bm] = bm_method_1_2_automask(total)
% Define a groupbm based only on the EPI time series of all subjs
global opt
reverseStr = '';
count = 0;
if ~isempty(opt.bm.automask_session_selector)
    session = opt.bm.automask_session_selector;
else
    session = 1:opt.session_number;
end
for g = 1:opt.group_number
    for i=1:opt.subject_number(g)
        for k=session
            count = count +1;
            [status, string] = system(['3dAutomask -prefix _temp_AUTOMASK.nii -clfrac 0.1 ',c_b(opt.DATA.FUNCTIONALS{g,k}(i,:))]);control(status,string);
            temp = spm_read_vols(spm_vol('_temp_AUTOMASK.nii'));
            if count == 1
                s = size(temp);
                group_bm = ones(s(1),s(2),s(3));
            end
            group_bm = group_bm.*temp;
            [status, string] = system('rm _temp_AUTOMASK.nii');control(status,string);
            percentDone = 100 * count / total;
            msg = sprintf('Percent done: %3.1f (%s)', percentDone,['3dAUTOMASK']);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end
end
return
end

function brainmask_ROI_intersection
global opt
cw = pwd;
cd (opt.folders.seed)
total = size(opt.DATA.ROIS,1);
reverseStr  = '';
err_indx = zeros(size(opt.DATA.ROIS,1),1);
war_indx = err_indx;
v_count = zeros(size(opt.DATA.ROIS,1),2);
for l=1:size(opt.DATA.ROIS,1)   %ciclo rois
    roi_str = c_b(opt.AUX.roi_names(l,:));
    roi_str = roi_str(1:end-4);
    [status,string] =system(['3dcalc -a ',opt.DATA.ROIS(l,:),' -b ',opt.AUX.bmask_path{1,1},' -prefix _bmi_',roi_str,'.nii -expr ''a*b''']);control(status,string);
    bef = spm_read_vols(spm_vol(c_b(opt.DATA.ROIS(l,:))));
    aft = spm_read_vols(spm_vol(['_bmi_',roi_str,'.nii']));
    v_bef = int64(nansum(bef(:)));   %can be NaN if produced by SPM
    v_aft = int64(nansum(aft(:)));
    if v_aft == 0
        err_indx(l,1) = 1;
    elseif v_bef ~= v_aft
        war_indx(l,1) = 1;
        v_count(l,1) = v_bef;
        v_count(l,2) = v_aft;
    end
    percentDone = 100 * l / total;
    msg = sprintf('Percent done: %3.1f (%s)', percentDone,roi_str);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');
warning off backtrace
if sum(err_indx) == total
    opt.AUX.seedmeasure = 0;
    opt.AUX.steps_total = opt.AUX.steps_total -opt.AUX.TODO_StoV - opt.AUX.TODO_DEL - opt.AUX.TODO_COH;
    warning('ATTENTION! No ROIs survived after brain mask intersection. Skipping seed-based calculation.');
elseif sum(err_indx)~=0
    for l=1:length(err_indx)
        if err_indx(l,1)
            str_ = c_b(opt.AUX.roi_names(l,:));
            str_ = str_(1:end-4);
            warning('ATTENTION! ROI %s is empty after intersection with brain mask. Skipping this ROI.',str_);
            [status,string] = system(['rm _bmi_',str_,'.nii']);control(status,string);
        end
    end
end
if sum(war_indx)~=0
%     fprintf('\n');
    for l=1:length(war_indx)
        if war_indx(l,1)
            str_ = c_b(opt.AUX.roi_names(l,:));
            str_ = str_(1:end-4);
            warning('ROI %s has been reshaped due to intersection with brain mask (from %.0d to %.0d voxels).',str_,v_count(l,1) ,v_count(l,2))
        end
    end
end
warning on backtrace
[R_bmi,list_bmi]= conn_dir('_bmi_*.nii');
opt.DATA.ROIS = R_bmi;
roi_names = [];
for k = 1:length(list_bmi)
    roi_names = strvcat(roi_names,list_bmi(k).name(6:end)); 
end
opt.AUX.roi_names = roi_names;  % potrebbe essere ridondante, ma per sicurezza riassocio i nomi a ciascuna roi
cd (cw);
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
return
end

function brainmask_ROI_intersection_SS
% this intersection produes no _bmi_ files. It is only for knowing the
% number of voxels surviving for each subjects (AND for removing ROIs with
% at least a subject empty seed)
global opt
total = size(opt.DATA.ROIS,1);
reverseStr  = '';
err_indx = zeros(size(opt.DATA.ROIS,1),1);
v_count = cell(total);
for l=1:size(opt.DATA.ROIS,1)   %ciclo rois
    roi_str = c_b(opt.AUX.roi_names(l,:));
    roi_str = roi_str(1:end-4);
    for g =1:opt.group_number
            for k=1:opt.subject_number(g)
                bef = spm_read_vols(spm_vol(c_b(opt.DATA.ROIS(l,:))));
                aft = opt.AUX.bmask_img{g,k}.*bef;
                %v_bef = int64(sum(bef(:)));
                v_aft = int64(nansum(aft(:)));
                if v_aft == 0
                    err_indx(l,1) = err_indx(l,1) + 1;
                end
                v_count{l}(g,k) = v_aft;
            end
    end
    percentDone = 100 * l / total;
    msg = sprintf('Percent done: %3.1f (%s)', percentDone,roi_str);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');
warning off backtrace
if sum(logical(err_indx)) == total
    opt.AUX.seedmeasure = 0;
    opt.AUX.steps_total = opt.AUX.steps_total -opt.AUX.TODO_StoV - opt.AUX.TODO_DEL - opt.AUX.TODO_COH;
    warning('ATTENTION! No ROIs survived after brain mask intersection. Skipping seed-based calculation.');
else
    %PLOTTING RESULS
    f = figure('Name','ROIs intersection subject-specific maks','NumberTitle','off');
    for g =1:opt.group_number
        subplot(opt.group_number,1,g)
        hold on
        colors = hsv(opt.subject_number(g));
        for l=1:size(opt.DATA.ROIS,1)
            leg_indx = [];
            %subj_names = [];
            tmp = c_b(opt.AUX.roi_names(l,:));
            roi_n{l} = tmp(1:end-4);
            for k = 1:opt.subject_number(g)
                %subj_names{k} = ['subj ',num2str(k)];
                h{k} = plot(l,v_count{l}(g,k),'.','Color',colors(k,:),'MarkerFaceColor',colors(k,:),'MarkerSize',7);
                leg_indx = [leg_indx, h{k}];
            end
        end
        set(gca,'XTick',[1:1:size(opt.DATA.ROIS,1)])
        set(gca,'XTickLabel',roi_n)
        xlim([0 (size(opt.DATA.ROIS,1) + 1)]);
        legend(leg_indx,opt.AUX.subject_names);
        ylabel('Voxel number');
        title(opt.AUX.group_names_prefix{g});
        hold off
    end
    saveas(f,[opt.folders.results,'/ER_',opt.folders.prject_name,'_ROIs_SSmask_intersection.fig'],'fig')
    %---clean empty rois
    if sum(err_indx)~=0
        for l=1:length(err_indx)
            if err_indx(l,1) > 0
                str_ = c_b(opt.AUX.roi_names(l,:));
                str_ = str_(1:end-4);
                warning('ATTENTION! ROI %s in %d subject(s) is empty after intersection with brain mask. Skipping this ROI.',str_,err_indx(l,1));
                opt.DATA.ROIS(l,:) = [];
                opt.AUX.roi_names(l,:) = [];
            end
        end
    end
end


warning on backtrace
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
return
end

function brainmask_ROI_intersection_RtoR
global opt
cw = pwd;
mkdir(opt.folders.firstlevel,'RtoR');
opt.RtoR.destination_dir = [opt.folders.firstlevel,'/RtoR'];
mkdir(opt.RtoR.destination_dir,opt.RtoR.name);
opt.RtoR.destination_dir = [opt.RtoR.destination_dir,'/',opt.RtoR.name];
cd (opt.RtoR.destination_dir)
[~,~] = system('rm _bmi*.nii *.txt');
total = size(opt.RtoR.ROIS,1);
reverseStr  = '';
err_indx = zeros(size(opt.RtoR.ROIS,1),1);
war_indx = err_indx;
v_count = zeros(size(opt.RtoR.ROIS,1),2);
for l=1:size(opt.RtoR.ROIS,1)   %ciclo rois
    roi_str = c_b(opt.RtoR.roi_names(l,:));
    %roi_str = roi_str(1:end-4);
    [status,string] =system(['3dcalc -a ',opt.RtoR.ROIS(l,:),' -b ',opt.AUX.bmask_path{1,1},' -prefix _bmi_',roi_str,'.nii -expr ''a*b''']);control(status,string);
    bef = spm_read_vols(spm_vol(c_b(opt.RtoR.ROIS(l,:))));
    aft = spm_read_vols(spm_vol(['_bmi_',roi_str,'.nii']));
    v_bef = int64(nansum(bef(:)));   %can be NaN if produced by SPM
    v_aft = int64(nansum(aft(:)));
    if v_aft == 0
        err_indx(l,1) = 1;
    elseif v_bef ~= v_aft
        war_indx(l,1) = 1;
        v_count(l,1) = v_bef;
        v_count(l,2) = v_aft;
    end
    percentDone = 100 * l / total;
    msg = sprintf('Percent done: %3.1f (%s)', percentDone,roi_str);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n');
warning off backtrace
if sum(err_indx) == total
    opt.AUX.TODO_RtoR = 0;
    opt.AUX.steps_total = opt.AUX.steps_total -1;
    warning('ATTENTION! No ROIs for ROI to ROI connectivity survived after brain mask intersection. Skipping ROI to ROI connectivity calculation.');
elseif sum(err_indx)~=0
    for l=1:length(err_indx)
        if err_indx(l,1)
            str_ = c_b(opt.RtoR.roi_names(l,:));
            %str_ = str_(1:end-4);
            warning('ATTENTION! ROI %s is empty after intersection with brain mask. Skipping this ROI.',str_);
            [status,string] = system(['rm _bmi_',str_,'.nii']);control(status,string);
        end
    end
end
if sum(war_indx)~=0
%     fprintf('\n');
    for l=1:length(war_indx)
        if war_indx(l,1)
            str_ = c_b(opt.RtoR.roi_names(l,:));
            %str_ = str_(1:end-4);
            warning('ROI %s has been reshaped due to intersection with brain mask (from %.0d to %.0d voxels).',str_,v_count(l,1) ,v_count(l,2))
        end
    end
end
warning on backtrace
[R_bmi,list_bmi]= conn_dir('_bmi_*.nii');
opt.RtoR.ROIS = R_bmi;
opt.RtoR.n_rois = length(list_bmi);
roi_names = [];
for k = 1:opt.RtoR.n_rois
    roi_names = strvcat(roi_names,list_bmi(k).name(6:end-4)); 
end
opt.RtoR.roi_names = roi_names;  % potrebbe essere ridondante, ma per sicurezza riassocio i nomi a ciascuna roi

%Now create parcellation_brain
if opt.RtoR.n_ovelapping == 0       %no everlapping voxels found before bmi.
    out_name = ['masks_',opt.RtoR.name,'.nii'];
    for l = 1:opt.RtoR.n_rois
        if l == 1
            hdr_o = spm_vol(list_bmi(l).name);
            mask1 = spm_read_vols(hdr_o);
            mask_sum = mask1;
            fid_rtor = fopen([out_name(1:end-3),'txt'],'w');
            fprintf(fid_rtor,'%s\n',c_b(opt.RtoR.roi_names(l,:)));
            continue
        end
        mask = spm_read_vols(spm_vol(list_bmi(l).name));
        mask = mask*l;
        mask_sum = mask_sum + mask;
        fprintf(fid_rtor,'%s\n',c_b(opt.RtoR.roi_names(l,:)));
    end
    hdr_o.pinfo(1) = 1; hdr_o.pinfo(2) = 0;
    hdr_o.fname = out_name;
    hdr_o.private.dat.fname = out_name;
    spm_write_vol(hdr_o,mask_sum);
    fclose(fid_rtor);
else
    warning off backtrace
    fprintf('\nChecking for overlapping ROIs after brain mask intersection...\n');
    out_name = ['RtoR_',opt.RtoR.name,'_masks.nii'];
    count = 0;
    for l = 1:opt.RtoR.n_rois
        if l == 1
            hdr_o = spm_vol(list_bmi(l).name);
            mask1 = spm_read_vols(hdr_o);
            tmp = (opt.RtoR.ovelapping == 1 & mask1 == 1);
            check = length(tmp(tmp > 0));
            if check > 0
                count = count +1;
                str_ = c_b(opt.RtoR.roi_names(l,:));
                %str_ = str_(1:end-4);
                v_count_before = sum(mask1(:));
                v_count_after = v_count_before -check;
                if v_count_after == 0
                    error('You should not use overlapping ROIs for ROI to ROI connectivity. Trying to remove overlapping voxels reults in an empty ROI.');
                else
                    warning('ROI %s has %.0d overlapping voxels. These voxel will be removed (from %.0d to %.0d voxels).',str_,check ,v_count_before, v_count_after);
                end
                mask1(opt.RtoR.ovelapping == 1) = 0;
            end
            mask_sum = mask1;
            fid_rtor = fopen([out_name(1:end-3),'txt'],'w');
            fprintf(fid_rtor,'%s\n',c_b(opt.RtoR.roi_names(l,:)));
            continue
        end
        mask = spm_read_vols(spm_vol(list_bmi(l).name));
        tmp = (opt.RtoR.ovelapping == 1 & mask == 1);
        check = length(tmp(tmp > 0));
        if check > 0
            count = count +1;
            str_ = c_b(opt.RtoR.roi_names(l,:));
            %str_ = str_(1:end-4);
            v_count_before = sum(mask(:));
            v_count_after = v_count_before -check;
            if v_count_after == 0
                error('You should not use overlapping ROIs for ROI to ROI connectivity. Trying to remove overlapping voxels reults in an empty ROI.');
            else
                warning('ROI %s has %.0d overlapping voxels. These voxel will be removed (from %.0d to %.0d voxels).',str_,check ,v_count_before, v_count_after);
            end
            mask(opt.RtoR.ovelapping == 1) = 0;
        end
        mask = mask*l;
        mask_sum = mask_sum + mask;
        fprintf(fid_rtor,'%s\n',c_b(opt.RtoR.roi_names(l,:)));
    end
    if count == 0
        fprintf('no problem found!\n');
    end
    hdr_o.pinfo(1) = 1; hdr_o.pinfo(2) = 0;
    hdr_o.fname = out_name;
    hdr_o.private.dat.fname = out_name;
    spm_write_vol(hdr_o,mask_sum);
    fclose(fid_rtor);
    opt.RtoR.ovelapping = [];       %save space
    warning on backtrace
end
opt.AUX.RtoR.masks_path = [opt.RtoR.destination_dir,'/',out_name];
cd (cw);
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
return
end

function brainmask_ROI_intersection_RtoR_mf
global opt
cw = pwd;
mkdir(opt.folders.firstlevel,'RtoR');
opt.RtoR.destination_dir = [opt.folders.firstlevel,'/RtoR'];
mkdir(opt.RtoR.destination_dir,opt.RtoR.name);
opt.RtoR.destination_dir = [opt.RtoR.destination_dir,'/',opt.RtoR.name];
cd (opt.RtoR.destination_dir)
[~,~] = system('rm RtoR*');   
output_name = ['RtoR_',opt.RtoR.name,'_masks.nii'];
[status,string] =system(['3dcalc -a ',opt.RtoR.roi_folder_or_file,' -b ',opt.AUX.bmask_path{1,1},' -prefix ',output_name,' -expr ''a*b''']);control(status,string);  % no need for bm selection. This modality works only in no SS case
bef = spm_read_vols(spm_vol(opt.RtoR.roi_folder_or_file));
hdr_o = spm_vol(output_name);
aft = spm_read_vols(hdr_o);
n_roi_aft_bmi = length(unique(aft(:))) - 1;
checker = opt.RtoR.mf_max_roi - n_roi_aft_bmi;
warning off backtrace
if checker == 0    % no Rois has been deleted, return
    for l = 1:opt.RtoR.mf_max_roi
        v_b = length(bef(bef == l));
        v_a = length(aft(aft == l));
        if v_b ~= v_a
            warning('ROI %s has been reshaped due to intersection with brain mask (from %.0d to %.0d voxels).',remove_endofline(opt.RtoR.roi_names(l,:)),v_b ,v_a)
        end
    end
    if opt.RtoR.mf_txt_exist
        [status,string] =system(['cp ',opt.RtoR.mf_txt_file,' ./',output_name(1:end-3),'txt']);control(status,string);
    end
    opt.RtoR.n_rois = opt.RtoR.mf_max_roi;
    opt.AUX.RtoR.masks_path = [opt.RtoR.destination_dir,'/',output_name];
    cd (cw);
    save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
    return
end
%------now handle the more dirty case of removed rois------------------
err_indx = zeros(opt.RtoR.mf_max_roi,1);
war_indx = err_indx;
v_count = zeros(opt.RtoR.mf_max_roi,2);
for l = 1:opt.RtoR.mf_max_roi
    v_b = length(bef(bef == l));
    v_a = length(aft(aft == l));
    if v_a == 0
        err_indx(l,1) = 1;
    elseif v_b ~= v_a
        war_indx(l,1) = 1;
        v_count(l,1) = v_b;
        v_count(l,2) = v_a;
    end
end
if sum(war_indx)~=0
    for l=1:length(war_indx)
        if war_indx(l,1)
            str_ = c_b(opt.RtoR.roi_names(l,:));
            warning('ROI %s has been reshaped due to intersection with brain mask (from %.0d to %.0d voxels).',remove_endofline(str_),v_count(l,1) ,v_count(l,2))
        end
    end
end
if sum(err_indx) == opt.RtoR.mf_max_roi
    opt.AUX.TODO_RtoR = 0;
    opt.AUX.steps_total = opt.AUX.steps_total -1;
    warning('ATTENTION! No ROIs survived after brain mask intersection. Skipping ROI to ROI calculation.');
elseif sum(err_indx)~=0
    for l=1:length(err_indx)
        if err_indx(l,1)
            str_ = c_b(opt.RtoR.roi_names(l,:));
            warning('ATTENTION! ROI %s is empty after intersection with brain mask. Skipping this ROI.',remove_endofline(str_));
        end
    end
end
warning on backtrace
fprintf('Re-writing mask file to account for ROIs lost...');
% now that the user is aware of the situation let's fix the mask
count = 0;
size_mask = size(aft);
new_mask = zeros(size_mask(1),size_mask(2),size_mask(3));
% fix names----------
if opt.RtoR.mf_txt_exist
    fid = fopen([output_name(1:end-3),'txt'],'w');
    opt.RtoR.roi_names(logical(err_indx),:) = [];
    for l = 1:size(opt.RtoR.roi_names,1)
        fprintf(fid,'%s',c_b(opt.RtoR.roi_names(l,:)));
    end
else
    indexs = 1:1:n_roi_aft_bmi;
    indexs = indexs';
    opt.RtoR.roi_names = num2str(indexs);
end
% --------------------
% fix rois
for l = 1:opt.RtoR.mf_max_roi
    if err_indx(l) == 0
        count = count +1;
        new_mask(aft == l) = count;
    end
end
hdr_o.pinfo(1) = 1; hdr_o.pinfo(2) = 0;
hdr_o.fname = output_name;
hdr_o.private.dat.fname = output_name;
spm_write_vol(hdr_o,new_mask);        
opt.RtoR.n_rois = n_roi_aft_bmi;
opt.AUX.RtoR.masks_path = [opt.RtoR.destination_dir,'/',output_name];
cd (cw);
save([opt.folders.results,'/ER_',opt.folders.prject_name,'.mat'],'opt');
fprintf('done!\n');
return
end

function do_3dRSFC(g,k,j) %g group, k subj, j session
global opt

% if opt.aCompCor.do
%     rp_str = [' -ort ',opt.aCompCor.X1D{g,j,k},opt.AUX.v_s_rp_str{g,j}]; %it already includes rp if they are present
% elseif opt.rp_regression.do == 1 && ~opt.aCompCor.do
%     rp_str = [' -ort ',c_b(opt.DATA.RP_1D{g,j}(k,:)),opt.AUX.v_s_rp_str{g,j}];
% else
%     rp_str = '';
% end

if opt.aCompCor.do || opt.rp_regression.do
    %rp_str = [' -ort ',opt.X1D{g,j,k},opt.AUX.v_s_rp_str{g,j}]; 
    rp_str = [' -ort ',opt.X1D{g,j,k}]; 
else
    rp_str = '';
end

if opt.AUX.TODO_RSFA
    rsfa_str = '';
else
    rsfa_str = ' -no_rsfa';
end

bmask_path = brainmask_selector(g,k,1);
bmask_string = brainmask_selector(g,k,2);
if opt.psc
    [status,string] = system(['3dTcat ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),' -rlt+ -prefix _detrended.nii']);control(status,string);
    [status,string] = system(['3dTstat -mask ',bmask_path,' -mean -prefix _mean.nii _detrended.nii']);control(status,string);
    
    media = spm_read_vols(spm_vol('_mean.nii'));
    media_ = mean(media(opt.AUX.bmask_tmp == 1));                % media in bmask
    v_e = media(opt.AUX.bmask_tmp == 1 & media < (media_*0.05)); %considero voxel errati quelli in cui la media temporale ?? < 5% della media su tutta la maschera
    if ~isempty(v_e)
        v_e = length(v_e);
        v_b = sum(opt.AUX.bmask_tmp(:));
        percent_voxel = v_e/v_b * 100;
        if percent_voxel >= 50
            error('You cannot transform data in percent signal change (psc) if data have zero-mean.');
        elseif percent_voxel >= 5
            error('Your group brain mask is beadly defined. Try another method for brain mask definition.');
        end
        opt.AUX.bmask_vxl_err = [opt.AUX.bmask_vxl_err;v_e];
        opt.AUX.bmask_tmp(opt.AUX.bmask_tmp == 1 & media < (media_*0.05)) = 0;
        [status,string] = system(['rm ',opt.AUX.bmask_path]);control(status,string);
        cd ([opt.folders.results,'/brainmask_m',num2str(opt.bm_method)]);
        spm_write_vol(opt.AUX.bmask_tmp_hdr,opt.AUX.bmask_tmp);
        cd (opt.folders.preprocessing)
        [status,string] = system('rm _mean.nii');control(status,string);
        [status,string] = system(['3dTstat -mask ',opt.AUX.bmask_path,' -mean -prefix _mean.nii _detrended.nii']);control(status,string);
    end
    
    [status,string] = system(['3dcalc -a ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),' -b _mean.nii -expr ''step(-(100*a/b)+200)*(100*a/b)'' -prefix ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string);      %also trucate values avoce 100%
    [status,string] = system('rm _detrended.nii');control(status,string);
    [status,string] = system('rm _mean.nii');control(status,string);
    [~,string] = system(['3dinfo -tr ',opt.AUX.output_name{g,j,k},'_psc.nii']);
    if str2num(string) ~= opt.tr{g,j}
        [status,string] = system(['3drefit -TR ',num2str(opt.tr{g,j}),' ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string);  
    end
    [status,string] = system(['3dRSFC -dt ',num2str(opt.tr{g,j}),rp_str,opt.tredrsfc_extraoption,rsfa_str,' ',bmask_string,' -prefix ',opt.AUX.output_name{g,j,k},' ',num2str(opt.filter_band(1)),' ',num2str(opt.filter_band(2)),' ',opt.AUX.output_name{g,j,k},'_psc.nii',opt.AUX.v_s_str{g,j}]);control(status,string);
    if ~opt.AUX.psc_save
        [status,string] = system(['rm ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string); 
    end
else
    [~,string] = system(['3dinfo -tr ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))]);
    if str2num(string) ~= opt.tr{g,j}
        [status,string] = system(['3drefit -TR ',num2str(opt.tr{g,j}),' ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))]);control(status,string);
    end
    [status,string] = system(['3dRSFC -dt ',num2str(opt.tr{g,j}),rp_str,opt.tredrsfc_extraoption,rsfa_str,' ',bmask_string,' -prefix ',opt.AUX.output_name{g,j,k},' ',num2str(opt.filter_band(1)),' ',num2str(opt.filter_band(2)),' ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),opt.AUX.v_s_str{g,j}]);control(status,string);
end
return
end

function do_3dTproject(g,k,j) %g group, k subj, j session
global opt

% if opt.aCompCor.do
%     rp_str = [' -ort ',opt.aCompCor.X1D{g,j,k},opt.AUX.v_s_rp_str{g,j}]; %it already includes rp if they are present
% elseif opt.rp_regression.do == 1 && ~opt.aCompCor.do
%     rp_str = [' -ort ',c_b(opt.DATA.RP_1D{g,j}(k,:)),opt.AUX.v_s_rp_str{g,j}];
% else
%     rp_str = '';
% end

if opt.aCompCor.do || opt.rp_regression.do
    %rp_str = [' -ort ',opt.X1D{g,j,k},opt.AUX.v_s_rp_str{g,j}]; 
    rp_str = [' -ort ',opt.X1D{g,j,k}]; 
else
    rp_str = '';
end

bmask_path = brainmask_selector(g,k,1);
bmask_string = brainmask_selector(g,k,2);
if opt.psc
    [status,string] = system(['3dTcat ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),' -rlt+ -prefix _detrended.nii']);control(status,string);
    [status,string] = system(['3dTstat -mask ',bmask_path,' -mean -prefix _mean.nii _detrended.nii']);control(status,string);
    
    media = spm_read_vols(spm_vol('_mean.nii'));
    media_ = mean(media(opt.AUX.bmask_tmp == 1));                % media in bmask
    v_e = media(opt.AUX.bmask_tmp == 1 & media < (media_*0.05)); %considero voxel errati quelli in cui la media temporale ?? < 5% della media su tutta la maschera
    if ~isempty(v_e)
        v_e = length(v_e);
        v_b = sum(opt.AUX.bmask_tmp(:));
        percent_voxel = v_e/v_b * 100;
        if percent_voxel >= 50
            error('You cannot transform data in percent signal change (psc) if data have zero-mean.');
        elseif percent_voxel >= 5
            error('Your group brain mask is beadly defined. Try another method for brain mask definition.');
        end
        opt.AUX.bmask_vxl_err = [opt.AUX.bmask_vxl_err;v_e];
        opt.AUX.bmask_tmp(opt.AUX.bmask_tmp == 1 & media < (media_*0.05)) = 0;
        [status,string] = system(['rm ',opt.AUX.bmask_path]);control(status,string);
        cd ([opt.folders.results,'/brainmask_m',num2str(opt.bm_method)]);
        spm_write_vol(opt.AUX.bmask_tmp_hdr,opt.AUX.bmask_tmp);
        cd (opt.folders.preprocessing)
        [status,string] = system('rm _mean.nii');control(status,string);
        [status,string] = system(['3dTstat -mask ',opt.AUX.bmask_path,' -mean -prefix _mean.nii _detrended.nii']);control(status,string);
    end
    
    [status,string] = system(['3dcalc -a ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),' -b _mean.nii -expr ''step(-(100*a/b)+200)*(100*a/b)'' -prefix ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string);      %also trucate values avoce 100%
    [status,string] = system('rm _detrended.nii');control(status,string);
    [status,string] = system('rm _mean.nii');control(status,string);
    [~,string] = system(['3dinfo -tr ',opt.AUX.output_name{g,j,k},'_psc.nii']);
    if str2num(string) ~= opt.tr{g,j}
        [status,string] = system(['3drefit -TR ',num2str(opt.tr{g,j}),' ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string);  
    end
    [status,string] = system(['3dTproject -input ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),'_psc.nii',opt.AUX.v_s_str{g,j},' -prefix ',opt.AUX.output_name{g,j,k},'_LFF.nii -passband ',num2str(opt.filter_band(1)),' ',num2str(opt.filter_band(2)),' -dt ',num2str(opt.tr{g,j}),rp_str,opt.tredrsfc_extraoption,' ',bmask_string,' -polort 2 -censor ',opt.censoring.censor1D{g,j,k},' -cenmode ',opt.censoring.mode]);control(status,string);
    opt.censoring.DOF{g,j}(k,1) = parse_3dTproject(string);
    if ~opt.AUX.psc_save
        [status,string] = system(['rm ',opt.AUX.output_name{g,j,k},'_psc.nii']);control(status,string); 
    end
else
    [~,string] = system(['3dinfo -tr ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))]);
    if str2num(string) ~= opt.tr{g,j}
        [status,string] = system(['3drefit -TR ',num2str(opt.tr{g,j}),' ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))]);control(status,string);
    end
    [status,string] = system(['3dTproject -input ',c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)),opt.AUX.v_s_str{g,j},' -prefix ',opt.AUX.output_name{g,j,k},'_LFF.nii -passband ',num2str(opt.filter_band(1)),' ',num2str(opt.filter_band(2)),' -dt ',num2str(opt.tr{g,j}),rp_str,opt.tredrsfc_extraoption,' ',bmask_string,' -polort 2 -censor ',opt.censoring.censor1D{g,j,k},' -cenmode ',opt.censoring.mode]);control(status,string);
    opt.censoring.DOF{g,j}(k,1) = parse_3dTproject(string);
end
return
end

function [dof] = parse_3dTproject(str)
indx_up = strfind(str,'==>');
indx_down = strfind(str, 'D.O.F.');
dof = str2num(str((indx_up+3):(indx_down-1)));
return
end

function psd_initialization
global opt 
opt.AUX.psd_nfft = 2^nextpow2(max([opt.N{:}]));
opt.AUX.psd_N = max([opt.N{:}]);
warning('For PSD calculation the NFFT is choosen based on the maximum number of volumes among your data set. This allows for PDS estimates to have the same number of spectral points. In this way you can average/compare them. However, you should remember that for N < max(N) the PSD will be smoothed.')
return
end

function gFC_initialization(ot)
global opt 
%set default
if ~isfield(ot,'gFC')
    %ot.gFC.thr = 0.25;
    ot.gFC.Vthr = [0.25 0.5 0.05];
else
%     if ~isfield(ot.gFC,'thr');
%         ot.gFC.thr = 0.25;
%     else
%         if numel(ot.gFC.thr) > 1
%             error('Too many input in ER.others.gFC.thr . Specify only a number');
%         end
%         if ot.gFC.thr >= 1
%             error('ER.others.gFC.thr must be < 1 [default is 0.25].');
%         end
%     end
    if ~isfield(ot.gFC,'Vthr'); 
        ot.gFC.Vthr = [0.25 0.5 0.05];
    else
        if numel(ot.gFC.Vthr) ~= 3
            error('Bad ER.others.gFC.Vthr spefication. Specify 3 numbers: t0 t1 and dt (default is [0.25 0.5 0.05]).');
        end
        if sum(ot.gFC.Vthr <= 1) ~= 3 
            error('ER.others.gFC.Vthr must be < 1 (default is [0.25 0.5 0.05]).');
        elseif ot.gFC.Vthr(1) > ot.gFC.Vthr(2)
            error('Bad definition of ER.others.gFC.Vthr (default is [0.25 0.5 0.05]).');
        end
    end
end
ot.gFC.hist = 1;
% opt.AUX.gFC.thr = ot.gFC.thr;
opt.AUX.gFC.Vthr = ot.gFC.Vthr;
opt.AUX.gFC.hist = ot.gFC.hist;
%opt.AUX.gFC.num_thr = floor((opt.AUX.gFC.Vthr(2) -
%opt.AUX.gFC.Vthr(1))./opt.AUX.gFC.Vthr(3) + 1); eliminato perch?? non
%sempre produce la stesso valore di afni, ?? un problema di gestione degli
%interi differente tra c e matlab. risolto prendendo il valore dirattamente
%dal numero di volumi
return
end

function DeC_initialization(ot)
global opt 
%set default
if ~isfield(ot,'DeC')
    ot.DeC.Thrs = [0.25 0.5 0.05];
else
    if ~isfield(ot.DeC,'Thrs') 
        ot.DeC.Thrs = [0.25 0.5 0.05];
    else
        if numel(ot.DeC.Thrs) ~= 3
            error('Bad ER.others.DeC.Thrs spefication. Specify 3 numbers: t0 t1 and dt (default is [0.25 0.5 0.05]).');
        end
        if sum(ot.DeC.Thrs <= 1) ~= 3 
            error('ER.others.DeC.Thrs must be < 1 (default is [0.25 0.5 0.05]).');
        elseif ot.DeC.Thrs(1) > ot.DeC.Thrs(2)
            error('Bad definition of ER.others.DeC.Thrs (default is [0.25 0.5 0.05]).');
        end
    end
end
opt.AUX.DeC.Thrs = ot.DeC.Thrs;
opt.AUX.DeC.num_thr = floor((opt.AUX.DeC.Thrs(2) - opt.AUX.DeC.Thrs(1))./opt.AUX.DeC.Thrs(3) + 1);
N = 200;        %IMPORTANT: change this accordingly to gFC HIST output.
D = 2/N;
rr = -1:D:(1-D); % building correlation index. 
% Looking for quality on float is problematic... 
% Results wouldn't be really equally spaced! Scaling data to int
rr= uint8(rr*100);
Thrs_vector = nan(opt.AUX.DeC.num_thr,1);
Thrs_afni_indx = Thrs_vector;
for l = 1:opt.AUX.DeC.num_thr
    Thrs_vector(l) = opt.AUX.DeC.Thrs(1) + (l-1)*opt.AUX.DeC.Thrs(3);
    tmp = find(rr >= uint8((Thrs_vector(l))*100));       
    Thrs_afni_indx(l) = tmp(1)-1;
end
opt.AUX.DeC.Thrs_vector = Thrs_vector;
opt.AUX.DeC.Thrs_afni_indx = Thrs_afni_indx;
return
end

function zscore_conversion(old_string,postfixZ,postfixM)
% zscore (and m*-1) computation
global opt
old_string = [old_string,'.nii'];
%list_file = dir(['*',old_string]);
if nargin == 2
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)
                bmask = uint8(brainmask_selector(g,k,3));
                base_name = opt.AUX.output_name{g,j,k};
                filename = [base_name,old_string];
                if j == 1 && g==1 && k==1; hdr = spm_vol(filename);hdr = rmfield(hdr,'pinfo'); if length(hdr) > 1; hdr(2:end) = []; hdr.private.dat.dim(end) = []; end;end;
                volume = spm_read_vols(spm_vol(filename));
                if length(size(volume)) > 3; volume = volume(:,:,:,1); end
                media = mean(volume(bmask==1));
                sd = std(volume(bmask == 1));
                zvolume = (volume - media)./sd;
                zvolume(bmask~=1) = 0;
                %save
                output_nameZ = [base_name,postfixZ,'.nii'];
                hdr.fname = output_nameZ;
                hdr.private.dat.fname = output_nameZ;
                spm_write_vol(hdr,zvolume);
            end
        end
    end
elseif nargin == 3      %also zero mean computation
    for g =1:opt.group_number
        for j=1:opt.session_number
            for k=1:opt.subject_number(g)  
                bmask = uint8(brainmask_selector(g,k,3));
                base_name = opt.AUX.output_name{g,j,k};
                filename = [base_name,old_string];
                if j == 1 && g==1 && k==1; hdr = spm_vol(filename);hdr = rmfield(hdr,'pinfo'); if length(hdr) > 1; hdr(2:end) = []; hdr.private.dat.dim(end) = []; end;end;
                volume = spm_read_vols(spm_vol(filename));
                if length(size(volume)) > 3; volume = volume(:,:,:,1); end
                media = mean(volume(bmask==1));
                sd = std(volume(bmask == 1));
                zvolume = (volume - media)./sd;
                mvolume = (volume./media)-double(bmask);
                zvolume(bmask~=1) = 0;
                mvolume(bmask~=1) = 0;
                %save
                output_nameZ = [base_name,postfixZ,'.nii'];
                hdr.fname = output_nameZ;
                hdr.private.dat.fname = output_nameZ;
                spm_write_vol(hdr,zvolume);
                output_nameM = [base_name,postfixM,'.nii'];
                hdr.fname = output_nameM;
                hdr.private.dat.fname = output_nameM;
                spm_write_vol(hdr,mvolume);
            end
        end
    end
end
return
end

function zscore_conversion_list(old_string,postfixZ,bmask,postfixM)
% zscore (and m*-1) computation
old_string = [old_string,'.nii'];
list_file = dir(['*',old_string]);
if nargin == 3
    for l = 1:length(list_file)
        if l == 1; hdr = spm_vol(list_file(l).name);hdr = rmfield(hdr,'pinfo');end
        volume = spm_read_vols(spm_vol(list_file(l).name));
        media = mean(volume(bmask==1));
        sd = std(volume(bmask == 1));
        zvolume = (volume - media)./sd;
        zvolume(bmask~=1) = 0;
        %save
        output_nameZ = [list_file(l).name(1:(end-length(old_string))),postfixZ,'.nii'];
        hdr.fname = output_nameZ;
        hdr.private.dat.fname = output_nameZ;
        spm_write_vol(hdr,zvolume);
    end
elseif nargin == 4      %also zero mean computation
    for l = 1:length(list_file)
        if l == 1; hdr = spm_vol(list_file(l).name);hdr = rmfield(hdr,'pinfo');end
        volume = spm_read_vols(spm_vol(list_file(l).name));
        media = mean(volume(bmask==1));
        sd = std(volume(bmask == 1));
        zvolume = (volume - media)./sd;
        mvolume = (volume./media)-double(bmask);
        zvolume(bmask~=1) = 0;
        mvolume(bmask~=1) = 0;
        %save
        output_nameZ = [list_file(l).name(1:(end-length(old_string))),postfixZ,'.nii'];
        hdr.fname = output_nameZ;
        hdr.private.dat.fname = output_nameZ;
        spm_write_vol(hdr,zvolume);
        output_nameM = [list_file(l).name(1:(end-length(old_string))),postfixM,'.nii'];
        hdr.fname = output_nameM;
        hdr.private.dat.fname = output_nameM;
        spm_write_vol(hdr,mvolume);
    end
end
return
end

function smooth_map(file_wildcard,g,k)
global opt
list_file = dir(['*',file_wildcard]);
for l = 1:length(list_file)
    [status,string] = system(['3dBlurInMask -input ',list_file(l).name,' -float -FWHM ',num2str(opt.AUX.special_smoothing_FWHM),' -prefix ',list_file(l).name,' ',brainmask_selector(g,k,2),' -overwrite -preserve']);control(status,string);
end           
return
end

function smooth_map_ALFF(file_wildcard,g,k)
%the only difference is the FWHM, not as special option but equal to the
%FWHM applied to the timeseries
global opt
list_file = dir(['*',file_wildcard]);
for l = 1:length(list_file)
    [status,string] = system(['3dBlurInMask -input ',list_file(l).name,' -float -FWHM ',num2str(opt.FWHM),' -prefix ',list_file(l).name,' ',brainmask_selector(g,k,2),' -overwrite -preserve']);control(status,string);
end           
return
end

function spectral_entropy(input)
output_name = [input(1:end-7),'SpEn.nii'];
[status,string] = system(['3dcalc -a ',input,' -prefix _sqrt.nii -expr ''sqrt(a)''']);control(status,string);
[status,string] = system('3dTstat -sum -prefix _sum.nii _sqrt.nii');control(status,string);
[status,string] = system('3dcalc -a _sqrt.nii -b _sum.nii -prefix _p.nii -expr ''a/b*ispositive(b-1)''');control(status,string);
[status,string] = system('3dcalc -a _p.nii -prefix _en.nii -expr ''(-1)*a*log(a)/log(2)''');control(status,string);
[status,string] = system(['3dTstat -sum -prefix ',output_name,' _en.nii']);control(status,string);
[status,string] = system('rm _sqrt.nii _sum.nii _en.nii _p.nii');control(status,string);
return
end

function correlation_entropy(input)
output_name = [input(1:end-12),'CoEn.nii'];
[status,string] = system(['3dcalc -a ',input,' -b _sum.nii -prefix _p.nii -expr ''a/b*ispositive(b-1)''']);control(status,string);
[status,string] = system('3dcalc -a _p.nii -prefix _en.nii -expr ''(-1)*a*log(a)/log(2)''');control(status,string);
[status,string] = system(['3dTstat -sum -prefix ',output_name,' _en.nii']);control(status,string);
[status,string] = system('rm _en.nii _p.nii');control(status,string);
return
end

function fastconn(epi,hdr_output,ss,tr,filter,nreg,output_name,group,subj,sess)
global opt
s = size(epi);

%OLD METHOD (memory consuming and slower). New method is around 75% faster
% epi_2d = reshape(epi,[s(1)*s(2)*s(3),s(4)]);
% ss_2d = permute(repmat(ss',[s(1) 1 s(2) s(3)]),[1 3 4 2]);
% ss_2d = reshape(ss_2d,[s(1)*s(2)*s(3),s(4)]);
% corr_1d = correlation(epi_2d',ss_2d');
%-------------------------------------------

%------------NEW METHOD---------------------
%it's simply a GLM using normalized data
epi_2d = reshape(epi,[s(1)*s(2)*s(3),s(4)]);
epi_2d = zscore(epi_2d');
ss = zscore(ss);
corr_1d = ss\epi_2d;
%-------------------------------------------

corr_3d = reshape(corr_1d,[s(1),s(2),s(3)]);

if strcmp(opt.prepro_mode,'3dTproject') 
    DOF = opt.censoring.DOF{group,sess}(subj,1);
else
    DOF = estimate_dof(s(4),nreg,tr,filter(1),filter(2));
    %DOF=max(0,s(4)*(min(1/(2*tr),filter(2))-max(0,filter(1)))/(1/(2*tr))-nreg);
end

z = atanh(corr_3d);

p=spm_Ncdf(z*sqrt(max(0,DOF-3)));
p=2*min(p,1-p);

%%%%%%%%%%%%%%%%
% save map .nii
%%%%%%%%%%%%%%%%
%TODO: add selector to avoid print all these maps
hdr_output.fname = [output_name,'_r.nii'];
hdr_output.private.dat.fname = [output_name,'_r.nii'];
spm_write_vol(hdr_output,corr_3d);

hdr_output.fname = [output_name,'_z.nii'];
hdr_output.private.dat.fname = [output_name,'_z.nii'];
spm_write_vol(hdr_output,z);

if strcmp(opt.prepro_mode,'3dTproject') 
    hdr_output.fname = [output_name,'_Z.nii'];
    hdr_output.private.dat.fname = [output_name,'_Z.nii'];
    spm_write_vol(hdr_output,z*sqrt(DOF-3));
end

hdr_output.fname = [output_name,'_p.nii'];
hdr_output.private.dat.fname = [output_name,'_p.nii'];
spm_write_vol(hdr_output,p);

return
end

function [C]=correlation(A,B)
%OLD FUNCTION used by fastconn
%remember: A and B must be a n-x-m matrix in which n are time points and m are voxels
An=bsxfun(@minus,A,mean(A,1));
Bn=bsxfun(@minus,B,mean(B,1));
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
C=sum(An.*Bn,1);
return
end

function do_3ddelay(epi,seed,output_name)
global opt
[status,string] = system(['3dROIstats -mask ',seed,' -quiet ',epi,' > _roi.1D']);control(status,string)
[status,string] = system(['3ddelay -input ',epi,' -ideal_file _roi.1D -fs ',num2str(opt.AUX.fs),' -T 0 -prefix ',output_name,'.nii ',opt.AUX.bmask_string,' -nodtrnd -uS']);control(status,string)
[status,string] = system('rm _roi.1D');control(status,string)
return
end

function coherency_initialization(ot)
% coherency prepare parameters and check options
global opt
if exist('cpsd_R2014b')~=2
    error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''./ER''))');
end
% if isempty(opt.N)     YOU CAN DELATE IT!
%     [~,N] = system(['3dinfo -nv ',c_b(opt.DATA.FUNCTIONALS{1,1}(1,:))]);        %attention: assuming that all scans have the same number of time points
%     opt.N = str2num(N);
% end

opt.coherency.bands = [opt.filter_band];
opt.coherency.HW_length = [];
if isfield(ot,'coh')
    if isfield(ot.coh,'band_division') && ~isempty(ot.coh.band_division) 
        if size(ot.coh.band_division,2) ~= 2
            error('Bad definition of coh.band_division. A nx2 matrix is needed.');
        end
        c_diff = ot.coh.band_division(:,2) -ot.coh.band_division(:,1);
        cc = c_diff(c_diff < 0);
        if ~isempty(cc)
            error('Bad definition of coh.band_division. Top limit is lower than the bottom one.');
        end 
        bot = opt.filter_band(1);
        top = opt.filter_band(2);
        c_bot = ot.coh.band_division(:,1) < bot | ot.coh.band_division(:,1) > top;
        c_top = ot.coh.band_division(:,2) > top | ot.coh.band_division(:,2) < bot;
        cc = sum(c_bot + c_top);
        if sum(cc) > 0
            error('Bad definition of coh.band_division. The division must be in the range of the filter band.');
        end         
        opt.coherency.bands = [opt.coherency.bands;ot.coh.band_division];
    end
    if isfield(ot.coh,'window_length') && ~isempty(ot.coh.window_length) 
       opt.coherency.HW_length = ot.coh.window_length;
       if opt.coherency.HW_length > opt.N % TO FIX now N is a cell
           error('Window for Welch''s method segmentation can''t be longer than epi time serie.');
       end
    end
end

if isempty(opt.coherency.HW_length)
    opt.coherency.HW_length = fix(1/(opt.tr*(opt.filter_band(1)))); %MUST BE FIXED, SINCE NOW TR IS A CELL    
    if opt.coherency.HW_length > opt.N % TO FIX now N is a cell
        error('The default value of window''s length for Welch''s method segmentation is too longer. Please define it manually via coh.window_length.');
    end
end

opt.coherency.NO = fix(opt.coherency.HW_length/2);      % 50% overlap
opt.coherency.N_seg = (opt.N-opt.coherency.NO)./(opt.coherency.HW_length-opt.coherency.NO); % TO FIX now N is a cell

war_integer = 0; if fix(opt.coherency.N_seg) ~= opt.coherency.N_seg; war_integer = 1; end
opt.coherency.N_seg = fix(opt.coherency.N_seg);

%--NFFT and frequency points
opt.coherency.NFFT=max(256,2^nextpow2(opt.coherency.HW_length)); % Points for FFT (N)
% if rem(opt.coherency.NFFT,2)            %mi assicuro sia pari per non avere rogne con le frequenze
%     opt.coherency.NFFT = opt.coherency.NFFT + 1;
% else
%     opt.coherency.NFFT = opt.coherency.NFFT;
% end
opt.coherency.freq = opt.AUX.fs/2*linspace(0,1,opt.coherency.NFFT/2+1);   
%lets cut the out-filter frequency (only the top)
ind = find(opt.coherency.freq<=opt.filter_band(2));
b = ind(end);
opt.coherency.freq = opt.coherency.freq(1:b);
%--------------------------

fprintf('\n_______________________________________________\n');
fprintf('\tCoherency parameters\n');
fprintf('Segment length: \t%.0d\n',opt.coherency.HW_length);
fprintf('Overlap length: \t%.0d\n',opt.coherency.NO);
fprintf('Averaged segments: \t%.0d\n',opt.coherency.N_seg);
fprintf('NFFT segment: \t\t%.0d\n',opt.coherency.NFFT);
fprintf('Freq. interval: \t0 %0.3f\n',opt.coherency.freq(end))
warning off backtrace
if opt.coherency.N_seg < 3
    warning('The number of averaged segments is low.');
end
if war_integer
    warning('signal:welchparse:MustBeInteger','The number of segments was not an integer, data will be truncated.');
end
warning on backtrace
fprintf('_______________________________________________\n');
return
end

function coherency(epi,hdr_epi,seed,output_name)
% Coherence evaluation
%
% CoherencY is a complex function of frequency. CoherencE is a real
% adimensional number between 0,1 . Delay is a real number (time).
% Sun et al. Neuroimage 2005

global opt

hdr_seed = spm_vol(seed);
seed = spm_read_vols(hdr_seed);

s = size(epi);  % Image Dimensions (Voxels)

reg = nan(s(4),1);

for l=1:s(4) 
    tmp = epi(:,:,:,l);
    reg(l,1) = mean(tmp(seed == 1));  %r stands for regressor
end

reg = zscore(reg');

WS = (1/opt.tr); %MUST BE FIXED, SINCE NOW TR IS A CELL

% Cpsd params


% Evaluates 'auto'cpsd for the regressor
[frr,W]=cpsd_R2014b(reg,reg,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);

FFT_L = length(W);  %acquiring FFT length from cpsd

epi = reshape(epi,[s(1)*s(2)*s(3),s(4)]);

epi = zscore(epi');          

reg_2d = repmat(reg',[1, s(1)*s(2)*s(3)]);

% cpsd of regressor with epi' timecourses
fre = cpsd_R2014b(reg_2d,epi,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);
% 'auto'cpsd of the epi
fee = cpsd_R2014b(epi,epi,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);

% Evaluation of the coherencY
fun = @(A,B) sqrt(A.*B);
frr = repmat(frr,[1, s(1)*s(2)*s(3)]);
normalization = bsxfun(fun,fee,frr);
coherency = fre./normalization;

% reshape
frr = reshape(frr',[s(1),s(2),s(3),FFT_L]);
fre = reshape(fre',[s(1),s(2),s(3),FFT_L]);
fee = reshape(fee',[s(1),s(2),s(3),FFT_L]);
coherency = reshape(coherency',[s(1),s(2),s(3),FFT_L]);

% compute coherence from coherency
coherence = coherency.*conj(coherency);
% compute delay from coherency
phase_spectrum = angle(coherency);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%               WRITE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Density files: (coherence(f) e phase(f) (ie. phase spectrum))
hdr_epi = rmfield(hdr_epi,'pinfo');

%lets cut the out-filter frequency (only the top)
ind = find(W<=opt.filter_band(2));
b = ind(end);

hdr_cohd = hdr_epi(1:b);
hdr_phad = hdr_cohd;
cohDens_output = [output_name,'_MA.nii'];   %MA, MAgnitude
phaDens_output = [output_name,'_PH.nii'];   %PH, Phase spectrum

for l = 1:b
    hdr_cohd(l).fname = cohDens_output;
    hdr_cohd(l).private.dat.fname = cohDens_output;
    spm_write_vol(hdr_cohd(l),coherence(:,:,:,l));
    hdr_phad(l).fname = phaDens_output;
    hdr_phad(l).private.dat.fname = phaDens_output;
    spm_write_vol(hdr_phad(l),phase_spectrum(:,:,:,l));
end
%--------------------------------------------------------

% Band-averaged chorence and delay

bands = opt.coherency.bands;    % first band is the preprocessing filter (i.e, all band)

hdr_ba_coh = hdr_epi(1:size(bands,1));
hdr_ba_del = hdr_ba_coh;
cohBAr_output = [output_name,'_MABA_r.nii'];     %MABA, MAgnitude Band Averanged
cohBAz_output = [output_name,'_MABA_z.nii'];     %fisher transformed
delBA_output = [output_name,'_DEBA.nii'];        %DEBA, Delay Band Averanged


for l = 1:size(bands,1)
    f1 = bands(l,1);
    f2 = bands(l,2);
    ind = find(W>=f1 & W<f2);
    a = ind(1);
    b = ind(end); 
    tau_den = 2*pi*sum(W(a:b)); % Sum of the frequency vector
    BA_coherence = (abs(sum(fre(:,:,:,a:b),4)).^2)./((sum(frr(:,:,:,a:b),4).*sum(fee(:,:,:,a:b),4)));
    BA_delay = sum(angle(fre(:,:,:,a:b)),4)./tau_den;
    hdr_ba_coh(l).fname=cohBAr_output;
    hdr_ba_coh(l).private.dat.fname = cohBAr_output;
    spm_write_vol(hdr_ba_coh(l),BA_coherence);
    hdr_ba_coh(l).fname=cohBAz_output;
    hdr_ba_coh(l).private.dat.fname = cohBAz_output;
    spm_write_vol(hdr_ba_coh(l),atanh(BA_coherence));
    hdr_ba_del(l).fname=delBA_output;
    hdr_ba_del(l).private.dat.fname = delBA_output;
    spm_write_vol(hdr_ba_del(l),BA_delay);
end    
    
return
end

function [r,z,p] = RtoR_coherency_matrix(data)

n = size(data,2);  %number of ROIs
np = size(data,1); %number of time points

r = ones(n);

for i = 1:(n-1)
    for j = 1+i:n
        r(i,j) = RtoR_coherency_calculation(data(:,i),data(:,j));
        r(j,i) = r(i,j);
    end
end

% lowerhalf = (tril(ones(n),-1)>0);
% r(lowerhalf) = r(lowerhalf');

z = atanh(r);
% % p-calculation is temporary. Maybe the DOF for the calculation depend on
% % the number of frequency points in the band. For the moment let's use the
% % n time points.
% 
% % Tstat = +/-Inf and p = 0 if abs(r) == 1, NaN if r == NaN.
% Tstat = r .* sqrt((n-2) ./ (1 - r.^2));
% p = zeros(n);
% p(lowerhalf) = 2*tpvalue(-abs(Tstat),nv-2);
% p = p + p' + diag(diag(r)); % Preserve NaNs on diag.

p = zeros(n);       % per il momento no p.

return
end

function [coherence] = RtoR_coherency_calculation(x,y)
global opt
WS = (1/opt.tr); %MUST BE FIXED, SINCE NOW TR IS A CELL
[fxx,W]=cpsd_R2014b(x,x,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);
fxy = cpsd_R2014b(x,y,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);
fyy = cpsd_R2014b(y,y,hanning(opt.coherency.HW_length,'periodic'),opt.coherency.NO,opt.coherency.NFFT,WS);

% coherency = fxy./sqrt(fxx.*fxy);
% coherence = coherency.*conj(coherency);

bands = opt.coherency.bands;   
f1 = bands(1,1);
f2 = bands(1,2);
ind = find(W>=f1 & W<f2);
a = ind(1);
b = ind(end); 
coherence = abs(sum(fxy(a:b)).^2)./(sum(fxx(a:b)).*sum(fyy(a:b)));
return
end

function [output] = brainmask_selector(group,subj,type)
%type 
% 1 bmask_path
% 2 bmask_string (inlcude -mask and path)
% 3 bmask_img
global opt
if opt.bm.use_ss
    switch type
        case 1
            output = opt.AUX.bmask_path{group,subj};
        case 2
            output = opt.AUX.bmask_string{group,subj};
        case 3
            output = opt.AUX.bmask_img{group,subj};
    end
else
    switch type
        case 1
            output = opt.AUX.bmask_path{1,1};
        case 2
            output = opt.AUX.bmask_string{1,1};
        case 3
            output = opt.AUX.bmask_img{1,1};
    end    
end
return
end

function acompcor_mask_definition(g,k)
% creation of masks for acompcor. Pay attention: voxels with zero-singal
% could be present in these masks, however, they will be removed in the
% actual acompcor function. 

% todo: it will be better to have the masks already corrected for this
% problem.

global opt

tissues = [opt.aCompCor.ROI.tissue];
if ~isempty(find(tissues == 4,1)) % if GSR, load all masks
    hdr = spm_vol( c_b(opt.DATA.WC1{g}(k,:)) ); wc1 = spm_read_vols(hdr);
    hdr = spm_vol( c_b(opt.DATA.WC2{g}(k,:)) ); wc2 = spm_read_vols(hdr);
    hdr = spm_vol( c_b(opt.DATA.WC3{g}(k,:)) ); wc3 = spm_read_vols(hdr);
else
    if ~isempty(find(tissues == 1,1)); hdr = spm_vol( c_b(opt.DATA.WC1{g}(k,:)) ); wc1 = spm_read_vols(hdr); end
    if ~isempty(find(tissues == 2,1)); hdr = spm_vol( c_b(opt.DATA.WC2{g}(k,:)) ); wc2 = spm_read_vols(hdr); end
    if ~isempty(find(tissues == 3,1)); hdr = spm_vol( c_b(opt.DATA.WC3{g}(k,:)) ); wc3 = spm_read_vols(hdr); end
end

for r = 1:opt.aCompCor.masknumb
    THR   = opt.aCompCor.ROI(r).thr;  
    ERODE = opt.aCompCor.ROI(r).erode;
    DIME  = opt.aCompCor.ROI(r).dime;
    switch opt.aCompCor.ROI(r).tissue    %not the best way....
        case 1
            roi = wc1;
            output_str = 'GM';
            opt.aCompCor.ROI(r).name = 'GM';
        case 2
            roi = wc2;
            output_str = 'WM';
            opt.aCompCor.ROI(r).name = 'WM';
        case 3
            roi = wc3;
            output_str = 'CSF';
            opt.aCompCor.ROI(r).name = 'CSF';
        case 4
            roi = wc1+wc2+wc3;
            roi(roi > 0) = 1;
            output_str = 'GS_MASK';
            opt.aCompCor.ROI(r).name = 'GS_MASK';
            THR = 0; % 0verride thr in this case. But you still can erode the final mask
    end
    %-----erode mask-------
    idx1=find(roi(:)>THR);
    if isempty(idx1); error(sprintf('No suprathreshold voxels in ROI file %s.',tmp));end
    [idxx,idxy,idxz]=ind2sub(size(roi),idx1);
    idxt=find(idxx>ERODE&idxx<size(roi,1)+1-ERODE&idxy>ERODE&idxy<size(roi,2)+1-ERODE&idxz>ERODE&idxz<size(roi,3)+1-ERODE);
    for n1=1:length(idxt), if (sum(sum(sum(roi(idxx(idxt(n1))+(-ERODE:ERODE),idxy(idxt(n1))+(-ERODE:ERODE),idxz(idxt(n1))+(-ERODE:ERODE))<THR,3),2),1))>1, idxt(n1)=0; end; end
    idxt=idxt(idxt>0);
    idx1=idx1(idxt);
    if length(idx1) < DIME + 2; error(sprintf('After erosion, too few suprathreshold voxels in ROI file %s.',[opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,k},'_',output_str]));end        
    roi=zeros(size(roi));roi(idx1)=1;
    if 1
        %check mask. Can be left 0
        hdr.fname = ['_aCompCor_',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,k},'_',output_str,'.nii'];
        hdr.private.dat.fname = ['_aCompCor_',opt.AUX.group_names_prefix{g},opt.AUX.subject_names{g,k},'_',output_str,'.nii'];
        spm_write_vol(hdr,roi);
    end
    opt.aCompCor.ROI(r).img{g,k} = roi;
end
return
end

function acompcor(g,k,j)  %g group, k subj, j session
% see neuroimage 37(2007) 90-101
%there are two differences respect to the original paper of CompCor. These
%differences are also present in CONN toolbox (see related papeper). 1) PCA is
%done separately in different noise ROIs (eg WM and CSF, Originally was a
%unique roi). 2) Before PCA data is orhogolized respect to head movements.

%TODO: 1) move rp regression out of the loop 2) extract roi witout loop
global opt
data = spm_read_vols(spm_vol(c_b(opt.DATA.FUNCTIONALS{g,j}(k,:))));
% Volume selector
Norig = size(data,4);
indx = volume_selector_matlab(Norig,g,j);  %the function that cares of no selection case
data(:,:,:,indx) = [];
%----------------
s_d = size(data);
Xall = [];
for r = 1:opt.aCompCor.masknumb
    roi = opt.aCompCor.ROI(r).img{g,k};
    vn = sum(roi(:));
    V = nan(s_d(4),vn); %prealloc, first dim = time, second dim = voxels
    for l = 1:s_d(4)
        tmp = data(:,:,:,l);
        V(l,:) = tmp(roi == 1);
    end
    % removal of costant and linear trend
    V = detrend(V);
    % lets remove voxels whose variance is equal zero (no signal in those voxels)
    % and Nan values
    stdv = std(V);
    a = zeros(1,length(stdv));
    a(stdv == 0) = 1;
    a(isnan(stdv)) = 1;
    V(:,logical(a)) = [];
    %------------------Orthogonalize V-------------------------------------
    COV = [];
    if opt.aCompCor.asCONN %CONN first extract the mean signal (mS), then compute PCA over mS regressed out
        mS = detrend(mean(V,2));
        COV = [COV,mS];
    end
    if opt.rp_regression.do  % extract PCA over data already ortogonalized to rp.
                          % In this way the model is maximally predictive. 
        if opt.aCompCor.rpOrtogonalize
            %load regressor
            rp = load(c_b(opt.DATA.RP_1D{g,j,k}));
%             if opt.AUX.v_s_do  
%                 rp(indx,:) = [];
%             end
            COV = [COV,detrend(rp)];
        end
    end
    if ~isempty(COV)
        V = V-COV*(COV\V);
    end
    %----------------------------------------------------------------------
    if opt.aCompCor.ROI(r).dime > 0
        if opt.aCompCor.TvarianceNormalise
            %tvariance normalization
            V = bsxfun(@rdivide,V,std(V));
        end
        if 0 %same results   verified (2018)
            [~,U,~] = pca(V);
            comp = U(:,1:opt.aCompCor.ROI(r).dime);
         else
%           [U,~,~] = svd(V);
%           comp = U(:,1:opt.aCompCor.ROI(r).dime);
            %weight each component (doesn't change the correlation
            %structure)
            %-------A bit of math--------------------------------------------------------------------------------
            % The principal components T of a matrix A are obtained by projecting A on the space
            % of the eigenvectors of the convariance matrix A'A, that we call V (also known as principal axes). 
            % T = AV
            % A can be written is SVD as
            % A = USV' where U are eigenvectors of AA'
            %                V are eigenvectors of A'A
            % Thus, 
            % T = USV'V = US
            % If we take the svd 
            % A'A = VS^2V'
            % AA' = US^2U
            % Thus, these are equivalent:
            %          [U1,U2] = svd((V*V'));
            %          comp1 = [mV,U1(:,1:4)*diag(sqrt(diag(U2(1:4,1:4))))];
            %          [U1,U2] = svd(V);
            %          comp2 = [mV,U1(:,1:4)*diag(diag(U2(1:4,1:4)))];
            %          [coeff score latent] = pca(V);        
            %          comp2 = [mV,score(:,1:4)];
            %----------------------------------------------------------------------------------------------------
            [U,P] = svd(V);
            if opt.aCompCor.asCONN %add the mean signal and remove one dimension
                comp = [mS,U(:,1:opt.aCompCor.ROI(r).dime-1)*diag(diag(P(1:opt.aCompCor.ROI(r).dime-1,1:opt.aCompCor.ROI(r).dime-1)))];           
            else
                comp = U(:,1:opt.aCompCor.ROI(r).dime)*diag(diag(P(1:opt.aCompCor.ROI(r).dime,1:opt.aCompCor.ROI(r).dime)));           
            end
        end
    else  %if dime == 0 simply compute the straight average
        comp = mean(V,2);
    end
    % first derivatives computation, if requested
    if opt.aCompCor.ROI(r).deri > 0
        d = [];
        for l = 1:opt.aCompCor.ROI(r).deri
            d = [d, [zeros(l,size(comp,2));diff(comp,l)] ]; % I have to add l zeros as first rows...
        end
    else 
        d = [];
    end
    X = [comp,d];
    % variance normalize the extracted componets
    X = X./std(X,0,1);
    if 0 % add a costant regressor. 
        %There is no need to add a costant term. 
        %EPI time series is detrended before regression by 3drsfc
        X = [X,ones(size(X,1),1)];
    end
    opt.aCompCor.ROI(r).X{g,j,k} = X;
    Xall = [Xall,X];
end
opt.aCompCor.X{g,j,k} = Xall;
% Xall_afni = zeros(Norig,size(Xall,2));  %now I need to restore the original shape for afni (the cutting of volume is made by it).
% Xall_afni(~indx,:) = Xall;
%save('_temp.txt','Xall_afni','-ascii');
save('_temp.txt','Xall','-ascii');
str_out = ['_aCompCor_',opt.AUX.output_name{g,j,k},'_regressors.1D'];
[~,~] = system(['1dcat _temp.txt > ',str_out]);
% if opt.rp_regression.do
%     tmp = c_b(opt.DATA.RP_1D{g,j}(k,:));
%     [~,~] = system(['1dcat ',tmp,' _temp.txt > ',str_out]);
%     clear tmp
% else
%     [~,~] = system(['1dcat _temp.txt > ',str_out]);
% end
[~,~] = system('rm _temp.txt');
opt.aCompCor.aCC1D{g,j,k} = [opt.folders.preprocessing,'/',str_out];
% opt.X1D{g,j,k} = [opt.folders.preprocessing,'/',str_out];
return
end

function rp_initialization(g,j,k)
% it creates 1D files, computes derivative if required, computes FDpower.
% TODO: add FDafni (maybe not worthed, it is just a different norm (L1 vs L2) 
% FD is computed on cut data
global opt
rp_path = c_b(opt.DATA.RP{g,j}(k,:));
rp = load(rp_path);
rp_s = size(rp,2);
% check if # of columns are different across subjs 
if ~isfield(opt.DATA,'RP_raw_columns')
    opt.DATA.RP_raw_columns = rp_s;
else
    if opt.DATA.RP_raw_columns ~= rp_s
        error('RP files have different number of columns across subjects.');
    end
end
%--------------------------------------------------
% compute FD (I don't need derivatives for this)
compute_FD(rp,g,j,k); %group, session, subj
%--------------------------------------------------

output_name = ['_rp_',opt.AUX.output_name{g,j,k},'.1D'];
% if 0 %test the effect of detrending
%     rp = detrend(rp);
%     save('_temp.txt','rp','-ascii');
%     [~,~] = system(['1dcat _temp.txt > ',output_name,'.1D']);
%     system('rm _temp.txt');
% end
if opt.rp_regression.derivate && rp_s == 6   % derivative
    [status,string] = system(['1d_tool.py -infile ',rp_path,' -set_nruns 1 -derivative -write _temp.1D']); control(status,string);
    [~,~] = system(['1dcat ',rp_path,opt.AUX.v_s_rp_str{g,j},' _temp.1D',opt.AUX.v_s_rp_str{g,j},' > ',output_name]);
    [~,~] = system('rm _temp.1D');
else
    [~,~] = system(['1dcat ',rp_path,opt.AUX.v_s_rp_str{g,j},' > ',output_name]);
end
% save path
opt.DATA.RP_1D{g,j,k} = [opt.folders.preprocessing,'/',output_name];


return
end


function censoring_initialization(g,j,k)
global opt
% FD is computed on cut data (in case of volume selector). The censoring
% file has the lenght of cut data (both in opt and in 1D file)
tmask = er_censoring_mask(opt.QC.FD{g,j,k},opt.censoring.value,opt.censoring.pre_TR,opt.censoring.post_TR);
% save info and vector in opt:
opt.censoring.censor{g,j,k} = tmask;
n_cen = sum(not(tmask));
opt.censoring.censor_number{g,j}(k,1) = n_cen;
% save on .1D for 3dTproject
cen_name = ['_censoring_',opt.AUX.output_name{g,j,k},'_censor.1D'];
opt.censoring.censor1D{g,j,k} = [opt.folders.preprocessing,'/',cen_name];
dlmwrite(cen_name,tmask);
%[~,~] = system(['1dcat _temp.txt > ',cen_name]);
%[~,~] = system('rm _temp.txt');
return
end


% function censoring_initialization_old(cen,g,j,k)
% global opt
% % Warning: at the moment censoring is based on FDafni, not on FDpower.
% 
% % name_file = c_b(cen.data{g,j}(k,:));
% name_file = c_b(cen.data{g,j,k});
% output_name = ['_censoring_',opt.AUX.output_name{g,j,k}];
% [status,string] = system(['1d_tool.py -infile ''',name_file,'[1..6]'' -set_nruns 1 -show_censor_count',opt.censoring.extraoption,' -censor_motion ',num2str(opt.censoring.value),' ',output_name,' -overwrite']);control(status,string);
% %[status,string] = system(['1d_tool.py -infile ''',name_file,'[1..6]'' -set_nruns 1 -show_censor_count',opt.censoring.extraoption,' -censor_motion ',num2str(opt.censoring.value),' ',output_name,' -collapse_cols weighted_enorm -weight_vec 1 1 1 1 1 1 -overwrite']);control(status,string);
% %[status,string] = system(['1d_tool.py -infile ''',name_file,'[1..6]'' -set_nruns 1 -show_censor_count',opt.censoring.extraoption,' -derivative -collapse_cols weighted_enorm -weight_vec 1 1 1 1 1 1 -moderate_mask -',num2str(opt.censoring.value),' ',num2str(opt.censoring.value),' -write_censor ',output_name,'_censor.1D -write_CENSORTR ',output_name,'_CENSORTR.txt -overwrite']);control(status,string);
% cen_name = [output_name,'_censor.1D'];
% tmp = load(cen_name);
% if opt.AUX.v_s_do
%     tmp = tmp(opt.AUX.v_s_matlab{g,j}(1):opt.AUX.v_s_matlab{g,j}(2));
% end
% n_cen = sum(not(tmp));
% opt.censoring.censor1D{g,j,k} = [opt.folders.preprocessing,'/',cen_name];
% opt.censoring.censor_number{g,j}(k,1) = n_cen;
% return
% end

function prepare_atlas
global opt
radius = 5;
cw = pwd;
[~,~]= system('rm -r ./atlas');
mkdir(cw,'atlas'); % we are in denoising folder
dest_dir = [cw,'/atlas'];
str_f = 'atlas_modify.m';
p = which(str_f);  %horrible way...but it works
folder = [p(1:(end-length(str_f))),'rois'];
cd (folder);
list_exist_atlas = dir(['power_atlas_',num2str(radius),'mm_*.nii']);
found = 0;
[~,string] = system(['3dinfo -ni -nj -nk ',opt.DATA.FUNCTIONALS{1,1}(1,:)]);%assuming data have all the same size
nxyz_v = str2num(string);
for l = 1:length(list_exist_atlas)
    [~, string] = system(['3dinfo -ni -nj -nk ',list_exist_atlas(l).name]);
    nxyz_a = str2num(string);
    if nxyz_a(1) == nxyz_v(1) && nxyz_a(2) == nxyz_v(2) && nxyz_a(3) == nxyz_v(3)
        found = 1;
        break
    end
end
if found
    % copy atlas to denoing directory
    [~,~]=system(['cp ',list_exist_atlas(l).name,' ',dest_dir,'/atlas.nii']);
    [~,~]=system(['cp ',list_exist_atlas(l).name(1:end-4),'_distance.mat ',dest_dir,'/atlas_distance.mat']);
else
    % create a new atlas
    cd (dest_dir)
    atlas_create([folder,'/Power2011_rois.txt'],radius,[c_b(opt.DATA.FUNCTIONALS{1,1}(1,:)),'''[0]'''],'atlas.nii',1);
    try  %ER auto "learning". Possibility of permission required
        [~,~]=system(['3dcopy atlas.nii ',folder,'/power_atlas_',num2str(radius),'mm_',num2str(nxyz_v(1)),'_',num2str(nxyz_v(2)),'_',num2str(nxyz_v(3)),'.nii']);
        [~,~]=system(['cp atlas_distance.mat ',folder,'/power_atlas_',num2str(radius),'mm_',num2str(nxyz_v(1)),'_',num2str(nxyz_v(2)),'_',num2str(nxyz_v(3)),'_distance.mat']);
    catch 
    end
end
%-------LOAD AND PREPARE VARIABLES-----------------------------------------
opt.DA.atlas.img_path = [dest_dir,'/atlas.nii'];
opt.DA.atlas.distance_path = [dest_dir,'/atlas_distance.mat'];
%let's load them
opt.DA.atlas.img = spm_read_vols(spm_vol(opt.DA.atlas.img_path));
opt.DA.atlas.distance = load(opt.DA.atlas.distance_path);
% create some variable that are common for each loop
indx = ones(size(opt.DA.atlas.distance.d));
indx = triu(indx,1); % upper diagonal (no diag)
opt.DA.atlas.indx = indx; 
opt.DA.atlas.distance.d_vect = opt.DA.atlas.distance.d(indx == 1); clear indx;
[opt.DA.atlas.distance.d_vect_sorted,opt.DA.atlas.distance.d_vect_sorted_indx] = sort(opt.DA.atlas.distance.d_vect);
cd (cw);
return
end

function project_summary
% print summary of current project.
%todo. Plot all possible information. Also this function should work
%separately,i.e, given the mat of the project it print all the info
global opt
c = fix(clock);
c_str = [num2str(c(4)),':',num2str(c(5)),':',num2str(c(6)),' ',num2str(c(2)),'/',num2str(c(3)),'/',num2str(c(1))];
fprintf('\n------------------------------------------------------------');
fprintf('\nProject name:\n\t%s',opt.folders.prject_name);
fprintf('\nProject folder:\n\t%s',opt.folders.results);
% fprintf('\nOptions:');
% fprintf('\n\tTR(s): %s',opt.folders.results);
% fprintf('\n\tFilter: %s',opt.folders.results);
% fprintf('\n\tRP regression: %s',opt.folders.results);
% fprintf('\n\taCompCor: %s',opt.folders.results);
% fprintf('\n\tPSC: %s',opt.folders.results);
fprintf('\nStarting time:\n\t%s',c_str);
fprintf('\n------------------------------------------------------------\n');
return
end

function str = remove_endofline(str)
% even if it doesn't find the end of line it returns the correct string
str = regexp(str,'[\f\n\r]','split');
str = str{1};
end

function str = r_ext(str)  %r_ext stands for remove extension 
%remove any chracters after the last dot
indx = find(str == '.', 1, 'last');
if isempty(indx)
    return
else
    str = str(1:(indx-1));
end
return
end

function [str] = c_b(str)    %c_b stands for cut blanks
%remove blank space at end of the string (issue due to conn_dir)
str(str==' ') = '';
return
end

function [str] = r_b(str)    %c_b stands for replace blanks
%replace blank space with '_' 
str(str==' ') = '_';
return
end

function str = fix_float_errors(str)
% remove AFNI warning "corrected tot float errors 
% apparently can't be turned off via environment variables

mark = 'float errors';

indx = strfind(str,mark);

if not(isempty(indx))
    str(1:(indx + length(mark))) = [];
end

return
end