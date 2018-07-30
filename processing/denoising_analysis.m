function denoising_analysis(g,k,j)  %g group, k subj, j session
% This function can be called only if rp has been performed
% it works separately for each subject

%BE AWERE. This works in the single subject GM space. However, at the
%begging, an intersection between the SS GM and the group brain mask is
%performed. In this way all subjs and sessions are are tested on a
%"common" space 
%However, when SS mode is ON, this is not true. Yet, the SS gm mask are
%still insteresected among sessions. So, session comparison is perfectly
%fine. (this last step is done becouse remotion of zerovoxel is performed
%also on the pro_data)

global opt
global min_voxel_mode
output_name_mat = ['ER_DA_',opt.AUX.output_name{g,j,k},'.mat'];

%-------------%avoid rerunning-------------
if opt.AUX.DA_skip_SS    % for debbuging purpose to use if ALL the data have already been extracted (save time).
    return               % to apply this, modify the variable insied easyrest
end

% this handle the case that some data have been already extracted (not
% all). This can happen in case of unexpected errors
if exist(output_name_mat,'file') == 2 && ~opt.DA.force_overwrite
    return
end

%-----------------------------------------

wc1_thres = 0.80;   % thr for gm voxels selection. In some patholgical population this value should be reduced
                    % this value is not used when GM SS modality is
                    % selected. The thr of the GM is used.
nresampled = 15000; % only for voxel to voxel correlation. Adjust this value according to the available memory
save_data = 1;   %usually you don't need to save extracted data.
flag_norm = 1;   % 4d Intensity normalization. 
                 % 0 none
                 % 1 mode
                 % 2 mean
min_voxel_mode = 1;     % 1/0; In this modality ROI for a given specificity will have the same number of voxels.
                        % cons: selected voxel are chosen randomly, the
                        % might be not adjacent
                        
                        
%------------------load raw and processed data
raw_data = single(spm_read_vols(spm_vol(c_b(opt.DATA.FUNCTIONALS{g,j}(k,:)))));
lff_path = [opt.folders.preprocessing,'/',opt.AUX.output_name{g,j,k},'_LFF.nii'];    % TODO, separate LFF and sLFF ? 
pro_data = single(spm_read_vols(spm_vol(lff_path)));
%------------------------------------------------------------

%------------------ Volume selector, only for raw_data
if opt.AUX.v_s_do  % the real cut is done by afni so at the end I need to restore the original size
    indx = volume_selector_matlab(size(raw_data,4),g,j);
    raw_data(:,:,:,indx) = [];
end
%------------------------------------------------------------


if opt.bm.use_ss && opt.bm_method == 2  % if using ss mask of gm avoid recomputing
    gm = opt.AUX.bmask_img{g,k};
else
    wc1 = spm_read_vols(spm_vol(opt.DATA.WC1{g}(k,:)));
    wc1(wc1<=wc1_thres) = 0;
    wc1(wc1>0) = 1;
    gm = wc1;
    if ~opt.bm.use_ss
        gm = gm.*opt.AUX.bmask_img{1,1};  
    end
    clear wc1;
end

s = size(raw_data);
N = s(4);
raw_data = reshape(raw_data,[s(1)*s(2)*s(3),s(4)]);
pro_data = reshape(pro_data,[s(1)*s(2)*s(3),s(4)]);
gm = reshape(gm,[s(1)*s(2)*s(3),1]);
for ss = 1:length(opt.DA.roi_specificity)
    for l = 1:length(opt.DA.roi_specificity(ss).roi)
        spe(ss).rois(l).img = opt.DA.roi_specificity(ss).roi(l).img;
        spe(ss).rois(l).img = reshape(spe(ss).rois(l).img,[s(1)*s(2)*s(3),1]);
        spe(ss).rois(l).name = opt.DA.roi_specificity(ss).roi(l).name_contract;
    end
end
%Atlas load and reshape
atlas = opt.DA.atlas;
atlas.img = reshape(opt.DA.atlas.img,[s(1)*s(2)*s(3),1]);
%atlas.distance = opt.DA.atlas.distance;

%remove zero voxels. It works also on the GM mask, ROIs and atlas.
%ATTENTION: for most BM modalities this step is not needed. Since, it is a
%time consuming step...considered to bypass it in specific case
[raw_data,pro_data,gm,spe,atlas] = remove_zero_voxels(raw_data,pro_data,gm,spe,atlas);
%-----------------------------------------------

% intensity normalization. Following Power let's normalize to a mode value of 1000;
% why the mode instead of the mean?? With the mode you don't need to define a
% brain mask.
switch flag_norm
    case 0
        % I don't know why don't normalize...anyway let's consider also
        % this case
        opt.DA.normMode = 'none'; normfactor = [];
    case 1 %mode normalization
        [raw_data,pro_data,normfactor] = mode_normalization(raw_data,pro_data);    % TIME CONSUMING STEP!!
        opt.DA.normMode = 'mode';        
    case 2 % mean normalization
        [raw_data,pro_data,normfactor] = mean_normalization(raw_data,pro_data);
end     
opt.DA.normFactor(g,k,j) = normfactor;
%-----------------------------------------------

%Computing DVARS (see Power et al, NeuroImage, 59(3), 2012) on the entire brain (no gm intersection)
% and SD computation over the all brain (hence SDgs)
% adn GS
[dvars,SD,GS] = dvars_SD_GS_computation(raw_data); 
dvars = [dvars(1);dvars];  % due to the differenziation I need to add a fake value. Power use 0 instead 
%-----------------------------------------------

%------------- GM intersection
raw_data(not(logical(gm)),:) = [];
pro_data(not(logical(gm)),:) = [];
gmv = sum(gm(:));
for ss = 1:length(spe)
    for l = 1:length(spe(ss).rois)
        spe(ss).rois(l).img(not(logical(gm))) = [];
        spe(ss).rois(l).vcount = sum(spe(ss).rois(l).img(:));
    end
end
atlas.img(not(logical(gm))) =[];
%voxel count
fprintf(opt.DA.fid,'\n%s',opt.AUX.output_name{g,j,k});
fprintf(opt.DA.fid,'\n%d\t\t',gmv);
for ss = 1:length(spe)
    for l = 1:length(spe(ss).rois)
        fprintf(opt.DA.fid,'%d\t',spe(ss).rois(l).vcount);
    end
    fprintf(opt.DA.fid,'\t');
end
%-----------------------------------------------

%it is convinient to have time in rows and voxels in columns:
raw_data = raw_data';
pro_data = pro_data';

% compute t-mean and store for tSNR computation.
mean_raw_data = mean(raw_data,1);

% load compcor
if opt.aCompCor.do
    comp = opt.aCompCor.X{g,j,k};
    % create string for compcor
    if ~isfield(opt.DA,'str_comp')
        str_comp = [];
        for r = 1:opt.aCompCor.masknumb
            str_comp = [str_comp,opt.aCompCor.ROI(r).name,' '];
            if opt.aCompCor.ROI(r).deri > 0
                der_str = '''';
                for l = 1:opt.aCompCor.ROI(r).deri
                    str_comp = [str_comp,opt.aCompCor.ROI(r).name,der_str];
                end
                str_comp = [str_comp,' '];
            end
        end
        str_comp = str_comp(1:end-1);
        opt.DA.str_comp = str_comp;
    end
end

% load rp.
%rp = load(c_b(opt.DATA.RP{g,j}(k,:)));
rp = load(c_b(opt.DATA.RP_1D{g,j,k}));
%------------------ Volume selector, rp----------------------
if opt.AUX.v_s_do  
    rp(indx,:) = [];
end
%------------------------------------------------------------
nrp = size(rp,2);
%there is a problem, spm put first traslation, fsl rotations. Now I try to
%guess what type is it considering that traslation should be aroud 50 times
%greater than rotations.
% it's not alwasys working (i.e, HCP!) This info must be provided
%             tmp_f3 = mean(mean(rp(:,1:3),2));
%             tmp_l3 = mean(mean(rp(:,4:6),2));
%             if tmp_f3 > tmp_l3
%                 tra = rp(:,1:3);
%                 rot = rp(:,4:6);
%                 rp_type_control_str = 'Frist columns: traslations';
%             else
%                 tra = rp(:,4:6);
%                 rot = rp(:,1:3);
%                 rp_type_control_str = 'Frist columns: rotations';
%             end

tra = rp(:,opt.AUX.rp_indx.tra);
rot = rp(:,opt.AUX.rp_indx.rot);

%FD Frame Displacement (see Power et al, NeuroImage, 59(3), 2012)
fd_tra = sum( abs( diff(tra)  ) ,2);
%rotaions must be converted from rad to mm. (assuming a radius of 50 mm)
switch lower(opt.AUX.rp_rot_unit)
    case {'rad'}
        rot_mm = rot*opt.AUX.rp_average_brain_radius; 
    case {'mm'} 
        %do nothing
        rot_mm = rot;
    case {'deg'}
        rot_mm = rot*opt.AUX.rp_average_brain_radius*2*pi/360;         
end
fd_rot = sum( abs( diff(rot_mm)  ) ,2);
fd = sum([fd_tra,fd_rot],2);
fd = [0; fd]; % due to the differenziation I need to add a 0.
fd_TRnorm = fd./opt.tr{g,j};  %lets consider also a normalized version of FD

%aboslute values and RMS.
%Pay attention: since the reference volume is arbitrary we can remove the
%mean from the rp before calculating abs and rms. However, the choice of the ref
%may change the results. 
rp_mm = [tra,rot_mm];
rp_mm_dm = bsxfun(@minus, rp_mm, mean(rp_mm));
tra_abs = sum(abs(rp_mm_dm(:,1:3)),2);
rot_abs = sum(abs(rp_mm_dm(:,4:6)),2);
rms_ = rms(rp_mm_dm);
rms_m = mean(rms_);

% define RP for later saving
RP.rp = rp;
RP.tra = tra;
RP.tra_abs = tra_abs;
RP.rot = rot;
RP.rot_mm = rot_mm;
RP.rot_abs = rot_abs;
RP.fd = fd;
RP.fd_mean = mean(fd);
RP.fd_TRnorm = fd_TRnorm;
RP.fd_TRnorm_mean = mean(fd_TRnorm);
RP.rms = rms_;
RP.rms_mean = rms_m;

%Now let's work on the raw gm series and remove one by one different
%effects
% no reg 
reg(1).name = 'raw data';
reg(1).name_contract = 'RAW';
reg(1).X = [];
reg(1).data = raw_data;
reg(1).specificity = [];
reg(1).dof = estimate_dof(N,0,opt.tr{g,j},0,inf);
% costant, linear, quadratic trend
[N,P] = size(raw_data);  % N time points, P voxels
reg(2).name = 'detrend (cost + lin + quad)';
reg(2).name_contract = 'DETREND'; 
reg(2).X = [ones(N,1),(1:N)'/N, ((1:N)'/N).^2];
reg(2).dof = estimate_dof(N,size(reg(2).X,2),opt.tr{g,j},0,inf);
% detrend + rp
reg(3).name = 'detrend + rp';
reg(3).name_contract = ['RP',num2str(size(rp,2))]; 
reg(3).X = [reg(2).X,rp];
reg(3).dof = estimate_dof(N,size(reg(3).X,2),opt.tr{g,j},0,inf);
%execute reg
reg(2).data = er_glmfit_RES(raw_data,reg(2).X); 
reg(3).data = er_glmfit_RES(raw_data,reg(3).X); 
if opt.aCompCor.do
    % detrend + compcor
    reg(4).name = 'detrend + compcor';
    reg(4).name_contract = opt.DA.str_comp; 
    reg(4).X = [reg(2).X,comp];
    reg(4).dof = estimate_dof(N,size(reg(4).X,2),opt.tr{g,j},0,inf);
    % detrend + rp + compcor
    reg(5).name = 'detrend + rp + compcor';
    reg(5).name_contract = [opt.DA.str_comp,' RP',num2str(size(rp,2))]; 
    reg(5).X = [reg(2).X,rp,comp];
    reg(5).dof = estimate_dof(N,size(reg(5).X,2),opt.tr{g,j},0,inf);
    %execute reg
    reg(4).data = er_glmfit_RES(raw_data,reg(4).X); 
    reg(5).data = er_glmfit_RES(raw_data,reg(5).X); 
end
% add here GS?

n_trials = length(reg);
reg(n_trials + 1).name = '3dRSFC';
reg(n_trials + 1).name_contract = 'FULL DENOISING'; 
reg(n_trials + 1).X = [];
reg(n_trials + 1).data = pro_data;
if opt.aCompCor.do; M = size(reg(5).X,2); else M = size(reg(3).X,2); end  % the case no rp is not expected
reg(n_trials + 1).dof = estimate_dof(N,M,opt.tr{g,j},opt.filter_band(1),opt.filter_band(2));
n_trials = length(reg);

reg = compute_explained_var_tSNR(reg,mean_raw_data);

%-------------------calculate specificity-------------
reg(1).specificity = [];
reg(2).specificity = calculate_specificity(reg(2).data,spe,reg(2).dof);
reg(3).specificity = calculate_specificity(reg(3).data,spe,reg(3).dof);
if opt.aCompCor.do
    reg(4).specificity = calculate_specificity(reg(4).data,spe,reg(4).dof);
    reg(5).specificity = calculate_specificity(reg(5).data,spe,reg(5).dof);
end
reg(n_trials).specificity = calculate_specificity(reg(n_trials).data,spe,reg(n_trials).dof);
%-----------------------------------------------------

%-------------------BOLD:MOTION correlation-----------
reg(1).bold_motion = [];
reg(2).bold_motion = bold_motion_correlation(reg(2).data,fd,dvars,reg(2).dof);
reg(3).bold_motion = bold_motion_correlation(reg(3).data,fd,dvars,reg(3).dof);
if opt.aCompCor.do
    reg(4).bold_motion = bold_motion_correlation(reg(4).data,fd,dvars,reg(4).dof);
    reg(5).bold_motion = bold_motion_correlation(reg(5).data,fd,dvars,reg(5).dof);
end
reg(n_trials).bold_motion = bold_motion_correlation(reg(n_trials).data,fd,dvars,reg(n_trials).dof);
%-----------------------------------------------------

%-------------------ROI-to-ROI/distance--------------
reg(1).roitoroi = do_roitoroi(reg(1).data,atlas,reg(1).dof,[]);
reg(2).roitoroi = do_roitoroi(reg(2).data,atlas,reg(2).dof,reg(1).roitoroi.r.matrix);
reg(3).roitoroi = do_roitoroi(reg(3).data,atlas,reg(3).dof,reg(1).roitoroi.r.matrix);
if opt.aCompCor.do
    reg(4).roitoroi = do_roitoroi(reg(4).data,atlas,reg(4).dof,reg(1).roitoroi.r.matrix);
    reg(5).roitoroi = do_roitoroi(reg(5).data,atlas,reg(5).dof,reg(1).roitoroi.r.matrix);
end
reg(n_trials).roitoroi = do_roitoroi(reg(n_trials).data,atlas,reg(n_trials).dof,reg(1).roitoroi.r.matrix);
%-----------------------------------------------------



% now compute voxel to voxel 
% 2ways 
       %1) full: effective computation (average across gm) (slow, memory problem)
       %2) sampled: sample around 10000 points and run corrcoef (approx., fast)
       
% % 0) correlation with average gm signal (test only) 
% for r = 2:n_trials
%     avg = mean(reg(r).data,2);
%     z_avg = zscore(avg);
%     z_epi = zscore(reg(r).data);
%     c = (z_avg\z_epi);  %check that the result is the same that fastconn
%     reg(r).vTv_avg = histo_related(c); 
% end
%1) full: effective computation (average across gm) (slow, memory problem)
% matlabpool open
% for r = 2:n_trials
%     An=bsxfun(@minus,reg(r).data,mean(reg(r).data,1));
%     An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
%     C = [];
%     
%         parfor l = 1:(P-1)
%             Ann = An(:,(l+1:P));
%             Bn = repmat(An(:,l),[1 P-l]);
%             C = [C,sum(Ann.*Bn)];
%             %C = sum(Ann.*Bn);
%             %C = sum( An(:,(l+1:P)).*repmat(An(:,l),[1 P-l]) );
%         end
%         
%         clear An Ann Bn
%     reg(r).vTv_full = histo_related(C); 
%     clear C;
% end
% matlabpool close

%2) sampled: sample around 10000 points and run corrcoef (approx., fast)
indx = zeros(1,P,'single');
if P > nresampled
    for l = 1:round(P/nresampled):P
        indx(1,l) = 1;
    end
end
for r = 2:n_trials
    reduced = reg(r).data;
    reduced(:,not(logical(indx))) = [];
    c = corrcoef(reduced); 
    id = triu(ones(size(c),'uint8'),1);
    c = c(id == 1);
    reg(r).vTv_sampled = histo_related(c); 
    clear c id;
end
clear reduced;

%------------dvars fd correlation---------------
% would be nice to look at this quantity at different regression stages.
% However, computing dvars again is tricky
DA.dvars_fd_corr.r = corrcoef(dvars,fd);
DA.dvars_fd_corr.zf = atanh(DA.dvars_fd_corr.r);
%-----------------------------------------------

%-------------------store variables-------------
DA.subj_name = opt.AUX.output_name{g,j,k};
DA.intensity_based.SD = SD;
DA.intensity_based.DVARS = dvars;
DA.intensity_based.GS = GS;
DA.rp_based = RP;
DA.reg = reg;
%-----------------------------------------------

%-----------do Single Subject plotting----------
denoising_analysis_plot(DA);
for ss = 1:length(spe)
    denoising_analysis_plot_specificity_timeseries(DA,ss);
end
denoising_analysis_plot_extra(DA,1);    %unit r
denoising_analysis_plot_extra(DA,0);    %unit z-score
%-----------------------------------------------


%-------------------save data-------------------
if ~save_data
    for r =1:length(DA.reg)
        DA.reg(r).data = [];
    end
end
save(output_name_mat,'DA');
%-----------------------------------------------

end

function [V,V2,gm,spe,atlas] = remove_zero_voxels(V,V2,gm,spe,atlas)
if nargin == 1
    stdv = std(V,[],2);
    a = zeros(size(stdv));
    a(stdv == 0) = 1;
    a(isnan(stdv)) = 1;
    V(logical(a),:) = [];
else
    stdv1 = std(V,[],2);
    stdv2 = std(V2,[],2);  %This is needed when SSBM is used but bm isnot 2
    a = zeros(size(stdv1));
    a(stdv1 == 0) = 1; a(isnan(stdv1)) = 1;
    a(stdv2 == 0) = 1; a(isnan(stdv2)) = 1;
    V(logical(a),:) = [];
    V2(logical(a),:) = [];
    gm(logical(a)) = [];
    for s = 1:length(spe)
        for l = 1:length(spe(s).rois)
            spe(s).rois(l).img(logical(a)) = [];
        end
    end
    atlas.img(logical(a)) = [];
end
return
end

function [V,V2,mode_] = mode_normalization(V,V2)
% it does not always converge..it seems to work better if I remove the
% extreme values
if 1
    tmp = V(:);
    percntiles = prctile(tmp,[5 99]);
    indx = tmp < percntiles(1) | tmp > percntiles(2);
    tmp(indx) = [];
    h = histfit(tmp,[],'kernel');
    x = get(h(2),'XData');
    y = get(h(2),'YData');
    mode_= x(find(y==max(y)));
else
    warning('You are in debugging mode, mode is aribirarily computed');
    mode_ = 10000;
end
V = (V/mode_)*1000;
V2 = (V2/mode_)*1000;
return
end

function [V,V2,mean_] = mean_normalization(V,V2)
% Fast normalization compared to the mode. It works perfectly when the
% images are already brain extracted. Otherwise mode is preferred
mean_ = mean(V(:));   % No zero voxel should be present at this stage.
V = (V/mean_)*1000;
V2 = (V2/mean_)*1000;
return
end

function reg = compute_explained_var_tSNR(reg,mean_)
% mean_ is the mean of raw data before any processing. Actually, data is
% rescaled before computing the mean, but this doesn't influence tsnr

ref = std(reg(1).data,[],1);
reg(1).explained_var.mean = 0;
reg(1).explained_var.std = 0;
reg(1).explained_var.points = [];
reg(1).tSNR.vec = mean_./ref;
reg(1).tSNR.mean = mean(reg(1).tSNR.vec);
reg(1).tSNR.std = std(reg(1).tSNR.vec);

for r = 2:length(reg)
    tmp = std(reg(r).data,[],1);
    frac = 100*(1-tmp./ref);
    reg(r).explained_var.mean = mean(frac);
    reg(r).explained_var.std = std(frac);
    reg(r).explained_var.points = length(frac);
    reg(r).tSNR.vec = mean_./tmp;
    reg(r).tSNR.mean = mean(reg(r).tSNR.vec);
    reg(r).tSNR.std = std(reg(r).tSNR.vec);
end
return
end

function [stats] = histo_related(V)
% create histograms from voxel to voxel measures.
stats.x = -0.99:0.01:0.99;
y = hist(V(:),stats.x);
%normalize to 1
my = sum(y); 
stats.y = y/my;
%compute mean and mode
stats.mean = mean(V(:));
stats.std = std(V(:));
return
end

function [RES] = er_glmfit_RES(V,X)
RES = V-X*(X\V);
% this is equivalent to 
% inversa = inv(X'*X);
% B = inversa*X'*Y;
% RES2 = Y -X*B;
% Both method work column wise, no need for cylces (adjust ER??)
return
end

function [dvars,SD,GS] = dvars_SD_GS_computation(V)
%DVARS (see Power et al, NeuroImage, 59(3), 2012)
V=detrend(V')';
dvars = sqrt(mean( (diff(V,1,2).^2) ,1))';
%differently form Power, now I am going to remove the mean otherwise the
%zero point would create problems.
%dvars = dvars-mean(dvars);

%----SD
SD = std(V,[],1);

%----GS
GS = mean(V,1);

return
end

function Z = calculate_specificity(X,spe,dof)
global min_voxel_mode

for s  = 1:length(spe)
    clear Y;
    Y.tmseries = [];
    if min_voxel_mode
        v_count = [];
        for l = 1:length(spe(s).rois)
            v_count = [v_count,spe(s).rois(l).vcount]; % no problem with the variable number of rois
        end
        Y.min_voxel = min(v_count);
    else
        Y.min_voxel = [];
    end
    for l = 1: length(spe(s).rois)
        tmp = X;
        tmp(:,not(logical(spe(s).rois(l).img))) = [];
        if min_voxel_mode
            tmp = tmp(:,1:Y.min_voxel);
        end
        mtmp = mean(tmp,2);
        if isempty(tmp)
            mtmp = nan(size(tmp,1),1);
            warning(['ROI ',spe(s).rois(l).name,' is empty']);
        end
        Y.tmseries = [Y.tmseries,mtmp];
    end
    Y.corr = corr(Y.tmseries);
    Y.corr_zfisher = atanh(Y.corr);
    Y.corr_zscore = Y.corr_zfisher.*sqrt(dof-3);
    
    % specificity----------------------------------------------------------
    % as defined in Chai2012 Neuroimage 59(2012) 1420-1428
    % S = (|Zt| -|Zr|)/(|Zt| +|Zr|)
    % Zt = correlation between PCC and MPFC
    % Zr = correlazion between non connected area eg. PCC and VIS
    % NB in this case there is no need to adjust for dof since it is an already
    % normalized measure
    Y.zf.t = Y.corr_zfisher(1,2);
    if length(spe(s).rois) == 4
        Y.zf.r = ( Y.corr_zfisher(1,3) + Y.corr_zfisher(1,4) )/2;
    elseif length(spe(s).rois) == 3
        Y.zf.r = Y.corr_zfisher(1,3);
    else
        error('I am not able to handle more than 4 rois for each specificity test.');
    end
    Y.zf.s = (abs(Y.zf.t) - abs(Y.zf.r) )./ (abs(Y.zf.t) + abs(Y.zf.r) );
    % also create zt zr and for z-score
    Y.zs.t = Y.corr_zscore(1,2);
    if length(spe(s).rois) == 4
        Y.zs.r = ( Y.corr_zscore(1,3) + Y.corr_zscore(1,4) )/2;
    elseif length(spe(s).rois) == 3
        Y.zs.r = Y.corr_zscore(1,3);
    end
    %----------------------------------------------------------------------
    %---names
    Y.names.t = [spe(s).rois(1).name,':',spe(s).rois(2).name];
    if length(spe(s).rois) == 4
        Y.names.r = [spe(s).rois(1).name,':<',spe(s).rois(3).name,',',spe(s).rois(4).name,'>'];
    elseif length(spe(s).rois) == 3
        Y.names.r = [spe(s).rois(1).name,':<',spe(s).rois(3).name];
    end
    for z = 1:length(spe(s).rois)
        Y.names.all{z} = spe(s).rois(z).name;
    end
    %--------
    Z(s) = Y;
end
return
end


function Z = bold_motion_correlation(Y,fd_v,dvars_v,dof)
%zscore trasformation of data and motion measures
Y = zscore(Y);
fd_v = zscore(fd_v);
dvars_v = zscore(dvars_v);
N = size(Y,2);  % NUMBER OF VOXELS
%check if they are all column vector
fd.corr = fd_v\Y;
fd.zf = atanh(fd.corr);
fd.zf_mean_p = sum(fd.zf(fd.zf > 0))/N;
fd.zf_mean_n = sum(fd.zf(fd.zf < 0))/N;
fd.zs = fd.zf.*sqrt(dof-3);
fd.zs_mean_p = sum(fd.zs(fd.zs > 0))/N;
fd.zs_mean_n = sum(fd.zs(fd.zs < 0))/N;
dvars.corr = dvars_v\Y;
dvars.zf = atanh(dvars.corr);
dvars.zf_mean_p = sum(dvars.zf(dvars.zf > 0))/N;
dvars.zf_mean_n = sum(dvars.zf(dvars.zf < 0))/N;
dvars.zs = dvars.zf.*sqrt(dof-3);
dvars.zs_mean_p = sum(dvars.zs(dvars.zs > 0))/N;
dvars.zs_mean_n = sum(dvars.zs(dvars.zs < 0))/N;
Z.fd = fd;
Z.dvars = dvars;
Z.N = N;
return
end

function Z = do_roitoroi(X,atlas,dof,C_ref)

img = atlas.img;
d = atlas.distance;
M = size(d.roi_mni,1);   %ORIGINAL NUMBER OF ROIS
N = size(X,1);           %NUMBER OF TIME POINTS

T = nan(N,M);
for l = 1:M
    tmp = X(:,img == l);
    T(:,l) = mean(tmp,2); %if the rois is empty it will be NaN

end
nan_indx = isnan(mean(T,1));
nan_number = sum(nan_indx);
roi_number = M-nan_number;

C = corr(T);
C_vect = C(atlas.indx == 1);
C_vect_zf = atanh(C_vect);
D_vect = d.d_vect;
%compute average curve
resolution = 6; %mm
l = min(D_vect):resolution:max(D_vect);
C_avg_zf = [];
for z = 2:(length(l))     
    C_avg_zf(z-1) = nanmean(C_vect_zf(D_vect <= l(z) & D_vect > l(z-1)));
end
C_avg_d = l(1:end-1) + diff(l)/2;


%DELTA R
if isempty(C_ref)
    C_ref = zeros*C;
end
DELTA_C = atanh(C)-atanh(C_ref);
DELTA_C_vect = DELTA_C(atlas.indx == 1);
%compute average curve
DELTA_C_avg = [];
for z = 2:(length(l))     
    DELTA_C_avg(z-1) = nanmean(DELTA_C_vect(D_vect <= l(z) & D_vect > l(z-1)));
end


% EXPORT (plus some minor computations)
Z.roi_number = roi_number;
Z.roi_nan = nan_indx;
Z.d_vect = D_vect;
Z.d_vect_avg = C_avg_d;
Z.d_vect_avg_resolution = resolution;

Z.r.matrix = C;
Z.r.vect = C_vect;
Z.r.vect_avg = tanh(C_avg_zf);

%z-fisher
Z.zf.matrix = atanh(C);
Z.zf.vect = C_vect_zf;
Z.zf.vect_avg = C_avg_zf;
%z-score
Z.zs.matrix = Z.zf.matrix * sqrt(dof-3);
Z.zs.vect   = Z.zf.vect * sqrt(dof-3);
Z.zs.vect_avg = C_avg_zf * sqrt(dof-3);
% DELTA R
Z.r.delta_matrix = tanh(DELTA_C);   
Z.r.delta_vect = tanh(DELTA_C_vect); 
Z.r.delta_vect_avg = tanh(DELTA_C_avg); 
Z.zf.delta_matrix = DELTA_C;
Z.zf.delta_vect = DELTA_C_vect;
Z.zf.delta_vect_avg = DELTA_C_avg;
Z.zs.delta_matrix = Z.zf.delta_matrix * sqrt(dof-3);
Z.zs.delta_vect   = Z.zf.delta_vect * sqrt(dof-3);
Z.zs.delta_vect_avg = DELTA_C_avg * sqrt(dof-3);
return
end

