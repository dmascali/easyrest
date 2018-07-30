function denoising_analysis_extractSS(g,k,j) %g group, k subj, j session
%this function could have been in the denoising function. However, in this
%manner it is possible to extract already produced mat file. Usefull in
%case of debbugging. 
global opt DA_extracted

%Fill DA_extracted for group/session comparison

load (['ER_DA_',c_b(opt.AUX.output_name{g,j,k}),'.mat']);



DA_extracted.fd{g,j,k} = DA.rp_based.fd_mean;
DA_extracted.fd_rms{g,j,k} = rms(detrend(DA.rp_based.fd));
DA_extracted.fd_TRnorm_rms{g,j,k} = rms(detrend(DA.rp_based.fd_TRnorm));
DA_extracted.fd_TRnorm{g,j,k} = DA.rp_based.fd_TRnorm_mean;
DA_extracted.rms{g,j,k} = DA.rp_based.rms_mean;
DA_extracted.rms_dvars{g,j,k} = rms(detrend(DA.intensity_based.DVARS));
DA_extracted.mean_dvars{g,j,k} = mean(DA.intensity_based.DVARS);
DA_extracted.rms_sd{g,j,k} = rms(detrend(DA.intensity_based.SD));
count = 0;
for r = 2:length(DA.reg)
    count = count + 1;
    DA_extracted.vTv{g,j,k,count} = DA.reg(r).vTv_sampled.y;
    DA_extracted.vTv_mean{g,j,k,count} = DA.reg(r).vTv_sampled.mean;
    DA_extracted.var_exp{g,j,k,count} = DA.reg(r).explained_var.mean;
    DA_extracted.tSNR_mean{g,j,k,count} = DA.reg(r).tSNR.mean;
    for s = 1:length(DA.reg(r).specificity)
        DA_extracted.specifcity{g,j,k,count,s} = DA.reg(r).specificity(s).zf.s;
        % common for all
        if g == 1 && k == 1 && j == 1
            DA_extracted.specifcity_names{s} = opt.DA.roi_specificity(s).name;
        end
    end
    %----boldmotion
    %DA_extracted.bold_motion.fd_zf{g,j,k,count} = DA.reg(r).bold_motion.fd.zf;
    DA_extracted.bold_motion.fd_zs{g,j,k,count} = DA.reg(r).bold_motion.fd.zs;
    %DA_extracted.bold_motion.fd_zf_mp{g,j,k,count} = DA.reg(r).bold_motion.fd.zf_mean_p;
    DA_extracted.bold_motion.fd_zs_mp{g,j,k,count} = DA.reg(r).bold_motion.fd.zs_mean_p;
    %DA_extracted.bold_motion.fd_zf_mn{g,j,k,count} = DA.reg(r).bold_motion.fd.zf_mean_n;
    DA_extracted.bold_motion.fd_zs_mn{g,j,k,count} = DA.reg(r).bold_motion.fd.zs_mean_n;
    %DA_extracted.bold_motion.dvars_zf{g,j,k,count} = DA.reg(r).bold_motion.dvars.zf;
    DA_extracted.bold_motion.dvars_zs{g,j,k,count} = DA.reg(r).bold_motion.dvars.zs;
    %DA_extracted.bold_motion.dvars_zf_mp{g,j,k,count} = DA.reg(r).bold_motion.dvars.zf_mean_p;
    DA_extracted.bold_motion.dvars_zs_mp{g,j,k,count} = DA.reg(r).bold_motion.dvars.zs_mean_p;
    %DA_extracted.bold_motion.dvars_zf_mn{g,j,k,count} = DA.reg(r).bold_motion.dvars.zf_mean_n;
    DA_extracted.bold_motion.dvars_zs_mn{g,j,k,count} = DA.reg(r).bold_motion.dvars.zs_mean_n;
    %--------------------
    
    %----ROItoROI
    DA_extracted.roitoroi.roi_number{g,j,k,count} = DA.reg(r).roitoroi.roi_number;
    DA_extracted.roitoroi.roi_nan{g,j,k,count} = DA.reg(r).roitoroi.roi_nan;
    DA_extracted.roitoroi.d_vect{g,j,k,count} = DA.reg(r).roitoroi.d_vect;
    DA_extracted.roitoroi.d_vect_avg{g,j,k,count} = DA.reg(r).roitoroi.d_vect_avg;
    %DA_extracted.roitoroi.zf.matrix{g,j,k,count} = DA.reg(r).roitoroi.zf.matrix;
    DA_extracted.roitoroi.zs.matrix{g,j,k,count} = DA.reg(r).roitoroi.zs.matrix;
    DA_extracted.roitoroi.zf.vect{g,j,k,count} = DA.reg(r).roitoroi.zf.vect;
    DA_extracted.roitoroi.zs.vect{g,j,k,count} = DA.reg(r).roitoroi.zs.vect;
    
    %DA_extracted.roitoroi.zf.delta_matrix{g,j,k,count} = DA.reg(r).roitoroi.zf.delta_matrix;
    DA_extracted.roitoroi.zs.delta_matrix{g,j,k,count} = DA.reg(r).roitoroi.zs.delta_matrix;
    DA_extracted.roitoroi.zf.delta_vect{g,j,k,count} = DA.reg(r).roitoroi.zf.delta_vect;
    DA_extracted.roitoroi.zs.delta_vect{g,j,k,count} = DA.reg(r).roitoroi.zs.delta_vect;
    %--------------------
    
    % common for all
    if g == 1 && k == 1 && j == 1
        DA_extracted.name_contract{count} = DA.reg(r).name_contract;
        DA_extracted.colors = opt.DA.colors;
    end
end


return
end

