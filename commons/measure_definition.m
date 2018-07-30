function [measure] = measure_definition
global opt

% for retro-compatibility
if ~isfield(opt.MEASURES,'RSFA')
    opt.MEASURES.RSFA = 1;
end
if ~isfield(opt.MEASURES,'ALFF')
    opt.MEASURES.ALFF = 1;
end
if ~isfield(opt,'prepro_mode')
    opt.prepro_mode = '3dRSFC';
end
%-------------------------------
measure(1).name_selector = 'ALFF';          %for selection of the measure
measure(1).name = 'ALFF';                   %for selection of the file
measure(1).dir_name = 'ALFF';               %for selection of the firstlevel folder
measure(1).output_dir = 'A-ALFF';             %output dir of the second level analysis
measure(1).seed_name = '';                  %for selection of the seed in the case of seedmeasure
measure(1).subbrick = 1;
measure(1).do = opt.MEASURES.ALFF;
measure(2).name_selector = 'mALFF';
measure(2).name = 'mALFF';
measure(2).dir_name = 'ALFF';
measure(2).output_dir = 'A-mALFF';         
measure(2).seed_name = '';       
measure(2).subbrick = 1;
measure(2).do = opt.MEASURES.ALFF;
measure(3).name_selector = 'fALFF';
measure(3).name = 'fALFF';
measure(3).dir_name = 'ALFF';
measure(3).output_dir = 'A-fALFF';         
measure(3).seed_name = '';
measure(3).subbrick = 1;
measure(3).do = opt.MEASURES.ALFF;
measure(4).name_selector = 'RSFA';
measure(4).name = 'RSFA';
measure(4).dir_name = 'ALFF';
measure(4).output_dir = 'A-RSFA';         
measure(4).seed_name = '';
measure(4).subbrick = 1;
measure(4).do = opt.MEASURES.RSFA;
measure(5).name_selector = 'mRSFA';
measure(5).name = 'mRSFA';
measure(5).dir_name = 'ALFF';
measure(5).output_dir = 'A-mRSFA';         
measure(5).seed_name = '';
measure(5).subbrick = 1;
measure(5).do = opt.MEASURES.RSFA;
measure(6).name_selector = 'fRSFA';
measure(6).name = 'fRSFA';
measure(6).dir_name = 'ALFF';
measure(6).output_dir = 'A-fRSFA';         
measure(6).seed_name = '';
measure(6).subbrick = 1;
measure(6).do = opt.MEASURES.RSFA;
measure(7).name_selector = 'zALFF';          %for selection of the measure
measure(7).name = 'zALFF';                   %for selection of the file
measure(7).dir_name = 'ALFF';               %for selection of the firstlevel folder
measure(7).output_dir = 'A-zALFF';             %output dir of the second level analysis
measure(7).seed_name = '';                  %for selection of the seed in the case of seedmeasure
measure(7).subbrick = 1;
measure(7).do = opt.MEASURES.ALFF;
measure(8).name_selector = 'mALFF-1';
measure(8).name = 'mALFF-1';
measure(8).dir_name = 'ALFF';
measure(8).output_dir = 'A-mALFF-1';         
measure(8).seed_name = '';         
measure(8).subbrick = 1;
measure(8).do = opt.MEASURES.ALFF;
measure(9).name_selector = 'zfALFF';
measure(9).name = 'zfALFF';
measure(9).dir_name = 'ALFF';
measure(9).output_dir = 'A-zfALFF';         
measure(9).seed_name = '';
measure(9).subbrick = 1;
measure(9).do = opt.MEASURES.ALFF;
measure(10).name_selector = 'zRSFA';
measure(10).name = 'zRSFA';
measure(10).dir_name = 'ALFF';
measure(10).output_dir = 'A-zRSFA';         
measure(10).seed_name = '';
measure(10).subbrick = 1;
measure(10).do = opt.MEASURES.RSFA;
measure(11).name_selector = 'mRSFA-1';
measure(11).name = 'mRSFA-1';
measure(11).dir_name = 'ALFF';
measure(11).output_dir = 'A-mRSFA-1';         
measure(11).seed_name = ''; 
measure(11).subbrick = 1;
measure(11).do = opt.MEASURES.RSFA;
measure(12).name_selector = 'zfRSFA';
measure(12).name = 'zfRSFA';
measure(12).dir_name = 'ALFF';
measure(12).output_dir = 'A-zfRSFA';         
measure(12).seed_name = '';
measure(12).subbrick = 1;
measure(12).do = opt.MEASURES.RSFA;
%REHO
measure(13).name_selector = 'REHO7';
measure(13).name = 'REHO7';
measure(13).dir_name = 'ReHo';
measure(13).output_dir = 'REHO7';         
measure(13).seed_name = '';       
measure(13).subbrick = 1;
measure(13).do = 1 * opt.MEASURES.REHO;
measure(14).name_selector = 'REHO19';
measure(14).name = 'REHO19';
measure(14).dir_name = 'ReHo';
measure(14).output_dir = 'REHO19';         
measure(14).seed_name = '';        
measure(14).subbrick = 1;
measure(14).do = 1 * opt.MEASURES.REHO;
measure(15).name_selector = 'REHO27';
measure(15).name = 'REHO27';
measure(15).dir_name = 'ReHo';
measure(15).output_dir = 'REHO27';         
measure(15).seed_name = '';         
measure(15).subbrick = 1;
measure(15).do = 1 * opt.MEASURES.REHO;
%SpEn
measure(16).name_selector = 'SpEn';
measure(16).name = 'SpEn';
measure(16).dir_name = 'SpEn';
measure(16).output_dir = 'SpEn';         
measure(16).seed_name = '';         
measure(16).subbrick = 1;
measure(16).do = 1 * opt.MEASURES.SpE;
% CoEn
measure(17).name_selector = 'zCoEn';
measure(17).name = 'zCoEn';
measure(17).dir_name = 'CoEn';
measure(17).output_dir = 'zCoEn';         
measure(17).seed_name = '';         
measure(17).subbrick = 1;
measure(17).do = 1 * opt.MEASURES.CoE;
measure(18).name_selector = 'CoEn';
measure(18).name = 'CoEn';
measure(18).dir_name = 'CoEn';
measure(18).output_dir = 'CoEn';         
measure(18).seed_name = '';         
measure(18).subbrick = 1;
measure(18).do = 1 * opt.MEASURES.CoE;
% special ALFF welch
measure(19).name_selector = 'mALFFw';         
measure(19).name = 'mALFFw';             
measure(19).dir_name = 'ALFFW';       
measure(19).output_dir = 'Aw-mALFFw';        
measure(19).seed_name = '';          
measure(19).subbrick = 1;
measure(19).do = 1 * opt.MEASURES.ALFFW;
measure(20).name_selector = 'ALFFw';           
measure(20).name = 'ALFFw';             
measure(20).dir_name = 'ALFFW';       
measure(20).output_dir = 'Aw-ALFFw';      
measure(20).seed_name = '';          
measure(20).subbrick = 1;
measure(20).do = 1 * opt.MEASURES.ALFFW;
measure(21).name_selector = 'mALFFw-1';         
measure(21).name = 'mALFFw-1';             
measure(21).dir_name = 'ALFFW';       
measure(21).output_dir = 'Aw-mALFFw-1';        
measure(21).seed_name = '';     
measure(21).subbrick = 1;
measure(21).do = 1 * opt.MEASURES.ALFFW;
measure(22).name_selector = 'zALFFw';           
measure(22).name = 'zALFFw';             
measure(22).dir_name = 'ALFFW';       
measure(22).output_dir = 'Aw-zALFFw';      
measure(22).seed_name = '';     
measure(22).subbrick = 1;
measure(22).do = 1 * opt.MEASURES.ALFFW;
% special ALFF multitaper
measure(23).name_selector = 'mALFFm';      
measure(23).name = 'mALFFm';           
measure(23).dir_name = 'ALFFM';      
measure(23).output_dir = 'Am-mALFFm';  
measure(23).seed_name = '';         
measure(23).subbrick = 1;
measure(23).do = 1 * opt.MEASURES.ALFFM;
measure(24).name_selector = 'ALFFm';           
measure(24).name = 'ALFFm';              
measure(24).dir_name = 'ALFFM';       
measure(24).output_dir = 'Am-ALFFm';         
measure(24).seed_name = '';          
measure(24).subbrick = 1;
measure(24).do = 1 * opt.MEASURES.ALFFM;
measure(25).name_selector = 'mALFFm-1';      
measure(25).name = 'mALFFm-1';           
measure(25).dir_name = 'ALFFM';      
measure(25).output_dir = 'Am-mALFFm-1';  
measure(25).seed_name = '';         
measure(25).subbrick = 1;
measure(25).do = 1 * opt.MEASURES.ALFFM;
measure(26).name_selector = 'zALFFm';           
measure(26).name = 'zALFFm';              
measure(26).dir_name = 'ALFFM';       
measure(26).output_dir = 'Am-zALFFm';         
measure(26).seed_name = '';        
measure(26).subbrick = 1;
measure(26).do = 1 * opt.MEASURES.ALFFM;
%REHO_normalized
measure(27).name_selector = 'zREHO7';
measure(27).name = 'zREHO7';
measure(27).dir_name = 'ReHo';
measure(27).output_dir = 'zREHO7';         
measure(27).seed_name = '';       
measure(27).subbrick = 1;
measure(27).do = 1 * opt.MEASURES.REHO;
measure(28).name_selector = 'zREHO19';
measure(28).name = 'zREHO19';
measure(28).dir_name = 'ReHo';
measure(28).output_dir = 'zREHO19';         
measure(28).seed_name = '';        
measure(28).subbrick = 1;
measure(28).do = 1 * opt.MEASURES.REHO;
measure(29).name_selector = 'zREHO27';
measure(29).name = 'zREHO27';
measure(29).dir_name = 'ReHo';
measure(29).output_dir = 'zREHO27';         
measure(29).seed_name = '';         
measure(29).subbrick = 1;
measure(29).do = 1 * opt.MEASURES.REHO;

last_m = length(measure);

if opt.MEASURES.DEL
    indx = 0;
    roi_n = size(opt.AUX.DEL.roi_names,1);
    for k = (last_m+1): (last_m + roi_n)
        indx = indx + 1;
        roi_name = c_b(opt.AUX.DEL.roi_names(indx,:));
        roi_name = roi_name(1:end-4);
        str = ['DEL_',roi_name];
        measure(k).name_selector = 'DELAY';
        measure(k).name = str;
        measure(k).dir_name = 'DELAY';
        measure(k).output_dir = str; 
        measure(k).seed_name = roi_name;    
        measure(k).subbrick = 1;
        measure(k).do = 1 * opt.MEASURES.DEL;
    end
end
if opt.MEASURES.StoV
    last_m = length(measure);
    indx = 0;
    roi_n = size(opt.AUX.StoV.roi_names,1);
    for k = (last_m+1): (last_m + roi_n)
        indx = indx + 1;
        roi_name = c_b(opt.AUX.StoV.roi_names(indx,:));
        roi_name = roi_name(1:end-4);
        str = ['StV_',roi_name];
        measure(k).name_selector = 'StV';
        measure(k).name = [str,'_z'];
        measure(k).dir_name = 'StV';
        measure(k).output_dir = str; 
        measure(k).seed_name = roi_name;   
        measure(k).subbrick = 1;
        measure(k).do = 1 * opt.MEASURES.StoV;
    end
end
if opt.MEASURES.StoV && strcmp(opt.prepro_mode,'3dTproject') % only in this modality is produced the z-score of StV
    last_m = length(measure);
    indx = 0;
    roi_n = size(opt.AUX.StoV.roi_names,1);
    for k = (last_m+1): (last_m + roi_n)
        indx = indx + 1;
        roi_name = c_b(opt.AUX.StoV.roi_names(indx,:));
        roi_name = roi_name(1:end-4);
        str = ['StV_',roi_name];
        measure(k).name_selector = 'ZStV';
        measure(k).name = [str,'_Z'];
        measure(k).dir_name = 'StV';
        measure(k).output_dir = ['Zscore_',str]; 
        measure(k).seed_name = roi_name;   
        measure(k).subbrick = 1;
        measure(k).do = 1 * opt.MEASURES.StoV;
    end
end
if opt.MEASURES.COH
    %MABA
    n_subband = size(opt.coherency.bands,1);
    roi_n = size(opt.AUX.COH.roi_names,1);
    for z = 1:n_subband
        last_m = length(measure);
        indx = 0;
        for k = (last_m+1): (last_m + roi_n)
            indx = indx + 1;
            roi_name = c_b(opt.AUX.COH.roi_names(indx,:));
            roi_name = roi_name(1:end-4);
            str = ['COH_',roi_name,'_MABA'];
            measure(k).name_selector = 'COH_M';
            measure(k).name = [str,'_z'];
            measure(k).dir_name = 'COH';
            measure(k).output_dir = [str,'_',num2str(opt.coherency.bands(z,1)),'_',num2str(opt.coherency.bands(z,2))]; 
            measure(k).seed_name = roi_name;   
            measure(k).subbrick = z;
            measure(k).do = 1 * opt.MEASURES.COH;
        end
    end
    %DEBA
    n_subband = size(opt.coherency.bands,1);
    roi_n = size(opt.AUX.COH.roi_names,1);
    for z = 1:n_subband
        last_m = length(measure);
        indx = 0;
        for k = (last_m+1): (last_m + roi_n)
            indx = indx + 1;
            roi_name = c_b(opt.AUX.COH.roi_names(indx,:));
            roi_name = roi_name(1:end-4);
            str = ['COH_',roi_name,'_DEBA'];
            measure(k).name_selector = 'COH_D';
            measure(k).name = str;
            measure(k).dir_name = 'COH';
            measure(k).output_dir = [str,'_',num2str(opt.coherency.bands(z,1)),'_',num2str(opt.coherency.bands(z,2))]; 
            measure(k).seed_name = roi_name;   
            measure(k).subbrick = z;
            measure(k).do = 1 * opt.MEASURES.COH;
        end
    end
end
if opt.MEASURES.gFC
    last_m = length(measure);
    measure(last_m+1).name_selector = 'gFC_Zz';         
    measure(last_m+1).name = 'gFC_Zz';                   
    measure(last_m+1).dir_name = 'gFC';               
    measure(last_m+1).output_dir = 'gFC_Zz';            
    measure(last_m+1).seed_name = '';                 
    measure(last_m+1).subbrick = 1;
    measure(last_m+1).do = 1 * opt.MEASURES.gFC;
    measure(last_m+2).name_selector = 'gFC_Zr';         
    measure(last_m+2).name = 'gFC_Zr';                   
    measure(last_m+2).dir_name = 'gFC';               
    measure(last_m+2).output_dir = 'gFC_Zr';            
    measure(last_m+2).seed_name = '';                 
    measure(last_m+2).subbrick = 1;
    measure(last_m+2).do = 1 * opt.MEASURES.gFC;
    measure(last_m+3).name_selector = 'gFC_P';         
    measure(last_m+3).name = 'gFC_P';                   
    measure(last_m+3).dir_name = 'gFC';               
    measure(last_m+3).output_dir = 'gFC_P';            
    measure(last_m+3).seed_name = '';                 
    measure(last_m+3).subbrick = 1;
    measure(last_m+3).do = 1 * opt.MEASURES.gFC;
    last_m = length(measure);
    for k = (last_m+1): (last_m + opt.AUX.gFC.num_thr)
        measure(k).name_selector = 'gFC_VT'; 
        measure(k).name = 'gFC_VT';
        measure(k).dir_name = 'gFC'; 
        measure(k).output_dir = ['gFC_VT_',num2str((opt.AUX.gFC.Vthr(1)+opt.AUX.gFC.Vthr(3)*(k-last_m -1)),2)];
        measure(k).seed_name = '';   
        measure(k).subbrick = k-last_m;
        measure(k).do = 1 * opt.MEASURES.gFC;
    end
end
if opt.MEASURES.DeC
    last_m = length(measure);
    for k = (last_m+1): (last_m + opt.AUX.gFC.num_thr)
        measure(k).name_selector = 'DeC'; 
        measure(k).name = 'DeC';
        measure(k).dir_name = 'DeC'; 
        measure(k).output_dir = ['DeC_',num2str(opt.AUX.DeC.Thrs_vector(k-last_m),2)];
        measure(k).seed_name = '';   
        measure(k).subbrick = k-last_m;
        measure(k).do = 1 * opt.MEASURES.DeC;
    end
    last_m = length(measure);
    for k = (last_m+1): (last_m + opt.AUX.gFC.num_thr)
        measure(k).name_selector = 'zDeC'; 
        measure(k).name = 'zDeC';
        measure(k).dir_name = 'DeC'; 
        measure(k).output_dir = ['zDeC_',num2str(opt.AUX.DeC.Thrs_vector(k-last_m),2)];
        measure(k).seed_name = '';   
        measure(k).subbrick = k-last_m;
        measure(k).do = 1 * opt.MEASURES.DeC;
    end
end

return
end

