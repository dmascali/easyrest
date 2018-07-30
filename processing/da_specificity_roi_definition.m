function da_specificity_roi_definition

global opt
radius = 6; % mm
base_name = 'ER_DA_specificity_';

[~,~]= system('rm -r ER_DA_specificity_*.nii');

%{
 %OLD ROIS
spe(1).roi(1).center = [0,-53,26];    %definition of Liu2013 (pnas,110(11):4392-4397) (see also: Mao2015, Amico2014). Alteratively we can use the old definition from fox2005: Tal=[-5,-49,40] in MNI is around = [-6,-52,40]
spe(1).roi(1).name = [base_name,'PCC.nii'];
spe(1).roi(1).name_contract = 'PCC';
spe(1).roi(2).center = [-1,49,-2];    %from Fox2005 converted in MNI by Chai2012, Neuroimage 1420-1428
spe(1).roi(2).name = [base_name,'MPFC.nii'];
spe(1).roi(2).name_contract = 'MPFC';
spe(1).roi(3).center = [-30,-88,0];   %visual cortex, from Chai2012, Neuroimage 1420-1428 taken from Van Dijk 2010
spe(1).roi(3).name = [base_name,'VIS_L.nii'];
spe(1).roi(3).name_contract = 'VIS_L';
spe(1).roi(4).center = [30,-88,0];    %visual cortex, from Chai2012, Neuroimage 1420-1428 taken from Van Dijk 2010
spe(1).roi(4).name = [base_name,'VIS_R.nii'];
spe(1).roi(4).name_contract = 'VIS_R';
%}

% I try to decenter the ROIs to avoid CSF signal
% DMN test     
spe(1).name = 'DMN';
spe(1).roi(1).center = [0,-52,27];    %Raichel 2011
spe(1).roi(1).name = [base_name,spe(1).name,'_PCC.nii'];
spe(1).roi(1).name_contract = 'PCC';
spe(1).roi(2).center = [0,54,6]; %[-1,54,27];    %Raichel 2011
spe(1).roi(2).name = [base_name,spe(1).name,'_MPFC.nii'];
spe(1).roi(2).name_contract = 'MPFC';
spe(1).roi(3).center = [30,-88,0]; % before was[30,-88,0]; % before was [7, -83, 2];   %Raichel 2011But there is too much signal in V1
spe(1).roi(3).name = [base_name,spe(1).name,'_VIS.nii'];
spe(1).roi(3).name_contract = 'VIS';

% AUD test     
spe(2).name = 'AUD';
spe(2).roi(1).center = [-38,-28, 12];    %HESCHL's gyrus
spe(2).roi(1).name = [base_name,spe(2).name,'_LAU.nii'];
spe(2).roi(1).name_contract = 'LHesc';
spe(2).roi(2).center = [40,-26, 12];    
spe(2).roi(2).name = [base_name,spe(2).name,'_RAU.nii'];
spe(2).roi(2).name_contract = 'RHesc';
spe(2).roi(3).center = [30,-88,0];   
spe(2).roi(3).name = [base_name,spe(2).name,'_VIS.nii'];
spe(2).roi(3).name_contract = 'VIS';

% SMN test     
spe(3).name = 'SMN';
spe(3).roi(1).center = [-48,-26,48];    %Postecentral gyrus
spe(3).roi(1).name = [base_name,spe(3).name,'_LMC.nii'];
spe(3).roi(1).name_contract = 'LPcg';
spe(3).roi(2).center = [48,-26,48];    %Raichel 2011
spe(3).roi(2).name = [base_name,spe(3).name,'_RMC.nii'];
spe(3).roi(2).name_contract = 'RPcg';
spe(3).roi(3).center = [30,-88,0];   
spe(3).roi(3).name = [base_name,spe(3).name,'_VIS.nii'];
spe(3).roi(3).name_contract = 'VIS';


% % LOW SIGNAL test     
% spe(4).name = 'ITG';
% spe(4).roi(1).center = [-54,-24,-26];    
% spe(4).roi(1).name = [base_name,spe(4).name,'_LITG.nii'];
% spe(4).roi(1).name_contract = 'LITG';
% spe(4).roi(2).center = [54,-24,-26];  
% spe(4).roi(2).name = [base_name,spe(4).name,'_RITG.nii'];
% spe(4).roi(2).name_contract = 'RITG';
% spe(4).roi(3).center = [30,-88,0];   
% spe(4).roi(3).name = [base_name,spe(4).name,'_VIS.nii'];
% spe(4).roi(3).name_contract = 'VIS';
% 
% % LOW SIGNAL test     
% spe(5).name = 'FMC';
% spe(5).roi(1).center = [-8,40,-18];    
% spe(5).roi(1).name = [base_name,spe(5).name,'_LFMC.nii'];
% spe(5).roi(1).name_contract = 'LFMC';
% spe(5).roi(2).center = [8,40,-18];  
% spe(5).roi(2).name = [base_name,spe(5).name,'_RFMC.nii'];
% spe(5).roi(2).name_contract = 'RFMC';
% spe(5).roi(3).center = [30,-88,0];   
% spe(5).roi(3).name = [base_name,spe(5).name,'_VIS.nii'];
% spe(5).roi(3).name_contract = 'VIS';
% 
% % Low-T2* test    
% spe(6).name = 'Putamen';
% spe(6).roi(1).center = [-28,-4,0];    %Raichel 2011
% spe(6).roi(1).name = [base_name,spe(6).name,'_L.nii'];
% spe(6).roi(1).name_contract = 'LPut';
% spe(6).roi(2).center = [26,-2,4];    %Raichel 2011
% spe(6).roi(2).name = [base_name,spe(6).name,'_R.nii'];
% spe(6).roi(2).name_contract = 'RPut';
% spe(6).roi(3).center = [30,-88,0];   
% spe(6).roi(3).name = [base_name,spe(6).name,'_VIS.nii'];
% spe(6).roi(3).name_contract = 'VIS';

spe(4).name = 'DMNold';
spe(4).roi(1).center = [0,-52,27];    %Raichel 2011
spe(4).roi(1).name = [base_name,spe(4).name,'_PCC.nii'];
spe(4).roi(1).name_contract = 'PCC';
spe(4).roi(2).center = [-1,54,27]; %[-1,54,27];    %Raichel 2011
spe(4).roi(2).name = [base_name,spe(4).name,'_MPFC.nii'];
spe(4).roi(2).name_contract = 'MPFC';
spe(4).roi(3).center = [30,-88,0]; % before was[30,-88,0]; % before was [7, -83, 2];   %Raichel 2011But there is too much signal in V1
spe(4).roi(3).name = [base_name,spe(4).name,'_VIS.nii'];
spe(4).roi(3).name_contract = 'VIS';

% % LOW SIGNAL test     
% spe(8).name = 'DMN_FMC';
% spe(8).roi(1).center = [-8,40,-18];    
% spe(8).roi(1).name = [base_name,spe(8).name,'_LFMC.nii'];
% spe(8).roi(1).name_contract = 'LFMC';
% spe(8).roi(2).center = [0,-52,27];  
% spe(8).roi(2).name = [base_name,spe(8).name,'_PCC.nii'];
% spe(8).roi(2).name_contract = 'RFMC';
% spe(8).roi(3).center = [30,-88,0];   
% spe(8).roi(3).name = [base_name,spe(8).name,'_VIS.nii'];
% spe(8).roi(3).name_contract = 'VIS';




for s = 1:length(spe)
    for l = 1:length(spe(s).roi)
        c = sphere_eq_string(spe(s).roi(l).center);
        %Data must all have the same spatial resolution
        [status,string] = system(['3dcalc -a ',c_b(opt.DATA.FUNCTIONALS{1,1}(1,:)),'''[1]'' -LPI -expr ''step(',num2str(radius*radius),'-(x',c{1},')*(x',c{1},')-(y',c{2},')*(y',c{2},')-(z',c{3},')*(z',c{3},'))'' -prefix ',spe(s).roi(l).name]); control(status,string);
        spe(s).roi(l).img = spm_read_vols(spm_vol(spe(s).roi(l).name));
    end
end

opt.DA.roi_specificity = spe;

return
end

function [c_str] =  sphere_eq_string(c)

c_str = cell(3,1);
for l = 1:3
    if c(l) >= 0
        c_str{l} = ['-',num2str(c(l))];
    else
        c_str{l} = ['+',num2str(abs(c(l)))];
    end
end
return
end

