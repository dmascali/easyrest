function rp_define_ordering_and_unit
global opt
switch lower(opt.AUX.rp_tye)
    case {'fsl'} %ok
        opt.AUX.rp_indx.rot = [1 2 3];
        opt.AUX.rp_indx.tra = [4 5 6];  
        opt.AUX.rp_rot_unit = 'rad';
    case {'spm'} %should be ok. Indeed they are plotted by SPM in degrees but they are saved in radians
        opt.AUX.rp_indx.rot = [4 5 6];
        opt.AUX.rp_indx.tra = [1 2 3];         
        opt.AUX.rp_rot_unit = 'rad';
    case {'afni'} %ok
        opt.AUX.rp_indx.rot = [1 2 3];
        opt.AUX.rp_indx.tra = [4 5 6];         
        opt.AUX.rp_rot_unit = 'deg';
    case {'hcp'} %verified
        opt.AUX.rp_indx.rot = [4 5 6];
        opt.AUX.rp_indx.tra = [1 2 3];  
        opt.AUX.rp_rot_unit = 'deg';
end
opt.AUX.rp_average_brain_radius = 50; %mm; 50 according to Power, 60 according to AFNI (see 1d_tool.py)
return
end

