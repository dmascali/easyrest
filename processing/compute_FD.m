function compute_FD(rp,g,j,k) %group, session, subj
%FD Frame Displacement (see Power et al, NeuroImage, 59(3), 2012)
global opt
%rp = load(c_b(opt.DATA.RP_1D{g,j,k}));

tra = rp(:,opt.AUX.rp_indx.tra);
rot = rp(:,opt.AUX.rp_indx.rot);

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
%------------------ Volume selector, rp----------------------
nrp = size(rp,1);
indx = volume_selector_matlab(nrp,g,j)';
if opt.AUX.v_s_do  
    fd(indx) = [];
end
%------------------------------------------------------------
opt.QC.FD{g,j,k} = fd;
opt.QC.FD_mean{g,j}(k,1) = mean(fd);
return
end