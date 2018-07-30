function ALFF_welch(psdw_name,bmask,output_name)
global opt
hdr = spm_vol(psdw_name);
psdw = spm_read_vols(hdr);
psdw = sqrt(2*psdw*(1/opt.tr)/(opt.psdw.ENBW*opt.N));           % se la finestra è rettangolare (ie, no window) c = 1/N -> cN = 1; Ricorda 2Amp = sqrt(psdw*2*Fs*N).  Ma ALFF = sqrt(2Amp.^2/N)
alffw = sum(psdw(:,:,:,opt.AUX.alffw.lowcutoff:opt.AUX.alffw.highcutoff),4);        %afni somma su tutto, tanto è filtrato! Se si usa windowing conviene limitarsi alla banda
normalization = mean(alffw(bmask == 1));
malffw = alffw./normalization;
%save maps
hdr = hdr(1);
hdr = rmfield(hdr,'pinfo');
alffw_name = [output_name,'ALFFw.nii'];
malffw_name = [output_name,'mALFFw.nii'];
hdr.fname = alffw_name;
hdr.private.dat.fname = alffw_name;
spm_write_vol(hdr,alffw);
hdr.fname = malffw_name;
hdr.private.dat.fname = malffw_name;
spm_write_vol(hdr,malffw);
%produce zscore a m-1 maps
Zalffw_name = [output_name,'zALFFw.nii'];
Mmalffw_name = [output_name,'mALFFw-1.nii'];
sd = std(alffw(bmask == 1));
Zalffw = (alffw-normalization)./sd;
Malffw = malffw-double(bmask);
hdr.fname = Zalffw_name;
hdr.private.dat.fname = Zalffw_name;
spm_write_vol(hdr,Zalffw);
hdr.fname = Mmalffw_name;
hdr.private.dat.fname = Mmalffw_name;
spm_write_vol(hdr,Malffw);
return
end