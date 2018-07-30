function ALFF_multitaper(psdm_name,bmask,output_name)
global opt
hdr = spm_vol(psdm_name);
psdm = spm_read_vols(hdr);
psdm = sqrt(2*psdm*(1/opt.tr)/(opt.psdm.ENBW*opt.N));           % se la finestra è rettangolare (ie, no window) c = 1/N -> cN = 1; Ricorda 2Amp = sqrt(psdw*2*Fs*N).  Ma ALFF = sqrt(2Amp.^2/N)
alffm = sum(psdm(:,:,:,opt.AUX.alffm.lowcutoff:opt.AUX.alffm.highcutoff),4);        %afni somma su tutto, tanto è filtrato! Se si usa windowing conviene limitarsi alla banda
normalization = mean(alffm(bmask == 1));
malffm = alffm./normalization;
%save maps
hdr = hdr(1);
hdr = rmfield(hdr,'pinfo');
alffm_name = [output_name,'ALFFm.nii'];
malffm_name = [output_name,'mALFFm.nii'];
hdr.fname = alffm_name;
hdr.private.dat.fname = alffm_name;
spm_write_vol(hdr,alffm);
hdr.fname = malffm_name;
hdr.private.dat.fname = malffm_name;
spm_write_vol(hdr,malffm);
%produce zscore a m-1 maps
Zalffm_name = [output_name,'zALFFm.nii'];
Mmalffm_name = [output_name,'mALFFm-1.nii'];
sd = std(alffm(bmask == 1));
Zalffm = (alffm-normalization)./sd;
Malffm = malffm-double(bmask);
hdr.fname = Zalffm_name;
hdr.private.dat.fname = Zalffm_name;
spm_write_vol(hdr,Zalffm);
hdr.fname = Mmalffm_name;
hdr.private.dat.fname = Mmalffm_name;
spm_write_vol(hdr,Malffm);
return
end