function psd_welch(epi_name,output_name)
global opt

hdr_epi = spm_vol(epi_name);
epi = spm_read_vols(hdr_epi);

s = size(epi);  % Image Dimensions (Voxels)

WS = (1/opt.tr);

epi = reshape(epi,[s(1)*s(2)*s(3),s(4)])';

[fee,W] = cpsd_R2014b(epi,epi,opt.psdw.win,opt.psdw.NO,opt.psdw.NFFT,WS); 

FFT_L = length(W);  %acquiring FFT length from cpsd

% reshape
fee = reshape(fee',[s(1),s(2),s(3),FFT_L]);

hdr_epi = rmfield(hdr_epi,'pinfo');
hdr_out = hdr_epi(1:FFT_L);

for l = 1:size(fee,4)
    hdr_out(l).fname = output_name;
    hdr_out(l).private.dat.fname = output_name;
    spm_write_vol(hdr_out(l),fee(:,:,:,l));
end

return
end