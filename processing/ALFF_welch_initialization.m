function ALFF_welch_initialization
global opt
freq_max = opt.filter_band(2);
freq_min = opt.filter_band(1);
opt.AUX.alffw.lowcutoff = find(opt.psdw.freq >(freq_min - opt.psdw.deltafreq/2) & opt.psdw.freq < (freq_min + opt.psdw.deltafreq/2));   
opt.AUX.alffw.highcutoff= find(opt.psdw.freq >(freq_max - opt.psdw.deltafreq/2) & opt.psdw.freq < (freq_max + opt.psdw.deltafreq/2));
return
end