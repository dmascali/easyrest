function ALFF_multitaper_initialization
global opt
freq_max = opt.filter_band(2);
freq_min = opt.filter_band(1);
opt.AUX.alffm.lowcutoff = find(opt.psdm.freq >(freq_min - opt.psdm.deltafreq/2) & opt.psdm.freq < (freq_min + opt.psdm.deltafreq/2));   
opt.AUX.alffm.highcutoff= find(opt.psdm.freq >(freq_max - opt.psdm.deltafreq/2) & opt.psdm.freq < (freq_max + opt.psdm.deltafreq/2));
return
end