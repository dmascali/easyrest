function psdw_initialization(ot)
global opt

if exist('cpsd_R2014b')~=2
    error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''./ER''))');
end

% if isempty(opt.N)   DELETE AS U CAN
%     [~,N] = system(['3dinfo -nv ',c_b(opt.DATA.FUNCTIONALS{1,1}(1,:))]);        %attention: assuming that all scans have the same number of time points
%     opt.N = str2num(N);
% end


opt.psdw.HW_length = [];
if isfield(ot,'psdw')
    if isfield(ot.psdw,'window_length') && ~isempty(ot.psdw.window_length) 
       opt.psdw.HW_length = ot.psdw.window_length;
       if opt.psdw.HW_length > opt.N
           error('Window for Welch''s method segmentation can''t be longer than epi time serie.');
       end
    end
end

if isempty(opt.psdw.HW_length)
    opt.psdw.HW_length = fix(1/(opt.tr*(opt.filter_band(1))));    
    if opt.psdw.HW_length > opt.N
        error('The default value of window''s length for Welch''s method segmentation is too longer. Please define it manually via psdw.window_length.');
    end
end

opt.psdw.NO = fix(opt.psdw.HW_length/2);      % 50% overlap
opt.psdw.N_seg = (opt.N-opt.psdw.NO)./(opt.psdw.HW_length-opt.psdw.NO);

war_integer = 0; if fix(opt.psdw.N_seg) ~= opt.psdw.N_seg; war_integer = 1; end
opt.psdw.N_seg = fix(opt.psdw.N_seg);

opt.psdw.NFFT=max(256,2^nextpow2(opt.psdw.HW_length));          %uso il default di cpsd. Praticamente sempre verr� usato nfft = 256; 
% if rem(opt.psdw.HW_length,2)
%     opt.psdw.NFFT = opt.psdw.HW_length+1;
% else
%     opt.psdw.NFFT = opt.psdw.HW_length;
% end

fprintf('\n_______________________________________________\n');
fprintf('\tPSD parameters (Welch''s method)\n');
fprintf('Segment length: \t%.0d\n',opt.psdw.HW_length);
fprintf('Overlap length: \t%.0d\n',opt.psdw.NO);
fprintf('Averaged segments: \t%.0d\n',opt.psdw.N_seg);
fprintf('NFFT segment: \t\t%.0d\n',opt.psdw.NFFT);
warning off backtrace
if opt.psdw.N_seg < 3
    warning('The number of averaged segments is low.');
end
if war_integer
    warning('signal:welchparse:MustBeInteger','The number of segments was not an integer, data will be truncated.');
end
warning on backtrace
fprintf('_______________________________________________\n');

opt.psdw.freq = opt.AUX.fs/2*linspace(0,1,opt.psdw.NFFT/2+1);      
opt.psdw.deltafreq = opt.psdw.freq(2) - opt.psdw.freq(1);
opt.psdw.win = hanning(opt.psdw.HW_length,'periodic');          %change here the window type. You could also use a rectwin (but what will happens at the edges?)
s1 = sum(opt.psdw.win);
s2 = sum(opt.psdw.win.^2);
c = s2/s1^2;
opt.psdw.ENBW = c;      %equivalent noise bandwith

return
end