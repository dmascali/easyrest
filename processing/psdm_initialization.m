function psdm_initialization(ot)
global opt
% problema con opt.N , qui richiediamo che tutti i dati abbiamo lo stesso
% numero di volumi

nw_def = 3;     %change here the default nw

if exist('pmtm_R2014b')~=2
    error('Improper ER installation. Add ER to your matlab path with all subfolders, i.e.: addpath(genpath(''./ER''))');
end

% if isempty(opt.N)   DELETE AS U CAN
%     [~,N] = system(['3dinfo -nv ',c_b(opt.DATA.FUNCTIONALS{1,1}(1,:))]);        %attention: assuming that all scans have the same number of time points
%     opt.N = str2num(N);
% end


opt.psdm.dpss_nw = [];
if isfield(ot,'psdm')
    if isfield(ot.psdm,'dpss_nw') && ~isempty(ot.psdm.dpss_nw) 
       opt.psdm.dpss_nw = ot.psdm.dpss_nw;
        if opt.psdm.dpss_nw < 0.75
           error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '0.75', 'Droplasttaper', 'false'));    
        end
    end
end

if isempty(opt.psdm.dpss_nw)    %use default
    opt.psdm.dpss_nw = nw_def;    
end

opt.psdm.dpss_k = [];
if isfield(ot,'psdm')
    if isfield(ot.psdm,'dpss_k') && ~isempty(ot.psdm.dpss_k) 
       opt.psdm.dpss_k = ot.psdm.dpss_k;
    end
end

if isempty(opt.psdm.dpss_k)     
    [E,V] = dpss(opt.N,opt.psdm.dpss_nw); %use default k
    % now lets drop the last component (k = 2nw-1)
    numvec = length(V);
    if numvec > 2
       E = E(:,1:numvec-1);
       V = V(1:numvec-1);
    else
       error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
    end
else
    [E,V] = dpss(opt.N,opt.psdm.dpss_nw,opt.psdm.dpss_k);
end

opt.psdm.dpss_k = length(V);
opt.psdm.dpss_E = E;
opt.psdm.dpss_V = V;

opt.psdm.NFFT=max(256,2^nextpow2(opt.N));          %uso il default di pmtm
% if rem(opt.N,2)
%     opt.psdm.NFFT = opt.N+1;
% else
%     opt.psdm.NFFT = opt.N;
% end

fprintf('\n_______________________________________________\n');
fprintf('\tPSD parameters (multitaper method)\n');
fprintf('NW: \t\t\t%.0d\n',opt.psdm.dpss_nw);
fprintf('K (# dpss): \t\t%.0d\n',opt.psdm.dpss_k);
fprintf('NFFT segment: \t\t%.0d\n',opt.psdm.NFFT);
fprintf('_______________________________________________\n');


opt.psdm.freq = opt.AUX.fs/2*linspace(0,1,opt.psdm.NFFT/2+1);      
opt.psdm.deltafreq = opt.psdm.freq(2) - opt.psdm.freq(1);
c = 1/mean((sum(E).^2));
opt.psdm.ENBW = c;      %equivalent noise bandwith

return
end