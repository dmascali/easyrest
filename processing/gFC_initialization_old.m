function gFC_initialization_old(ot,CoE)
%%%%% NO USED FUNCTION
%this is an old version of gFC_initialization. In this version is possible
%to choose to output HIST from gFC
global opt 
%set default
opt.AUX.gFC.deleteHIST = 0;
if ~isfield(ot,'gFC')
    %ot.gFC.thr = 0.25;
    ot.gFC.Vthr = [0.25 0.5 0.05];
    if CoE
        ot.gFC.hist = 1;
    else
        ot.gFC.hist = 0;
    end
else
%     if ~isfield(ot.gFC,'thr');
%         ot.gFC.thr = 0.25;
%     else
%         if numel(ot.gFC.thr) > 1
%             error('Too many input in ER.others.gFC.thr . Specify only a number');
%         end
%         if ot.gFC.thr >= 1
%             error('ER.others.gFC.thr must be < 1 [default is 0.25].');
%         end
%     end
    if ~isfield(ot.gFC,'Vthr'); 
        ot.gFC.Vthr = [0.25 0.5 0.05];
    else
        if numel(ot.gFC.Vthr) ~= 3
            error('Bad ER.others.gFC.Vthr spefication. Specify 3 numbers: t0 t1 and dt (default is [0.25 0.5 0.05]).');
        end
        if sum(ot.gFC.Vthr <= 1) ~= 3 
            error('ER.others.gFC.Vthr must be < 1 (default is [0.25 0.5 0.05]).');
        elseif ot.gFC.Vthr(1) > ot.gFC.Vthr(2)
            error('Bad definition of ER.others.gFC.Vthr (default is [0.25 0.5 0.05]).');
        end
    end
    if ~isfield(ot.gFC,'hist');
        if CoE
            ot.gFC.hist = 1;
        else
            ot.gFC.hist = 0;
        end
    else
        if CoE 
            if ~ot.gFC.hist
                ot.gFC.hist = 1;
                opt.AUX.gFC.deleteHIST = 1;
            end
        end
        if ~(ot.gFC.hist == 1 || ot.gFC.hist == 0)
            error('ER.others.gFC.hist is a logical variable. [default is 0].');
        end
    end
end
% opt.AUX.gFC.thr = ot.gFC.thr;
opt.AUX.gFC.Vthr = ot.gFC.Vthr;
opt.AUX.gFC.hist = ot.gFC.hist;
opt.AUX.gFC.num_thr = floor((opt.AUX.gFC.Vthr(2) - opt.AUX.gFC.Vthr(1))./opt.AUX.gFC.Vthr(3) + 1);
return
end