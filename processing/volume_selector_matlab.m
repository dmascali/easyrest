function indx = volume_selector_matlab(N,g,j) %g group, j session
% compute indx for removing unwanted volumes.  This function is needed only when you need to load preprocessed
% data, not already cutted by afni.
% indx is 1 where data must be removed.
global opt
if isempty(opt.AUX.v_s_matlab{g,j})
    indx = false(1,N);
else
    indx = ones(1,N);
    M1 = opt.AUX.v_s_matlab{g,j}(1); % lower bound
    M2 = opt.AUX.v_s_matlab{g,j}(2); % upper bound
    indx(1,M1:M2) = 0;
    indx = logical(indx);
end
return
end
