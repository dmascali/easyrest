function dof = estimate_dof(N,M,TR,f1,f2)
%assuming stationary signal and aproximately white (not colored spectra)
%within the frequency band of interest
if nargin <= 3 % no filter
    f1 = 0;
    f2 = 1/(2*TR);
else
    %handle some actual and future situations
    if f2 == 99999    %afni code for no filtering
        f2 = 1/(2*TR);
    end
    % in case I will switch to a future different notation
    if f2 == inf
        f2 = 1/(2*TR);
    end
end
dof = max(0,N*2*TR*(f2-f1)-M);
return
end