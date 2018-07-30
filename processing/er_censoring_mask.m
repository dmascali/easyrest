function [tmask,n_cens] = er_censoring_mask(vect,thr,pre_TR,post_TR)
%given a vector VECT the script creates a temporal (binary) mask which is 
% 0 where VECT  > THR
% 1 where VECT <= THR
% this is the same convention of ANFI 3dTproject (censored points are indicated with 0s)
%
% Additionaly, you can censor any number of previous TR (pre_TR) or post TR (post_TR).

if nargin == 2
    pre_TR = [];
    post_TR = [];
elseif nargin == 3
    post_TR = [];
end

N = length(vect);

tmask = (vect > thr);

indx = find(tmask);

if ~isempty(pre_TR) && pre_TR > 0
   for l = 1:pre_TR
       indxtemp = indx; indxtemp(indxtemp <= l) = [];
       tmask(indxtemp-l) = 1;
   end
end

if ~isempty(post_TR) && post_TR > 0
   for l = 1:post_TR
       indxtemp = indx; indxtemp(indxtemp >= (N-l+1)) = [];
       tmask(indxtemp+l) = 1;
   end
end

n_cens = sum(tmask);

% from logic to double
% plus invert 0s with 1s (to follow ANFI's convention)
tmask = uint8(~tmask);

return
end

