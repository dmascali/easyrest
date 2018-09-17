function [tmask,n_cens] = er_censoring_mask(vect,thr,varargin)
%given a vector VECT the script creates a temporal (binary) mask which is 
%   0 where VECT  > THR  (to be censored)
%   1 where VECT <= THR  (to be kept)
% this is the same convention of ANFI 3dTproject (censored points are 
% indicated with 0s)
% If an array of THR is provided, multiple masks will be created.

%Additional options can be specified using the following parameters (each 
% parameter must be followed by its value ie,'param1',value1,'param2',value2)
%
%Besides the standard mode (based on absolute THR values) you can create
% masks based on different algorithm:
%   - 'mode' = ['standard'/'random'/'top'/'bottom'] {default = 'standard'}
%       -standard- works as described above.
%       -random  - create random masks. In this modality, the THR values
%                  specify the number of volumes to be removed (or the
%                  percentage of volumes, depending on the'Units' property).
%       -top     - remove the top N (or percent, depending on the'Units' 
%                  property) points based on the VECT values. 
%       -bottom  - remove the bottom N (or percent, depending on the'Units' 
%                  property) points based on the VECT values. 
%   - 'Units' = ['N'/'percent']  {default = 'N'}
%      specify the units (number of points or percentage) for modality 
%      'random','top' and 'bottom'  
%
% In the standard mode only, you can specify the following properties:
%   - 'preTR' = [integer value],
%      to censor any number of previous TR {default = []}
%   - 'postTR' = [integer value],
%      to censor any number of post TR {default = []}

%--------------VARARGIN----------------------
params  =  {'preTR','postTR',    'mode',  'units'};
defParms = {     [],      [],'standard',      'N'};
legalValues{1} = [];
legalValues{2} = [];
legalValues{3} = {'standard','random','top','bottom'};
legalValues{4} = {'percent','N'};
[pre_TR,post_TR,Mode,Units] = parse_varargin(params,defParms,legalValues,varargin);
% %------------------------------------------

%vect must be a row vector
if not(isrow(vect))
    vect = vect';
end

N_masks = length(thr);
if strcmpi(Mode,'standard'); Units = 'vector unit';end
fprintf('\nCreating %d temporal mask(s) in modality ''%s'' (units: ''%s'')...',N_masks,Mode,upper(Units));
N = length(vect);
tmask = zeros(N_masks,N);

switch Mode
    case {'standard'}
        for l = 1:N_masks
            tmask(l,:) = cens_standard(vect,thr(l),pre_TR,post_TR,N);
        end
    case {'random','top','bottom'}
        %------------------------------------------------------------------
        switch Units %this section is common to all
            case {'n'}
                if ~isempty(find(thr> N,1))
                    error('In the current modality THR represents the number of volumes to be randomly censored, thus, thr must be <= N.');
                end                 
            case {'percent'}
                if ~isempty(find(thr> 100 | thr < 0,1))
                    error('In the current modality THR represents the percent of volumes to be randomly censored, thus, thr must be <= 100.');
                end 
                % convert thr from percent to N
                thr = round(thr.*N/100);
        end
        %------------------------------------------------------------------
        switch Mode
            case {'random'}
                for l = 1:N_masks
                    indx = randsample(N,thr(l),'false'); %random sample without replacements
                    tmask(l,indx) = 1;
                end
            case {'top','bottom'}
                if strcmpi(Mode,'top')
                    sort_str = 'descend';
                else
                    sort_str = 'ascend';
                end
                for l = 1:N_masks
                    [~,indx] = sort(vect,2,sort_str);
                    indx = indx(1:thr(l));
                    tmask(l,indx) = 1;
                end
        end
        %------------------------------------------------------------------
end

n_cens = sum(tmask,2)';

% from logic to double
% plus invert 0s with 1s (to follow ANFI's convention)
tmask = uint8(~tmask);

fprintf('done!\n');

return
end

function tmask = cens_standard(vect,thr,pre_TR,post_TR,N)

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

return
end

