function X = er_acompcor(data,rois,dime,confounds,varargin)
%ER_ACOMPCOR(DATA,ROIS,DIME) extracts signals from DATA using ROIS as masks.
% -DIME specifies the number/type of signals. If DIME = 0 only the mean
%  signal will be extracted, if DIME > 0 the first n=DIME principal components
%  will be extracted (following the aCompCor approach).
% -DATA can be a matrix or a path to a nifti file, if DATA is a matrix, the 
%  first dimension must be time. 
% -ROIS is a cell array containg either matrices or paths to nifti files.
%  ROIS must be binary.
%
%ER_ACOMPCOR(DATA,ROIS,DIME,CONFOUNDS) if a matrix of confounds variable is
% provided, these will be regressed out of data before extracting  any signal. 
%
%Additional options can be specified using the following parameters (each 
% parameter must be followed by its value ie,'param1',value1,'param2',value2):
%
%   'firstmean' :   ['on'/'off'] If 'on', the first extracted component is 
%                   the mean signal, then the PCA is performed on data 
%                   ortogonalised to the mean signal {default='off'}
%   'deirvatives' : is a vector of roi length that specifies the degree of
%                   derivatives to be computed on the extracted singals
%                   {default=[],which is the same as zeros(1,length(rois))
%   'TvarNormalise':['on'/'off'], if set to 'on' the data is normalised by
%                   its temporal variance {default='off'}
%   'DataCleared' : ['true'/'false'] if 'true' the function avoids to seek
%                   for NaNs or voxels with zero signals (to save time)
%                   {default='false'}
%
% see aCompCor paper for details: neuroimage 37(2007) 90-101
%__________________________________________________________________________
% First version: 2016        at CMRR
%                August 2018 at MARBILab
% Author:
%   Daniele Mascali
%   Enrico Fermi Center, MARBILab, Rome
%   August, 2018
%   danielemascali@gmail.com


if nargin < 3
    error('Not enough input.');
end

%--------------VARARGIN----------------------
params  =  {'firstmean','derivatives','TvarNormalise','DataCleared'};
defParms = {      'off',           [],          'off',      'false'};
legalValues{1} = {'on','off'};
legalValues{2} = {[]};
legalValues{3} = {'on','off'};
legalValues{4} = {'true','false'};
[firstmean,deri,TvarNormalise,DataCleared] = parse_varargin(params,defParms,legalValues,varargin);
%---convert chars to logical variables-----
firstmean     = char2logical(firstmean);
TvarNormalise = char2logical(TvarNormalise);
DataCleared   = char2logical(DataCleared);
%--------------------------------------------

if ~iscell(rois)
    error('Please provide rois as cell, i.e., rois = {''path1'',''path2.nii''} or rois = {matrix1,matrix2}');
end
n_rois = length(rois);
if length(dime)~=n_rois
    error('Please specify a dime for each rois e.g., dime = [5 5]');
end
if ~isempty(deri)
    %check if there is one value for each roi
    if length(deri) ~= n_rois
        error('Please specify a derivative value for each roi e.g., ...,''der'',[1 0]');
    end
else
    deri = zeros(1,n_rois);
end


%--------------------------------------------------------------------------
%------LOADING DATA  and reshape-------------------------------------------
if ischar(data)  %in case data is a path to a nifti file
    data = spm_read_vols(spm_vol(data));
    s = size(data);
    data = reshape(data,[s(1)*s(2)*s(3),s(4)])';
else
    % In this case the first dimension must be time!
    s = size(data);
    n_dimension = length(siz);
    if n_dimension > 2
        data = reshape(data,[s(1), prod(s(2:end))]);
    end
end
%--------------------------------------------------------------------------

for r = 1:n_rois
    %----------------------------------------------------------------------
    %------LOADING ROI, Checking compatibility with data and reshape-------
    if ischar(rois{r})
        ROI = spm_read_vols(spm_vol(rois{r}));
        sr = size(ROI);
        %check if s and sr are identical in the first 3 dimensions
        if ~logical(sr == s(1:3))
            error(sprintf('ROI %d does not have the same dimension of data',r));
        end
        ROI = reshape(ROI,s(1)*s(2)*s(3));
    else
        ROI = rois{r};
        sr = size(ROI);
        if length(sr) > 2  %it's a 3d volume
            if ~isequal(s(2:end),sr) 
                error(sprintf('ROI %d does not have the same dimension of data',r));
            end
            % ok, reshape
            ROI = reshape(ROI,s(1)*s(2)*s(3));
        elseif isvector(ROI)  %in this case the data is assumed to be 2D
            if numel(ROI)~= s(2)
                error(sprintf('ROI %d does not have the same dimension of data',r));
            end
            % no need to reshape
        else
            error(sprintf('ROI %d does not have the same dimension of data',r));
        end
    end
    %----------------------------------------------------------------------
    %--------------------------------
    %check if ROI is binary
    un = unique(ROI(:));
    if length(un) > 2 || ~logical(un == [0 1])
        error(sprintf('ROI %d is not binary',r));
    end
    %--------------------------------
    % data extraction
    indx = find(ROI);
    V = data(:,indx);
    %--------------------------------
    % removal of costant and linear trend
    V = detrend(V);
    % lets remove voxels whose variance is equal zero (no signal in those voxels)
    % and Nan values
    if ~DataCleared
        stdv = std(V);
        a = zeros(1,length(stdv));
        a(stdv == 0) = 1;
        a(isnan(stdv)) = 1;
        V(:,logical(a)) = [];
    end
    %------------------Orthogonalise V-------------------------------------
    COV = [];
    if firstmean % as done in CONN: first extract the mean signal (mS), then compute PCA over data ortogonalised to mS
        mS = detrend(mean(V,2));
        COV = [COV,mS];
    end
    if ~isempty(confounds)  % extract PCA over data already ortogonalised to confounds.
        COV = [COV,detrend(confounds)];
    end
    if ~isempty(COV)
        V = V-COV*(COV\V);
    end
    %----------------------------------------------------------------------
    if dime(r) > 0
        if TvarNormalise
            %tvariance normalization
            V = bsxfun(@rdivide,V,std(V));
        end
        if 0 %same results   verified (2018)
            [~,U,~] = pca(V);
            comp = U(:,1:opt.aCompCor.ROI(r).dime);
        else
            %-------A bit of math--------------------------------------------------------------------------------
            % The principal components T of a matrix A are obtained by projecting A on the space
            % of the eigenvectors of the convariance matrix A'A, that we call V (also known as principal axes).
            % T = AV
            % A can be written is SVD as
            % A = USV' where U are eigenvectors of AA'
            %                V are eigenvectors of A'A
            % Thus,
            % T = USV'V = US
            % If we take the svd
            % A'A = VS^2V'
            % AA' = US^2U
            % Thus, these are equivalent:
            %          [U1,U2] = svd((V*V'));
            %          comp1 = [mS,U1(:,1:4)*diag(sqrt(diag(U2(1:4,1:4))))];
            %          [U1,U2] = svd(V);
            %          comp2 = [mS,U1(:,1:4)*diag(diag(U2(1:4,1:4)))];
            %          [coeff score latent] = pca(V);
            %          comp2 = [mS,score(:,1:4)];
            %----------------------------------------------------------------------------------------------------
            [U,P] = svd(V);
            if asCONN %add the mean signal and remove one dimension
                comp = [mS,U(:,1:dime(r)-1)*diag(diag(P(1:dime(r)-1,1:dime(r)-1)))];
            else
                comp = U(:,1:dime(r))*diag(diag(P(1:dime(r),1:dime(r))));
            end
        end
    else  %if dime == 0 simply compute the straight average
        comp = mean(V,2);
    end
    % derivatives computation, if requested
    if deri(r) > 0
        d = [];
        for l = 1:deri(r)
            d = [d, [zeros(l,size(comp,2));diff(comp,l)] ]; % I have to add l zeros as first rows...
        end
    else 
        d = [];
    end
    X = [comp,d];
    % variance normalize the extracted componets
    X = X./std(X,0,1);
end
    

return
end