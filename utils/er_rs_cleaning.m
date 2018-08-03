function [res,des_info,X] = er_rs_cleaning(Y,concat_index,TR,polort,pass_band,ort,cens,varargin)
% ER_RE_CLEANING mimics AFNI 3dTproject
%
% Y must be NxV (points x voxels)
% in case of concatenated runs use concat_index [1 240 480] to specify the starting index of each run. Otherwise put
% it as empty.
%
% polort defined as in afni NB: polort > 2 is not supported as it is a waste of DOF if low frequences are removed
% by the pass band filter.
%
% ort are the condounds timeseries against which the data will be ortoganalized.
% You can construct ort so to have one regressor across all runs (not
% optimal but it can be useful sometimes) or separate each run regressrs
% (eg, you can concatenate rp parameters across runs = 6 regressors or
% have 6*n_run regressors). 
% the mean will be removed from ort.
%
% cens is a binary vector. Zeros indicate the point to be removed
% 
% optional arguments:
% 'cenmode' can be 'KILL' or 'ZERO' (def: 'ZERO') (see 3dTproject for
% info).




params  = {'cenmode','ortdemean'}; %NB: In case of multiple sessions with one separate set of ort for each
                                        %session, deaming all X or separately each ort session is the same. 
                                        % So, ortdemean should be always
                                        % set to 1
defparms = {'zero', 1};
legalvalues{1} = {'zero','kill'};
legalvalues{2} = [0 1];
[cenmode,ortdemean] = parse_varargin(params,defparms,legalvalues,varargin);


if polort == -1
    warning('The intercept is not modeled. Be sure Y has been demeaned');
end

if ~isempty(concat_index)
    run_number = length(concat_index);
    if concat_index(1) ~= 1
        error('concat_index(1) must be 1.')
    end
    if run_number < 2 % no concat case
        N = size(Y,1);
        run_number = 1;
    else
        N = zeros(run_number,1);
        for l = 2:run_number
            N(l-1) = concat_index(l) - concat_index(l-1);  
        end
        N(end) = size(Y,1) - concat_index(end) + 1;
    end
else
    N = size(Y,1);
    run_number = 1;
end

X_dim = zeros(run_number,1);
for r = 1:run_number
    %add Legendre pol
    X_diag{r} = er_leg_pol(N(r),polort,0);
    if ~isempty(pass_band)
        % sines and cosines to apply a band_pass filter
        X_diag{r} = [X_diag{r},BandPassOrt(N(r),TR,pass_band(1),pass_band(2),1)];
    end
    X_dim(r) = size(X_diag{r},2);
end

X = [];
for l_row = 1:run_number
    X_row = [];
    for l_col = 1:run_number
        if l_row == l_col
            X_row = [X_row,X_diag{l_row}];
        else
            X_row = [X_row,zeros(N(l_row),X_dim(l_col))];
        end
    end
    X = [X;X_row];
end

%add ort regressors
if ~isempty(ort)
    s_ort = size(ort);
    if s_ort(1) ~= size(X,1) 
        error('ort regressors don''t have the correct size');
    end
    %demean ort
    if ortdemean
        ort = ort -mean(ort,1);
    end
    M_ort = size(ort,2);
    X = [X,ort];
else
    M_ort = 0;
end

%get now some basic info
N_tot = size(Y,1);
M = size(X,2);
M_band = size(X,2) - M_ort -(polort+1)*run_number;

%censoring initialization
if ~isempty(cens)
    if length(cens) ~= N_tot
        error('cens length doesn''t match the time points in Y.');
    end
    N_cens = sum(cens == 0);
else 
    N_cens = 0;
end

%calculate dof
dof = N_tot - M - N_cens;

% Print Model info
fprintf('\nNumber of time points \t\t%d',N_tot);
fprintf('\nNumber of total regressors \t%d',M);
fprintf('\nNumber of Legendre pol \t\t%d',(polort+1)*run_number);
fprintf('\nNumber of passband regressors \t%d',M_band);
fprintf('\nNumber of ort regressors \t%d',M_ort);
fprintf('\nNumber of censored points \t%d',N_cens);
fprintf('\nNumber of Residual DoF \t\t%d\n',dof);

if dof <= 0
    error('Not enough dof to proceed.');
elseif dof < 10
    warning('Low nuber of DoF. Results might be not reliable.');
end

%do censoring
if ~isempty(cens)
    if strcmpi(cenmode,'kill')
        Y(cens == 0,:) = [];
        X(cens == 0,:) = [];
    elseif strcmpi(cenmode,'zero')
        Y(cens == 0,:) = 0;
        X(cens == 0,:) = 0;
    end
end

%let's do the glm
[~,res] = er_glm_fit(Y,X,[],0);

des_info.N = N_tot;
des_info.Mtot = M;
des_info.Mpol = (polort+1)*run_number;
des_info.Mband = M_band;
des_info.Mort = M_ort;
des_info.Mcens = N_cens;
des_info.dof = dof;
des_info.X = X;
if ~isempty(cens)
    des_info.cenmode = cenmode;
else
    des_info.cenmode = [];
end

return
end


function X = er_leg_pol(N,order,exclude_costant)
%return leg pol for regression till order 2
% use exclude_costant flag to exlude constant term
%
% order -1 return empty matrix

if order == -1
    X = [];
    return
elseif order == 0
    if ~exclude_costant
        X = ones(N,1);
    else
        X = [];
    end
elseif order == 1
    if ~exclude_costant
        X = [ones(N,1), linspace(-1,1,N)'];    
    else
        X = [linspace(-1,1,N)'];    
    end
elseif order == 2
    if ~exclude_costant
        X = [ones(N,1), linspace(-1,1,N)',(linspace(-1,1,N).^2)'];    
    else
        X = [linspace(-1,1,N)',(linspace(-1,1,N).^2)'];   
    end  
end

return
end