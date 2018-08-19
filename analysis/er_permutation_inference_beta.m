function [STATS] = er_permutation_inference_beta(Y,X,C,exchange,perm_number,varargin)
%ER_PERMUTATION_INFERENCE(Y,X,C,EXCHANGE,PERM_NUMBER) performs permuations
% based inference within a GLM framework using the predictor matrix X and
% response variables in Y. 
% - Y can be a vector of observations or an N-dimensional matrix with the 
% first dimension containing observations and the other dimensions running 
% on different response variables (eg, time x voxels; or time x variable1 x
% variable 2 ...).
% - C is a row vector containing the contrast for statistical inference 
% (only one contrast at a time is testable).  
% - EXCHANGE is a vector idicating the exchange blocks in case of repeated
% meaures designs (eg., paired t-test). EXCHANGE must have as many rows/col
% as there are observarions/subjects e.g., EXCHANGE = [1 2 3 4 1 2 3 4].
% Set EXCHANGE as an empty vector if it is not required.
% - PERM_NUMBER specifies the number of permutations.
%
%Property list (details below):
%           'TestType'   = 'onesample'/'ttest'(default ttest)
%           'ynan'       = 'skip' / 'remove'  (default skip)
%           'constant'   = 'on'   / 'off'     (default on)
%           'SqrMaType'  = [0,1,2]            (default 0)
%           'tail'       = 'both' / 'single'  (default both)
%           'zscore'     = 'on'   / 'off'     (default off)
%           'showplot'   = 'on'   / 'off'     (default on)
%           'showaitbar' = 'on'   / 'off'     (default on)
%
%Property details:
% TesType   - ['ttest','onesample'] Only t-tests are possible at the moment (no
%             F-test). For two-sample t-test or t-test on not-constant
%             responding variables use the 'ttest' option. In case of
%             one-sample t-test (i.e., testing a constant regressor like 
%             X = [ones(n,1)]) use the option 'onesample'. In this case for
%             each permutation the sign of each data point is randomly
%             flipped. {default = ttest}
% ynan      - ['skip'/'remove'] NaNs in Y are controlled with this property
%             'skip': no fit will be computed on responding variables 
%             containing NaNs. 
%             'remove': for those repsonding variables containing NaNs their
%             lenght will be adjusted removing NaNs. Thus, the fit will be 
%             calculated for every variable in Y with enough dof. In this
%             modality the FWE correction is deprecated.
%             NB: NaNs in X will be automatically removed, no matter of this
%             property. {default = 'skip'}
% constant  - ['on'/'off'] Adds the constant term (intercept) into the model.
%             This option is ignored in case 'zscore' is set to 'on'.
%             {default = 'on'}
% SqrMaType - [0,1,2] In case dependent variables form a square matrix
%             (e.g., Y = Nx10x10) specificies how to deal with repeated
%             variables or irrelevant elements. This option is critical for
%             correct FDR computation (parametric and not).
%               0 = consider all elements (no repetitions) 
%               1 = simmetric square matrix & diagonal NOT relevant. This is
%                   the case of common connectivity matrix.  
%               2 = simmetric square matrix & diagonal relevant. {default = 0}
% tail      - ['both'/'single'] Two or one side alternative hypothesis test.
%             This property has effect only for parametric p-FDR (i.e.,
%             it doesn't affect permutations. Permutations tests will be 
%             output always in single and both tails mode).{default = 'both'}
% zscore    - ['on'/'off'] Standardizes X and Y (i.e., it removes the mean 
%             and divides by std both X and Y; mean = 0; std = 1). It is a
%             usefull trick to compute Pearson's correlation via GLM (b == r) 
%             {default = 'off'}
% ShowPlot  - ['on'/'off'] Plots results (only available if ndims(Y) < 4)
%             and Tmax distributions. {default = 'on'}
% ShoWaitBar- ['on'/'off'] Displays progression bar (it slows computation)
%             {default = 'on'}
%
%Usage example:
% -General linear regression test:
%   Y = randn(30,10,5); X = [randn(30,3)];
%   stats = er_permutation_inference_beta(Y,X,[0 1],[],500);
% -Two-sample t-test + 1 covariate:
%   Y = randn(30,10,5); 
%   X = [[ones(20,1);zeros(10,1)],[zeros(20,1);ones(10,1)],randn(30,1)];
%   stats = er_permutation_inference_beta(Y,X,[1 -1],[],500,'constant','off');
% -Pearson's correlation
%   Y = randn(30,10,5); X = [randn(30,1)];
%   stats = er_permutation_inference_beta(Y,X,[1],[],500,'zscore','on');
% -One sample t-test
%   Y = randn(30,10,5); Y(1:10,:,:) = 1; X = [ones(30,1)];
%   stats = er_permutation_inference_beta(Y,X,[1],[],500,'testype','onesample');
%__________________________________________________________________________
%
%   Version 1.0
%   Copyright (C) 2018   danielemascali@gmail.com

global SquareMatrix

if nargin == 0
    help er_permutation_inference_beta
    return
end

%--------------VARARGIN----------------------
params  =  {'tail','zscore','ynan','constant','alpha','showplot','ShoWaitBar','SqrMaType','NetNames','TesType'};
defParms = {'both',   'off','skip',      'on',   0.05,      'on',        'on',         0 ,        [],  'ttest'};
legalValues{1} = {'both','single'};     %tail
legalValues{2} = {'on','off'};          %zscore
legalValues{3} = {'skip','remove'};     %ynan
legalValues{4} = {'on','off'};          %constant
legalValues{5} = [];                    %alpha: only affects Critic_t
legalValues{6} = {'on','off'};          %showplot
legalValues{7} = {'on','off'};          %ShoWaitBar it really slows computation
legalValues{8} = [0 1 2];               %SqrMaType
legalValues{9} = [];                    %NetNames works only for square matrix (xlabel == ylabel)
legalValues{10}= {'onesample','ttest'}; %TesType
[tail,zscore_flag,ynan,constant,alpha_level,ShowPlot,ShoWaitBar,SqrMaType,NetNames,testype] = parse_varargin(params,defParms,legalValues,varargin);
%---convert chars to logical variables-----
% zscore and ynan don't require char2logical convertions, they
% are directly pass though er_glm
ShowPlot = char2logical(ShowPlot);
ShoWaitBar = char2logical(ShoWaitBar);
constant = char2logical(constant);
%------------------------------------------

siz = size(Y);
n = siz(1);

if constant && not(strcmpi(testype,'onesample') || zscore) 
    X = [X,ones(n,1)];
end
    
n_dimension = length(siz);
if n_dimension == 3 && siz(2)== siz(3)
    SquareMatrix = 1;
    m = siz(2);
    warning off backtrace; warning(sprintf(['Square matrix detected, SqrMaType = ',num2str(SqrMAType),'.'])); warning on backtrace;
    switch SqrMaType
        case 1 %simmetric & diagonal not relevant
            index = find(tril(ones(m,m),-1));
            Y = Y(:,index);
            M = m*(m-1)/2;
        case 2 %simmetric & diagonal relevant
            index = find(tril(ones(m,m),0));
            Y = Y(:,index); 
            M = m*(m-1)/2 + m;
        case 0 % NOT simmetric matrix, all elements will be considered
            index = [];
            Y = reshape(Y,[n,m*m]);
            M = m*m;
    end
elseif n_dimension >= 2  %this contain also the case of an Nx1 vector
    SquareMatrix = 0;
    Y = reshape(Y,[n, prod(siz(2:end))]);
    M = prod(siz(2:end));
end

% check for NaNs in X
nanXindex = logical(sum(isnan(X),2));
if sum(nanXindex) > 0
    X(nanXindex,:) = [];
    Y(nanXindex,:) = [];
    n = size(Y,1);
    warning off backtrace; warning(sprintf(['There are NaNs in the predictors (i.e.,X). They will be removed from both X and Y.\nNew observation number = ',num2str(n),'.'])); warning on backtrace;
end

%original stats, not permuted data
[STATS_orig] = er_glm_fit_beta(Y,X,C,'zscore',zscore_flag,'ynan',ynan,'constant','off','tail',tail,'showplot','off'); T = STATS_orig.T;
Pfdr_noPermutation = mafdr(STATS_orig.P,'BHFDR','true','Showplot','false');

if not(isempty(exchange))
    %Set up exchange blocks
    blks=unique(exchange); 
    %Number of blocks
    n_blks=length(blks);
    %Number of observations per block
    sz_blk=n/n_blks; 
    blk_ind=zeros(n_blks,sz_blk);
    for i=1:length(blks)
        blk_ind(i,:)=find(blks(i)==exchange);
    end
end

%add missing zeros in Contrast vector
C = [C, zeros(1,(size(X,2)-size(C,2)))];  

ind_nuisance=find(~C);
if isempty(ind_nuisance)
    %No nuisance predictors
else
    %Regress out nuisance predictors and compute residual
    b=X(:,ind_nuisance)\Y;
    resid_y=Y-X(:,ind_nuisance)*b;
end

max_permutation = factorial(n);
%TODO adjust to consider the paired case
if perm_number > max_permutation
    warning('Adjusting the number of permutations.');
    perm_number = max_permutation;
end

Tdistr = nan(perm_number,1);  %for fwe
Tdistr_pos = Tdistr; Tdistr_neg = Tdistr;
pvals_pos = zeros(1,M);       %for fdr
if strcmp(ynan,'skip')
    NaNindex = logical(sum(isnan(Y),1));
    pvals_pos(NaNindex) = NaN;
end
pvals_neg = pvals_pos;
pvals_ttails = pvals_pos;

if ShoWaitBar; h = waitbar(0,['Running ',num2str(perm_number),' permutations...']); end
for mm = 1:perm_number
    if mm == 1 %original data (must be included so that p >= 1/perm_number)
            STATS.T = T; 
    else  
        if isempty(ind_nuisance)
            if exist('blk_ind','var')
                for j=1:n_blks
                    y_perm(blk_ind(j,:),:)=...
                    Y(blk_ind(j,randperm(sz_blk)),:);  
                end
            else
                y_perm = Y(randperm(n),:);
            end
        else
            if exist('blk_ind','var')
                for j=1:n_blks
                    resid_y(blk_ind(j,:),:)=...
                    resid_y(blk_ind(j,randperm(sz_blk)),:);
                end        
            else
                resid_y = resid_y(randperm(n),:);
            end
            %Freedman & Lane
            %Add permuted residual back to nuisance signal, giving a realisation 
            %of the data under the null hypothesis.
            y_perm = resid_y+X(:,ind_nuisance)*b; %this is not stricktly required see winkler2014 neuroimage
        end
        if strcmpi(testype,'onesample')
            y_perm = y_perm.*repmat(sign(rand(n,1)-0.5),1,M);
            [STATS] = er_glm_fit_beta(y_perm,X,C,'zscore',zscore_flag,'ynan',ynan,'constant','off','PermMode','on');
        elseif strcmpi(testype,'ttest')
            [STATS] = er_glm_fit_beta(y_perm,X,C,'zscore',zscore_flag,'ynan',ynan,'constant','off','PermMode','on');
        end
    end
    %----for fwe---------
    [maxval,maxind] = max(abs(STATS.T(:)));
    Tdistr(mm) = STATS.T(maxind);
    Tdistr_pos(mm) = maxval;
    Tdistr_neg(mm) = min(STATS.T(:));
    %----for fdr---------
    pvals_pos = pvals_pos + (T <= STATS.T);
    pvals_neg = pvals_neg + (T >= STATS.T);
    pvals_ttails = pvals_ttails + (abs(T) <= abs(STATS.T));
    if ShoWaitBar; waitbar(mm/perm_number,h); end
end
if ShoWaitBar; close(h); end

%----for fdr---------
PPermBased.P(1).P = pvals_ttails/perm_number;   PPermBased.fdr(1).P = mafdr(PPermBased.P(1).P,'BHFDR','true','Showplot','false');
PPermBased.P(2).P = pvals_pos/perm_number;      PPermBased.fdr(2).P = mafdr(PPermBased.P(2).P,'BHFDR','true','Showplot','false');
PPermBased.P(3).P = pvals_neg/perm_number;      PPermBased.fdr(3).P = mafdr(PPermBased.P(3).P,'BHFDR','true','Showplot','false');
PPermBased.P(1).name = 'p Two-Tails';
PPermBased.P(2).name = 'p One-Tail Positive';
PPermBased.P(3).name = 'p One-Tail Negative';
PPermBased.fdr(1).name = 'p-FDR Two-Tails';
PPermBased.fdr(2).name = 'p-FDR One-Tail Positive';
PPermBased.fdr(3).name = 'p-FDR One-Tail Negative';
%---------------------

% %As NBS
% %Sort p-values
% [pvals_neg,ind_srt]=sort(pvals_neg);
% 
% %Simes procedure
% pvals_neg=(pvals_neg<=(1:M)/M*alpha_level);
% ind=find(pvals_neg); 
% 
% if ~isempty(ind)
%     %Maximum index
%     ind=ind(length(ind));
%     PPfdrasNBS = zeros(M,1);
%     PPfdrasNBS(ind_srt(1:ind))=1; 
% else
%     PPfdrasNBS = zeros(M,1);
% end


%----for fwe---------
if ShowPlot
    figure; 
    set(gcf,'name','Permutation Based Inference, max(T) distribution');
    subplot(3,1,1);
    hist(Tdistr(:),sqrt(perm_number));
    xlabel('T-score');
    ylabel('Count');
    title(['Distribution of max(T,across GLM) TT for ',num2str(perm_number),' permutations']);

    subplot(3,1,2);
    hist(Tdistr_pos(:),sqrt(perm_number));
    xlabel('T-score');
    ylabel('Count');
    title(['Distribution of max(T,across GLM) Positive Tail for ',num2str(perm_number),' permutations']);

    subplot(3,1,3);
    hist(Tdistr_neg(:),sqrt(perm_number));
    xlabel('T-score');
    ylabel('Count');
    title(['Distribution of max(T,across GLM) Negative Tail for ',num2str(perm_number),' permutations']);

end

PPermBased.fwe(1).P = nan(1,M); PPermBased.fwe(1).name = 'p-FWE Two-Tails';
Critic_t(1,2) = prctile(abs(Tdistr),100-100*alpha_level);
Critic_t(1,1) = -Critic_t(1,2);
PPermBased.fwe(1).Critic_t = Critic_t;
PPermBased.fwe(2).P = nan(1,M); PPermBased.fwe(2).name = 'p-FWE One-Tail Positive';
Critic_t(1,2) = prctile(Tdistr,100-100*alpha_level);
Critic_t(1,1) = NaN;
PPermBased.fwe(2).Critic_t = Critic_t;
PPermBased.fwe(3).P = nan(1,M); PPermBased.fwe(3).name = 'p-FWE One-Tail Negative';
Critic_t(1,1) = prctile(Tdistr,100*alpha_level);
Critic_t(1,2) = NaN;
PPermBased.fwe(3).Critic_t = Critic_t;
PPermBased.Tdistr_fwe = Tdistr;
% for l = 1:M
%     PPermBased.fwe(1).P(1,l) = mean((abs(Tdistr) >= abs(T(l))));
%     PPermBased.fwe(2).P(1,l) = mean(Tdistr_pos >= T(l));  %BEFORE WAS Tdistr >=
%     PPermBased.fwe(3).P(1,l) = mean(Tdistr_neg <= T(l));  %BEFORE WAS Tdistr <=
% end
PPermBased.fwe(1).P = mean(bsxfun(@ge, abs(Tdistr),     abs(T)),1);
PPermBased.fwe(2).P = mean(bsxfun(@ge, Tdistr_pos,          T ),1);% ge:>=  %BEFORE WAS Tdistr >=
PPermBased.fwe(3).P = mean(bsxfun(@le, Tdistr_neg,          T ),1);% le:<=  %BEFORE WAS Tdistr <=
%---------------------



%-----reshaping to oringal form--------
if SquareMatrix  %else it is already ready
    for l = 1:3
        PPermBased.fwe(l).P = reshape_to_SquareMatrix(PPermBased.fwe(l).P,SqrMaType,m,index);
        PPermBased.fdr(l).P = reshape_to_SquareMatrix(PPermBased.fdr(l).P,SqrMaType,m,index);
        PPermBased.P(l).P = reshape_to_SquareMatrix(PPermBased.P(l).P,SqrMaType,m,index);
    end
%    PPermBased.fdr(4).P = reshape_to_matrix(PPfdrasNBS,matrix_type,m,index);
    Pfdr_noPermutation =  reshape_to_SquareMatrix(Pfdr_noPermutation ,SqrMaType,m,index);
    STATS_orig.CB = reshape_to_SquareMatrix(STATS_orig.CB,SqrMaType,m,index);
    STATS_orig.T  = reshape_to_SquareMatrix(STATS_orig.T, SqrMaType,m,index);
    STATS_orig.P  = reshape_to_SquareMatrix(STATS_orig.P, SqrMaType,m,index);
else
    for l = 1:3
        PPermBased.fwe(l).P = squeeze(reshape(PPermBased.fwe(l).P,[1,siz(2:end)]));
        PPermBased.fdr(l).P = squeeze(reshape(PPermBased.fdr(l).P,[1,siz(2:end)]));
        PPermBased.P(l).P   = squeeze(reshape(PPermBased.P(l).P,[1,siz(2:end)]));
    end
%    PPermBased.fdr(4).P = reshape_to_matrix(PPfdrasNBS,matrix_type,m,index);
    Pfdr_noPermutation = squeeze(reshape(Pfdr_noPermutation,[1,siz(2:end)]));
    STATS_orig.CB      = squeeze(reshape(STATS_orig.CB,     [1,siz(2:end)]));
    STATS_orig.T       = squeeze(reshape(STATS_orig.T,      [1,siz(2:end)]));
    STATS_orig.P       = squeeze(reshape(STATS_orig.P,      [1,siz(2:end)]));
    
end
%--------------------------------------

STATS = [];
STATS.model = STATS_orig.model;
STATS.model.Exchange_blocks = exchange;
STATS.CB = STATS_orig.CB;
STATS.T =  STATS_orig.T;
STATS.P =  STATS_orig.P;
STATS.fdrP = Pfdr_noPermutation;
STATS.PermBasedP = PPermBased;
STATS.PermBasedP.perm_number = perm_number;

if ShowPlot && n_dimension < 4
    er_permutation_plot(STATS,NetNames)
end

return
end


function er_permutation_plot(STATS,names)
figure;
set(gcf,'name','Permutation Based Inference')
pos = [50 50 1100 1200];
set(gcf,'Position',pos);

n_row = 4;
n_col = 3;

if STATS.model.zscore
    Blim = [-1,1];
    B_name = 'Pearson''s r';
else
    Blim = [];
    B_name = 'Contrast*Beta';
    %TODO
end

%t limits
maxx = max(STATS.T(STATS.T ~= inf));
minn = min(STATS.T(STATS.T ~= inf));
% for simmetric colorbar
if abs(maxx) > abs(minn)
    minn = -abs(maxx);
    maxx = abs(maxx);
else
    minn = -abs(minn);
    maxx = abs(minn);
end

plot_matrix(n_row,n_col,1,STATS.CB,Blim,names,B_name);
plot_matrix(n_row,n_col,2,STATS.T,[minn,maxx],names,'T-values');
plot_matrix(n_row,n_col,3,STATS.P,[0 0.05],names,'Parametric P-values (Two-Tails) + FDR');

%fdr on parametric P-values
P = STATS.fdrP;
[I,J] = find(P <= 0.05 & P > 0.01);
hstar1 = text(J,I,'*','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);
[I,J] = find(P <= 0.01 & P > 0.001);
hstar2 = text(J,I,'**','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);
[I,J] = find(P <= 0.001);
hstar3 =text(J,I,'***','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);

plot_matrix(n_row,n_col,6,STATS.PermBasedP.P(1).P,[0 0.05],names,STATS.PermBasedP.P(1).name);
plot_matrix(n_row,n_col,4,STATS.PermBasedP.P(2).P,[0 0.05],names,STATS.PermBasedP.P(2).name);
plot_matrix(n_row,n_col,5,STATS.PermBasedP.P(3).P,[0 0.05],names,STATS.PermBasedP.P(3).name);

plot_matrix(n_row,n_col,9,STATS.PermBasedP.fdr(1).P,[0 0.05],names,STATS.PermBasedP.fdr(1).name);
plot_matrix(n_row,n_col,7,STATS.PermBasedP.fdr(2).P,[0 0.05],names,STATS.PermBasedP.fdr(2).name);
plot_matrix(n_row,n_col,8,STATS.PermBasedP.fdr(3).P,[0 0.05],names,STATS.PermBasedP.fdr(3).name);
% P = STATS.PPermBased.fdr(4).P;
% [I,J] = find(P);
% hstar1 = text(I,J,'*','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);

plot_matrix(n_row,n_col,12,STATS.PermBasedP.fwe(1).P,[0 0.05],names,STATS.PermBasedP.fwe(1).name);
plot_matrix(n_row,n_col,10,STATS.PermBasedP.fwe(2).P,[0 0.05],names,STATS.PermBasedP.fwe(2).name);
plot_matrix(n_row,n_col,11,STATS.PermBasedP.fwe(3).P,[0 0.05],names,STATS.PermBasedP.fwe(3).name);


return
end

function plot_matrix(n_row,n_col,index,matrix,limits,names,title_)
global SquareMatrix
s = subplot(n_row,n_col, index);
if isempty(limits)
    imagesc(matrix);
else
   imagesc(matrix,limits);
end
hold on
if SquareMatrix
    pbaspect([1 1 1])
end

pos = get(s,'position');
hb = colorbar('location','eastoutside');
set(s,'position',pos);

set(gca,'ytick',[1:length(names)])
set(gca,'yTickLabel',names,'fontsize',8);
set(gca,'xtick',[1:length(names)])
set(gca,'xTickLabel',names,'fontsize',8);

title([title_]);
return
end


function  P = reshape_to_SquareMatrix(Y,matrix_type,m,index)
switch matrix_type
    case 1 %simmetric, diagonal not relevant       
        index_diag = find(diag(ones(1,m)));
        P = zeros(m,m);
        P(index) = Y;
        P = P'; P(index) = Y;
        P(index_diag) = nan; %putting diagonal values to nan;      
    case 2 %simmetric, diagonal relevant
        index_diag = find(diag(ones(1,m)));
        P = zeros(m,m);
        P(index) = Y;
        P = P'; P(index_diag) = 0; P(index) = Y;
    case 0 %not simmetric case
        P = reshape(Y,[m,m]);
end
return
end

function out = char2logical(inp)
%convert char varargin to logical
switch inp
    case {'on','true'}
        out = 1;
    case {'off','false'}
        out = 0;
end
out = logical(out);
return
end

