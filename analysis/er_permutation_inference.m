function [STATS] = er_permutation_inference(Y,X,exchange,C,perm_number,alpha_level,matrix_type,zscore_flag,names)
%not working for one-sample t-test
disp_waitbar = 0;  %it really slows the computation

siz = size(Y);
m = size(Y,2);  %m is not necessarily the number of tests. Number of tests is M
n = size(Y,1);

n_dimension = length(siz);
if n_dimension == 3
    if siz(2)~= siz(3)
        error('It is not an ROI to ROI matrix');
    end
    matrix_mode = 1;
    switch matrix_type
        case 0 %simmetric & diagonal not relevant
            index = tril(ones(m,m),-1);
            index = find(index);
            M = m*(m-1)/2;
            YY = zeros(n,M);
            YY = Y(:,index);
            Y = YY; clear YY;
        case 1 %simmetric & diagonal relevant
            index = tril(ones(m,m),0);
            index = find(index);
            M = m*(m-1)/2 + m;
            YY = zeros(n,M);
            YY = Y(:,index);
            Y = YY; clear YY;            
        case 2 % NOT simmetric matrix, all elements will be considered
            Y = reshape(Y,[n,m*m]);
            M = m*m;
    end
elseif n_dimension == 2
    matrix_mode = 0;
    M = m;
else
    error('Not supported Y format');
end

%original stats, not permuted data
[STATS_orig] = er_glm_fit(Y,X,C,zscore_flag); T = STATS_orig.T;
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

ind_nuisance=find(~C);
if isempty(ind_nuisance)
    %No nuisance predictors
else
    %Regress out nuisance predictors and compute residual
    b=zeros(length(ind_nuisance),M);
    resid_y=zeros(n,M);
    b=X(:,ind_nuisance)\Y;
    resid_y=Y-X(:,ind_nuisance)*b;
end

max_permutation = factorial(n);
if perm_number > max_permutation
    warning('Adjusting the number of permutation');
    perm_number = max_permutation;
end

Tdistr = nan(perm_number,1);  %for fwe
pvals_pos = zeros(1,M);       %for fdr
pvals_neg = pvals_pos;
pvals_ttails = pvals_pos;

if disp_waitbar
    h = waitbar(0,['Running ',num2str(perm_number),' permutations...']);
end
for mm = 1:perm_number
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
            y_perm = resid_y+[X(:,ind_nuisance)]*b;
        else
            resid_y = resid_y(randperm(n),:);
            y_perm = resid_y+[X(:,ind_nuisance)]*b;
        end
    end
    [STATS] = er_glm_fit(y_perm,X,C,zscore_flag);
    %----for fwe---------
    [~,maxind] = max(abs(STATS.T(:)));
    Tdistr(mm) = STATS.T(maxind);
    %----for fdr---------
    pvals_pos = pvals_pos + (STATS_orig.T <= STATS.T);
    pvals_neg = pvals_neg + (STATS_orig.T >= STATS.T);
    pvals_ttails = pvals_ttails + (abs(STATS_orig.T) <= abs(STATS.T));
    if disp_waitbar; waitbar(mm/perm_number,h); end
end
if disp_waitbar; close(h); end;

%----for fdr---------
pvals_ttails = pvals_ttails/perm_number; PPermBased.fdr(1).P = mafdr(pvals_ttails,'BHFDR','true','Showplot','false');
pvals_pos = pvals_pos/perm_number;       PPermBased.fdr(2).P = mafdr(pvals_pos,'BHFDR','true','Showplot','false');
pvals_neg = pvals_neg/perm_number;       PPermBased.fdr(3).P = mafdr(pvals_neg,'BHFDR','true','Showplot','false');
PPermBased.fdr(1).name = 'FDR Two-Tails';
PPermBased.fdr(2).name = 'FDR One-Tail Positive';
PPermBased.fdr(3).name = 'FDR One-Tail Negative';
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
%figure; hist(Tdistr(:),sqrt(perm_number))

PPermBased.fwe(1).P = nan(1,M); PPermBased.fwe(1).name = 'FWE Two-Tails';
Critic_t(1,2) = prctile(abs(Tdistr),100-100*alpha_level);
Critic_t(1,1) = -Critic_t(1,2);
PPermBased.fwe(1).Critic_t = Critic_t;
PPermBased.fwe(2).P = nan(1,M); PPermBased.fwe(2).name = 'FWE One-Tail Positive';
Critic_t(1,2) = prctile(Tdistr,100-100*alpha_level);
Critic_t(1,1) = NaN;
PPermBased.fwe(2).Critic_t = Critic_t;
PPermBased.fwe(3).P = nan(1,M); PPermBased.fwe(3).name = 'FWE One-Tail Negative';
Critic_t(1,1) = prctile(Tdistr,100*alpha_level);
Critic_t(1,2) = NaN;
PPermBased.fwe(3).Critic_t = Critic_t;
PPermBased.Tdistr_fwe = Tdistr;
for l = 1:M
    PPermBased.fwe(1).P(1,l) = mean((abs(Tdistr) >= abs(T(l))));
    PPermBased.fwe(2).P(1,l) = mean(Tdistr >= T(l));
    PPermBased.fwe(3).P(1,l) = mean(Tdistr <= T(l));
end
%---------------------


%-----reshaping to oringal form--------
if matrix_mode  %else it is already ready
    for l = 1:3
        PPermBased.fwe(l).P = reshape_to_matrix(PPermBased.fwe(l).P,matrix_type,m,index);
        PPermBased.fdr(l).P = reshape_to_matrix(PPermBased.fdr(l).P,matrix_type,m,index);
    end
%    PPermBased.fdr(4).P = reshape_to_matrix(PPfdrasNBS,matrix_type,m,index);
    Pfdr_noPermutation =  reshape_to_matrix(Pfdr_noPermutation ,matrix_type,m,index);
    STATS_orig.B = reshape_to_matrix(STATS_orig.B,matrix_type,m,index);
    STATS_orig.T = reshape_to_matrix(STATS_orig.T,matrix_type,m,index);
    STATS_orig.P = reshape_to_matrix(STATS_orig.P,matrix_type,m,index);
end
%--------------------------------------

STATS = [];
STATS.model.X = X;
STATS.model.Contrast = C;
STATS.model.Exchange_blocks = exchange;
STATS.model.zscore_transformation = zscore_flag;
STATS.B = STATS_orig.B;
STATS.T = STATS_orig.T;
STATS.P = STATS_orig.P;
STATS.Pfdr = Pfdr_noPermutation;
STATS.PPermBased = PPermBased;
STATS.PPermBased.perm_number = perm_number;

%er_permutation_plot(STATS,names)

return
end


function er_permutation_plot(STATS,names)
figure;
set(gcf,'name','Permutation Based Inference')
pos = [35 35 1039 746];
set(gcf,'Position',pos);

n_row = 3;
n_col = 3;

if STATS.model.zscore_transformation == 1
    Blim = [-1,1];
    B_name = 'Pearson''s r';
else
    Blim = [];
    B_name = 'Beta values';
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

plot_matrix(n_row,n_col,1,STATS.B,Blim,names,B_name);
plot_matrix(n_row,n_col,2,STATS.T,[minn,maxx],names,'T-values');
plot_matrix(n_row,n_col,3,STATS.P,[0 0.05],names,'Parametric P-values (Two-Tails) + FDR');

%fdr on parametric P-values
P = STATS.Pfdr;
[I,J] = find(P <= 0.05 & P > 0.01);
hstar1 = text(I,J,'*','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);
[I,J] = find(P <= 0.01 & P > 0.001);
hstar2 = text(I,J,'**','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);
[I,J] = find(P <= 0.001);
hstar3 =text(I,J,'***','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);

plot_matrix(n_row,n_col,6,STATS.PPermBased.fdr(1).P,[0 0.05],names,STATS.PPermBased.fdr(1).name);
plot_matrix(n_row,n_col,4,STATS.PPermBased.fdr(2).P,[0 0.05],names,STATS.PPermBased.fdr(2).name);
plot_matrix(n_row,n_col,5,STATS.PPermBased.fdr(3).P,[0 0.05],names,STATS.PPermBased.fdr(3).name);
% P = STATS.PPermBased.fdr(4).P;
% [I,J] = find(P);
% hstar1 = text(I,J,'*','HorizontalAlignment','center','fontsize',15,'FontWeight','Bold','Color',[1 1 1]);

plot_matrix(n_row,n_col,9,STATS.PPermBased.fwe(1).P,[0 0.05],names,STATS.PPermBased.fwe(1).name);
plot_matrix(n_row,n_col,7,STATS.PPermBased.fwe(2).P,[0 0.05],names,STATS.PPermBased.fwe(2).name);
plot_matrix(n_row,n_col,8,STATS.PPermBased.fwe(3).P,[0 0.05],names,STATS.PPermBased.fwe(3).name);


return
end

function plot_matrix(n_row,n_col,index,matrix,limits,names,title_)

s = subplot(n_row,n_col, index);
if isempty(limits)
    imagesc(matrix);
else
   imagesc(matrix,limits);
end
hold on
pbaspect([1 1 1])

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






function  P = reshape_to_matrix(Y,matrix_type,m,index)

switch matrix_type
    case 0 %simmetric, diagonal not relevant
        
        index_diag = find(diag(ones(1,m)));
        P = zeros(m,m);
        P(index) = Y;
        P = P'; P(index) = Y;
        P(index_diag) = nan; %putting diagonal values to nan;
      
    case 1 %simmetric, diagonal relevant
        
        index_diag = find(diag(ones(1,m)));
        P = zeros(m,m);
        P(index) = Y;
        P = P'; P(index_diag) = 0; P(index) = Y;
       
    case 2 %not simmetric case
        
        P = reshape(Y,[m,m]);
end
    
return
end

function [STATS] = er_glm_fit(Y,X,C,zscore_flag)
siz = size(Y);
n = siz(1);
M = siz(2);
n_dimension = length(siz);
if n_dimension == 3
    if siz(2)~= siz(3)
        error('It is not a ROI to ROI matrix');
    end
    Y = reshape(Y,[n, M*M]);
elseif n_dimension == 2
    %
else
    error('Not supported Y format');
end
if zscore_flag
    Y = zscore(Y);
    X = zscore(X);
end
rang = rank(X);
warning off; C = [C, zeros(1,(size(X,2)-size(C,2)))];warning on;  %add missing zeros
% inversa = inv(X'*X);
% for l = 1:size(Y,2)
%     B(:,l) = inversa*X'*Y(:,l);
%     res = Y(:,l)- X*B(:,l);
%     sigma2(:,l) = res'*res/(n - rang);
%     for j=1:size(C,1)
%         T(j,l) = C(j,:)*B(:,l)/sqrt(sigma2(:,l)*C(j,:)*inversa*C(j,:)');
%     end
% end
B = (X\Y);
RES = Y-X*B;
sigma2 = sum(RES.^2,1)/(n-rang);
T = C*B./sqrt(sigma2*(C*inv(X'*X)*C'));
P =  2 * tcdf(-abs(T), n-rang);
%reshape
if n_dimension == 3
    STATS.T = nan(size(C,1),M,M);
    STATS.P = nan(size(C,1),M,M);
    for j=1:size(C,1)
        STATS.B(j,:,:) = reshape(B(j,:),M,M);
        STATS.T(j,:,:) = reshape(T(j,:),M,M);
        STATS.P(j,:,:) = reshape(P(j,:),M,M);
    end
elseif n_dimension == 2
    STATS.T = nan(size(C,1),M);
    STATS.P = nan(size(C,1),M);
    for j=1:size(C,1)
        STATS.B(j,:) = B(j,:);
        STATS.T(j,:) = T(j,:);
        STATS.P(j,:) = P(j,:);
        
    end
end
return
end