function [STATS,RES] = er_glm_fit_beta(Y,X,C,varargin)
%ER_GLM(Y,X,C) fits a general linear model using the predictor matrix X and
% response variables in Y. 
% C is a matrix containing contrast vectors for inference (t-test; one 
% contrast for each row as in SPM). Use C = eye(size(X,2)) to test all predictors.
% Put C = [], if you just need residuals and don't care for statistical
% inference (usefull if you just want to regress out confounds from data).
%
%       [STATS]     = ER_GLM(..) returns statitical inference.
%       [STATS,RES] = ER_GLM(..) returns statitical inference and
%                                residuals.
% 
% Y can be a vector of observations or an N-dimensional matrix with the 
% first dimension containing observations and the other dimensions running 
% on different response variables (eg, time x voxels; or time x variable1 x
% variable 2 ...).
%
% Property list (details below):
%         "ynan"     = "skip" / "remove" (default skip)
%         "constant" = "on"   / "off"    (default on)
%         "rmvobs"   = [4 6 34]           (default [])
%         "tail"     = "both" / "right" / "left" / "single" (default both)
%         "zscore"   = "on"   / "off"    (default off)
%         "r2t"      = "on"   / "off"    (default off)
%         "Z-stat"   = "on"   / "off"    (default off)
%         "fdr"      = "on"   / "off"    (default off)
%         "PermMode" = "on"   / "off"    (default off) 
%         'showplot' = 'on'   / 'off'    (default on)
% 
% NaNs in Y are controlled with the property "ynan". NB: NaNs in X will be 
% automatically removed, no matter of the following properties.
%         "ynan" = "skip" (default) No fit will be computed on responding
%                variables containing NaNs.
%         "ynan" = 'remove' For those repsonding variables containing NaNs
%                their lenght will be adjusted removing NaNs. Thus, the fit will
%                be calculated for every variable in Y with enough dof.
%
% The constant term (intercept) is added to the model by default (unless
% the "zscore" option is turned on). You can disalbe the constant term via
% the property "constant"
%         "constant" = "on"/"off" (default on) Add a constant term to X.
%
% "rmvobs" is a vector of integer values that specifies observations that 
% you want to remove from the model. (default [])   
%
% Two-tails p-value are output by default. To swich to single tail test use
% the property "tail"
%         "tail" = "both" / "right" / "left" / "single" (default both)
% The "single" mode is a simultaneous left and right test (generally
% deprecated).
%
% You can standardize X and Y via the property "zscore" which removes the
% mean and divides by std both X and Y (mean = 0; std = 1). It is a usefull
% trick to compute Pearson's correlation via GLM. In this case you can also
% apply the property "r2t" to change the method to compute T and thus 
% P-value, ie, not from the model residual variance but from the r2t equation.
%         "zscore" = "on"/"off" (default off). Standardize X and Y. In this
%                  case B = r. 
%         "r2t"    = "on"/"off" (default off). Compute P from T calculated
%                  via eq: T = sqrt( (n-2)./(1-r^2) ). The default is to
%                  use GLM residual error to compute T.
%
% The default beahviour is to output (in STATS): B, T and P. You can
% use the property "Z-stat" to also output z-score. 
%         "Z-stat" = "on"/"off" (default off). Output also Z-score values.
%
% The property "fdr" enables FDR correction for multiple comparisons. The
% corrected p-values will be stored in STATS.Pfdr
%         "fdr"    = "on"/"off" (default off). Returns also FDR p-values.
%
% "PermMode" is a non verbose modality with no P-value computation (It's handy
% when this function is called several times by er_permutation_inference.)
%        "PermMode"= "on"/"off" (default off)
%__________________________________________________________________________
%
%   Version 1.0
%   Copyright (C) 2018   danielemascali@gmail.com

if nargin == 0
    help er_glm_fit_beta
    return
end

%--------------VARARGIN----------------------
params  =  {'tail','zscore','Z-stat','r2t', 'ynan','constant','fdr','permmode','showplot','rmvobs','spearman'};
defParms = {'both',   'off',   'off','off', 'skip',      'on','off',     'off',      'on',      [],     'off'};
legalValues{1} = {'both','right','left','single'};
legalValues{2} = {'on','off'};
legalValues{3} = {'on','off'};
legalValues{4} = {'on','off'};
legalValues{5} = {'skip','remove'};
legalValues{6} = {'on','off'};
legalValues{7} = {'on','off'};
legalValues{8} = {'on','off'};
legalValues{9} = {'on','off'};
legalValues{10} = [];
legalValues{11} = {'on','off'};
[tail,zscore_flag,zStat_flag,r2t,ynan,constant,fdr,PermMode,showplot,rmvobs,spearman] = parse_varargin(params,defParms,legalValues,varargin,1);
% %---convert chars to logical variables-----
% zscore_flag = char2logical(zscore_flag);
% zStat_flag = char2logical(zStat_flag);
% r2t = char2logical(r2t);
% constant = char2logical(constant);
% fdr = char2logical(fdr);
% PermMode = char2logical(PermMode);
% showplot = char2logical(showplot);
% %------------------------------------------

siz = size(Y);
n = siz(1);
n_dimension = length(siz);

if n_dimension > 2
    Y = reshape(Y,[n, prod(siz(2:end))]);
end

%remove observations if desired
if ~isempty(rmvobs)
    Y(rmvobs,:) = [];
    X(rmvobs,:) = [];
    n = size(Y,1);
end 

if spearman
    %force zscore on
    zscore_flag = 1;
end

warning off backtrace;
% check for NaNs in X
nanXindex = logical(sum(isnan(X),2));
if sum(nanXindex) > 0
    X(nanXindex,:)   = [];
    Y(nanXindex,:) = [];
    n = size(Y,1);
    if ~PermMode
        warning(sprintf(['There are NaNs in the predictors (i.e.,X). They will be removed from both X and Y.\nNew observation number = ',num2str(n),'.']));
    end
end

%check for NaNs in Y
remove_ynan = 0;
if logical(sum(sum(isnan(Y),1)))  %check if thereis at least one NaN
    switch ynan
        case {'skip'}
            if ~PermMode
                warning('There are NaNs in some of the responding variables (Y). For those variables no fit will be performed. You can use the option "ynan" to change this behaviour');
            end
        case {'remove'}
            remove_ynan = 1;
            if ~PermMode
                warning('There are NaNs in some of the responding variables (Y). NaNs will be removed and the fit will be performed separately for each reponding variable.');
            end
    end
end
warning on backtrace;


    
if ~remove_ynan  %skip case (default)  
    if spearman
        X = tiedrank(X);
        Y = tiedrank(Y);
    end
    if zscore_flag
        Y = zscore(Y);
        X = zscore(X);
    else %add constant term. There is no need to add the constant term after zscore convertion.   
        if constant; X = [X,ones(n,1)]; end  
    end
    
    rang = rank(X);
    
    B = (X\Y);
    RES = Y-X*B;
    %R2 and R2 adj computation
    SSres = sum(RES.^2,1);
    SStot = sum(bsxfun(@minus,Y,mean(Y)).^2,1);
    R2    = 1-SSres./SStot;
    R2adj = 1-(1-R2).*( (n-1)./(n-rang) ); % assuming rang == size(X,2). Non singular desing. 
    
    %-------------------------
    sigma2 = SSres/(n-rang);
else  %remove yNaNs
    %add constant term. There is no need to add the constant term after zscore convertion.   
    if constant && ~zscore_flag; X = [X,ones(n,1)]; end  
    
    NaNcount = sum(isnan(Y),1);   
    Ynan.n = n - NaNcount;
    rang = rank(X); 
    
    %find indexes
    NaNindex = logical(NaNcount);   %where thare are at least one NaN in Y
    Ynan.index{1} = find(not(NaNindex));                     %NOT NaN index  
    Ynan.index{2} = find( (Ynan.n > (rang+1)) & NaNindex);   %NaN index & valid (enough dof)
    
    % initialize to NaN
    Ynan.B       = nan(rang,size(Y,2));
    Ynan.sigma2  = nan(1,size(Y,2));
    Ynan.notNaNx = nan(n,length(Ynan.index{2})); %this variable is needed for index{2} ie., when there are NaNs

    [Ynan.B(:,Ynan.index{1}),Ynan.sigma2(1,Ynan.index{1}),~] = do_fit_remove_ynan(Y(:,Ynan.index{1}),X,rang,zscore_flag);
    
    for l = 1:length(Ynan.index{2})
       [Ynan.B(:,Ynan.index{2}(l)),Ynan.sigma2(1,Ynan.index{2}(l)),Ynan.notNaNx(:,l)] = do_fit_remove_ynan(Y(:,Ynan.index{2}(l)),X,rang,zscore_flag); 
    end
    RES = [];
    Ynan.notNaNx = logical(Ynan.notNaNx);
end

% C adjustment must be here since size(X,2) can change in case the
% constant term has been added
if ~isempty(C)
    warning off; C = [C, zeros(size(C,1),(size(X,2)-size(C,2)))];warning on;  %add missing zeros
    C_number = size(C,1);
end

if ~isempty(C)
    %--------------------------T----------------------------------
    if ~remove_ynan  %skip case (default)
        if C_number == 1
            if zscore_flag && r2t   %case Pearson (T not from model)
                T = B .* sqrt((n-2) ./ (1 - B.^2)); 
                CB = C*B;  %it's always just B
            else
                CB = C*B;
                T = CB./sqrt(sigma2*(C*inv(X'*X)*C'));
            end
        else
            CB = C*B;
            T = CB./sqrt( repmat(sigma2,[C_number,1]) .* diag((C*inv(X'*X)*C')) );
        end
    else
        [CB,T] = t_remove_ynan(Ynan,X,C,C_number,r2t,zscore_flag);
    end
    %------------------------------------------------------------

    %-------------------------P----------------------------------
    if ~PermMode
        T = double(T);
        if ~remove_ynan  %skip case (default)
            switch tail
                case {'both'}
                    P =  2 * tcdf(-abs(T), (n-rang-zscore_flag));   %zscoring is equivalent to add a constant vector in the model
                case {'right'}
                    P = tcdf(-T,(n-rang-zscore_flag));
                case {'left'}
                    P = tcdf(T,(n-rang-zscore_flag));
                case {'single'}
                    P =  tcdf(-abs(T), (n-rang-zscore_flag));
            end
        else
            P = p_remove_ynan(T,Ynan,rang,tail,zscore_flag);
            B = Ynan.B;
        end
        if fdr || showplot
            Pfdr = nan(size(P));
            for l = 1:C_number
                Pfdr(l,:) = mafdr(P(l,:),'BHFDR','true','Showplot','false');
            end
        end
    end
    %------------------------------------------------------------
    
    %-------------------------Z----------------------------------
    if zStat_flag
        % to obtain the Z-score the equations would be
        %    Z(T>0) = -norminv(tcdf(-T(T>0),dof));  %-sqrt(2)*erfcinv(2*P(T<0)); 
        %    Z(T<0) =  norminv(tcdf( T(T<0),dof));  %sqrt(2).*erfcinv(2*P(T>0))
        % to avoid recomputing the P values, I adjust the above calculated
        % P-values according to the chosen tail. However, that doesn't work
        % for left and right tail (not enough precision to compute 1-P when
        % P is close to 1)
        Z = zeros(size(T));
        switch tail
            case {'both'}
                Z(T>0) = -norminv(P(T>0)./2);  
                Z(T<0) =  norminv(P(T<0)./2);   
            case {'left', 'right'}
                if ~remove_ynan  %skip case (default)
                    Z(T>0) = -norminv(tcdf(-T(T>0),(n-rang-zscore_flag)));
                    Z(T<0) =  norminv(tcdf( T(T<0),(n-rang-zscore_flag)));
                else
                    Ptmp = p_remove_ynan(T,Ynan,rang,'single',zscore_flag);
                    Z(T>0) = -norminv(Ptmp(T>0));  
                    Z(T<0) =  norminv(Ptmp(T<0)); 
                end
             case {'single'}  
                Z(T>0) = -norminv(P(T>0));  
                Z(T<0) =  norminv(P(T<0));    
        end
        Z(Z==Inf)=38.5;  %z-score are limited to such value
        Z(Z==-Inf)=-38.5;
    end
    %------------------------------------------------------------
    
    
    %--------------------reshape and output----------------------
    if ~PermMode
        STATS.model.X = X;
        STATS.model.C = C;
        if exist('R2','var');   STATS.model.R2    = R2;     end;   %R2 is not available in ynan remove mode.
        if exist('R2adj','var');STATS.model.R2adj = R2adj;  end;
        STATS.model.constant = constant;
        STATS.model.ynan = ynan;
        STATS.model.zscore = zscore_flag;
        STATS.model.r2t = r2t;
        if n_dimension > 2
            STATS.B  = reshape(B,[rang,siz(2:end)]); %reshaping checked (with an epi). It works
            STATS.CB = reshape(CB,[C_number,siz(2:end)]); 
            STATS.T  = reshape(T, [C_number,siz(2:end)]);
            STATS.P  = reshape(P, [C_number,siz(2:end)]);
            if zStat_flag;      STATS.Z    = reshape(Z,     [C_number,siz(2:end)]); end
            if fdr || showplot; STATS.Pfdr = reshape(Pfdr,  [C_number,siz(2:end)]); end
            if remove_ynan;     STATS.N    = reshape(Ynan.n,siz(2:end));        else STATS.N = n; end
        elseif n_dimension == 2
            STATS.B  = B;
            STATS.CB = CB; 
            STATS.T  = T;
            STATS.P  = P;
            if zStat_flag;  STATS.Z = Z;      end
            if remove_ynan; STATS.N = Ynan.n; else STATS.N = n; end
            if fdr;         STATS.Pfdr = Pfdr;end
        end
        if showplot
            Y = reshape(Y,[n, siz(2:end)]);
            er_UI_GLMplot(STATS,Y)
        end
    else
        if n_dimension > 2
            STATS.T = reshape(T,[C_number,siz(2:end)]);
        elseif n_dimension == 2
            STATS.T = T;
        end
    end

    %------------------------------------------------------------
else
    STATS = [];
end


return
end

function [B,sigma2,notNaNx] = do_fit_remove_ynan(Y,X,rang,zscore_flag)

%check for NaNs in Y
nanindex = logical(sum(isnan(Y),2));
if sum(nanindex) > 0
    Y(nanindex,:) = [];
    X(nanindex,:) = [];
    notNaNx = not(nanindex);
else
    notNaNx = [];
end

if zscore_flag  % the constant hasn't been added
    Y = zscore(Y);
    X = zscore(X);
end

B = (X\Y);
RES = Y-X*B;
sigma2 = sum(RES.^2,1)/(size(Y,1)-rang);
return
end

function [CB,T] = t_remove_ynan(Ynan,X,C,C_number,r2t,zscore_flag)

T  = nan(C_number,size(Ynan.B,2));
CB = nan(C_number,size(Ynan.B,2));

if C_number == 1 && zscore_flag && r2t   %case Pearson (T not from model)
    T = Ynan.B .* sqrt((Ynan.n-2) ./ (1 - Ynan.B.^2)); 
    CB= Ynan.B;
else
    CB(:,Ynan.index{1}) = C*Ynan.B(:,Ynan.index{1});
    %in case of z-score convertion, X needs to be converted
    if zscore_flag; Xtmp = zscore(X); else Xtmp = X; end
    T (:,Ynan.index{1}) = CB(:,Ynan.index{1})./sqrt( repmat(Ynan.sigma2(1,Ynan.index{1}),[C_number,1]) .* diag((C*inv(Xtmp'*Xtmp)*C')) );
    for l = 1:length(Ynan.index{2})
        %in case of z-score convertion, X needs to be converted
        Xtmp = X(Ynan.notNaNx(:,l),:);
        if zscore_flag;  Xtmp = zscore(X); end 
        CB(:,Ynan.index{2}(l)) = C*Ynan.B(:,Ynan.index{2}(l));
        T (:,Ynan.index{2}(l)) = CB(:,Ynan.index{2}(l))./sqrt( repmat(Ynan.sigma2(1,Ynan.index{2}(l)),[C_number,1]) .* diag((C*inv(Xtmp'*Xtmp)*C')));
    end
end

return
end

function P = p_remove_ynan(T,Ynan,rang,tail,zscore_flag)

P = nan(size(T));

switch tail
    case {'both'}
        if ~isempty(Ynan.index{1})
            P(:,Ynan.index{1}) =  2 * tcdf(-abs(T(:,Ynan.index{1})), (Ynan.n(Ynan.index{1}(1))-rang-zscore_flag));
        end
        for l = 1:length(Ynan.index{2})
            P(:,Ynan.index{2}(l)) =  2 * tcdf(-abs(T(:,Ynan.index{2}(l))), (Ynan.n(Ynan.index{2}(l))-rang-zscore_flag));
        end
    case {'right'}
        if ~isempty(Ynan.index{1})
            P(:,Ynan.index{1}) =  tcdf(-T(:,Ynan.index{1}), (Ynan.n(Ynan.index{1}(1))-rang-zscore_flag));
        end
        for l = 1:length(Ynan.index{2})
            P(:,Ynan.index{2}(l)) =  tcdf(-T(:,Ynan.index{2}(l)), (Ynan.n(Ynan.index{2}(l))-rang-zscore_flag));
        end
    case {'left'}
        if ~isempty(Ynan.index{1})
            P(:,Ynan.index{1}) =  tcdf(T(:,Ynan.index{1}), (Ynan.n(Ynan.index{1}(1))-rang-zscore_flag));
        end
        for l = 1:length(Ynan.index{2})
            P(:,Ynan.index{2}(l)) =  tcdf(T(:,Ynan.index{2}(l)), (Ynan.n(Ynan.index{2}(l))-rang-zscore_flag));
        end
    case {'single'}
        if ~isempty(Ynan.index{1})
            P(:,Ynan.index{1}) =  tcdf(-abs(T(:,Ynan.index{1})), (Ynan.n(Ynan.index{1}(1))-rang-zscore_flag));
        end
        for l = 1:length(Ynan.index{2})
            P(:,Ynan.index{2}(l)) =  tcdf(-abs(T(:,Ynan.index{2}(l))), (Ynan.n(Ynan.index{2}(l))-rang-zscore_flag));
        end
end

return
end