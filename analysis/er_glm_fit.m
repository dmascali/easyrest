function [STATS,RES] = er_glm_fit(Y,X,C,zscore_flag)

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
if ~isempty(C)
    warning off; C = [C, zeros(size(C,1),(size(X,2)-size(C,2)))];warning on;  %add missing zeros
    C_number = size(C,1);
end
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

if ~isempty(C)
    if C_number == 1
        T = C*B./sqrt(sigma2*(C*inv(X'*X)*C'));
    else
        T = C*B./sqrt( repmat(sigma2,[C_number,1]) .* diag((C*inv(X'*X)*C')) );
    end

    % P
    P = zeros(size(T));
    P(T>0) = tcdf(-T(T>0),n-rang);
    P(T<0) = tcdf(T(T<0),n-rang);
    %P =  2 * tcdf(-abs(T), n-rang); %which one? use a flag to swith

    % Calculatin Z-score. tcdf has a poor precision so z-score is limited to
    % something like 8. 
    Z = zeros(size(T));
    Z(T>0) = sqrt(2).*erfcinv(2*P(T>0));
    Z(T<0) = -sqrt(2)*erfcinv(2*P(T<0));

    %reshape
    if n_dimension == 3
        STATS.B = nan(C_number,M,M);
        STATS.T = nan(C_number,M,M);
        STATS.P = nan(C_number,M,M);
        STATS.Z = nan(C_number,M,M);
        for j=1:C_number
            STATS.B(j,:,:) = reshape(B(j,:),M,M);
            STATS.T(j,:,:) = reshape(T(j,:),M,M);
            STATS.P(j,:,:) = reshape(P(j,:),M,M);
            STATS.Z(j,:,:) = reshape(Z(j,:),M,M);
        end
    elseif n_dimension == 2
        STATS.B = B;
        STATS.T = T;
        STATS.P = P;
        STATS.Z = Z;
    end
else
    STATS = [];
end


return
end