function [idx, C] = clustering(X, thresold)
    %%% Clustering
    nUE = size(X, 2);
    idx =  zeros(nUE, 1);
    
    %%% First layer
    X_norm = zeros(size(X));
    for k = 1:nUE
        X_norm(:,k) = X(:,k) / norm(X(:,k));
    end
    [idx1, C1] = kmeansth(X_norm, 0.4, nUE);
    
    %%% Second layer
    for k = 1:nUE
        X_A(1,k) =  norm(X(:,k));
    end
    X_B=10*log10(X_A./min(X_A));
    for i = min(idx1):max(idx1)
        if sum(idx1==i) > 1
            [idx2, CA] = kmeansth(X_B(:,idx1==i), thresold, sum(idx1==i));
        else
            idx2 = 1;
        end
        idx(idx1==i) = max(idx) + idx2;
    end
    
    if any(idx <= 0)
        error('Clustering failed');
    end
    
    %%% Center point
    nGroups = numel(unique(idx));
    C = zeros(nGroups, size(X, 1));
    for iter_group = 1 : nGroups
        C(iter_group,:) = sum(X(:,idx==iter_group), 2) / sum(idx==iter_group);
    end
end
