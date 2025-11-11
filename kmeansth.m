function [idxglobal, Cglobal] = kmeansth(X, thresold, K)

% % %     X = X.';
% % %     idx = [1*ones(3,1); 2*ones(2,1); 3*ones(3,1); 4*ones(2,1)];
% % %     nGroups = numel(unique(idx));
% % %     C = zeros(nGroups, size(X, 2));
% % %     for iter_group = 1 : nGroups
% % %         C(iter_group,:) = sum(X(idx==iter_group,:)) / sum(idx==iter_group);
% % % % %         X_ = X(idx==iter_group, :);
% % % % %         count = size(X_, 1);
% % % % %         amplitude = 0;
% % % % %         for i = 1 : count
% % % % %             amplitude = amplitude + norm(X_(i,:));
% % % % %             C(iter_group, :) = C(iter_group, :) + ...
% % % % %                 X_(i,:) / norm(X_(i,:));
% % % % %         end
% % % % %         C(iter_group, :) = amplitude/count * C(iter_group, :) / norm(C(iter_group, :));
% % %     end
%%% kmeans(X,mid, "Replicates",2);
    Nkmean=10;
    X = X.';
    left = 1;
    % M -- allusers
    right = K;

    %thresold = thresold * mean((sum(X.^2, 2)));
    distglobal=10^10;
    numgroup=K;
    %C---K*P的矩阵，存储的是K个聚类质心位置
    while left < right
        mid = max(fix(left + (right - left) / 2), 1);
        distopt=10^10;
        for nkm=1:Nkmean
            [idxcur,Ccur,~,Dcur] = kmeans(X,mid);
            if mid==1
                distcur=max(Dcur);
            else
                distcur=max(min(Dcur.'));
            end
            if distcur<distopt
                idx=idxcur;
                C=Ccur;
                D=Dcur;
                distopt=distcur;
            end
        end
        if (distopt < thresold || distopt < distglobal) && (numgroup>=mid || distglobal>thresold)
            Cglobal=C;
            Dglobal=D;
            idxglobal=idx;
            distglobal=distopt;
            numgroup=mid;
        end
        if distopt <= thresold
            right = mid; 
        else
            % more group
            left = mid;
        end
        if right-left <= 1
            distopt=10^10;
            for nkm=1:Nkmean
                [idxcur,Ccur,~,Dcur] = kmeans(X,left);
                if left==1
                    distcur=max(Dcur);
                else
                    distcur=max(min(Dcur.'));
                end
                if distcur<distopt
                    idx=idxcur;
                    C=Ccur;
                    D=Dcur;
                    distopt=distcur;
                end
            end
            if (distopt < thresold || distopt < distglobal) && (numgroup>=left || distglobal>thresold)
                Cglobal=C;
                Dglobal=D;
                idxglobal=idx;
                distglobal=distopt;
                numgroup=left;
            end
            if distopt <= thresold
                break;
            else
                distopt=10^10;
                for nkm=1:Nkmean
                    [idxcur,Ccur,~,Dcur] = kmeans(X,right);
                    if right==1
                        distcur=max(Dcur);
                    else
                        distcur=max(min(Dcur.'));
                    end
                    if distcur<distopt
                        idx=idxcur;
                        C=Ccur;
                        D=Dcur;
                        distopt=distcur;
                    end
                end
                if (distopt < thresold || distopt < distglobal) && (numgroup>=right || distglobal>thresold)
                    Cglobal=C;
                    Dglobal=D;
                    idxglobal=idx;
                    distglobal=distopt;
                end
                break;
            end
        end
    end
    
    nGroups = length(unique(idxglobal));
    if nGroups==1
                        dtest=max(Dglobal);
                    else
                        dtest=max(min(Dglobal.'));
    end
             
   
    if dtest >  thresold
        disp('why larger than threshold');
    end

