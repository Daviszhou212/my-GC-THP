function [V_BBs, W_svds, HG_est, Sg_real, r_RZF, statistics_RZF] = calG2(H_eff, R_smoothed, idx, sorted_group_idxs, Sg, nUE, nRB, S, Pt_By_N0, miu)
    
    nGroups = length(unique(idx));
    groups = cell(1, nGroups);
    for iter_group = 1:nGroups
        groups{iter_group} = find(idx == iter_group);
    end
    
    Sg_real = zeros(nGroups, nRB);
    V_BBs = zeros(S, S, nRB);
    HG_est = zeros(2, S, nUE, nRB);
    r_RZF = zeros(S, nRB);
    statistics_RZF = zeros(nUE, 2, nRB);
    W_svds = zeros(2, 2, nUE, nRB);

    for idxRB = 1:nRB
        isBadGroups = zeros(1, nGroups);
        for iter_UE = 1:nUE
            [w_svd,~,~] = svd(H_eff(:, :, iter_UE, idxRB));
            W_svds(:, :, iter_UE, idxRB) = w_svd;
            H_eff(:, :, iter_UE, idxRB) = w_svd' * H_eff(:, :, iter_UE, idxRB);
        end
        H_n = H_eff(:, :, :, idxRB);

        RxPorts = zeros(nUE, 2);   % 2*UE ports (0:unscheduled / 1:otherwise)
        Hs = zeros(S, S);          % Effective channel of scheduled ports
        pfs_Factor = zeros(1, S);  % PFS-factor of scheduled ports

        idxTxStream = 1;
        UserIdxsOfStreams = zeros(1, S);
        
        pre_sum_rate = 0;
        for iter_Group = 1:nGroups
            idxGroup = sorted_group_idxs(iter_Group);
            
            if Sg(idxGroup) == 0
                continue;
            end

            for iter_Stream = 1 : Sg(idxGroup)

                [userIdxs_, portsIdxs_] = find(RxPorts == 0); % Find ports that have not been scheduled
                pfs_sum_rate = zeros(1, numel(userIdxs_));

                for i = 1:numel(userIdxs_)

                    if ~ismember(userIdxs_(i)+(portsIdxs_(i)==2)*nUE, groups{idxGroup})
                        continue % 该用户需要在当前分组
                    end

                    Hs(idxTxStream, :) = H_n(portsIdxs_(i), :, userIdxs_(i));
                    H_ = Hs(1:idxTxStream, :); % 当前已调度端口信道矩阵
                    pfs_Factor(idxTxStream) = 1 / (R_smoothed(userIdxs_(i)) + eps);

                    % Weighted sum rate
                    pfs_sum_rate(i) = calWSR(H_, Pt_By_N0, miu, pfs_Factor(1:idxTxStream));

                end

                [cur_sum_rate, i_] = max(pfs_sum_rate); % 贪婪调度新一端口，使得当前加权和速率最大
                if cur_sum_rate < pre_sum_rate
                    isBadGroups(idxGroup) = 1;
                    break
                end        
                
                RxPorts(userIdxs_(i_), portsIdxs_(i_)) = 1;
                Hs(idxTxStream, :) = H_n(portsIdxs_(i_), :, userIdxs_(i_));
                % pfs_Factor(idxTxStream) = Queues(userIdxs_(i_)) / (R_smoothed(userIdxs_(i_)) + eps);
                pfs_Factor(idxTxStream) = 	1 / (R_smoothed(userIdxs_(i_)) + eps);
                statistics_RZF(userIdxs_(i_), 1, idxRB) = ...
                    1 + statistics_RZF(userIdxs_(i_), 1, idxRB);
                UserIdxsOfStreams(idxTxStream) = userIdxs_(i_);

                idxTxStream = idxTxStream + 1;
                pre_sum_rate = cur_sum_rate;
                Sg_real(idxGroup, idxRB) = Sg_real(idxGroup, idxRB) + 1;
            end

        end

        nRemainedStreams = S - sum(Sg_real(:, idxRB));
        for iter_Stream = 1 : nRemainedStreams
            [userIdxs_, portsIdxs_] = find(RxPorts == 0); % Find ports that have not been scheduled
            pfs_sum_rate = zeros(1, numel(userIdxs_));
            for i = 1:numel(userIdxs_)

                if isBadGroups(idx(userIdxs_(i)+(portsIdxs_(i)==2)*nUE)) == 1
                    continue
                end

                Hs(idxTxStream, :) = H_n(portsIdxs_(i), :, userIdxs_(i));
                H_ = Hs(1:idxTxStream, :); % 当前已调度端口信道矩阵
                pfs_Factor(idxTxStream) = 1 / (R_smoothed(userIdxs_(i)) + eps);

                % Weighted sum rate
                pfs_sum_rate(i) = calWSR(H_, Pt_By_N0, miu, pfs_Factor(1:idxTxStream));
            end

            [cur_sum_rate, i_] = max(pfs_sum_rate); % 贪婪调度新一端口，使得当前加权和速率最大
            if cur_sum_rate < pre_sum_rate
                break
            end

            RxPorts(userIdxs_(i_), portsIdxs_(i_)) = 1;
            Hs(idxTxStream, :) = H_n(portsIdxs_(i_), :, userIdxs_(i_));
            pfs_Factor(idxTxStream) = 	1 / (R_smoothed(userIdxs_(i_)) + eps);
            statistics_RZF(userIdxs_(i_), 1, idxRB) = ...
                1 + statistics_RZF(userIdxs_(i_), 1, idxRB);
            UserIdxsOfStreams(idxTxStream) = userIdxs_(i_);

            idxTxStream = idxTxStream + 1;
            pre_sum_rate = cur_sum_rate;
            Sg_real(idx(userIdxs_(i_)+(portsIdxs_(i_)==2)*nUE), idxRB) = ...
                Sg_real(idx(userIdxs_(i_)+(portsIdxs_(i_)==2)*nUE), idxRB) + 1;
        end

        % Digital precoding (Giving schedulded users)
        nStreamsReal = sum(Sg_real(:,idxRB));
        G = calRzfPrecoder(Hs(1:nStreamsReal,:), Pt_By_N0, miu);
        V_BBs(:, :, idxRB) = [G, zeros(S,S-nStreamsReal)];

        % H_eff * G
        for idxUE = 1:nUE
            HG_est(:, :, idxUE, idxRB) = [H_eff(:, :, idxUE, idxRB) * G, zeros(2,S-nStreamsReal)];
        end
        
        tmp = abs(Hs(1:nStreamsReal,:) * G) .^ 2;
        r_RZF(:, idxRB) = [log2(1 + diag(tmp) ./ (1 + sum(tmp, 2) - diag(tmp))); zeros(S-nStreamsReal,1)];
        for iter_Stream = 1 : nStreamsReal
            statistics_RZF(UserIdxsOfStreams(iter_Stream), 2, idxRB) = ...
                statistics_RZF(UserIdxsOfStreams(iter_Stream), 2, idxRB) + ...
                r_RZF(iter_Stream, idxRB);
        end

    end
end