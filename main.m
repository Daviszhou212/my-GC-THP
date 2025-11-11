clear; close all; clc;
tic;
rng(40);

%%% Parameters setting
run("parameters.m")

% %%% Variables for statistics
% R_All = zeros(nSlotPerFrame * nTTiPerSlot, nUE);
% SINR_p_s = [];
% 
% %% % Simulation
% %% -1 Codebook design DFT Codebook
% [Tx_CB, Tx_CB_size] = dftCodebook(Tx_M, Tx_N, 15, 60);
% % Rx DFT Codebook for Combining
% [Rx_CB, Rx_CB_size] = dftCodebook(Rx_M, Rx_N, 90, 90);
% 
% 
% %% -2 Calculating RSRP & Selecting Tx/Rx Beams
% load([channel_path, '1_1'], 'H'); % Load the Channel Data
% % load('1_1.mat', 'H');
% [selectedRxBeams, lamda] = designAnalogBeams(H, S, Tx_CB, Rx_CB);
% % selectedRxBeams => (nUE, 1)
% % lamda => (Tx_CB_size, 30)
% 
% %% -3 clustering
% thresold = inf;
% [idx, C] = clustering(lamda,thresold);
% % C代表每个用户组的平均RSRP向量=>(J, Tx_CB_size)
% % idx表示各个用户对应的组号
% % 由于每个用户最多两流，所以扩展为nUE*2个虚拟单天线用户
% idx = repmat(reshape(idx, [], 1), [2,1]);
% 
% %% -4 Calculating the channel correlation matrix
% nGroups = length(unique(idx));
% rankRcron = zeros(nGroups, 1);
% 
% for ki = 1:nGroups
%     % 对于每一个用户组，计算平均信道空域相关矩阵
%     Rcron(:, :, ki) = calRcron(C(ki, :).', Tx_CB);
% % %     rankRcron(ki) = rank(Rcron(:, :, ki));
%     rankRcron(ki) = 2;
%     % 强制规定每个用户最多2流
% end
% 
% %% 5 Optimizing analog precoding and scheduling schemes
% run("calFS.m")
% save('single_lamda_scatter.mat')
load("single_lamda_scatter.mat")


%% 6 Digital precoding & TTI Scheduling given optimization result
r_imagUser = zeros(size(r_tlt));
for kj=1:nGroups
    r_imagUser(:,idx==kj) = r_tlt(:, sum(Kv(1:kj-1))+1:sum(Kv(1:kj)));
end
% 真实的用户DE速率
r_realUser = r_imagUser(:,1:nUE) + r_imagUser(:,nUE+1:2*nUE);
% % R_smoothed = zeros(1, nUE); % The exponentially smoothed average data rate over the past Tc TTIs
% % R_smoothed = qo' * r_realUser * nRB * nCarrierPerRB * nSymbolsPerTTI / LenOfTTI;
R_smoothed = 3.1364e+07 * ones(1, nUE);
Queues = pSize * poissrnd(10 * LenOfTTI / pInterval, 1, nUE);


for iter_slot = 1:nSlotPerFrame
    %% 6-1 Digital precoding

    Sg_real = zeros(nGroups, nRB, kL);
    V_RFs = zeros(nTx, S, kL);
    V_BBs = zeros(S, S, nRB, kL);
    HG_est = zeros(2, S, nUE, nRB, kL); % Estimated H_eff * G
    W_svd = zeros(2, 2, nUE, nRB, kL);
    r_RZF = zeros(S, nRB, kL);
    statistics_RZF = zeros(nUE, 2, nRB, kL);

    for timeSharing_idx = 1:kL
        % Sorting the groups
        r_weighted_group = zeros(1, nGroups);
        for iter_UE = 1:2*nUE
            r_weighted_group(idx(iter_UE)) = r_weighted_group(idx(iter_UE)) + ...
                r_imagUser(timeSharing_idx, iter_UE) / (R_smoothed(mod(iter_UE-1, nUE)+1) + eps);
        end
% %         for iter_UE = 1:nUE
% %             r_weighted_group(idx(iter_UE)) = r_weighted_group(idx(iter_UE)) + ...
% %                 r_realUser(timeSharing_idx, iter_UE) / (R_smoothed(iter_UE) + eps);
% %         end
        [~, sorted_group_idxs] = sort(r_weighted_group, 'descend');
        % sorted_group_idxs对不同用户组赋予了优先级
        % Analog precoder 部分连接
        F = [];
        for i = 1:S
            F = blkdiag(F, F_tlt(:, i, timeSharing_idx));
        end
        V_RFs(:,:,timeSharing_idx) = F;

        % 有效信道计算
        load([channel_path, num2str(mod(iter_slot-1,4)+1), '_1'], 'H'); % Load the Channel Data
        % load('1_1.mat', 'H');
        H_eff = zeros(2, S, nUE, nRB);
        for idxUE = 1:nUE
            W_RF = blkdiag(Rx_CB(:, selectedRxBeams(idxUE)), Rx_CB(:, selectedRxBeams(idxUE)));
            for idxRB = 1:nRB
                H_eff(:, :, idxUE, idxRB) = W_RF' * H(:, :, idxRB, idxUE) * F;
            end
        end

        % Digital precoding
        % 调用calG2函数，为每一个方案提前计算短时调度
        [V_BBs(:,:,:,timeSharing_idx), W_svd(:,:,:,:,timeSharing_idx), HG_est(:,:,:,:,timeSharing_idx), ...
            Sg_real(:,:,timeSharing_idx), r_RZF(:,:,timeSharing_idx), statistics_RZF(:,:,:,timeSharing_idx)] = ...
                calG2(H_eff, R_smoothed, idx, sorted_group_idxs, Sv_tlt(:, timeSharing_idx), nUE, nRB, S, Pt_By_N0, Ms*miu2);
    end


    %% 6-2 TTI scheduling
    timeSharing_count = zeros(kL, 1);
    r_rzf = squeeze(sum(squeeze(statistics_RZF(:,2,:,:)), 2))';
% %     qo = fmincon(@(q)CalUq(q, r_rzf), ones(kL, 1) / kL, [], [], ones(1, kL), 1, zeros(kL, 1), ones(kL, 1), []);
    for iter_TTI = 1:nTTiPerSlot
% %         timeSharing_idx = randsrc(1, 1, [1:kL; reshape(qo, [1, kL])]);
% %         while (timeSharing_count(timeSharing_idx) / (sum(timeSharing_count) + eps)) > qo(timeSharing_idx)
% %             timeSharing_idx = randsrc(1, 1, [1:kL; reshape(qo, [1, kL])]);
% %         end
        wsr = sum(r_rzf./(R_smoothed+eps), 2);
        [~,timeSharing_idx] = max(wsr);
        timeSharing_count(timeSharing_idx) = timeSharing_count(timeSharing_idx) + 1;

        F = V_RFs(:,:,timeSharing_idx);

        load([channel_path, num2str(mod(iter_slot-1,4)+1), '_', num2str(iter_TTI)], 'H'); % Load the Channel Data
        % load('1_1.mat', 'H');

        Queues = Queues + pSize * poissrnd(LenOfTTI / pInterval, 1, nUE);
        pfs_Factor = 1 ./ (R_smoothed + eps)';
        Data_currentTTI = zeros(1, nUE); % The achievable rate of current TTI
        
        % 逐个RB调度
        for idxRB = 1:nRB
            % 获取指定{F,Sv}方案在该RB上的总流数
            nStreamsReal = sum(Sg_real(:,idxRB,timeSharing_idx));

            %%% User scheduling (Using estimated channel)
            RxPorts = zeros(nUE, 2); % 2*UE ports (0:unscheduled / 1:otherwise)
            % 被调度了就设置为1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Different from the analog [RxPorts] above, the [RxPorts] here are digital ports
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 根据各个射频链在各RB上的RZF速率排序，按照链路逐个贪婪选取
            [~, TxStream_idxs] = sort(r_RZF(1:nStreamsReal, idxRB, timeSharing_idx), 'descend');
            for iter_TxStream = 1:nStreamsReal
                idxTxStream = TxStream_idxs(iter_TxStream);
                
                [userIdxs_, portsIdxs_] = find(RxPorts == 0); % Find ports that have not been scheduled

                data_est = zeros(1, numel(userIdxs_)); % Use the measured aging channel to estimate the rate
                W_BBs = zeros(1, 2, numel(userIdxs_));
                
                % 在该链路上，利用slot的第一个信道遍历预估所有用户的速率 
                for i = 1:numel(userIdxs_)
                    HG_ = HG_est(:, 1:nStreamsReal, userIdxs_(i), idxRB, timeSharing_idx);
                    W_BB = HG_(:, idxTxStream)' * pinv(HG_ * HG_' + eye(2));
                    % W_BB = (portsIdxs_(i)==1)*[1, 0] + (portsIdxs_(i)==2)*[0, 1];

                    tmp_ = W_BB * HG_;
                    tmp_ = abs(tmp_).^2;
                    noise = trace(W_BB * W_BB');
                    data_est(i) = log2(1 + tmp_(idxTxStream) / (noise + sum(tmp_) - tmp_(idxTxStream)));
                    W_BBs(:, :, i) = W_BB;
                end

                % data_est = min(nSymbolsPerTTI * nCarrierPerRB * data_est, Queues(userIdxs_(i)));
                
                % 将该流贪婪分配给PF效用函数最高的用户
                [~, i_] = max(data_est ./ (R_smoothed(userIdxs_) + eps)); % Greedy scheduling according to PFS criterion
                if data_est(i_) < r_RZF(idxTxStream, idxRB, timeSharing_idx)
                    [~, i_] = max(data_est);
                end

                % 标记已调用
                RxPorts(userIdxs_(i_), portsIdxs_(i_)) = 1;

                %%% Calculating the rate (Using realistic channel)
                k_ = userIdxs_(i_);
                W_RF = blkdiag(Rx_CB(:, selectedRxBeams(k_)), Rx_CB(:, selectedRxBeams(k_)));
                HG_ = W_svd(:,:,k_,idxRB,timeSharing_idx)' * W_RF' * H(:, :, idxRB, k_) * F * V_BBs(:, 1:nStreamsReal, idxRB, timeSharing_idx);
                % % W_BB = W_BBs(:, :, i_);
                W_BB = HG_(:, idxTxStream)' * pinv(HG_ * HG_' + eye(2));
                tmp_ = W_BB * HG_;
                tmp_ = abs(tmp_).^2;
                noise = trace(W_BB * W_BB');
                Data_currentTTI(k_) = Data_currentTTI(k_) + ...
                    log2(1 + tmp_(idxTxStream) / (noise + sum(tmp_) - tmp_(idxTxStream)));
                SINR_p_s(end+1) = tmp_(idxTxStream) / (noise + sum(tmp_) - tmp_(idxTxStream));
            end

        end
        
% %         disp([r_rzf(timeSharing_idx,:); Data_currentTTI]');
        r_rzf(timeSharing_idx,:) = 0.5 * r_rzf(timeSharing_idx,:) + 0.5 * Data_currentTTI;

        Data_currentTTI = nSymbolsPerTTI * nCarrierPerRB * Data_currentTTI;
        % Data_currentTTI = min(Data_currentTTI, Queues);
        R_currentTTI = Data_currentTTI / LenOfTTI;
        Queues = Queues - Data_currentTTI;
        R_smoothed = (Tc - 1) / Tc * R_smoothed + 1 / Tc * R_currentTTI;

        R_All((iter_slot-1)*nTTiPerSlot+iter_TTI, :) = R_currentTTI;
    end

end


disp(sum(mean(R_All)) / 1e6);
disp(sum(log(mean(R_All))));
toc;

save('SINR_p_s.mat', 'SINR_p_s');
r_p_s = mean(R_All) / 1e6;
save('r_p_s.mat', 'r_p_s');