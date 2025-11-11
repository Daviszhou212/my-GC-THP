clear; close all; clc;
tic;
rng(40);

%% Parameters setting
run("parameters.m")

%% % Simulation
%% -1 Codebook design DFT Codebook
[Tx_CB, Tx_CB_size] = dftCodebook(Tx_M, Tx_N, 15, 60);
% Rx DFT Codebook for Combining
[Rx_CB, Rx_CB_size] = dftCodebook(Rx_M, Rx_N, 90, 90);


%% -2 Calculating RSRP & Selecting Tx/Rx Beams
load([channel_path, '1_1'], 'H'); % Load the Channel Data
% load('1_1.mat', 'H');
[selectedRxBeams, lamda] = designAnalogBeams(H, S, Tx_CB, Rx_CB);
% selectedRxBeams => (nUE, 1)
% lamda => (Tx_CB_size, 30)

%% -3 clustering
thresold = inf;
[idx, C] = clustering(lamda,thresold);
% C代表每个用户组的平均RSRP向量=>(J, Tx_CB_size)
% J = nGropus
% idx表示各个用户对应的组号=>(nUE,1)
% 由于每个用户最多两流，所以扩展为nUE*2个虚拟单天线用户[idx;idx]
idx = repmat(reshape(idx, [], 1), [2,1]);

%% -4 Calculating the channel correlation matrix
nGroups = length(unique(idx));
rankRcron = zeros(nGroups, 1);
for ki = 1:nGroups
    % 对于每一个用户组，计算平均信道空域相关矩阵
    Rcron(:, :, ki) = calRcron(C(ki, :).', Tx_CB);
    % %     rankRcron(ki) = rank(Rcron(:, :, ki));
    rankRcron(ki) = 2;
    % 强制规定每个用户最多2流
end
fprintf('开始长时优化\n');
%% 5 Optimizing analog precoding and scheduling schemes
run("calFS.m")
save('single_lamda_scatter.mat')
% load('single_lamda_scatter.mat')
%% 6 Digital precoding & TTI Scheduling given optimization result
r_imagUser = zeros(size(r_tlt));
for kj=1:nGroups
    r_imagUser(:,idx==kj) = r_tlt(:, sum(Kv(1:kj-1))+1:sum(Kv(1:kj)));
end
% 用户DE速率
R_DE = r_imagUser(:,1:nUE) + r_imagUser(:,nUE+1:2*nUE);
% % R_smoothed = zeros(1, nUE); % The exponentially smoothed average data rate over the past Tc TTIs
R_smoothed = qo' * r_realUser * nRB * nCarrierPerRB * nSymbolsPerTTI / LenOfTTI;
% % R_smoothed = 3.1364e+07 * ones(1, nUE);
Queues = pSize * poissrnd(10 * LenOfTTI / pInterval, 1, nUE);

fprintf('开始短时调度\n');
for iter_slot = 1:nSlotPerFrame
    %% 6 TTI scheduling
    timeSharing_count = zeros(kL, 1);
    for iter_TTI = 1:nTTiPerSlot
        timeSharing_idx = randsrc(1, 1, [1:kL; reshape(qo, [1, kL])]);
        while (timeSharing_count(timeSharing_idx) / (sum(timeSharing_count) + eps)) > qo(timeSharing_idx)
            timeSharing_idx = randsrc(1, 1, [1:kL; reshape(qo, [1, kL])]);
        end
        % wsr = sum(r_rzf./(R_smoothed+eps), 2);
        % [~,timeSharing_idx] = max(wsr);
        timeSharing_count(timeSharing_idx) = timeSharing_count(timeSharing_idx) + 1;
        timeSharing_record(iter_slot, iter_TTI) = timeSharing_idx;
        F = [];
        for i = 1:S
            F = blkdiag(F, F_tlt(:, i, timeSharing_idx));
        end

        % load([channel_path, num2str(mod(iter_slot-1,4)+1), '_', num2str(iter_TTI)], 'H'); % Load the Channel Data
        load('../channel/1_1.mat', 'H');
        % 计算有效信道
        H_eff = zeros(2, S, nUE, nRB);
        for idxUE = 1:nUE
            W_RF = blkdiag(Rx_CB(:, selectedRxBeams(idxUE)), Rx_CB(:, selectedRxBeams(idxUE)));
            for idxRB = 1:nRB
                % H_eff使用真实的H和F计算
                H_eff(:, :, idxUE, idxRB) = W_RF' * H(:, :, idxRB, idxUE) * F;
            end
        end
        r_weighted_group = zeros(1, nGroups);
        for iter_UE = 1:2*nUE
            r_weighted_group(idx(iter_UE)) = r_weighted_group(idx(iter_UE)) + ...
                r_imagUser(timeSharing_idx, iter_UE) / (R_smoothed(mod(iter_UE-1, nUE)+1) + eps);
        end
        [~, sorted_group_idxs] = sort(r_weighted_group, 'descend');
        % [~, sorted_group_idxs] = sort(ones(1, nGroups), 'descend'); % <--- 使用一个"无信息"的排序
        % sorted_group_idxs = randperm(nGroups);
        %[V_BBs, W_svd, HG_est, Sg_real, r_RZF, statistics_RZF] = ...
        %     calG2(H_eff, R_smoothed, idx, sorted_group_idxs, Sv_tlt(:, timeSharing_idx), nUE, nRB, S, Pt_By_N0, Ms*miu2);
        [V_BBs, W_svd, HG_est, Sg_real, r_RZF, statistics_RZF] = ...
           Random_calG2(H_eff, R_smoothed, idx, sorted_group_idxs, Sv_tlt(:, timeSharing_idx), nUE, nRB, S, Pt_By_N0, Ms*miu2);
        % 更新队列与PF因子
        Queues = Queues + pSize * poissrnd(LenOfTTI / pInterval, 1, nUE);
        pfs_Factor = 1 ./ (R_smoothed + eps)';

        % 从calG2函数的结果中获取R_real与R_DE比较
        R_real((iter_slot-1)*nTTiPerSlot+iter_TTI, :) = squeeze(mean(statistics_RZF(:,2,:), 3));

        % 转换为速率
        Data_currentTTI = squeeze(sum(statistics_RZF(:,2,:), 3)).';
        Data_currentTTI_bits = nSymbolsPerTTI * nCarrierPerRB * Data_currentTTI;
        R_currentTTI = Data_currentTTI_bits / LenOfTTI;
        Queues = Queues - Data_currentTTI_bits;
        R_smoothed = (Tc - 1) / Tc * R_smoothed + 1 / Tc * R_currentTTI;
    end

end
for j = 1:nGroups
    group_rate_DE(j) = qo' * mean(R_DE(:, find(idx(1:nUE)==j)),2);
    group_rate_DE_mat_J(j,:) = mean(R_DE(:, find(idx(1:nUE)==j)),2);
end

for j = 1:nGroups
    group_rate_real(j) = mean(mean(R_real(:, find(idx(1:nUE)==j)),2),1);
    group_rate_real_mat(j,:) = mean(R_real(:, find(idx(1:nUE)==j)),2);
end

for i=1:size(group_rate_real_mat, 2)
    group_rate_real_DE_mat(:, i) = group_rate_DE_mat_J(:,timeSharing_record(i));
end


r_DE_avg = (qo' * R_DE);

UserRealRates_avg = mean(R_real, 1);