function [V_BBs, W_svds, HG_est, Sg_real, r_RZF, statistics_RZF] = Random_calG2(H_eff, R_smoothed_unused, idx, sorted_group_idxs_unused, Sg, nUE, nRB, S, Pt_By_N0, miu)
% calG_random: DE模型的基础调度器（随机调度）


nGroups = length(unique(idx));
groups = cell(1, nGroups);
for iter_group = 1:nGroups
    groups{iter_group} = find(idx == iter_group);
end

% 初始化输出
Sg_real = zeros(nGroups, nRB);
V_BBs = zeros(S, S, nRB);
HG_est = zeros(2, S, nUE, nRB);
r_RZF = zeros(S, nRB);
statistics_RZF = zeros(nUE, 2, nRB); % (:,1,:) 存流计数, (:,2,:) 存速率
W_svds = zeros(2, 2, nUE, nRB);

for idxRB = 1:nRB
    % --- 1. SVD合并器 (与calG2相同) ---
    % 这步是用户侧的合并器，与BS调度无关，必须保留
    isBadGroups = zeros(1, nGroups);
    for iter_UE = 1:nUE
        [w_svd,~,~] = svd(H_eff(:, :, iter_UE, idxRB));
        W_svds(:, :, iter_UE, idxRB) = w_svd;
        H_eff(:, :, iter_UE, idxRB) = w_svd' * H_eff(:, :, iter_UE, idxRB);
    end
    H_n = H_eff(:, :, :, idxRB);

    % --- 2. [核心修改] 随机用户选择 (替换贪婪逻辑) ---
    idxTxStream = 1;
    UserIdxsOfStreams = zeros(1, S); % 存储被调度流对应的 *真实* 用户索引 (1..nUE)
    Hs = zeros(S, S);                % 存储被调度流的有效信道
    Sg_real_temp = zeros(nGroups, 1);

    % 遍历所有组，执行长时尺度的 'Sg' 计划
    for j = 1:nGroups
        if Sg(j) > 0
            % 1. 找出该组的所有虚拟用户
            virtual_users_in_group = groups{j};

            % 2. 确定要调度的流数 (不能超过组内虚拟用户数)
            num_to_schedule = min(Sg(j), numel(virtual_users_in_group));

            if num_to_schedule == 0
                continue;
            end

            % 3. 随机选择 'num_to_schedule' 个虚拟用户
            % randperm(N, K) 从 N 个元素中随机选 K 个
            selected_indices = randperm(numel(virtual_users_in_group), num_to_schedule);
            selected_virtual_users = virtual_users_in_group(selected_indices);

            % 4. 将这些用户添加到调度集
            for k = 1:numel(selected_virtual_users)
                if idxTxStream > S
                    break; % 总RF链已用完，停止调度
                end

                v_idx = selected_virtual_users(k); % 虚拟用户索引 (1..2*nUE)

                % 5. 解码虚拟用户 -> (真实用户, 端口)
                if v_idx > nUE
                    real_user_idx = v_idx - nUE;
                    port_idx = 2;
                else
                    real_user_idx = v_idx;
                    port_idx = 1;
                end

                % 6. 存储调度结果
                Hs(idxTxStream, :) = H_n(port_idx, :, real_user_idx);
                UserIdxsOfStreams(idxTxStream) = real_user_idx;
                Sg_real_temp(j) = Sg_real_temp(j) + 1;

                % 更新流计数
                statistics_RZF(real_user_idx, 1, idxRB) = statistics_RZF(real_user_idx, 1, idxRB) + 1;

                idxTxStream = idxTxStream + 1;
            end
        end
        if idxTxStream > S
            break; % 总RF链已用完
        end
    end

    nStreamsReal = idxTxStream - 1;
    Sg_real(:, idxRB) = Sg_real_temp; % 记录本RB实际调度的流

    % --- 3. RZF预编码和速率计算 (与calG2相同) ---
    if nStreamsReal > 0
        % 计算RZF数字预编码器
        G = calRzfPrecoder(Hs(1:nStreamsReal,:), Pt_By_N0, miu); %
        V_BBs(:, :, idxRB) = [G, zeros(S,S-nStreamsReal)];

        % H_eff * G (用于TTI调度，在DE验证中非必须，但保持输出一致)
        for idxUE = 1:nUE
            HG_est(:, :, idxUE, idxRB) = [H_eff(:, :, idxUE, idxRB) * G, zeros(2,S-nStreamsReal)];
        end

        % 计算RZF速率 (bps/Hz)
        tmp = abs(Hs(1:nStreamsReal,:) * G) .^ 2;
        r_RZF(1:nStreamsReal, idxRB) = log2(1 + diag(tmp) ./ (1 + sum(tmp, 2) - diag(tmp)));

        % 累加速率到 statistics_RZF
        for iter_Stream = 1 : nStreamsReal
            statistics_RZF(UserIdxsOfStreams(iter_Stream), 2, idxRB) = ...
                statistics_RZF(UserIdxsOfStreams(iter_Stream), 2, idxRB) + ...
                r_RZF(iter_Stream, idxRB);
        end
    else
        % 如果此RB上没有调度流
        V_BBs(:, :, idxRB) = zeros(S,S);
    end
end
end