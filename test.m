R_direct = zeros(Ms, Ms, nGroups);
H_from_1_1 = H; % 明确我们使用的是 1_1.mat 的 H

for j = 1:nGroups
    % 找到属于第j组的真实用户索引
    users_in_group_j = find(idx(1:nUE) == j); 
    if isempty(users_in_group_j)
        continue;
    end
    
    R_group_j_sum = zeros(Ms, Ms); % 用于累加该组所有用户的SCM
    
    for k_idx = 1:length(users_in_group_j)
        k = users_in_group_j(k_idx); % 真实用户索引 (1..nUE)
        
        % H_k 的维度是 [nRx, nTx, nRB]
        % nTx = S * Ms
        H_k = squeeze(H_from_1_1(:,:,:,k)); 
        
        % 将H_k重塑为 [nRx, Ms, S, nRB] 以便分离各个子阵列
        H_k_reshaped = reshape(H_k, [nRx, Ms, S, nRB]);
        
        R_user_k_sum = zeros(Ms, Ms);
        % 对所有接收天线(nRx)、所有子阵列(S)、所有资源块(nRB)进行平均
        for rb = 1:nRB
            for s = 1:S
                for rx = 1:nRx
                    % h_vec 是一个 [1 x Ms] 的行向量
                    h_vec = squeeze(H_k_reshaped(rx, :, s, rb)); 
                    % 计算外积 (h_vec' * h_vec) 得到 [Ms x Ms] 矩阵
                    R_user_k_sum = R_user_k_sum + (h_vec' * h_vec);
                end
            end
        end
        
        % 计算用户k的平均SCM
        R_user_k_avg = R_user_k_sum / (nRB * S * nRx);
        R_group_j_sum = R_group_j_sum + R_user_k_avg;
    end
    
    % 计算组j的平均SCM
    R_direct(:,:,j) = R_group_j_sum / length(users_in_group_j);
end

% 1. 比较总功率 (Trace)
trace_Rcron = squeeze(arrayfun(@(i) trace(Rcron(:,:,i)), 1:nGroups));
trace_R_direct = squeeze(arrayfun(@(i) trace(R_direct(:,:,i)), 1:nGroups));
disp('比较 Trace (Rcron vs R_direct):');
disp([trace_Rcron, trace_R_direct]);

% 2. 比较矩阵相似度 (Frobenius 范数)
diff_norm = squeeze(arrayfun(@(i) norm(Rcron(:,:,i) - R_direct(:,:,i), 'fro'), 1:nGroups));
disp('矩阵差异 (Frobenius 范数):');
disp(diff_norm);