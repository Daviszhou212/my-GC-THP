% 执行long time-scale optimization
% 计算每组包含的用户个数矢量 Kv
Kv = zeros(nGroups, 1);
for i = 1:numel(idx)
    Kv(idx(i)) = Kv(idx(i)) + 1;
end
% 这个文件中包含tlt的变量就是记录各个长时方案的
%% 初始化
% 是为了找到一个比较好的初始变量
weight = ones(sum(Kv), 1);
[Sv, F, r_propose] = CalS(idx, Kv, Ms, S, Tx_CB, Rcron, PT, miu1, weight, rankRcron);
Sv_tlt(:, 1) = Sv;
F_tlt(:, :, 1) = F;
r_tlt(1, :) = r_propose;
kL = 1;
weight_rate = (ones(kL, 1) / kL)' * r_tlt;
% 确保每个用户都有一个非零速率，也就是说每一个用户组至少被一个{F,Sv}调度了
while ~isempty(find(weight_rate == 0, 1))
    kL = kL + 1;
    weight(weight_rate > 0) = 0;
    [Sv, F, r_propose] = CalS(idx, Kv, Ms, S, Tx_CB, Rcron, PT, miu1, weight, rankRcron);
    Sv_tlt(:, kL) = Sv;
    F_tlt(:, :, kL) = F;
    r_tlt(kL, :) = r_propose;
    weight_rate = (ones(kL, 1) / kL)' * r_tlt;
end
iter = 1;
rateU_pre = 0;
%% 迭代优化速率最大化
while (1)
    % cal qo
    r_imagUser = zeros(size(r_tlt));
    for kj=1:nGroups
        r_imagUser(:,idx==kj) = r_tlt(:, sum(Kv(1:kj-1))+1:sum(Kv(1:kj)));
    end
    r_realUser = r_imagUser(:,1:nUE) + r_imagUser(:,nUE+1:2*nUE);
    % 基于当前的{F,Sv}的集合，使用fmincon优化效用函数
    % qo = fmincon(@(q)CalUq(q, r_realUser), ones(kL, 1) / kL, [], [], ones(1, kL), 1, zeros(kL, 1), ones(kL, 1), []);
    % 使用cvx替代：
    cvx_begin quiet % 'quiet' 模式可以抑制求解过程的输出
        cvx_solver Mosek_5
        variable qo(kL) % 定义 q 是一个 kL x 1 的变量
        rate_per_user = qo.' * r_realUser;
        maximize( geo_mean(rate_per_user + 10^-6) )
        subject to
            sum(qo) == 1;
            qo >= 0;
    cvx_end
    % 计算加权速率
    weight_rate = qo' * r_tlt;
    rateU = CalUq(qo, r_realUser);
    fprintf('PF效用函数增加%f\n',-rateU + rateU_pre);
    % 收敛检查
    if iter > 20 || abs(rateU - rateU_pre) < 10^ - 4
        break;
    end
    rateU_pre = rateU;
    % 剪枝，概率过低的{F,Sv}被删去
    [qval, qindx] = sort(qo, 'descend');
    kL = length(find(qo >= 10^ - 4));
    Sv_tlt = Sv_tlt(:, qindx(1:kL));
    F_tlt = F_tlt(:, :, qindx(1:kL));
    r_tlt = r_tlt(qindx(1:kL), :);

    % weighted
    % 得到新的PF权重，然后根据新权重，找到一个新的{F,Sv}
    weight = 1 ./ (weight_rate.' + 10^ - 4);
    [Sv, F, r_propose] = CalS(idx, Kv, Ms, S, Tx_CB, Rcron, PT, miu1, weight, rankRcron);
    kL = kL + 1;
    Sv_tlt(:, kL) = Sv;
    F_tlt(:, :, kL) = F;
    r_tlt(kL, :) = r_propose;
    iter = iter + 1;
end