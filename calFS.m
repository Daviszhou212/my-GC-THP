% 计算每组用户个数矢量 Kv
Kv = zeros(nGroups, 1);

for i = 1:numel(idx)
    Kv(idx(i)) = Kv(idx(i)) + 1;
end

%%%%%%%%%%%%%%initilization
weight = ones(sum(Kv), 1);
[Sv, F, r_propose] = CalS(idx, Kv, Ms, S, Tx_CB, Rcron, PT, miu1, weight, rankRcron);
Sv_tlt(:, 1) = Sv;
F_tlt(:, :, 1) = F;
r_tlt(1, :) = r_propose;
kL = 1;
weight_rate = (ones(kL, 1) / kL)' * r_tlt;

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
%%% 迭代优化速率最大化
while (1)
    % cal qo
    r_imagUser = zeros(size(r_tlt));
    for kj=1:nGroups
        r_imagUser(:,idx==kj) = r_tlt(:, sum(Kv(1:kj-1))+1:sum(Kv(1:kj)));
    end
    r_realUser = r_imagUser(:,1:nUE) + r_imagUser(:,nUE+1:2*nUE);
    qo = fmincon(@(q)CalUq(q, r_realUser), ones(kL, 1) / kL, [], [], ones(1, kL), 1, zeros(kL, 1), ones(kL, 1), []);
    % 计算加权速率
    weight_rate = qo' * r_tlt;
    rateU = CalUq(qo, r_realUser);
    disp(rateU - rateU_pre)

    if iter > 20 || abs(rateU - rateU_pre) < 10^ - 4
        break;
    end

    rateU_pre = rateU;
    [qval, qindx] = sort(qo, 'descend');
    kL = length(find(qo >= 10^ - 4));
    Sv_tlt = Sv_tlt(:, qindx(1:kL));
    F_tlt = F_tlt(:, :, qindx(1:kL));
    r_tlt = r_tlt(qindx(1:kL), :);

    % weighted
    weight = 1 ./ (weight_rate.' + 10^ - 4);
    [Sv, F, r_propose] = CalS(idx, Kv, Ms, S, Tx_CB, Rcron, PT, miu1, weight, rankRcron);
    kL = kL + 1;
    Sv_tlt(:, kL) = Sv;
    F_tlt(:, :, kL) = F;
    r_tlt(kL, :) = r_propose;
    iter = iter + 1;
end