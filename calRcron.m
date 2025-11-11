function Rcron = calRcron(lamda, Tx_CB)
% 根据BS码本与RSRP计算平均空间协方差矩阵
    m = size(lamda, 1);
    [n, ~] = size(Tx_CB);
    Rcron = zeros(n, n);

    for i = 1:m
        Rcron = Rcron + lamda(i) * Tx_CB(:, i) * Tx_CB(:, i)';
    end

end