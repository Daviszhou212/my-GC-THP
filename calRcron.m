function Rcron = calRcron(lamda, Tx_CB)
    m = size(lamda);
    [n, ~] = size(Tx_CB);
    Rcron = zeros(n, n);

    for i = 1:m
        Rcron = Rcron + lamda(i) * Tx_CB(:, i) * Tx_CB(:, i)';
    end

end