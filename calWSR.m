function wsr = calWSR(H, PT, miu, pfs_Factor)
    nStreams = size(H, 1);

    % RZF precoder
    G = H' * pinv(H * H' + miu * eye(nStreams));
    P = diag(PT / nStreams ./ diag(G' * G));
    G = G * sqrt(P);

    tmpMatrix = H * G;
    tmpMatrix = abs(tmpMatrix) .^ 2;

    % Weighted sum rate
    wsr = 0;
    for iter_Stream = 1 : nStreams
        tmp_ = tmpMatrix(iter_Stream, :);
        wsr = wsr + pfs_Factor(iter_Stream) * ...
            log2(1 + tmp_(iter_Stream) / (1 + sum(tmp_) - tmp_(iter_Stream)));
    end
end