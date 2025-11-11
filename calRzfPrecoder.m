function G = calRzfPrecoder(H, PT, miu)
    nStreams = size(H, 1);

    % RZF precoder
    G = H' * pinv(H * H' + miu * eye(nStreams));
    P = diag(PT / nStreams ./ diag(G' * G));
    G = G * sqrt(P);
end