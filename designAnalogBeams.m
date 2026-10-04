function [selectedRxBeams, lamda] = designAnalogBeams(H, S, Tx_CB, Rx_CB)
    [~, ~, nRB, nUE] = size(H);
    Tx_CB_size = size(Tx_CB, 2);
    Rx_CB_size = size(Rx_CB, 2);

    selectedRxBeams = zeros(nUE, 1);                                    % Selected Tx & Rx beams
    lamda = zeros(Tx_CB_size, nUE);                                     % RSRP energy (Given Rx Beam)

    for k = 1:nUE
        BeamPowers = zeros(Tx_CB_size, Rx_CB_size);                     % RSRP energy

        for t = 1:Tx_CB_size
            x = repmat(Tx_CB(:, t), S, 1); % Tx scanning beam (All sub-arrays adopt the same beam)

            for r = 1:Rx_CB_size
                W_RF = blkdiag(Rx_CB(:, r), Rx_CB(:, r)); % analog combiner

                % Summing power over different subcarriers (Only one sampling from each RB)
                for iter_RB = 1:nRB
                    H_k = H(:, :, iter_RB, k);
                    BeamPowers(t, r) = BeamPowers(t, r) + sum(abs(W_RF' * H_k * x).^2);
                end
            end
        end

        [~, selectedRxBeams(k)] = find(BeamPowers == max(max(BeamPowers)));
% %         [~, selectedRxBeams(k)] = max(sum(BeamPowers));
        lamda(:, k) = BeamPowers(:, selectedRxBeams(k)) / nRB / 2;
    end
end