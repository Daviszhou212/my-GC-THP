function [selectedRxBeams, lamda] = designAnalogBeams(H, S, Tx_CB, Rx_CB)
% --- 模拟长时尺度波束搜索 (对应论文 Section III-B) ---
    %
    % 功能:
    % 1. 模拟所有Tx-Rx波束对的扫描过程，计算RSRP (参考信号接收功率)。
    % 2. 为每个用户(UE)从Rx码本中选择一个功率最强的固定模拟接收波束(Analog Combiner)。
    % 3. 计算并输出RSRP向量(lamda)，即B-CSTI (波束域信道统计信息)，
    %    该信息将用于后续的用户分群(Clustering)
    [~, ~, nRB, nUE] = size(H);
    Tx_CB_size = size(Tx_CB, 2);
    Rx_CB_size = size(Rx_CB, 2);

    selectedRxBeams = zeros(nUE, 1);                                    % Selected Tx & Rx beams
    lamda = zeros(Tx_CB_size, nUE);                                     % RSRP energy (Given Rx Beam)

    for k = 1:nUE
        BeamPowers = zeros(Tx_CB_size, Rx_CB_size);                     % RSRP energy
    % 对于每个用户，测试所有的Tx与Rx波束的组合
        for t = 1:Tx_CB_size
            x = repmat(Tx_CB(:, t), S, 1); % Tx scanning beam (All sub-arrays adopt the same beam)

            for r = 1:Rx_CB_size
                W_RF = blkdiag(Rx_CB(:, r), Rx_CB(:, r)); % analog combiner

                % Summing power over different subcarriers (Only one sampling from each RB)
                for iter_RB = 1:nRB
                    H_k = H(:, :, iter_RB, k);
                    % 计算接收功率
                    BeamPowers(t, r) = BeamPowers(t, r) + sum(abs(W_RF' * H_k * x).^2);
                end
            end
        end
        % 寻找最大接收功率的波束组合，记录其Rx波束索引
        % 在后续过程中，将被用作该用户固定不变的模拟合并器
        [~, selectedRxBeams(k)] = find(BeamPowers == max(max(BeamPowers)));
% %         [~, selectedRxBeams(k)] = max(sum(BeamPowers));
        % 计算RSRP
        lamda(:, k) = BeamPowers(:, selectedRxBeams(k)) / nRB / 2;
    end
end