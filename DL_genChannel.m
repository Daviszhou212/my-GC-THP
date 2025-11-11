clear; close all; clc;

tic;
% rng(4);
rng(60);

% Load system parameters from parameters.m
parameters;


% Other Channel Parameters
Scenario = '3GPP_38.901_UMa_NLOS';   % 3GPP 38.901 UMa - 3GPP NR Urban Macro-Cell
                                    % For typical terrestrial base stations deployed above rooftop in densely populated urban
                                    % areas. Parameters cover LOS, NLOS and O2I. The maximum cell radius is about 1
                                    % km, the carrier frequency can be varied from 500 MHz to 100 GHz, the BS height is 25
                                    % meters and the MT height can vary from 1.5 to 22.5 m
                                    %
                                    % LOS
Scenario1 = '3GPP_38.901_UMa_LOS';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Failure: TOO LARGE !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Channel matrix (preallocating memory)
% H = repmat(1+1j, [nRx, nTx, nCarrierPerRB, nRB, nUE, nTTiPerSlot, nSlotPerFrame]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Failure: TOO LARGE !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel generation by QuaDRiGa
l = qd_layout;                                          % Create new QuaDRIGa layout
l.simpar.center_frequency = Fc;                         % Set center frequency
l.simpar.use_absolute_delays = 1;                       % Include delay of the LOS path
l.simpar.show_progress_bars = 0;                        % Disable progress bars

l.no_tx = 1;                                            % Set # of BSs
l.tx_array = qd_arrayant('3gpp-mmw', Tx_M, Tx_N, Fc, Tx_poInd, [], Tx_spacing, Tx_Mg, Tx_Ng);
l.tx_position = [0; 0; Tx_height];                      % Set BS's posittions
l.tx_track.orientation = [0; -Tx_downtilt; pi/6];       % Set BS's orientation
% h_build = qd_builder.gen_cdl_model('NR-CDL-A',Fc,UE_Speed/3.6);

l.no_rx = nUE;                                          % Set # of users
for iter_user = 1 : nUE
    l.rx_array(iter_user) = qd_arrayant('3gpp-mmw', Rx_M, Rx_N, Fc, Rx_poInd, [], Rx_spacing, Rx_Mg, Rx_Ng);

    t = qd_track('linear', (nSlotPerFrame*nTTiPerSlot-1)*LenOfTTI*UE_Speed/3.6, 2*pi*rand);
    t.name = ['Rx', num2str(iter_user)];

    % Generate random position
    % 5 cluster * 5 users (25 users in clusters)
    % and the other 5 users uniformly distributed in the cell
    if mod(iter_user-1, 5) == 0 || iter_user > 25
        x = [100, 100]; % Arbitrarily set, as long as the following condition is not met
        while (abs(x(1)) > 0.5*sqrt(3)) || (abs(x(1))/sqrt(3)+abs(x(2)) > 1)
            x = 2 * rand(1,2) - 1;      % Evenly distributed in the hexagonal cell
        end
    else
        x = x + 0.03 * randn(1, 2);
    end
    x1 = cellDiameter/2 * (x + [0.5*sqrt(3), 0.5]);
    
    t.initial_position = [x1(1); x1(2); Rx_height];

    t.set_speed(UE_Speed/3.6);
    t.interpolate('time', LenOfTTI, [], [], true);
    t.no_segments = t.no_snapshots;                     % Use spatial consisteny for mobility

    l.rx_track(iter_user) = t;
end


l.set_scenario(Scenario, [], [], 0);                     % Set scenario and SC_lambda
c = l.get_channels_seg(1, iter_user,1);
H = repmat(1+1j, [nRx, nTx, nRB, nUE]);   % Channel matrix of each single TTI (preallocating memory)
H_temp = repmat(1+1j, [nRx, nTx]);

% Progress tracking
total_files = nSlotPerFrame * nTTiPerSlot;
fprintf('开始生成信道数据...\n');
fprintf('总共需要生成 %d 个信道文件 (%d 时隙 × %d TTI)\n', total_files, nSlotPerFrame, nTTiPerSlot);
fprintf('用户数量: %d, 资源块数量: %d\n', nUE, nRB);
fprintf('========================================\n');

file_count = 0;
for iter_slot = (1 : nSlotPerFrame)
    fprintf('正在处理时隙 %d/%d...\n', iter_slot, nSlotPerFrame);
    for iter_TTI = 1 : nTTiPerSlot
        fprintf('正在处理TTI %d/%d...\n', iter_TTI, nTTiPerSlot);
        for iter_user = 1 : nUE

            c = l.get_channels_seg(1, iter_user, (iter_slot-1)*nTTiPerSlot+iter_TTI);  % Generate channel coefficients
            % h_ = c.fr(subcarrierSpacing*nCarrierPerRB*nRB, nRB*nCarrierPerRB);                       % Transforms the channel into frequency domain
                                                                                      
            h_ = c.fr(subcarrierSpacing*nCarrierPerRB*nRB, nRB);  %    (Only one sampling of each RB)                                                                          
            for iter_RB = 1 : nRB
 
                    h_n = h_(:,:,iter_RB);
                    
                    % Adjust the order
                    %
                    %               6  8     14 16      (H/V are different polalizations)
                    %               5  7     13 15      before: 1H 1V 2H 2V 3H 3V 4H 4V 5H...
                    %  Tx_array:
                    %               2  4     10 12      after : 1H 2H 3H 4H 1V 2V 3V 4V 5H...
                    %               1  3      9 11
                    %
                    
%                     for i_TxP = 1 : Tx_Mg*Tx_Ng  % index of panel at BS side
%                         for i_RxP = 1 : Rx_Mg*Rx_Ng  % index of panel at UE side
%                             
%                             H_temp((i_RxP-1)*2*Rs+1:(i_RxP-1)*2*Rs+Rs, (i_TxP-1)*2*Ms+1:(i_TxP-1)*2*Ms+Ms) = ...
%                                 h_n((i_RxP-1)*2*Rs+1:2:i_RxP*2*Rs, (i_TxP-1)*2*Ms+1:2:i_TxP*2*Ms);
%                             H_temp((i_RxP-1)*2*Rs+Rs+1:i_RxP*2*Rs, (i_TxP-1)*2*Ms+1:(i_TxP-1)*2*Ms+Ms) = ...
%                                 h_n((i_RxP-1)*2*Rs+2:2:i_RxP*2*Rs, (i_TxP-1)*2*Ms+1:2:i_TxP*2*Ms);
%                             H_temp((i_RxP-1)*2*Rs+1:(i_RxP-1)*2*Rs+Rs, (i_TxP-1)*2*Ms+Ms+1:i_TxP*2*Ms) = ...
%                                 h_n((i_RxP-1)*2*Rs+1:2:i_RxP*2*Rs, (i_TxP-1)*2*Ms+2:2:i_TxP*2*Ms);
%                             H_temp((i_RxP-1)*2*Rs+Rs+1:i_RxP*2*Rs, (i_TxP-1)*2*Ms+Ms+1:i_TxP*2*Ms) = ...
%                                 h_n((i_RxP-1)*2*Rs+2:2:i_RxP*2*Rs, (i_TxP-1)*2*Ms+2:2:i_TxP*2*Ms);
%                         end
%                     end
%                     H(:, :, iter_CP+nCarrierPerRB*(iter_RB-1), iter_user) = H_temp;
                    
                    H(:, :, iter_RB, iter_user) = h_n;
            end   
        end
        save([channel_path, num2str(iter_slot), '_', num2str(iter_TTI)], 'H');
        
        % Update progress
        file_count = file_count + 1;
        progress_percent = file_count / total_files * 100;
        fprintf('  TTI %d/%d 完成 - 文件: %s_%s.mat (总进度: %.1f%%)\n', ...
                iter_TTI, nTTiPerSlot, num2str(iter_slot), num2str(iter_TTI), progress_percent);
    end
    fprintf('时隙 %d 完成!\n', iter_slot);
end

fprintf('========================================\n');
fprintf('信道数据生成完成! 共生成 %d 个文件\n', file_count);
toc;