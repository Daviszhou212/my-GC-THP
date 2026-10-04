%channelData
channel_path = 'E:/linbaishuo/channel/';

AllUsers = 1:30;
% % pg_dB = randn(1,30) + [-112*ones(1,5), -108*ones(1,5), -109*ones(1,5), -106*ones(1,5), -115*ones(1,5), -105-10*rand(1,5)];
% % AllUsers = [1,2,3,9,10,12,13,15,16,17];
% % pg_dB = randn(1,10) + [-112*ones(1,3), -108*ones(1,2), -114*ones(1,3), -106*ones(1,2)];

% % factors = CalScalingFactor(channel_path, AllUsers, pg_dB);
% % save('factors_30.mat', 'factors');
% % load factors_30
% % warning("Using restored factors. Don't forget to check it.");

%%% System Parameters
% Package
pSize = 8 * 25000;    % The fixed package size in [bit]
pInterval = 0.01;     % Average time interval between packages' arrivals
% Power
Pt_By_N0_dB = 120;    % P_T/N_0
Pt_By_N0 = 10^(Pt_By_N0_dB / 10);
% System
all_len = 1;
miu1 = 1e-12;
miu2 = 1e-12;

PT = Pt_By_N0;

nUE = numel(AllUsers);             % # of users
nRB = 20;             % # of RB

nCarrierPerRB = 12;   % # of subcarriers of each RB
subcarrierSpacing = 120e3;
nSlotPerFrame = 4;    % # of slots in each frame
nTTiPerSlot = 160;    % # of TTIs in each slot
nSymbolsPerTTI = 14;  % # of symbols in each TTI
Fc = 28e9;            % Center frequency in [Hz]
LenOfTTI = 125e-6;    % Length of one TTI in [s]
% Cell
cellDiameter = 200;   % Diameter of one cell in [m] (BS's max coverage distance)
% Antenna at UE
Rx_M = 2;             % # of vertical elements on each panel
Rx_N = 2;             % # of horizontal elements on each panel
Rx_Mg = 1;            % # of nested panels in a column
Rx_Ng = 1;            % # of nested panels in a row
Rx_spacing = 0.5;     % Element spacing in [位]
Rx_poInd = 2;         % Polarization indicator = 2: H/V polarized elements
Rx_height = 1.5;      % The height of UE antenna in [meter]
UE_Speed = 3;         % Speed of UEs in [km/h]
Rs = Rx_M * Rx_N;     % # of antennas in each sub-array
nRx = 2 * Rs * Rx_Mg * Rx_Ng;  % # of total antennas at UE side
% Antenna at BS
Tx_M = 16;            % # of vertical elements on each panel
Tx_N = 16;            % # of horizontal elements on each panel
Tx_Mg = 2;            % # of nested panels in a column
Tx_Ng = 2;            % # of nested panels in a row
Tx_spacing = 0.5;     % Element spacing in [Lambda]
Tx_poInd = 3;         % Polarization indicator = 3: +/-45 degree polarized elements
Tx_height = 25;       % The height of BS antenna in [meter]
Tx_downtilt = atan((Tx_height - Rx_height) / cellDiameter);  % BS's downtilt angle in [rad]
Ms = Tx_M * Tx_N;     % # of antennas in each sub-array
nPanelsPair = Tx_Mg * Tx_Ng;  % # of antenna panels pair (each with 2 polalizations)
S = 2 * nPanelsPair;  % # of sub-arrays/RF-chains
nTx = S * Ms;         % # of total antennas at BS side
thresold = 0.2;       % Thresold for clustering

%%% Algorithm Parameters
Threshold = 15;       % Interference threshold (dB)
Tc = 320;             % Exponentially smoothed average data rate over the past
                      %  Tc TTIs to asymptotically approach PF objective