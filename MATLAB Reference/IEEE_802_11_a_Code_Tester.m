% Jai Mangal
% PhD22104
% IEEE 802.11a Based Wireless Communication System
% IIIT Delhi
% Date: 20 January 2025

%% Clear Workspace

close all; % Close all figure windows
clear; % Clear workspace variables
clc; % Clear command window

%% Signal Parameters

N_FFT = 64; % FFT Size
fc_hz = 5e9; % Carrier frequency (5 GHz)
fs_hz = 20e6; % Sampling frequency (20 MHz)
ts_sec = 1/fs_hz; % Time period (50 ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Transmitter Designing %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Short Preamble

S_k = sqrt(13/6)*[0,0,1+1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,0,0,0,0,-1-1i,0,0,0,-1-1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0,0,1+1i,0,0]; % Short preamble sequence [1x53] 
virtual_subcarrier = zeros(1,N_FFT-length(S_k)); % Virtual subcarriers [1x11]
Short_preamble_slot_Frequency = [virtual_subcarrier(1:6),S_k,virtual_subcarrier(7:11)]; % Frequency-domain preamble [1x64]
Short_preamble_slot_Time = ifft(ifftshift(Short_preamble_slot_Frequency)); % Time-domain preamble [1x64]
Short_preamble = repmat(Short_preamble_slot_Time(1:16),1,10); % Repeat short preamble [1x160]

%% Long Preamble

L_k = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1]; % Long preamble sequence [1x53] 
virtual_subcarrier = zeros(1,N_FFT-length(L_k)); % Virtual subcarriers [1x11]
Long_preamble_slot_Frequency = [virtual_subcarrier(1:6),L_k,virtual_subcarrier(7:11)]; % Frequency-domain preamble [1x64]
Long_preamble_slot_Time = ifft(ifftshift(Long_preamble_slot_Frequency)); % Time-domain preamble [1x64]
Long_preamble = [Long_preamble_slot_Time(33:64),Long_preamble_slot_Time,Long_preamble_slot_Time]; % Repeat long preamble [1x160]

%% Payload

M = 4; % QPSK modulation
bits_per_symb = log2(M); % Bits per symbol
rng('default'); % Set random seed for reproducibility
n_bits = 96; % Number of bits in one frame

% data_payload_1 = randi([0, 1], 1, n_bits); % Generate random payload 1 [1x96]
% data_payload_2 = randi([0, 1], 1, n_bits); % Generate random payload 2 [1x96]

data_payload_1 = [0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0];
data_payload_2 = [0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0];

data_payload_1_reshape = reshape(data_payload_1, 2, 0.5*n_bits)'; % Reshape data to make symbols [1x48]
data_payload_2_reshape = reshape(data_payload_2, 2, 0.5*n_bits)'; % Reshape data to make symbols [1x48]

data_payload_1_mod  = zeros(1,length(data_payload_1_reshape)); % Predefine array for QPSK Modulation [1x48]
data_payload_2_mod  = zeros(1,length(data_payload_2_reshape)); % Predefine array for QPSK Modulation [1x48]

% QPSK modulation for payload 1
for i = 1:length(data_payload_1_reshape)
    if data_payload_1_reshape(i,1) == 0 && data_payload_1_reshape(i,2) == 0
        data_payload_1_mod(i) = (1 + 1i)*(1/sqrt(2)); % Symbol [00]
    elseif data_payload_1_reshape(i,1) == 0 && data_payload_1_reshape(i,2) == 1
        data_payload_1_mod(i) = (-1 + 1i)*(1/sqrt(2)); % Symbol [01]
    elseif data_payload_1_reshape(i,1) == 1 && data_payload_1_reshape(i,2) == 0
        data_payload_1_mod(i) = (-1 - 1i)*(1/sqrt(2)); % Symbol [10]
    else
        data_payload_1_mod(i) =(1 - 1i)*(1/sqrt(2)); % Symbol [11]
    end
end

% QPSK modulation for payload 2
for i = 1:length(data_payload_2_reshape)
    if data_payload_2_reshape(i,1) == 0 && data_payload_2_reshape(i,2) == 0
        data_payload_2_mod(i) = (1 + 1i)*(1/sqrt(2)); % Symbol [00]
    elseif data_payload_2_reshape(i,1) == 0 && data_payload_2_reshape(i,2) == 1
        data_payload_2_mod(i) = (-1 + 1i)*(1/sqrt(2)); % Symbol [01]
    elseif data_payload_2_reshape(i,1) == 1 && data_payload_2_reshape(i,2) == 0
        data_payload_2_mod(i) = (-1 - 1i)*(1/sqrt(2)); % Symbol [10]
    else
        data_payload_2_mod(i) = (1 - 1i)*(1/sqrt(2)); % Symbol [11]
    end
end

pilot = [1,1,1,-1]; % Pilot symbols [1x4]

% Frame 1 construction with virtual subcarriers and pilots
data_frame_1 = [virtual_subcarrier(1:6),data_payload_1_mod(1:5),pilot(1),data_payload_1_mod(6:18),pilot(2),data_payload_1_mod(19:24),0,data_payload_1_mod(25:30),pilot(3),data_payload_1_mod(31:43),pilot(4),data_payload_1_mod(44:48),virtual_subcarrier(7:11)]; % [1x64]
% Frame 2 construction with virtual subcarriers and pilots
data_frame_2 = [virtual_subcarrier(1:6),data_payload_2_mod(1:5),pilot(1),data_payload_2_mod(6:18),pilot(2),data_payload_2_mod(19:24),0,data_payload_2_mod(25:30),pilot(3),data_payload_2_mod(31:43),pilot(4),data_payload_2_mod(44:48),virtual_subcarrier(7:11)]; % [1x64]

data_frame_1_ifft = ifft(ifftshift(data_frame_1)); % Inverse FFT for frame 1 [1x64]
data_frame_2_ifft = ifft(ifftshift(data_frame_2)); % Inverse FFT for frame 2 [1x64]

data_1_TX_payload = [data_frame_1_ifft(49:64),data_frame_1_ifft]; % Add cyclic prefix to frame 1 [1x80]
data_2_TX_payload = [data_frame_2_ifft(49:64),data_frame_2_ifft]; % Add cyclic prefix to frame 2 [1x80]

%% Frame Combination

Frame_Tx = [Short_preamble,Long_preamble,data_1_TX_payload,data_2_TX_payload]; % Combine short preamble, long preamble, and payloads into the transmission frame [160+160+80+80 = 1x480]

%% Oversampling

oversampling_rate_Tx = 2; % Oversampling rate
Frame_Tx_Oversamp = zeros(1,oversampling_rate_Tx*length(Frame_Tx));
Frame_Tx_Oversamp(1:oversampling_rate_Tx:end) = Frame_Tx; % Upsample the frame for transmission [1x960]

%% Root Raised Cosine Filter (Tx)

rolloff_factor_Tx = 0.5; % Rolloff factor for the RRC filter
length_rrc_Tx = 10; % Length of the RRC filter
rrc_filter_tx = rcosdesign(rolloff_factor_Tx,length_rrc_Tx,oversampling_rate_Tx,'sqrt'); % Coefficients of Root raised cosine filter [1x21]

%% Convolve the Tx frame with the RRC filter

len_inp_sig = length(Frame_Tx_Oversamp); % Length of input signal
len_coeff = length(rrc_filter_tx); % Length of coefficients of RRC filter
len_out_sig = len_inp_sig + len_coeff - 1; % Length of output signal [1x980]
inp_sig_padded = [zeros(1, len_coeff - 1), Frame_Tx_Oversamp, zeros(1, len_coeff - 1)]; % Padding x

% Initialize output array
Tx_signal = zeros(1, len_out_sig);  

% Perform convolution
for n = 1:len_out_sig
    for k = 1:len_coeff
        Tx_signal(n) = Tx_signal(n) + rrc_filter_tx(k) * inp_sig_padded(n + k - 1);  % Multiply and accumulate
    end
end

%% TX Signal

Tx_signal = repmat(Tx_signal,1,10); % Repeat the signal for transmission [1x9800]

%% %% Over The Air Transmission (Channel)

snr_dB_values = [100, 100]; % Example: SNR values from 0 dB to 20 dB in steps of 2
results = struct(); % Initialize a structure to store results for each SNR

for j = 1:length(snr_dB_values)
    snr_dB = snr_dB_values(j); % Current SNR value
    Tx_signal_power = mean(abs(Tx_signal).^2);  % Average power of the signal
    snr_linear = 10^(snr_dB / 10); % Convert SNR from dB to linear scale
    noise_power = Tx_signal_power / snr_linear; % Noise power for the desired SNR
    noise = sqrt(noise_power) * randn(1, length(Tx_signal));  % Gaussian noise with zero mean
    Tx_OTA_signal = Tx_signal; % Add noise to the signal

%% Start Reciver & Capture Packets

rng('default'); % Set random seed for reproducibility
np_packets_capture = 3000; % Number of packets to capture
rx_start = 1; % Random starting index for packet capture
Rx_signal = Tx_OTA_signal(rx_start:np_packets_capture+rx_start-1); % Capture the received signal [1x3000]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Receiver Designing %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Root Raised Cosine filter (Rx)

rolloff_factor_rx = 0.5; % Rolloff factor for the RRC filter
length_rrc_rx = 10; % Length of the RRC filter
oversampling_rate_rx = 2; % Oversampling rate for reception
rrc_filter_rx = rcosdesign(rolloff_factor_rx,length_rrc_rx,oversampling_rate_rx,'sqrt'); % Coefficients of Root raised cosine filter [1x21]

%% Convolve the Rx frame with the RRC filter

len_inp_sig = length(Rx_signal); % Length of input signal
len_coeff = length(rrc_filter_rx); % Length of coefficients of RRC filter
len_out_sig = len_inp_sig + len_coeff - 1; % Length of output signal [1x3020]
inp_sig_padded = [zeros(1, len_coeff - 1), Rx_signal, zeros(1, len_coeff - 1)]; % Padding x

% Initialize output array
Rx_filter_signal = zeros(1, len_out_sig);

% Convolve the received signal with the RRC filter [1x3020]
for n = 1:len_out_sig
    for k = 1:len_coeff
        Rx_filter_signal(n) = Rx_filter_signal(n) + rrc_filter_rx(k) * inp_sig_padded(n + k - 1);  % Multiply and accumulate
    end
end

%% Packet Detection

delay_param = 16; % Delay parameter for correlation
window_length = 32; % Length of the correlation window
corr_arr = zeros(1,length(Rx_signal)-delay_param+1-window_length); % Initialize correlation array
peak_arr = zeros(1,length(Rx_signal)-delay_param+1-window_length); % Initialize power array
corr_temp_arr = zeros(1,window_length); % Temporary array for correlation
peak_temp_arr = zeros(1,window_length); % Temporary array for power

% Perform correlation and power calculations
for n=1:length(Rx_signal)-delay_param+1-window_length
    for k=1:window_length
        corr_temp_arr(k) = Rx_signal(n+k-1)*Rx_signal(n+k-1+delay_param);  % Cross-correlation
        peak_temp_arr(k) = abs(Rx_signal(n+k-1+delay_param))^2;  % Power of signal
    end
    corr_arr(n) = sum(corr_temp_arr); % Sum of correlation
    peak_arr(n) = sum(peak_temp_arr); % Sum of power
end

corr_out = (abs(corr_arr).^2)./(peak_arr.^2); % Normalize correlation to detect peaks

%% Packet Selection

packet_threshold = 0.75; % Correlation threshold for detecting packets
packet_idx_arr = zeros(1, length(corr_out)); % Pre-allocate array for indices
idx_count = 0; % Counter for the number of valid indices

% Loop through corr_out and check for values greater than packet_threshold
for i = 1:length(corr_out)
    if corr_out(i) > packet_threshold
        idx_count = idx_count + 1; % Increment the counter
        packet_idx_arr(idx_count) = i; % Store the index in the array
    end
end
packet_idx_arr = packet_idx_arr(1:idx_count); % Trim unused elements from the array

% Create shifted arrays for comparison
temp_1 = [packet_idx_arr, 0]; % Shifted array for comparison
temp_2 = [0, packet_idx_arr]; % Shifted array for comparison
temp_3 = temp_1 - temp_2; % Difference between shifted arrays

packet_front = zeros(1, length(temp_3)); % Pre-allocate array for packet front indices
packet_front_count = 0; % Counter for the number of packet start indices

% Loop through temp_3 and check for values greater than 300
for i = 1:length(temp_3)
    if temp_3(i) > 300
        packet_front_count = packet_front_count + 1; % Increment the counter
        packet_front(packet_front_count) = i; % Store the index in the array
    end
end
packet_front = packet_front(1:packet_front_count); % Trim unused elements from the array

% Extract start indices of packets
packet_front_idx = packet_idx_arr(packet_front);

% Define threshold length
threshold_length = 230; % Minimum length of the packet

% packet_idx = 0;
% Select the packet start index
for x=1:length(packet_front_idx)-1
    if corr_out(packet_front_idx(x)+threshold_length) > packet_threshold
        packet_idx = packet_front_idx(x)+length_rrc_rx+1; % Adjust for RRC filter delay
        break; % Exit loop once a valid packet is found
    end
end

%% Downsampling

% Downsample the received signal using the specified oversampling rate (rx_frame)
rx_frame = Rx_filter_signal(packet_idx:oversampling_rate_rx:oversampling_rate_rx*length(Frame_Tx)+packet_idx-1); % [1x480] Frame length

%% Coarse CFO Estimation

% Define the length of the Short preamble slot for coarse CFO estimation
Short_preamble_slot_length = 16;

% Calculate the complex conjugate product between the two relevant sections of the received frame
prod_consq_frame_coarse = rx_frame(Short_preamble_slot_length*5+1:Short_preamble_slot_length*6)*rx_frame(Short_preamble_slot_length*6+1:Short_preamble_slot_length*7)'; % [1x16]*[16x1]

% Estimate the coarse frequency offset
freq_coarse_est = (-1/(2*pi*Short_preamble_slot_length*ts_sec))*atan2(imag(prod_consq_frame_coarse), real(prod_consq_frame_coarse));

% Apply the coarse frequency offset to the received frame
rx_frame_after_coarse = rx_frame.*exp(-1i*2*pi*freq_coarse_est*ts_sec*(0:480-1)); % [1x480]

%% Fine CFO Estimation

% Calculate the complex conjugate product between the relevant sections of the received frame for fine CFO estimation
prod_consq_frame_fine = rx_frame_after_coarse(Short_preamble_slot_length*12+1:Short_preamble_slot_length*16)*rx_frame_after_coarse(Short_preamble_slot_length*16+1:Short_preamble_slot_length*20)'; % [1x64]*[64x1]=[1x1]

% Estimate the fine frequency offset
freq_fine_est = (-1/(2*pi*64*ts_sec))*atan2(imag(prod_consq_frame_fine), real(prod_consq_frame_fine));

% Apply the fine frequency offset to the received frame
rx_frame_after_fine = rx_frame_after_coarse.*exp(-1i*2*pi*freq_fine_est*ts_sec*(0:480-1)); % [1x160]

%% Channel Estimation

% Extract sections of the fine frequency-offset received frame for channel estimation
Long_preamble_1 = rx_frame_after_fine(Short_preamble_slot_length*12+1:Short_preamble_slot_length*16); % [1x64]
Long_preamble_2 = rx_frame_after_fine(Short_preamble_slot_length*16+1:Short_preamble_slot_length*20); % [1x64]

% Perform FFT to transform the sections into the frequency domain
Long_preamble_1_After_FFT = fftshift(fft(Long_preamble_1)); % [1x64]
Long_preamble_2_After_FFT = fftshift(fft(Long_preamble_2)); % [1x64]

% Estimate the channel by averaging the FFT of both Long preamble sections and taking the conjugate of the Long preamble
H_est = 0.5*(Long_preamble_1_After_FFT+Long_preamble_2_After_FFT).*conj(Long_preamble_slot_Frequency); % [1x64]

% Perform IFFT to get the   channel estimate in the time domain
H_est_time = ifft(ifftshift(H_est)); % [1x64]

%% One tap Equalizer

% Equalize the first payload frame
RX_Payload_1_time = rx_frame_after_fine(320+1:400); % [1x80]
RX_Payload_1_no_CP = RX_Payload_1_time(17:end); % [1x64]
RX_Payload_1_Frequency = fftshift(fft(RX_Payload_1_no_CP)); % [1x64]
RX_Payload_1_Frequency_Equalizer = RX_Payload_1_Frequency./H_est; % [1x64]

% Equalize the second payload frame
RX_Payload_2_time = rx_frame_after_fine(400+1:480); % [1x80]
RX_Payload_2_no_CP = RX_Payload_2_time(17:end); % [1x64]
RX_Payload_2_Frequency = fftshift(fft(RX_Payload_2_no_CP)); % [1x64]
RX_Payload_2_Frequency_Equalizer = RX_Payload_2_Frequency./H_est; % [1x64]

%% De-Mapping

% Remove pilot symbols and extract data from the frequency domain for the first payload
RX_Payload_1_no_Equalizer = [RX_Payload_1_Frequency(7:11),RX_Payload_1_Frequency(13:25),RX_Payload_1_Frequency(27:32),RX_Payload_1_Frequency(34:39),RX_Payload_1_Frequency(41:53),RX_Payload_1_Frequency(55:59)]; % [1x48]

% Remove pilot symbols and extract data from the equalized frequency domain for the first payload
RX_Payload_1_no_pilot = [RX_Payload_1_Frequency_Equalizer(7:11),RX_Payload_1_Frequency_Equalizer(13:25),RX_Payload_1_Frequency_Equalizer(27:32),RX_Payload_1_Frequency_Equalizer(34:39),RX_Payload_1_Frequency_Equalizer(41:53),RX_Payload_1_Frequency_Equalizer(55:59)]; % [1x48]

% Remove pilot symbols and extract data from the frequency domain for the second payload
RX_Payload_2_no_Equalizer = [RX_Payload_2_Frequency(7:11),RX_Payload_2_Frequency(13:25),RX_Payload_2_Frequency(27:32),RX_Payload_2_Frequency(34:39),RX_Payload_2_Frequency(41:53),RX_Payload_2_Frequency(55:59)]; % [1x48]

% Remove pilot symbols and extract data from the equalized frequency domain for the second payload
RX_Payload_2_no_pilot = [RX_Payload_2_Frequency_Equalizer(7:11),RX_Payload_2_Frequency_Equalizer(13:25),RX_Payload_2_Frequency_Equalizer(27:32),RX_Payload_2_Frequency_Equalizer(34:39),RX_Payload_2_Frequency_Equalizer(41:53),RX_Payload_2_Frequency_Equalizer(55:59)]; % [1x48]

%% AGC For Rx Data Payload 1

% Initialize Mapped Points Array
RX_Payload_1_Final = zeros(size(RX_Payload_1_no_pilot)); %[1x48]

for idx = 1:length(RX_Payload_1_no_pilot)
    % Get the current point
    RX_Payload_1_no_pilot_current = RX_Payload_1_no_pilot(idx);
    
    % Initialize Mapped Real and Imaginary Parts
    RX_Payload_1_mapped_real = 0; % Real part of the mapped point
    RX_Payload_1_mapped_imag = 0; % Imaginary part of the mapped point
    
    % Decision Logic for Real Part
    if real(RX_Payload_1_no_pilot_current) > 0
        RX_Payload_1_mapped_real = 1 / sqrt(2); % Map to +1/sqrt(2) for positive real part
    elseif real(RX_Payload_1_no_pilot_current) < 0
        RX_Payload_1_mapped_real = -1 / sqrt(2); % Map to -1/sqrt(2) for negative real part
    end
    
    % Decision Logic for Imaginary Part
    if imag(RX_Payload_1_no_pilot_current) > 0
        RX_Payload_1_mapped_imag = 1 / sqrt(2); % Map to +1/sqrt(2) for positive imaginary part
    elseif imag(RX_Payload_1_no_pilot_current) < 0
        RX_Payload_1_mapped_imag = -1 / sqrt(2); % Map to -1/sqrt(2) for negative imaginary part
    end
    
    % Combine Mapped Real and Imaginary Parts
    RX_Payload_1_Final(idx) = RX_Payload_1_mapped_real + 1j * RX_Payload_1_mapped_imag;
end


%% AGC For Rx Data Payload 2

% Initialize Mapped Points Array
RX_Payload_2_Final = zeros(size(RX_Payload_2_no_pilot)); %[1x48]

for idx = 1:length(RX_Payload_2_no_pilot)
    % Get the current point
    RX_Payload_2_no_pilot_current = RX_Payload_2_no_pilot(idx);
    
    % Initialize Mapped Real and Imaginary Parts
    RX_Payload_2_mapped_real = 0; % Real part of the mapped point
    RX_Payload_2_mapped_imag = 0; % Imaginary part of the mapped point
    
    % Decision Logic for Real Part
    if real(RX_Payload_2_no_pilot_current) > 0
        RX_Payload_2_mapped_real = 1 / sqrt(2); % Map to +1/sqrt(2) for positive real part
    elseif real(RX_Payload_2_no_pilot_current) < 0
        RX_Payload_2_mapped_real = -1 / sqrt(2); % Map to -1/sqrt(2) for negative real part
    end
    
    % Decision Logic for Imaginary Part
    if imag(RX_Payload_2_no_pilot_current) > 0
        RX_Payload_2_mapped_imag = 1 / sqrt(2); % Map to +1/sqrt(2) for positive imaginary part
    elseif imag(RX_Payload_2_no_pilot_current) < 0
        RX_Payload_2_mapped_imag = -1 / sqrt(2); % Map to -1/sqrt(2) for negative imaginary part
    end
    
    % Combine Mapped Real and Imaginary Parts
    RX_Payload_2_Final(idx) = RX_Payload_2_mapped_real + 1j * RX_Payload_2_mapped_imag;
end

%% QPSK Demdoulation For Rx Data Payload 1

% Preallocate the demodulated bits array for better performance
RX_Payload_1_demod = zeros(1, 2 * length(RX_Payload_1_Final)); %[1x96]

% Loop through each symbol
for i = 1:length(RX_Payload_1_Final)
    % Extract real and imaginary parts of the current symbol
    RX_Payload_1_demod_real = real(RX_Payload_1_Final(i));
    RX_Payload_1_demod_imag = imag(RX_Payload_1_Final(i));
    
    % Determine bits based on QPSK mapping
    if RX_Payload_1_demod_real > 0 && RX_Payload_1_demod_imag > 0
        % Mapping for 00 -> 0.707 + i*0.707
        RX_Payload_1_demod(2*i-1:2*i) = [0, 0];
    elseif RX_Payload_1_demod_real < 0 && RX_Payload_1_demod_imag > 0
        % Mapping for 01 -> -0.707 + i*0.707
        RX_Payload_1_demod(2*i-1:2*i) = [0, 1];
    elseif RX_Payload_1_demod_real < 0 && RX_Payload_1_demod_imag < 0
        % Mapping for 10 -> -0.707 - i*0.707
        RX_Payload_1_demod(2*i-1:2*i) = [1, 0];
    elseif RX_Payload_1_demod_real > 0 && RX_Payload_1_demod_imag < 0
        % Mapping for 11 -> 0.707 - i*0.707
        RX_Payload_1_demod(2*i-1:2*i) = [1, 1];
    end
end

%% QPSK Demdoulation For Rx Data Payload 2

% Preallocate the demodulated bits array for better performance
RX_Payload_2_demod = zeros(1, 2 * length(RX_Payload_2_Final)); %[1x96]

% Loop through each symbol
for i = 1:length(RX_Payload_2_Final)
    % Extract real and imaginary parts of the current symbol
    RX_Payload_2_demod_real = real(RX_Payload_2_Final(i));
    RX_Payload_2_demod_imag = imag(RX_Payload_2_Final(i));
    
    % Determine bits based on QPSK mapping
    if RX_Payload_2_demod_real > 0 && RX_Payload_2_demod_imag > 0
        % Mapping for 00 -> 0.707 + i*0.707
        RX_Payload_2_demod(2*i-1:2*i) = [0, 0];
    elseif RX_Payload_2_demod_real < 0 && RX_Payload_2_demod_imag > 0
        % Mapping for 01 -> -0.707 + i*0.707
        RX_Payload_2_demod(2*i-1:2*i) = [0, 1];
    elseif RX_Payload_2_demod_real < 0 && RX_Payload_2_demod_imag < 0
        % Mapping for 10 -> -0.707 - i*0.707
        RX_Payload_2_demod(2*i-1:2*i) = [1, 0];
    elseif RX_Payload_2_demod_real > 0 && RX_Payload_2_demod_imag < 0
        % Mapping for 11 -> 0.707 - i*0.707
        RX_Payload_2_demod(2*i-1:2*i) = [1, 1];
    end
end

%% EVM Calculation (Before AGC)

% Calculate error vector (difference between transmitted and received symbols)
error_vector = [RX_Payload_1_no_pilot, RX_Payload_2_no_pilot] - [data_payload_1_mod, data_payload_2_mod];

% Calculate EVM as the RMS value of the error magnitude normalized by RMS value of transmitted symbols
evm = sqrt(mean(abs(error_vector).^2)) / sqrt(mean(abs([data_payload_1_mod, data_payload_2_mod]).^2));

% Convert EVM to dB
evm_dB = 20 * log10(evm);

%% EVM Calculation (After AGC)

% Calculate error vector (difference between transmitted and received symbols)
error_vector_AGC = [RX_Payload_1_Final, RX_Payload_2_Final] - [data_payload_1_mod, data_payload_2_mod];

% Calculate EVM as the RMS value of the error magnitude normalized by RMS value of transmitted symbols
evm_AGC = sqrt(mean(abs(error_vector_AGC).^2)) / sqrt(mean(abs([data_payload_1_mod, data_payload_2_mod]).^2));

% Convert EVM to dB
evm_AGC_dB = 20 * log10(evm_AGC);

%% BER Calculation

% Calculate the number of errors in the received payloads compared to the original data
Error_bits = sum([abs(sign(data_payload_1-RX_Payload_1_demod)),abs(sign(data_payload_2-RX_Payload_2_demod))]);

% Calculate the Bit Error Rate (BER)
BER = Error_bits/(length(data_payload_1)+length(data_payload_2));

results(j).evm = evm_dB; % Store EVM values for each SNR
results(j).evm_AGC = evm_AGC_dB; % Store EVM values for each SNR
results(j).ber = BER; % Store BER values for each SNR

end

%% Results

evm_dB_values = [results.evm];
evm_AGC_dB_values = [results.evm_AGC];
evm_AGC_dB_values(evm_AGC_dB_values == -inf) = -40;
ber_values = [results.ber];
ber_values(ber_values == 0) = 1e-6;

%% Plots & Figures

% EVM v/s SNR
figure;
plot(snr_dB_values, evm_dB_values, 'o-','LineWidth',2); hold on
plot(snr_dB_values, evm_AGC_dB_values, 'o-','LineWidth',2);
xlabel('SNR (dB)');
ylabel('EVM (dB)');
title('EVM v/s SNR');
legend('Before AGC','After AGC');
grid on;

% BER v/s SNR
figure;
semilogy(snr_dB_values, ber_values, 'o-','LineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER v/s SNR');
grid on;