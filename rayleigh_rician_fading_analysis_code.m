% Define parameters
N = 2000; % Length of the signal (covers range from -1000 to 1000)
snr_values = [0, 10, 20, 30]; % Different SNR levels
Tm_values = [0, 0.001, 0.01, 0.1]; % Different time variation levels (max Doppler shift)
fd = 10; % Doppler frequency (example value)

% Generate a random input signal
input_signal = randi([0, 1], N, 1) * 2 - 1; % Binary signal with values -1 and 1
impulse_signal = [zeros(N/2, 1); 1; zeros(N/2-1, 1)]; % Impulse signal at the center

% Define channel types
channel_types = {'Ray_TI_FF', 'Ray_TI_FS', 'Ray_TV_FS', 'Ric_TV_FS', 'Awgn', 'Ray_TV_FF'};

% Initialize results storage
results = struct();

% Simulate and analyze each channel type
for i = 1:length(channel_types)
    ch_type = channel_types{i};
    results.(ch_type) = struct();
    
    for snr_idx = 1:length(snr_values)
        snr = snr_values(snr_idx);
        for Tm_idx = 1:length(Tm_values)
            Tm = Tm_values(Tm_idx);
            
            % Simulate channel
            ch_output = simulate_channel(input_signal, ch_type, snr, Tm, fd);
            
            % Calculate performance metrics
            noise = ch_output - input_signal;
            signal_power = mean(abs(ch_output).^2);
            noise_power = mean(abs(noise).^2);
            ber = sum(input_signal ~= sign(real(ch_output))) / N;
            
            % Store results
            field_name = ['snr_', num2str(snr), '_Tm_', strrep(num2str(Tm), '.', '_')];
            results.(ch_type).(field_name) = struct(...
                'signal_power', signal_power, ...
                'noise_power', noise_power, ...
                'ber', ber);
        end
    end
end

% Plot impulse responses
figure;
for i = 1:length(channel_types)
    ch_type = channel_types{i};
    ch_output = simulate_channel(impulse_signal, ch_type, snr_values(end), Tm_values(end), fd);
    
    subplot(3, 2, i);
    plot(-1000:999, abs(ch_output), 'LineWidth', 1);
    title(['Impulse Response: ', ch_type]);
    xlabel('Sample Index');
    ylabel('Amplitude');
end

% Plot BER results
figure;
for i = 1:length(channel_types)
    ch_type = channel_types{i};
    subplot(3, 2, i);
    ber_values = [];
    x_labels = {};
    for snr_idx = 1:length(snr_values)
        for Tm_idx = 1:length(Tm_values)
            Tm = Tm_values(Tm_idx);
            field_name = ['snr_', num2str(snr_values(snr_idx)), '_Tm_', strrep(num2str(Tm), '.', '_')];
            ber_values = [ber_values, results.(ch_type).(field_name).ber];
            x_labels{end+1} = sprintf('SNR=%d, Tm=%0.3f', snr_values(snr_idx), Tm);
        end
    end
    bar(ber_values);
    title(['BER for ', ch_type]);
    xlabel('Conditions (SNR, Tm)');
    ylabel('Bit Error Rate');
    xticks(1:length(ber_values));
    xticklabels(x_labels);
    xtickangle(45);
end

% Plot BER results for varying SNR (Tm constant)
figure;
Tm_constant = Tm_values(2); % Select a constant Tm value
for i = 1:length(channel_types)
    ch_type = channel_types{i};
    subplot(3, 2, i);
    ber_values = [];
    x_labels = {};
    for snr_idx = 1:length(snr_values)
        snr = snr_values(snr_idx);
        field_name = ['snr_', num2str(snr), '_Tm_', strrep(num2str(Tm_constant), '.', '_')];
        ber_values = [ber_values, results.(ch_type).(field_name).ber];
        x_labels{end+1} = sprintf('SNR=%d', snr);
    end
    bar(ber_values);
    title(['BER vs SNR for ', ch_type, ' (Tm=', num2str(Tm_constant), ')']);
    xlabel('SNR');
    ylabel('Bit Error Rate');
    xticks(1:length(ber_values));
    xticklabels(x_labels);
    xtickangle(45);
end

% Plot BER results for varying Tm (SNR constant)
figure;
snr_constant = snr_values(2); % Select a constant SNR value
for i = 1:length(channel_types)
    ch_type = channel_types{i};
    subplot(3, 2, i);
    ber_values = [];
    x_labels = {};
    for Tm_idx = 1:length(Tm_values)
        Tm = Tm_values(Tm_idx);
        field_name = ['snr_', num2str(snr_constant), '_Tm_', strrep(num2str(Tm), '.', '_')];
        ber_values = [ber_values, results.(ch_type).(field_name).ber];
        x_labels{end+1} = sprintf('Tm=%0.3f', Tm);
    end
    bar(ber_values);
    title(['BER vs Tm for ', ch_type, ' (SNR=', num2str(snr_constant), ')']);
    xlabel('Tm');
    ylabel('Bit Error Rate');
    xticks(1:length(ber_values));
    xticklabels(x_labels);
    xtickangle(45);
end

% Function to simulate channel (from previous implementation)
function ch_output = simulate_channel(ch_input, ch_type, snr, Tm, fd)
    % Simulate the specified channel and add noise based on SNR
    N = length(ch_input);
    t = (0:N-1)';
    
    switch ch_type
        case 'Ray_TI_FF'
            % Time-invariant Rayleigh channel with flat fading
            h = (randn(N, 1) + 1i*randn(N, 1))/sqrt(2);
            ch_output = ch_input .* h;
            
        case 'Ray_TI_FS'
            % Time-invariant Rayleigh channel with frequency selective fading
            num_paths = 3; % Number of multipath components
            delays = [0, 2, 5]; % Path delays
            h = (randn(num_paths, 1) + 1i*randn(num_paths, 1))/sqrt(2);
            ch_output = zeros(N, 1);
            for k = 1:num_paths
                ch_output = ch_output + [zeros(delays(k), 1); ch_input(1:end-delays(k))] .* h(k);
            end
            
        case 'Ray_TV_FS'
            % Time-varying Rayleigh channel with frequency selective fading
            num_paths = 3; % Number of multipath components
            delays = [0, 2, 5]; % Path delays
            h = zeros(N, num_paths);
            for k = 1:num_paths
                h(:, k) = (randn(N, 1) + 1i*randn(N, 1))/sqrt(2) .* cos(2*pi*fd*t + randn);
            end
            ch_output = zeros(N, 1);
            for k = 1:num_paths
                ch_output = ch_output + [zeros(delays(k), 1); ch_input(1:end-delays(k))] .* h(:, k);
            end
            
        case 'Ric_TV_FS'
            % Time-varying Rician channel with frequency selective fading
            K = 3; % Rician K-factor
            num_paths = 3; % Number of multipath components
            delays = [0, 2, 5]; % Path delays
            h_los = sqrt(K/(K+1)) * exp(1i*2*pi*fd*t);
            h_nlos = sqrt(1/(K+1)) * (randn(N, num_paths) + 1i*randn(N, num_paths))/sqrt(2);
            h = h_los + h_nlos;
            ch_output = zeros(N, 1);
            for k = 1:num_paths
                ch_output = ch_output + [zeros(delays(k), 1); ch_input(1:end-delays(k))] .* h(:, k);
            end
            
        case 'Awgn'
            % Additive White Gaussian Noise channel
            ch_output = awgn(ch_input, snr, 'measured');
            
        case 'Ray_TV_FF'
            % Time-varying Rayleigh channel with flat fading
            h = (randn(N, 1) + 1i*randn(N, 1))/sqrt(2) .* cos(2*pi*fd*t + randn);
            ch_output = ch_input .* h;
            
        otherwise
            error('Unknown channel type');
    end
    
    % Add AWGN to the output signal
    ch_output = awgn(ch_output, snr, 'measured');
end
