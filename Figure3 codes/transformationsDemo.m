function RV = transformationsDemo(outPath) 



    % Generate original signal components and total
    fs = 1000;
    total_time = -0.5:1/fs:2.5;
    f_theta = 8;

    % Generate individual harmonics and total signal
    individual_harmonics = zeros(6, length(total_time));
    harmonic_weights = [1 0.5 0.3 0.2 0.15 0.1];
    signal = zeros(size(total_time));

    for i = 1:length(harmonic_weights)
        individual_harmonics(i,:) = harmonic_weights(i) * sin(2*pi*f_theta*i*total_time);
        signal = signal + individual_harmonics(i,:);
    end

    % Compute STFT with full frequency range for accurate reconstruction
    window_length = round(fs/2);
    hop_size = round(window_length/4);
    window = hann(window_length);
    [S, F, T] = spectrogram(signal, window, window_length-hop_size, [], fs);
    T = T + total_time(1);
    reconstructed_stft = istft(S, window_length, hop_size, length(window), fs);

    % Compute wavelet transforms
    freq_bins = 1:80;
    time_points = total_time;
    
    % classic CWT:
    [W_cwt, f_cwt] = cwt(signal, 'amor', fs, 'FrequencyLimits', freq_bins([1,end]), 'VoicesPerOctave', 24);   % amor: Morlet wavelet
    reconstructed_cwt = icwt(W_cwt, f_cwt, freq_bins([1,end]), 'amor', fs);

    % Long wavelet (50 cycles)
    W_long = zeros(length(freq_bins), length(time_points));
    for f_idx = 1:length(freq_bins)
        f = freq_bins(f_idx);
        s_long = (50/(2*pi*f))*fs;
        t_long = -4*s_long/fs:1/fs:4*s_long/fs;
        psi_long = (pi*s_long)^(-0.25) * exp(2*pi*1i*f*t_long) .* exp(-t_long.^2/(2*(s_long/fs)^2));
        W_long(f_idx,:) = conv(signal, psi_long, 'same') / sqrt(s_long);
    end

    % Short wavelet (5 cycles)
    W_short = zeros(length(freq_bins), length(time_points));
    for f_idx = 1:length(freq_bins)
        f = freq_bins(f_idx);
        s_short = (5/(2*pi*f))*fs;
        t_short = -4*s_short/fs:1/fs:4*s_short/fs;
        psi_short = (pi*s_short)^(-0.25) * exp(2*pi*1i*f*t_short) .* exp(-t_short.^2/(2*(s_short/fs)^2));
        W_short(f_idx,:) = conv(signal, psi_short, 'same') / sqrt(s_short);
    end

    % Z-score normalize short wavelet
    W_short_z = zscore(abs(W_short), 0, 2);
    W_short_z = W_short_z .* exp(1i * angle(W_short));

    % Compute reconstructions
    C_delta = 0.776;
    reconstructed_long = zeros(size(time_points));
    reconstructed_short = zeros(size(time_points));
    reconstructed_z = zeros(size(time_points));

    for f_idx = 1:length(freq_bins)
        scale_long = (50/(2*pi*freq_bins(f_idx)))*fs;
        scale_short = (5/(2*pi*freq_bins(f_idx)))*fs;

        reconstructed_long = reconstructed_long + real(W_long(f_idx,:)) / sqrt(scale_long);
        reconstructed_short = reconstructed_short + real(W_short(f_idx,:)) / sqrt(scale_short);
        reconstructed_z = reconstructed_z + real(W_short_z(f_idx,:)) / sqrt(scale_short);
    end

    
    
    
    
    % Create figure
    % figure('Position', [100 100 1200 1200]);
    figArr = class_figArr;  % a customized class defines arrangement of subplots.
    figArr.X = [1 1];  figArr.dx = [0.2, 0.2, 0.2];   
    figArr.Y = [ones(1,7)];     figArr.dy = [0.5, 0.2*ones(1,4), 0.6 0.7 1.3]; 
    figArr.GM = [1, 1; ... 
                 2, 3; ...
                 4, 5; ...
                 6, 7; ...
                 8, 9; ...
                 10, 11; ...
                 12, 13];
    figArr.width = 1.0;
    figArr.AsRatio = 13/10;
    figArr.screenOffset = [0.5, 0.8];
    figArr.hideInnerAxis = "off";
    [CF, CA, CAnno] = create(figArr);
    iCA = 0;
    
    % Modified Panel A code
    iCA = iCA +1;
    axes(CA(iCA));    % subplot(6,1,1);
    harmonic_colors = winter(6);
    harmonic_colors = harmonic_colors(end:-1:1,:);
    for i = 1:size(individual_harmonics,1)
        plot(total_time, individual_harmonics(i,:), 'Color', [harmonic_colors(i,:) 0.7]);
    end
    plot(total_time, signal, 'k', 'LineWidth', 1.5);
    ylabel('Amplitude', 'Color', 'k');   xlabel('Time');
    CA(iCA).YAxisLocation = 'left';
    ylim([-2 2]);   xlim([0, 1]);
    % Improved legend formatting
    CL = legend({'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'Sum'}, 'Orientation', 'horizontal');
    CL.Position(1) = sum(CA(1).Position([1,3])) - CL.Position(3);
    CL.Position(2) = sum(CA(1).Position([2,4]));
    CL.Box = 'off';
    CA(1).XAxisLocation = 'bottom';
    CT = title('Original Signal: Theta with Harmonics', 'Unit', 'normalized');
    CT.Position(2) = CT.Position(2) + CL.Position(4)/CA(1).Position(4);
    grid on;
    xlim([0, 1]);
    set(gca, 'XTick', 0:0.25:2);

    
    % Compute FFT and plot
    fft_signal = fft(signal);
    frequencies = linspace(0, fs, length(fft_signal));
    amplitudes = abs(fft_signal);
    phases = angle(fft_signal);

    % Only show up to Nyquist frequency
    nyquist = length(frequencies) / 2;
    frequencies = frequencies(1:nyquist);
    amplitudes = amplitudes(1:nyquist);
    phases = phases(1:nyquist);

    % Perform reconstruction using inverse FFT
    reconstructed_fft = real(ifft(fft_signal));

    % Panel B: FFT Magnitude and Phase
    iCA = iCA +1;
    axes(CA(iCA)); % subplot(6, 2, 3);
    yyaxis left;   set(gca, 'YColor', 'k');
    CP(1) = plot(frequencies, amplitudes, 'k', 'LineWidth', 1);
    ylabel('Amplitude');
    hold on;
    yyaxis right;   set(gca, 'YColor', 'b');
    CP(2) = plot(frequencies, rad2deg(phases), 'b.', 'LineWidth', 0.5);
    % ylabel('Phase (°)');
    xlabel('Frequency (Hz)');
    CL = legend(CP, {'FFT Amplitude', 'Phase (°)'}, 'Orientation', 'horizontal');
    CL.Position(1) = sum(CA(2).Position([1,3])) - CL.Position(3);
    CL.Position(2) = sum(CA(2).Position([2,4]));
    CL.Box = 'off';   
    CT = title('FFT Magnitude and Phase', 'Unit', 'normalized');
    CT.Position(2) = CT.Position(2) + CL.Position(4)/CA(2).Position(4);
    grid on; 
    a = axis; axis([0 80 a(3) a(4)])

    
    % Panel C: FFT Reconstruction
    iCA = iCA +1;
    axes(CA(iCA));    % subplot(6, 2, 4);
    yyaxis right;   set(gca, 'YColor', [0.7350 0.0780 0.1840]);
    CP(2) = plot(total_time, reconstructed_fft, 'color', [0.7350 0.0780 0.1840]);
    ylabel(['Amplitude'], 'color', 'k');
    yyaxis left;   set(gca, 'YColor', 'k');
    CP(1) = plot(total_time, signal, ':k', 'linewidth', 2);
    xlim([0, 1]);
    set(gca, 'XTick', 0:0.25:2);
    CL = legend(CP, {'Original', 'Reconstructed'}, 'Orientation', 'horizontal');
    CL.Position(1) = sum(CA(3).Position([1,3])) - CL.Position(3);
    CL.Position(2) = sum(CA(3).Position([2,4]));
    CL.Box = 'off';   
    xlabel('Time');
    CT = title('FFT Reconstruction', 'Unit', 'normalized');
    CT.Position(2) = CT.Position(2) + CL.Position(4)/CA(3).Position(4);
    

    % Modified STFT code to ensure proper amplitude preservation
    window_length = round(fs/2);
    hop_size = round(window_length/4);
    window = hann(window_length);
    % Normalize window for energy preservation
    window = window / sqrt(sum(window.^2));
    [S, F, T] = spectrogram(signal, window, window_length-hop_size, [], fs);
    T = T + total_time(1);
    reconstructed_stft = istft(S, window_length, hop_size, length(window), fs);
    
    
    % Magnitudes and reconstructions to be plotted:
    Titles =          ["STFT Magnitude",       "CWT",                 "Long Wavelet (50 cycels)",     "Short Wavelet (5 cycles)",      "Z-scored Short Wavelet";...
                       "STFT Reconstruction",  "CWT Reconstruction",  "Long Wavelet Reconstruction",  "Short Wavelet Reconstruction",  "Z-scored Reconstruction"];
    Frequencies =     {F,                      f_cwt,                 freq_bins,                      freq_bins,                       freq_bins};
    Ts =              {T,                      total_time,            total_time,                     total_time,                      total_time};
    Magnitudes =      {S,                      W_cwt,                 W_long,                         W_short,                         W_short_z};              
    Reconstructions = {reconstructed_stft,     reconstructed_cwt,     reconstructed_long,             reconstructed_short,             reconstructed_z};
    
    
    for iPanel = 1:size(Titles, 2)
        % imagesc Panels:
        iCA = iCA +1;
        axes(CA(iCA)); 
        pcolor(Ts{iPanel}, Frequencies{iPanel}, abs(Magnitudes{iPanel}));  
        shading flat;
        title(Titles(1, iPanel));
        ylabel('Frequency (Hz)');
        CB = colorbar('location', 'west');
        colormap('turbo');
        hold on;
        for i = 1:6
            yline(f_theta*i, '--w', ['H' num2str(i)], 'Alpha', 0.5);
        end
        ylim([1 80]);  xlim([0, 1])
        set(gca, 'XTick', 0:0.25:2);
        if iPanel ~= size(Titles, 2); CA(iCA).XTickLabel = []; end
        if iPanel == size(Titles, 2); CA(iCA).XLabel.String = 'Time'; end
        CA(iCA).YAxisLocation = 'right';
        CB.Label.String = 'Power (dB)';
        CB.Position(1) = CA(iCA).Position(1) - 0.02;
        CB.Position(2) = CA(iCA).Position(2);
        CB.Position(3) = 0.012;
        CB.Position(4) = CA(iCA).Position(4);
        CB.AxisLocation = 'out';


        % Reconstruction Panels:
        iCA = iCA +1;
        axes(CA(iCA));    
        yyaxis right;   set(gca, 'YColor', [0.7350 0.0780 0.1840]);
        CP(2) = plot(total_time(1:length(Reconstructions{iPanel})), Reconstructions{iPanel}, 'color', [0.7350 0.0780 0.1840]);
        ylabel(['Amplitudes'], 'color', 'k');
        yyaxis left;   set(gca, 'YColor', 'k');
        CP(1) = plot(total_time, signal, ':k', 'linewidth', 2);
        xlim([0, 1]);
        set(gca, 'XTick', 0:0.25:2);
        if iPanel ~= size(Titles, 2); CA(iCA).XTickLabel = []; end
        if iPanel == size(Titles, 2); CA(iCA).XLabel.String = 'Time'; end
        title(Titles(2, iPanel));
    end
    

    % Overall figure adjustments
    CSP = sgtitle('Signal Transformations: From Ground Truth to Increasing Abstraction', 'FontSize', 12);
    % colormap(parula);
    set(gcf, 'Color', [0.95 0.95 0.95]);
    linkaxes(CA([1,3:end]), 'x');
    

    % Modified ISTFT function for better amplitude preservation
    function x = istft(S, window_length, hop_size, nfft, fs)
        num_frames = size(S, 2);
        expected_length = (num_frames-1)*hop_size + window_length;
        x = zeros(expected_length, 1);
        w = hann(window_length);
        % Normalize window
        w = w / sqrt(sum(w.^2));
        norm_buffer = zeros(expected_length, 1);

        for i = 1:num_frames
            start_idx = (i-1)*hop_size + 1;
            end_idx = start_idx + window_length - 1;
            frame = real(ifft(S(:,i), nfft));
            frame = frame(1:window_length);
            x(start_idx:end_idx) = x(start_idx:end_idx) + frame.*w;
            norm_buffer(start_idx:end_idx) = norm_buffer(start_idx:end_idx) + w.^2;
        end
        % Careful normalization to preserve amplitude
        x = x ./ (norm_buffer + eps);
    end
  
    
    saveas(CF, [char(outPath), 'transformationsDemo.fig']);


    % close all;
    RV = 1;  return;
    
end


