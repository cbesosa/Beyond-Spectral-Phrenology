function [data, CF, CA] = cwtData(outPath, data, voicesPerOctave, cwtFrequencyLimit, varargin)

    % parameters:
    inP = inputParser;
    addParameter(inP, 'nSegPlot', 10);
    addParameter(inP, 'resoPlot', 1/512);    % resolution, will floor the step size.
    addParameter(inP, 'LPlot', 20);       % max length of each plot.
    
    inP.parse(varargin{:});
    nSegPlot = inP.Results.nSegPlot;
    resoPlot = inP.Results.resoPlot;
    LPlot = inP.Results.LPlot;
    
    cwtFrequency = 2.^(log(max(cwtFrequencyLimit))/log(2):-1/voicesPerOctave:log(min(cwtFrequencyLimit))/log(2));
    recordRate = data.read('recordRate');
    T = data.read('T');
    [AID, reAxis] = data.axisID("T");
    plotStepSize = floor(resoPlot./AID.mode);
    plotEndInd = ceil(LPlot./AID.mode);
        
    
    % preallocate:
    cwtCoefficients = nan(length(cwtFrequency), length(T));
    coi = nan(1, length(T));
    iF = 0;
    
    
    % non-isometric segments:
    [data2, dataArray] = data.segments("T");
    indT = 0;
    for iSeg = 1:length(dataArray)
        countf2 = fprintf(['cwt Seg# : ', num2str(iSeg),'/', num2str(length(dataArray))]);

        Tsub = dataArray(iSeg).read('T');
        LFPsub = dataArray(iSeg).read('LFP');
        
        fb = cwtfilterbank('SamplingFrequency', recordRate, ...
                           'SignalLength', length(LFPsub), ...
                           'FrequencyLimits', cwtFrequency([end,1]), ...
                           'VoicesPerOctave', voicesPerOctave);
        [cwtCoefficientsSub, frequencies, coiSub, fb2] = fb.wt(LFPsub);
        % use short-time Fourier transform:
%         [stftResult, f, t] = stft(signal, Fs, 'Window', hamming(256), 'OverlapLength', 128, 'FFTLength', 512);
%         amplitudeSpectrum = abs(stftResult);
%         phaseSpectrum = angle(stftResult);

        for iT = 1:length(Tsub)
            cwtCoefficientsSub(frequencies <= coiSub(iT)+0.001, iT) = nan;
        end
        if abs(frequencies(1) - cwtFrequency(1)) < 0.001
            cwtCoefficients(1:length(frequencies), indT+(1:length(Tsub))) = cwtCoefficientsSub;
            coi(indT+(1:length(Tsub))) = coiSub;
        else
            error("Yu: frequency missmatch")
        end
        indT = indT + length(Tsub);

        
        % plot cwt of segments:
        if iSeg <= nSegPlot
            iF = iF +1; CF(iF) = figure;
            CA{iF}(1) = subplot(3,1,1);
            ind2plot = 1:plotStepSize:min(length(Tsub),plotEndInd);
            imagesc(Tsub(ind2plot), frequencies, abs(cwtCoefficientsSub(:,ind2plot)));  
            axis tight;
            set(gca, 'YDir', 'normal', 'YScale', 'log');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title('CWT Magnitude Spectrum');
            colorbar('location', 'east');

            [~, ind8] = min(abs(frequencies - 8));

            CA{iF}(2) = subplot(3,1,2);
            plot(Tsub(ind2plot), angle(cwtCoefficientsSub(ind8,ind2plot)));  

            CA{iF}(3) = subplot(3,1,3);
            plot(Tsub(ind2plot), LFPsub(ind2plot));  hold on;
            plot(Tsub(ind2plot), real(cwtCoefficientsSub(ind8,ind2plot)));  hold on;

            linkaxes(CA{iF}, 'x');
            
            speed = mean(readState(dataArray(iSeg), 'speed'));
            saveas(CF(iF), [char(outPath), 'cwt_seg_', num2str(iSeg),'_speed_', num2str(round(speed)),'.fig'])
            1;
        end
        
        fprintf(1, repmat('\b',1,countf2));
        1;
    end
    
    
    % save data:
    data = data.write("cwtF", cwtFrequency);
    data = data.write("cwtCoeff(cwtF,T)", cwtCoefficients);
    

    % plot cwtPSD:
    iF = iF +1; CF(iF) = figure;
    meanCwt = nanmean(abs(cwtCoefficients'));
    psd = 2*meanCwt(2:end-1)./((cwtFrequency(1:end-2)-cwtFrequency(3:end)));
    subplot(1,2,1)
    loglog(cwtFrequency, meanCwt, '.-');
    subplot(1,2,2)
    loglog(cwtFrequency(2:end-1), psd, '.-');
    
    speed = mean(readState(data, 'speed'));
    saveas(CF(iF), [char(outPath), 'cwtPSD_','speed_', num2str(round(speed)),'.fig'])
    

    % plot all cwt:
    iF = iF +1; CF(iF) = figure;
    ind2plot = 1:plotStepSize:min(length(T),plotEndInd);
    imagesc(AID.mode*plotStepSize*(0:length(ind2plot)-1), cwtFrequency, abs(cwtCoefficients(:,ind2plot)));   % When you use imagesc with non-continuous or non-uniformly spaced vectors for the axes, you might not encounter an error directly, but the visualization could be misleading or not reflect the actual distribution of your data along those axes.
    axis tight;
    set(gca, 'YDir', 'normal', 'YScale', 'log');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('CWT Magnitude Spectrum');
    colorbar;

    % Overlay COI using a patch
    hold on; % Ensure the COI is added to the existing plot
    coiPatch = area(AID.mode*plotStepSize*(0:length(ind2plot)-1), coi(ind2plot), min(cwtFrequency));
    coiPatch.FaceColor = [1,1,1]*0.8;
    coiPatch.FaceAlpha = 0.5; % Make the COI patch semi-transparent
    coiPatch.EdgeColor = 'none'; % Remove the patch edge color
    ylim(cwtFrequency([end,1]));
    
    saveas(CF(iF), [char(outPath), 'cwtAll_','speed_', num2str(round(speed)),'.fig'])
    
        
    return
end