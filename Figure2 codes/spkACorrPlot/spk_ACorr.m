function RV = spk_ACorr(inPath, outPath, varargin)    
% spiking activity versus theta, for high and low speed.

    % parameters:
    IP = inputParser;
    addParameter(IP, 'cellCategory', ["region", "isInter"]);  % classfications to group spk together.
    IP.parse(varargin{:});
    P = IP.Results;
    
    CCatName = P.cellCategory;
    
    
    % load data and by speed:
    load([char(inPath), 'dataLite.mat']);
    data = dataLite;
    
    
    % read data:
    cellInd = data.read("cellInd");
    
    
    % obtain info for cell categories:
    CCat = classify_cells(data, CCatName);
    
    
    % by running speed:
    [dataBySpeed, dataSpeedSorted, DOF, meanSpeed, speedRange] = sortSpeed(data, 1.01, 0.5, 3);
    data = dataBySpeed;
    
    
    % calculated ACorr and PSD of ACorr for each cell:
    weight = zeros(1, length(meanSpeed));
    for iS = 1:length(data)
        for iD = 1:length(data{iS})
            countf = fprintf(['iS/iD: ', num2str(iS), '/', num2str(iD)]);
            
            % read sub data for all cells:
            dataSub = data{iS}(iD);
            spkAll = dataSub.read('spk');
            spkOffsetAll = dataSub.read('spkOffset');
            cellRegion = readType(dataSub, 'region');
            isIntern = readType(dataSub, 'isIntern');
            firingRate = readType(dataSub, 'firingRate');
            recordRate = dataSub.read('recordRate');
            T = dataSub.read('T');
            
            % calculate for individual cells:
            for iC = 1:length(cellInd)
                % read data for a single cell:
                spkSub = spkAll(iC, :);
                spkOffsetSub = spkOffsetAll(iC, :);
                spkT = T(full(spkSub)~=0)' + spkOffsetSub(full(spkSub)~=0);  % increase resolution.
                
                % calculate ACorr:
                if sum(spkSub) ~= 0
                    [ACsub(iC, :), lag] = xcorr(full(spkSub), full(spkSub), 'coeff');
                else 
                    [~, lag] = xcorr(full(spkSub), full(spkSub), 'coeff');
                    ACsub(iC, :) = zeros(1, length(spkSub)*2-1);
                end
                
                % calculate PSD of ACorr:
                dummyT = lag/recordRate;
                [frequency, PSDsub(iC, :), all_in_one] = psd1(movmean(ACsub(iC, :), 3), recordRate);
                
                % calculate PDF of 1/ISI:
                pdfISIsub(iC) = calculate_pdf_ISI(spkT);
%                     temp = 1./ISI{iS}{iC};
%                     tempAll = temp;
%                     for ik = 1:length(temp)-1
%                         temp = temp(1:end-1) + temp(2:end);    % this is wrong, but why this generate gamma band?!!!.
%                         tempAll = [tempAll , temp];
%                     end
%                     [numCounts, binEdges] = histcounts(tempAll, 0:1:512);
%                     binCenters = (binEdges(1:end-1) + binEdges(2:end))/2;
                1;
            end
            
            % calculate for cells categories:
            for pointerS = 1:numel(CCat.cellInd)
                % construct merged data for a category:
                [~, indC] = ismember(CCat.cellInd{pointerS}, cellInd);
                spkCatMerged = sum(spkAll(indC, :), 1);
                spkTCatMerged = nan(sum(spkCatMerged),1);
                    pointer = 1;
                    for iiC = 1:length(indC)
                        spkSub = spkAll(indC(iiC), :);
                        spkOffsetSub = spkOffsetAll(indC(iiC), :);
                        spkT = T(full(spkSub)~=0)' + spkOffsetSub(full(spkSub)~=0);
                        spkTCatMerged(pointer:pointer+sum(spkAll(indC(iiC), :))-1) = spkT;
                        pointer = pointer + sum(spkAll(indC(iiC), :));
                    end
                    spkTCatMerged = sort(spkTCatMerged);

                % calculate ACorr:
                if sum(spkCatMerged) ~= 0
                    [CatACsub(pointerS, :), lag] = xcorr(full(spkCatMerged), full(spkCatMerged), 'coeff');
                    if any(isnan(CatACsub(pointerS, :)))
                        error("Yu:?")
                    end
                else 
                    CatACsub(pointerS, :) = zeros(1, length(spkCatMerged)*2-1);
                end
                
                % calculate PSD of ACorr:
                dummyT = lag/recordRate;
                [frequency, CatPSDsub(pointerS, :), all_in_one] = psd1(movmean(CatACsub(pointerS, :), 3), recordRate);
                
                % calculate PDF of 1/ISI:
                CatPdfISIsub(pointerS) = calculate_pdf_ISI(spkTCatMerged);
            end
            
            
            % make average for each speed:
            if iD == 1
                AC{iS} = ACsub;
                PSD{iS} = PSDsub;
                pdfISI{iS} = pdfISIsub;
                
                CatAC{iS} = CatACsub;
                CatPSD{iS} = CatPSDsub;
                CatPdfISI{iS} = CatPdfISIsub;
            else
                AC{iS} = (AC{iS}*weight(iS) +ACsub)/(weight(iS)+1);
                PSD{iS} = (PSD{iS}*weight(iS) +PSDsub)/(weight(iS)+1);
                pdfISI{iS} = arrayfun(@(struct1, struct2)weightedAverageStructFields(struct1, struct2, weight(iS), 1), pdfISI{iS}, pdfISIsub);
                
                CatAC{iS} = (CatAC{iS}*weight(iS) +CatACsub)/(weight(iS)+1);
                CatPSD{iS} = (CatPSD{iS}*weight(iS) +CatPSDsub)/(weight(iS)+1);
                CatPdfISI{iS} = arrayfun(@(struct1, struct2)weightedAverageStructFields(struct1, struct2, weight(iS), 1), CatPdfISI{iS}, CatPdfISIsub);
            end
            weight(iS) = weight(iS) +1;
            
            fprintf(1, repmat('\b',1,countf));
            1;
        end
    end
    
    % save:
    dataSample = data{end}(end);
    save([char(outPath), 'workSpace.mat'], 'AC', 'dummyT', 'PSD', 'frequency', 'pdfISI', 'dataSample', 'DOF', 'meanSpeed', 'speedRange', ...
                                   'CCat', 'CatAC',     'CatPSD',           'CatPdfISI', '-v7.3');  % data by speed is too large.
    RV = 1;
end    