function RV = modbyHz(inPath, outPath, varargin)
% spike modulation by frequency.


    % parameters:
    inP = inputParser;
    addParameter(inP, 'cellCategory', ["region", "isInter"]);  % classfications to group spk together.
    
    % Fourier and wavelet parameters for phase and amp detection (cwtCoefficients):
    voicesPerOctave = 6;
    cwtFrequencyLimit = [1, 512];
    sigLevel = 1.0;    % threshold for (gamma) amplitude, only consider the top ratio sigLevel oscillations.
    
    inP.parse(varargin{:});
    CCatName = inP.Results.cellCategory;
    

    % load data:
    load([char(inPath), 'dataLite.mat']);
    data = dataLite;
    [channels, markers] = findMarker(data, 'maxTheta');
    data = data.interpAxis('channel', channels);
    
    
    % read data:
    cellInd = data.read("cellInd");
    for i = 1:length(CCatName)
        CCat{i} = readType(data, CCatName(i));
        CCatU{i} = unique(CCat{i});
        CatSzie(i) = length(CCatU{i});
    end
    firingRate = readType(data, 'firingRate');
    recordRate = data.read('recordRate');
    
    
    % recursive loop for arbitrary dimension of cellCategories (become a function since spk_ACorr):
    CatCellInd = cell([CatSzie,1]);
    pointerV = num2cell(nan(length(CatSzie), 1));    % vector pointer.
    pointerS = 0;                         % scalar pointer.
    while pointerS < prod(CatSzie)
        pointerS = pointerS +1;
        [pointerV{:}] = ind2sub(CatSzie, pointerS);   
        
        CatIndTemp = cellInd;
        for i = 1:length(pointerV)
            CTypeValue = readType(data, CCatName(i));
            CatIndTemp = intersect(CatIndTemp, cellInd(CTypeValue == CCatU{i}(pointerV{i})));
        end
        CatCellInd{pointerS} = CatIndTemp;    % it is the cellInd, not indCellInd.
    end
    
    
    % perform cwt for phase and amp:
    [data, CF, CA] = cwtData(outPath, data, voicesPerOctave, cwtFrequencyLimit, 'nSegPlot', 5);    
    cwtF = data.read('cwtF');
    
    
    % sort by speed:
    [dataBySpeed, dataSpeedSorted, DOF, meanSpeed] = sortSpeed(data, 1.01, 0.0, 3);
    data = dataBySpeed;
    for iS = 1:length(data)
        data{iS} = class_data.concatenate(data{iS}, 'T');
    end
    

    % core loop:
    for iS = 1:length(data)
        countf = fprintf(['iSpeed: ', num2str(iS), ' ; ']);
        T = data{iS}.read('T');
        spk = data{iS}.read('spk');
        cwtCoeff = data{iS}.read('cwtCoeff');
        spkOffsetAll = data{iS}.read('spkOffset');
        
        modStrength{iS} = nan(length(cellInd), length(cwtF));
        modStrength2{iS} = nan(length(cellInd), length(cwtF));
        modStrengthCat{iS} = nan(prod(CatSzie), length(cwtF));
        modStrength2Cat{iS} = nan(prod(CatSzie), length(cwtF));
        for iF = 1:length(cwtF)   % 37  %%%%
            countf3 = fprintf(['iF: ', num2str(iF), ' ; ']);
            NNans = sum(isnan(cwtCoeff(iF, :))); 
            NNoNans = length(T)-NNans; 
            [sortedValues, sortedIndices] = sort(abs(cwtCoeff(iF, :)), 'descend');
            indT = sortedIndices(NNans+(1:round(NNoNans*sigLevel)));
            indT = sort(indT);     % top percentage w/o nans.
            
            Tvalid = T(indT);
            cwtCoeffValid = cwtCoeff(iF, indT);
            
            % calculate modbyHz for individual cells:
            temp = class_data;
            temp = temp.write("T", Tvalid);
            temp = temp.write("phase(T)", angle(cwtCoeffValid));
            temp = temp.write("cellInd", cellInd);
            temp = temp.write("spkBool(cellInd, T)", spk(:, indT));
            PM = class_phaseModulation(temp, "T", "spkBool", "phase", 12, 0.5);
            
            modStrength{iS}(:, iF) = PM.read('concentration');
            modStrength2{iS}(:, iF) = PM.read('concentration2');
            
            % calculate for cells categories:
            spkCatMerged = nan(numel(CatCellInd), size(spk, 2));
            for pointerS = 1:numel(CatCellInd)
                [~, indC] = ismember(CatCellInd{pointerS}, cellInd);
                spkCatMerged(pointerS, :) = sum(spk(indC, :), 1);
            end
            temp = temp.write("cellCatInd", 1:numel(CatCellInd));
            temp = temp.write("spkBool(cellCatInd, T)", spkCatMerged(:,indT));
            CatPM = class_phaseModulation(temp, "T", "spkBool", "phase", 12, 0.5);
                
            modStrengthCat{iS}(:, iF) = CatPM.read('concentration');
            modStrength2Cat{iS}(:, iF) = CatPM.read('concentration2');
            
            fprintf(1, repmat('\b',1,countf3));
        end

        fprintf(1, repmat('\b',1,countf));
        1;
    end
    
    
    % save:
    dataSample = data{end}(end);
    save([char(outPath), 'modbyHz.mat'], 'modStrength', 'modStrength2', 'cwtF', 'PM', 'dataSample', 'DOF', 'meanSpeed',...
                                         'modStrengthCat', 'modStrength2Cat', 'CatPM', 'inP', 'CCatU', 'CatCellInd', '-v7.3');
    RV = 1;
    
end