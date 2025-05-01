

function RV = loadData(inPath, outPath, par, varargin)
% load all necessary data from raw.

    % parameters:
    IP = inputParser;
    addParameter(IP, "testL", 30);   % maxTimeEachStateValue for test data.
    addParameter(IP, "fullL", 3600);
    addParameter(IP, "state2exhaust", ["isRun", "isSleep", "hasPosition"]);
    parse(IP, varargin{:}); 
    P = IP.Results;
    
    if par == "test"
        maxTimeEachStateValue = P.testL;  
    elseif par == "full"
        maxTimeEachStateValue = P.fullL;  
    end

    load([char(inPath), 'data.mat']); 

    % dataset specific entries need to be filled:
    data = data.write("channel", [], "value");
    data = data.write("oriIndex(channel)", [], "value");
    data = data.write("shank(channel)", [], "value");
    data = data.write("region(channel)", [], "value");
    data = data.write("channelMarker(channel)", {}, "value");
    
    data = data.write("T", [], "value");
    data = data.write("LFP(channel, T)", [], "value");
    data = data.write("stateName", string([]), "value");
    data = data.write("state(stateName, T)", [], "value");
    
    data = data.write("cellInd", [], "value");
    data = data.write("spk(cellInd,T)", [], "value");
    data = data.write("spkOffset(cellInd,T)", [], "value");    % spk timing offset to the down sampled time grid.
    data = data.write("cellType(typeName,cellInd)", {}, "value");  % use 0/1 insteaad of string in case there are cells subjected to multiply types.
    data = data.write("typeName", string([]), "value");
    
    % fill the entries:
    if contains(inPath, 'IZ') && contains(inPath, 'sess')
        data = loadData_Zutshi(data);
    elseif contains(inPath, 'i01_maze')
        data = loadData_hc5(data);    %  data loaded must be continuous for velocity vector calcualtion.
    elseif contains(inPath, '_2') && length(char(data.read("animalNameSession"))) == 13
        data = loadData_LFPEEG(data);
    elseif all(isstrprop(char(data.read("animalNameSession")), 'digit')) && length(char(data.read("animalNameSession"))) == 4
        data = loadData_Burke_Nick(data);
    elseif data.read("datasetName") == "English"
        data = loadData_English(data);
    elseif data.read("datasetName") == "hc6"
        data = loadData_hc6(data);
    elseif data.read("datasetName") == "000552"
        data = loadData_000552(data);
    elseif data.read("animalName") == "539"
        data = loadData_539(data);
    else
        error('Yu:')
    end
    
    % shorten data base on whether is test or full:
    X = readState(data, "X");
    Y = readState(data, "Y");
    hasPosition = zeros(size(X));
    hasPosition(X ~= 0 & Y ~= 0 & ~isnan(X) & ~isnan(Y)) = 1;
    data = writeState(data, "hasPosition", hasPosition);
    data = shortenData(data, P.state2exhaust, maxTimeEachStateValue);    % length to keep!.
    
    
    % interp some nans in positions:
    Xori = readState(data, 'X');    Yori = readState(data, 'Y');
    T = data.read("T");
    nanInd = isnan(Xori) | isnan(Yori);
    Xori(nanInd) = [];  Yori(nanInd) = [];  
    T2 = T;   T2(nanInd) = [];  

    X = interp1(T2, Xori, T);   Y = interp1(T2, Yori, T);
    dT2 = T2(2:end)-T2(1:end-1);
    indGap = find(dT2 > 0.25);
    for i = 1:length(indGap)
        X(T > T2(indGap(i)) & T < T2(indGap(i)+1)) = nan;
        Y(T > T2(indGap(i)) & T < T2(indGap(i)+1)) = nan;
    end
    numberofnansInterpted = sum(nanInd) - sum(isnan(X) | isnan(Y));
    data = data.write('notes', [data.read('notes'), 10, 'numberofnansInterpted = ', num2str(numberofnansInterpted)]);
    data = writeState(data, 'X', X);    data = writeState(data, 'Y', Y);
    
    
    % add velocity vector:
    velocityX = nan*X;  
    velocityX(1) = (X(2)-X(1))*data.read('recordRate');
    velocityX(2:end-1) = (X(3:end)-X(1:end-2))*data.read('recordRate')/2;
    velocityX(end) = (X(end)-X(end-1))*data.read('recordRate');
    velocityY = nan*Y;
    velocityY(1) = (Y(2)-Y(1))*data.read('recordRate');
    velocityY(2:end-1) = (Y(3:end)-Y(1:end-2))*data.read('recordRate')/2;
    velocityY(end) = (Y(end)-Y(end-1))*data.read('recordRate');
    data = writeState(data, 'velocityX', velocityX);
    data = writeState(data, 'velocityY', velocityY);
    data = writeState(data, 'speed', sqrt(velocityX.^2 +velocityY.^2));    % check if it is the same for hc5 data.
    
    
    % dimension consistency for:
    channelMarker = data.read('channelMarker');
    if isempty(channelMarker)
        channel = data.read('channel');
        data = data.write("channelMarker(channel)", repmat({string([])},length(channel),1), "value");
    end
    data = data.sortAxis("channel", 'ascend');
    
    
    % add spiking rate:
    spk = data.read("spk");
    if ~isempty(spk)
        for i = 1:size(spk, 1)
            firingRate(i) = sum(spk(i,:))/length(spk(i,:))*data.read("recordRate");
        end
        data = writeType(data, "firingRate", firingRate);
    end
    
    
    % save:
    L = round(length(data.read("T"))/data.read("recordRate"));
    fileID = fopen([char(outPath),'Length = ', num2str(L), ' seconds'], 'wb');
    fclose(fileID);
    data = data.write('notes', [data.read('notes'), 10, 'total raw data Length = ', num2str(L), ' seconds']);
    t0 = tic;
    save([char(outPath), 'data.mat'], 'data','-v7.3');
    t1 = toc(t0);
    fileID = fopen([char(outPath),'time_2_save(data) = ', num2str(round(t1)), ' seconds'], 'wb');
    fclose(fileID);
    
    RV = 1; return;        
end


function data = loadData_hc6(data)


    dataPath = char(data.read("dataPath"));
    LFPfolderPath = [dataPath, 'EEG\'];
    LFPfiles = dir(LFPfolderPath);  LFPfiles = LFPfiles(~[LFPfiles.isdir]);
    for i = 1:length(LFPfiles)
        LFPfilesNameShort = strrep(strrep(LFPfiles(i).name, lower([char(data.read('animalName')),'eeg']), ''),'.mat','');
        LFPfilesNameShort = strrep(strrep(LFPfilesNameShort, [char(data.read('animalName')),'eeg'], ''),'.mat','');
        parts = split(LFPfilesNameShort, '-');
        day(i) = str2num(parts{1});
        session(i) = str2num(parts{2});
        tetrode(i) = str2num(parts{3});
        1;
    end
    

    %% initialization & parameters:
    load([LFPfolderPath, LFPfiles(1).name]);   % eeg loaded.
    samprate = eeg{day(1)}{session(1)}{tetrode(1)}.samprate;
    recordRate = samprate;
    data = data.write("recordRate", recordRate, "table");
    
               
    %% load channel indices:
    load([dataPath, lower(char(data.read('animalName'))), 'tetinfo.mat']);   % tetinfo loaded.
    oriIndex = [];   geoIndex = [];   shank = [];   region = [];
    day2read = str2num(data.read("animalSession")); 
    for uDay = day2read
        for daySession = unique(session(day == uDay))
            if isempty(oriIndex)
                [oriIndex, depth, numcellsMin, numcellsMax] = deal(nan(1, length(tetinfo{uDay}{daySession})));
                region = repmat("", [1, length(tetinfo{uDay}{daySession})]);
            end
            if daySession == min(unique(session(day == uDay)))
                [numcellsMinDay{uDay}, numcellsMaxDay{uDay}] = deal(nan(1, length(tetinfo{uDay}{daySession})));
            end
            for i = 1:length(tetinfo{uDay}{daySession})
                if isempty(tetinfo{uDay}{daySession}{i}) 
                    continue; 
                end
                if tetinfo{uDay}{daySession}{i}.depth{:} == 0
                    continue;
                end   
                if ~isfield(tetinfo{uDay}{daySession}{i}, 'area')
                    continue;
                end
                oriIndex(i) = i;
                depth(i) = tetinfo{uDay}{daySession}{i}.depth{:};
                region(i) = string(tetinfo{uDay}{daySession}{i}.area);
                
                numcellsMin(i) = min(numcellsMin(i), tetinfo{uDay}{daySession}{i}.numcells);
                numcellsMax(i) = max(numcellsMax(i), tetinfo{uDay}{daySession}{i}.numcells);
                
                numcellsMinDay{uDay}(i) = min(numcellsMinDay{uDay}(i), tetinfo{uDay}{daySession}{i}.numcells);
                numcellsMaxDay{uDay}(i) = max(numcellsMaxDay{uDay}(i), tetinfo{uDay}{daySession}{i}.numcells);  % the numcells are consistent within a day except few exceptions.
            end
            1;
        end
    end
    depth(isnan(oriIndex)) = [];
    region(isnan(oriIndex)) = [];
    oriIndex(isnan(oriIndex)) = [];
    [depth, ind] = sort(depth);
    region = region(ind);
    oriIndex = oriIndex(ind);
    
    geoIndex = 1:length(oriIndex);
    shank = oriIndex;
    
    data = data.write("oriIndex(channel)", oriIndex, "value");
    data = data.write("channel", geoIndex, "value");
    data = data.write("shank(channel)", shank, "value");
    data = data.write("region(channel)", region, "value");
    data = data.write("depth(channel)", depth*0.0265, "depth");   %  *0.0265 in mm.
    
    
    %% load LFP channels:
    Amplification = 1000;
    SampleRate_LFP = samprate;
    
    % load from the EEG folder  (day, session, tetrode data)
    LFP = [];
    T = [];
    dayState = [];
    sessionState = [];
    X = [];
    Y = [];
    Angle = [];
    isRun = [];
    isSleep = [];
    for uDay = day2read
        LFPday = [];
        Tday = [];
        for daySession = unique(session(day == uDay))
            countf = fprintf(['Loading LFP Day ', num2str(uDay), ', Session ', num2str(daySession)]);
            LFPsession = [];
            if exist([char(dataPath), char(lower(data.read('animalName'))), 'pos', num2str(uDay, '%02d'),'.mat'], 'file') == 2
                load([char(dataPath), char(lower(data.read('animalName'))), 'pos', num2str(uDay, '%02d'),'.mat']);    % pos loaded;
                load([char(dataPath), char(lower(data.read('animalName'))), 'task', num2str(uDay, '%02d'),'.mat']);    % task loaded;
            else
                load([char(dataPath), char((data.read('animalName'))), 'pos', num2str(uDay, '%02d'),'.mat']);    % pos loaded;
                load([char(dataPath), char((data.read('animalName'))), 'task', num2str(uDay, '%02d'),'.mat']);    % task loaded;
            end
            for iTetrode = oriIndex
                ind1 = find(day == uDay);
                ind2 = intersect(ind1, find(session == daySession));
                ind3 = intersect(ind2, find(tetrode == iTetrode));
                load([LFPfolderPath, LFPfiles(ind3).name]);    % eeg loaded;

                LFPsub = eeg{uDay}{daySession}{iTetrode}.data';
                Tsub = eeg{uDay}{daySession}{iTetrode}.starttime + ...
                          (0:length(LFPsub)-1)/samprate;
                
                % match LFP grids to pos grids:
                TposSub = pos{uDay}{daySession}.data(:,1);   
                indTemp = find(Tsub> TposSub(1) & Tsub< TposSub(end));
                Tsub = Tsub(indTemp);
                LFPsub = LFPsub(indTemp);
                
                if iTetrode ~= oriIndex(1)      % dealing with inconsistent number of records among tetrodes.
                    if ~isequal(Tsub, Tsession)
                        TsubOld = Tsub;
                        Tsub = max(Tsub(1), Tsession(1)):1/samprate:min(Tsub(end), Tsession(end));
                        LFPsub = interp1(TsubOld, LFPsub, Tsub);
                        LFPsession = interp1(Tsession, LFPsession', Tsub(:))';
                        1;
                    end
                end
                LFPsession = [LFPsession; LFPsub(:)'];
                Tsession  = Tsub;
                
                1;
            end
            Xsub = interp1(TposSub, pos{uDay}{daySession}.data(:,2), Tsession);
            Ysub = interp1(TposSub, pos{uDay}{daySession}.data(:,3), Tsession);
            AngleSub = interp1(TposSub, pos{uDay}{daySession}.data(:,4), Tsession);
            
            LFPday = [LFPday, LFPsession];
            Tday = [Tday, Tsession];
            dayState = [dayState, Tsession*0+uDay];
            sessionState = [sessionState,  Tsession*0+daySession];
            if length(task{uDay}) < daySession
                isRun = [isRun, Tsession*0-1];
                isSleep = [isSleep, Tsession*0-1];
            else
                if isfield(task{uDay}{daySession}, 'type')
                    if task{uDay}{daySession}.type == "run"
                        isRun = [isRun, Tsession*0+1];
                    else
                        isRun = [isRun, Tsession*0+0];
                    end
                    if task{uDay}{daySession}.type == "sleep"
                        isSleep = [isSleep, Tsession*0+1];
                    else
                        isSleep = [isSleep, Tsession*0+0];
                    end
                else
                    isRun = [isRun, Tsession*0-1];
                    isSleep = [isSleep, Tsession*0-1];
                end
            end
            X = [X, Xsub];
            Y = [Y, Ysub];
            Angle = [Angle, AngleSub];
            fprintf(1, repmat('\b',1,countf));
            1;
        end
        LFP = [LFP, LFPday];
        T = [T, Tday+uDay*3600*24];
        1;
    end
    LFP = LFP/Amplification;
    X(X == 0) = nan;
    Y(Y == 0) = nan;
    Angle(Angle == 0) = nan;
    
    data = data.write("T", T, "value");
    data = data.write("LFP(channel, T)", LFP, "value");
    data = data.write("LFP_std", std(LFP(:)), "table");
    
    data = writeState(data, 'day', dayState);
    data = writeState(data, 'session', sessionState);
    data = writeState(data, 'X', X);
    data = writeState(data, 'Y', Y);
    data = writeState(data, 'angle', Angle);
    data = writeState(data, 'isRun', isRun);
    data = writeState(data, 'isSleep', isSleep);
   
      
    %% load state channels: (time series that are functions of time, e.g. position, speed, stimulation, etc.)
    data = data.write('notes', [data.read('notes'), 10, 'Full LFP records with nans for missing position data']);
                   
    
    %% load events: (, can be converted to state channels by function)
    % bad performance when deviding into segments, abandon?
    
    
    
    %% load spike channels:
    cellInd = 0;
    spkT = {};
    load([char(dataPath), char(lower(data.read('animalName'))), 'cellinfo','.mat']);    % cellinfo loaded;
    for uDay = day2read

        Ncells2read = numcellsMinDay{uDay} *nan;
        indTemp = find(numcellsMinDay{uDay} - numcellsMaxDay{uDay} == 0 & numcellsMinDay{uDay} > 0);
        Ncells2read(indTemp) = numcellsMinDay{uDay}(indTemp);
        Ncells2read = Ncells2read(oriIndex);
        
        if exist([char(dataPath), char(lower(data.read('animalName'))), 'spikes', num2str(uDay, '%02d'),'.mat'], 'file') == 2
            load([char(dataPath), char(lower(data.read('animalName'))), 'spikes', num2str(uDay, '%02d'),'.mat']);    % spikes loaded;
        else
            load([char(dataPath), char((data.read('animalName'))), 'spikes', num2str(uDay, '%02d'),'.mat']);    % spikes loaded;
        end
        iChannel = 0;
        for iTetrode = oriIndex
            iChannel = iChannel +1;
            if isnan(Ncells2read(iChannel))
                continue;
            end
            countf = fprintf(['Loading spk Day ', num2str(uDay), ', iChannel ', num2str(iChannel)]);
            for iCellSub = 1:Ncells2read(iChannel)
                spkTsub = [];
                cellinfoStruct = [];
                for daySession = unique(session(day == uDay))
                    spkStruct = spikes{uDay}{daySession}{iTetrode}{iCellSub};
                    if ~isempty(spkStruct)
                        cellinfoStruct2 = cellinfo{uDay}{daySession}{iTetrode}{iCellSub};
                        if ~isempty(cellinfoStruct)   % check if I am looking at the same cell across sessions.
                            if abs(cellinfoStruct2.spikewidth - cellinfoStruct.spikewidth) > 0.01
                                error("Yu: please check"); 
                            end
                        end
                        cellinfoStruct = cellinfoStruct2;
                    end
                    if ~isempty(spkStruct)
                        if ~isempty(spkStruct.data)
                            spkTsub = [spkTsub, spkStruct.data(:,1)'];
                        end
                    end
                    1;
                end
                if isempty(spkTsub)
%                     fprintf(1, repmat('\b',1,countf));
                    continue; 
                end
                cellInd = cellInd +1;
                spkT{cellInd} = spkTsub +uDay*3600*24;
                nSpks(cellInd) = length(spkTsub);
                if isempty(cellinfoStruct.spikewidth)
                    spkWidth(cellInd) = 0;
                else
                    spkWidth(cellInd) = cellinfoStruct.spikewidth;
                end
                cellChannel(cellInd) = iChannel;
                cellRegion(cellInd) = string(cellinfoStruct.area);
                if isfield(cellinfoStruct, 'type')
                    if cellinfoStruct.type == "principal"
                        isPyd(cellInd) = 1;
                    else
                        isPyd(cellInd) = 0;
                    end
                    if cellinfoStruct.type == "inter"
                        isInter(cellInd) = 1;
                    else
                        isInter(cellInd) = 0;
                    end
                else
                    isPyd(cellInd) = -1;
                    isInter(cellInd) = -1;
                end
            end
            fprintf(1, repmat('\b',1,countf));
            1;
        end
    end

    
    
    boolSpikes = sparse(cellInd, length(T));
    spkOffset = sparse(cellInd, length(T));
    for iC = 1:cellInd
        countf = fprintf(['interp cells to LFP, Cell ', num2str(iC), '.']);
        indT = interp1(T, 1:numel(T), spkT{iC}, 'nearest', 'extrap');
        offset = spkT{iC} - T(indT);
        boolSpikes(iC, indT) = 1;
        spkOffset(iC, indT) = offset;
        fprintf(1, repmat('\b',1,countf));
    end
    
    
    data = data.write("cellInd", 1:cellInd, "value");   % cellInd is the totClu in the orginal data.
    data = data.write("spk(cellInd,T)", boolSpikes, "value");
    data = data.write("spkOffset(cellInd,T)", spkOffset, "value");
    data = writeType(data, "isInter", isInter);
    data = writeType(data, "isPyd", isPyd);
    data = writeType(data, "channel", cellChannel);
    data = writeType(data, "region", cellRegion);
    data = writeType(data, "nSpks", nSpks);
    data = writeType(data, "spkWidth", spkWidth);
    
    
    1;
end



