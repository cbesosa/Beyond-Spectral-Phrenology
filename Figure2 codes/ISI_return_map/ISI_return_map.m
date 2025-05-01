function RV = ISI_return_map(inPath, outPath)    
% spiking activity versus theta, for high and low speed.

    
    % load data and by speed:
    load([char(inPath), 'dataLite.mat']);
    data = dataLite;
    [dataBySpeed, dataSpeedSorted, nSegStr, meanSpeed] = sortSpeed(data, 2.01, 0.0, 2);
    data = dataBySpeed;
    for iS = 1:length(data)
        data{iS} = class_data.concatenate(data{iS}, 'T');
    end
    
    % calculated corr and PSD of corr for each cell:
    cellInd = data{1}(1).read("cellInd");
    weight = zeros(1, length(data));
    
    for iS = 1:length(data)
        for iD = 1:length(data{iS})
            countf = fprintf(['iS/iD: ', num2str(iS), '/', num2str(iD)]);
            
            dataSub = data{iS}(iD);
            spkOffsetAll = dataSub.read('spkOffset');
            spkAll = dataSub.read('spk');
            cellRegion = readType(dataSub, 'region');
            isInter = readType(dataSub, 'isInter');
            firingRate = readType(dataSub, 'firingRate');
            recordRate = dataSub.read('recordRate');
            T = dataSub.read('T');
            
            for iC = 1:length(cellInd)
                spkOffsetSub = spkOffsetAll(iC, :);
                spkSub = spkAll(iC, :);
                spkT = T(spkSub~=0)' + spkOffsetSub(spkSub~=0);
                
                if iD == 1  % pre-allocate:
                    ISI{iS}{iC} = nan(1,100*3600*10);
                    ISI_pair{iS}{iC} = nan(2,100*3600*10);
                    pointer_ISI(iS, iC) = 0;
                    pointer_ISI_pair(iS, iC) = 0;
                end
                ISI_sub = spkT(2:end)-spkT(1:end-1);
                ISI_pair_sub = [ISI_sub(1:end-1); ISI_sub(2:end)];  % first row is n, second row is n+1.
                
                ISI{iS}{iC}(pointer_ISI(iS, iC) +(1:size(ISI_sub,2))) = ISI_sub;
                ISI_pair{iS}{iC}(:, pointer_ISI_pair(iS, iC) +(1:size(ISI_pair_sub,2))) = ISI_pair_sub;
                
                pointer_ISI(iS, iC) = pointer_ISI(iS, iC) +size(ISI_sub,2);
                pointer_ISI_pair(iS, iC) = pointer_ISI_pair(iS, iC) +size(ISI_pair_sub,2);
                1;
            end
            
            fprintf(1, repmat('\b',1,countf));
            1;
        end
    end
    
    for iS = 1:length(data)
        for iC = 1:length(cellInd)
            ISI{iS}{iC} = ISI{iS}{iC}(1:pointer_ISI(iS, iC));
            ISI_pair{iS}{iC} = ISI_pair{iS}{iC}(:, 1:pointer_ISI_pair(iS, iC));
            
            [n, xedgesOut, yedgesOut] = histcounts2(ISI_pair{iS}{iC}(1,:), ISI_pair{iS}{iC}(2,:), ...
                                                    10.^(-3:0.1:2),  10.^(-3:0.1:2));
            n_normalized{iS}{iC} = n/max(n(:));
            if all(n(:) == 0)
                n_normalized{iS}{iC} = n;
            end
        end
    end
    
    
    %% plot by cells:
    figArr = class_figArr;
    figArr.X = ones(1,length(data));  figArr.dx = [0.25 0.08*ones(1,length(data)-1) 0.50];
    figArr.Y = [1 1 1];     figArr.dy = [0.22 0.19 0.07 0.55]; 
    figArr.GM = [1:length(data); length(data)+1:2*length(data); 2*length(data)+1:3*length(data)];
    figArr.width = 0.8;
    figArr.AsRatio = 13/10;
    figArr.screenOffset = [0.2, 0.2];
    figArr.figOut = 'ISI return map';
    figArr.hideInnerAxis = 'off';

    if 0
    for iC = 1:length(cellInd)
        [CF(iC), CA{iC}] = create(figArr);
        
        for iS = 1:length(data)
            
            CP1(iS) = plot(CA{iC}(iS), ISI_pair{iS}{iC}(1,:), ISI_pair{iS}{iC}(2,:),'.');
            CA{iC}(iS).XScale = 'log';
            CA{iC}(iS).YScale = 'log';
            CA{iC}(iS).YLim = 10.^[-3, 2];
            CA{iC}(iS).XLim = 10.^[-3, 2];
            CA{iC}(iS).Title.String = ['mean speed (cm/s)',10, '= ', num2str(meanSpeed(iS))];
            
            
            
            
            binCentersX = (xedgesOut(1:end-1) + xedgesOut(2:end))/2;
            binCentersY = (yedgesOut(1:end-1) + yedgesOut(2:end))/2;
            
            CP2(iS) = pcolor(CA{iC}(length(data)+iS), binCentersX, binCentersY, n_normalized{iS}{iC});
            shading(CA{iC}(length(data)+iS), 'flat');
            shading(CA{iC}(length(data)+iS), 'interp');
            CA{iC}(length(data)+iS).XScale = 'log';
            CA{iC}(length(data)+iS).YScale = 'log';
            CA{iC}(length(data)+iS).YLim = 10.^[-3, 2];
            CA{iC}(length(data)+iS).XLim = 10.^[-3, 2];
            CA{iC}(length(data)+iS).XLabel.String = 'ISI (second)';
 
            CA{iC}(length(data)+iS).Colormap = class_color.cut(class_color.map(304), 0.1, 1);
            CA{iC}(length(data)+iS).CLim(1) = 0;
            if iS > 1
                CA{iC}(length(data)+iS).CLim(2) = CA{iC}(length(data)+1).CLim(2);
            end
            
            CP3(iS) = pcolor(CA{iC}(2*length(data)+iS), 1./binCentersX, 1./binCentersY, n_normalized{iS}{iC});
            shading(CA{iC}(2*length(data)+iS), 'flat');
            shading(CA{iC}(2*length(data)+iS), 'interp');
            CA{iC}(2*length(data)+iS).XScale = 'log';
            CA{iC}(2*length(data)+iS).YScale = 'log';
            CA{iC}(2*length(data)+iS).YLim = 10.^[-1, 3];
            CA{iC}(2*length(data)+iS).XLim = 10.^[-1, 3];
            CA{iC}(2*length(data)+iS).XLabel.String = 'ISI (Hz)';
 
            CA{iC}(2*length(data)+iS).Colormap = class_color.cut(class_color.map(304), 0.1, 1);
            CA{iC}(2*length(data)+iS).CLim(1) = 0;
            if iS > 1
                CA{iC}(2*length(data)+iS).CLim(2) = CA{iC}(2*length(data)+1).CLim(2);
            end
                
            1;
        end
        
        linkaxes(CA{iC}(1:2*length(data)), 'x');
        linkaxes(CA{iC}(1:2*length(data)), 'y');
        CB = colorbar('peer', CA{iC}(end));
        CB.Location = 'east';
        CB.Position(1) = sum(CA{iC}(end).Position([1,3]))+0.1;
        CB.Position(3) = 0.02;
        CB.Label.String = 'Occurance Normalized';
        
        STH = suptitle( ['ISI return map', ...
                         'Cell ID = ',  num2str(cellInd(iC)), '; region = ', char(cellRegion(iC))...
                         10, 'firingRate = ', num2str(firingRate(iC)), '; isIntern = ', num2str(isInter(iC)), ...
                         ]);
        STH.FontSize = 10;
        
        saveas(CF(iC), [char(outPath), 'ISI_return_map', num2str(cellInd(iC)),'_region_', char(cellRegion(iC)), '_isIntern_', num2str(isInter(iC)),'.fig']);
        1;
    end
    end
    
    
    % to save:
    cellTypeRegion_all = [];
    frequency1 = [];
    frequency2 = [];
    n_normalized_merged_all = [];
    
    % plotting loop:
    cellType = repmat("", 1, length(cellInd));
    cellType(isInter == 0) = "pyd";
    cellType(isInter == 1) = "inter";
    cellType(isInter == -1) = "unknown";
    
    cellTypeAll = unique(cellType);
    cellTypeAll(cellTypeAll == "") = [];
    cellRegionAll = unique(cellRegion);
    cellRegionAll(cellRegionAll == "") = [];
    % plot region and cell type merged:
    for iR = 1:length(cellRegionAll)
        for iT = 1:length(cellTypeAll)
            [CF, CA] = create(figArr);

            for iS = 1:length(data)
                
                ind = intersect(find(cellType == cellTypeAll(iT)), find(cellRegion == cellRegionAll(iR)));  
                nCell(iR, iT) = length(ind);
                ISI_pair_merged = nan(2,0);
                weight = 0;
                n_normalized_merged = n_normalized{1}{1}*0;
                for iC = ind
                    countf = fprintf(['iC: ', num2str(iC), '/', num2str(ind(end))]);
                    ISI_pair_merged = [ISI_pair_merged, ISI_pair{iS}{iC}];
                    if firingRate(iC) < 1
                        dWeight = firingRate(iC);
                    else
                        dWeight = 1;
                    end
                    n_normalized_merged = n_normalized_merged*weight + n_normalized{iS}{iC}*dWeight;
                    weight = weight + dWeight;
                    n_normalized_merged = n_normalized_merged/weight;
                    if weight == 0
                        n_normalized_merged(:) = 0;
                    end
                    fprintf(1, repmat('\b',1,countf));
                end
                    
                % CP1(iS) = plot(CA(iS), ISI_pair_merged(1,:), ISI_pair_merged(2,:),'.', 'linewidth', 0.5);
                plot(CA(iS), ISI_pair_merged(1,:), ISI_pair_merged(2,:),'.', 'linewidth', 0.5);
                CA(iS).XScale = 'log';
                CA(iS).YScale = 'log';
                CA(iS).YLim = 10.^[-3, 2];
                CA(iS).XLim = 10.^[-3, 2];
                CA(iS).Title.String = ['mean speed (cm/s)',10, '= ', num2str(meanSpeed(iS))];


%                 [n, xedgesOut, yedgesOut] = histcounts2(ISI_pair_merged(1,:), ISI_pair_merged(2,:), ...
%                                                         10.^(-3:0.1:2),  10.^(-3:0.1:2));

                binCentersX = (xedgesOut(1:end-1) + xedgesOut(2:end))/2;
                binCentersY = (yedgesOut(1:end-1) + yedgesOut(2:end))/2;

                CP2(iS) = pcolor(CA(length(data)+iS), binCentersX, binCentersY, n_normalized_merged);
                shading(CA(length(data)+iS), 'flat');
                shading(CA(length(data)+iS), 'interp');
                CA(length(data)+iS).XScale = 'log';
                CA(length(data)+iS).YScale = 'log';
                CA(length(data)+iS).YLim = 10.^[-3, 2];
                CA(length(data)+iS).XLim = 10.^[-3, 2];
                CA(length(data)+iS).XLabel.String = 'ISI (second)';

                CA(length(data)+iS).Colormap = class_color.cut(class_color.map(304), 0.1, 1);
                CA(length(data)+iS).CLim(1) = 0;
                if iS > 1
                    CA(length(data)+iS).CLim(2) = CA(length(data)+1).CLim(2);
                end

                CP3(iS) = pcolor(CA(2*length(data)+iS), 1./binCentersX, 1./binCentersY, n_normalized_merged);
                shading(CA(2*length(data)+iS), 'flat');
                shading(CA(2*length(data)+iS), 'interp');
                CA(2*length(data)+iS).XScale = 'log';
                CA(2*length(data)+iS).YScale = 'log';
                CA(2*length(data)+iS).YLim = 10.^[-1, 3];
                CA(2*length(data)+iS).XLim = 10.^[-1, 3];
                CA(2*length(data)+iS).XLabel.String = 'ISI (Hz)';

                CA(2*length(data)+iS).Colormap = class_color.cut(class_color.map(304), 0.1, 1);
                CA(2*length(data)+iS).CLim(1) = 0;
                if iS > 1
                    CA(2*length(data)+iS).CLim(2) = CA(2*length(data)+1).CLim(2);
                end
                
                %  to save:
                if iS == 1
                    iS_all = [];
                    frequency1 = 1./binCentersX;
                    frequency2 = 1./binCentersY;
                end
                iS_all = [iS_all, iS];
                n_normalized_merged_all(iR, iT, iS,:,:) = n_normalized_merged;

                1;
            end

            linkaxes(CA(1:2*length(data)), 'x');
            linkaxes(CA(1:2*length(data)), 'y');
            CB = colorbar('peer', CA(end));
            CB.Location = 'east';
            CB.Position(1) = sum(CA(end).Position([1,3]))+0.1;
            CB.Position(3) = 0.02;
            CB.Label.String = 'Occurance Normalized';

            STH = suptitle( ['ISI return map', ...
                             'Cell ID = ',  'type merged', '; region = ', char(cellRegionAll(iR))...
                             10, '; cellType = ', char(cellTypeAll(iT)), ...
                             ]);
            STH.FontSize = 10;

            saveas(CF, [char(outPath), 'ISI_return_map', '_merged','_region_', char(cellRegionAll(iR)), '_type_', char(cellTypeAll(iT)),'.fig']);
            
            1;
        end
    end
    
    ISI_merged = class_data;
    ISI_merged = ISI_merged.write('iS', iS_all);
    ISI_merged = ISI_merged.write('meanSpeed(iS)', meanSpeed);
    ISI_merged = ISI_merged.write('F1', frequency1);
    ISI_merged = ISI_merged.write('F2', frequency2);
    ISI_merged = ISI_merged.write('n_normalized_merged(region, type, iS, F1, F2)', n_normalized_merged_all);
    ISI_merged = ISI_merged.write('region', cellRegionAll);
    ISI_merged = ISI_merged.write('type', cellTypeAll);
    ISI_merged = ISI_merged.write('nCell(region, type)',nCell);
    save([char(outPath), 'ISI_merged.mat'], 'ISI_merged', '-v7.3')
    
    close all;
    RV = 1; return;

end







