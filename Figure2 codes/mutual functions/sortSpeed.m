function [dataBySpeed, dataSpeedSorted, DOF, meanSpeed, speedRange] = sortSpeed(data, segL_, overlap_, N_, Vrange_)
% N_ : number of groups of speed with equal number of segments
% Vrange_ = [V1,V2; V3,V4; ...], Vrange overload N_.
    
    if nargin < 2;  segL_ = 1.05;   end
    if nargin < 3;  overlap_ = 0.5;  end
    if nargin < 4;  N_ = 3;  end
    if nargin < 5;  Vrange_ = [];  end
    
    
    %% segments sorted by speed:
    [dataArray, ~] = data.equalSegments("T", segL_, overlap_);
    speed = nan(length(dataArray), 1);
    for i = 1:length(dataArray)
        speed(i) = mean(readState(dataArray(i), 'speed'));
        dataArray(i) = dataArray(i).write('meanSpeed', speed(i), 'table');
        dataArray(i) = dataArray(i).write('segL', segL_, 'table');
%         % remove some records to save space:
%         dataArray(i) = dataArray(i).remove("iRipple");
%         dataArray(i) = dataArray(i).remove("tagRipple");
%         dataArray(i) = dataArray(i).remove("ripple");
%         dataArray(i) = dataArray(i).remove("speed_ori");
%         dataArray(i) = dataArray(i).remove("LFP_ori");
%         dataArray(i) = dataArray(i).remove("LFP_preFixed");
%         dataArray(i) = dataArray(i).remove("spkTimes_ori");
        
    end
    [speed, ind] = sort(speed, 'ascend');
    dataSpeedSorted = dataArray(ind);
    
    
    %% determine the speedRange:
    speedNoNan = speed(~isnan(speed));
    % ticksTemp = round(1:(length(speedNoNan)-1)/N_:length(speedNoNan));
    ticksTemp = 1:(length(speedNoNan)-1)/N_:length(speedNoNan);
    for i = 1:length(ticksTemp)-1
        ticks(i,1:2) = [ceil(ticksTemp(i)), floor(ticksTemp(i+1))];
    end
    if ~isempty(Vrange_)   % overload
        clear ticks
        for i = 1:size(Vrange_,1)
            [~, ind] = find(speedNoNan >= Vrange_(i,1));
            ticks(i,1) = min(ind);
            [~, ind] = find(speedNoNan <= Vrange_(i,2));
            ticks(i,2) = max(ind);
        end
    end
    for i = 1:size(ticks,1)
        speedRange(i,1:2) = [speedNoNan(ticks(i,1)), speedNoNan(ticks(i,2))];
        meanSpeed(i) = mean(speedNoNan(ticks(i,1):ticks(i,2)));
    end
    
    
    %% classify speed into ranges:
    DOFstr = '';
    for i = 1:size(speedRange,1)
        ind = find(speed >= speedRange(i,1) & speed <= speedRange(i,2));
        dataBySpeed{i} = dataSpeedSorted(ind);
        DOFstr = [DOFstr, num2str(length(ind)), ', '];
        DOF(i) = length(ind);
    end
    DOFstr(end) = [];
    
    return;
end