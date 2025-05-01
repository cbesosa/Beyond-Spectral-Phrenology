


function [ID, reAxis] = axisID(obj, name_)

% axis identity
% analyze and show statistics of the axis values
% reAxis is the corresponding isometric axis.

[nameShort, axisName] = obj.nameID;

for i = 1:length(nameShort)
    if isempty(axisName{i})
        
        axis = obj.read(nameShort(i));
        [ID(i).type, ID(i).mode] = classifySequence(axis);
        ID(i).length = estimateSequenceLength(axis, ID(i).type, ID(i).mode);
        
        if strcmp(ID(i).type, "isometric") 
            reAxis(i).mode = ID(i).mode;                         % grid size.
            reAxis(i).ticks = axis(1):reAxis(i).mode:axis(end);
            
        elseif strcmp(ID(i).type, "isometricWithExceptions") 
            Nintervals = ceil((axis(end)-axis(1))/ID(i).mode);
            reAxis(i).mode = (axis(end)-axis(1))/Nintervals;
            reAxis(i).ticks = axis(1):reAxis(i).mode:axis(end);
            
        elseif strcmp(ID(i).type, "discrete") 
            reAxis(i).mode = nan;
            reAxis(i).ticks = axis;
            
        elseif strcmp(ID(i).type, "single") 
            reAxis(i).mode = nan;
            reAxis(i).ticks = axis;
        else
            reAxis(i).mode = nan;
            reAxis(i).ticks = axis;            
        end
        
        % reAxis(i).accuCoeff = class_data.calculateAccuCoeff(axis, axis*0, reAxis(i).ticks);   % accuracy coefficient. 1st order data.
        % reAxis(i).accuCoeffMin = 1/2*sqrt(abs(ID(i).mode));  
        % the natural minimum accuracy coefficient between data points.

    end
end

if nargin == 1
    return;
elseif nargin == 2
    ind = find(nameShort == name_);
    ID = ID(ind);
    reAxis = reAxis(ind);
end
    
end


function [type, axisMode] = classifySequence(S)
% statistics of the difference of a sequence.

if isnumeric(S)
else
    type = "non-numeric";
    axisMode = nan;
    return
end
if isempty(S)
    type = "empty";
    axisMode = nan;
    return
end
if length(S) == 1
    type = "single";
    axisMode = nan;
    return
end

diff = S(2:end) -S(1:end-1);
obj = class_PDF(diff);

uniqueDataCounts = obj.read('uniqueDataCounts');
% merge almost same values:
i = 1;
while i < size(uniqueDataCounts ,2)
    if 2*abs(uniqueDataCounts(1,i) - uniqueDataCounts(1,i+1))/abs(uniqueDataCounts(1,i) + uniqueDataCounts(1,i+1)) < 0.0001
        uniqueDataCounts(1,i) = (uniqueDataCounts(1,i)*uniqueDataCounts(2,i) + uniqueDataCounts(1,i+1)*uniqueDataCounts(2,i+1))/...
                                (uniqueDataCounts(2,i) + uniqueDataCounts(2,i+1));
        uniqueDataCounts(2,i) = (uniqueDataCounts(2,i) + uniqueDataCounts(2,i+1));
        uniqueDataCounts(:,i+1) = [];
    else
        i = i+1;
    end
end

if size(uniqueDataCounts ,2) == 1
    type = "isometric";
    axisMode = uniqueDataCounts(1,1);
elseif uniqueDataCounts(2,1)/length(diff) >= 0.75
    type = "isometricWithExceptions";
    axisMode = uniqueDataCounts(1,1);
elseif size(uniqueDataCounts ,2)/length(diff) >= 0.5
    type = "discrete";
    axisMode = uniqueDataCounts(1,1);
else
    type = "irregular";
    axisMode = mode(diff);
    % error('Yu: Let me add more types.');
end


    1;
end


function L = estimateSequenceLength(S, type, mode)

    if strcmp(type, "isometric") 
        L = mode * length(S);
            
    elseif strcmp(type, "isometricWithExceptions") 
        L = mode * length(S);

    elseif strcmp(type, "discrete") 
        L = nan;
        
    else
        L = nan;
    end

end






