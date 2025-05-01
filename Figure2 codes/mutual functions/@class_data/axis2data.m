
% find the data related to the axisName:

function [nameShort, axes, Dim, indInValue] = axis2data(obj, axisName)

nameShort = string(0);
axes = {};
Dim = []; 
nFound = 0;

% find the data:
[nameShortAll, axesAll, ~, ind] = obj.nameID;
for i = 1:length(nameShortAll)
    if isempty(axesAll{i})
        if strcmp(nameShortAll(i), axisName)
            axisInd = i;  
            axisValue = obj.read(nameShortAll(i));
        end
    else
        for j = 1:length(axesAll{i}) 
            if strcmp(axesAll{i}(j), axisName)
                
                nFound = nFound +1;
                nameShort(nFound) = nameShortAll(i);
                axes{nFound} = axesAll{i};
                Dim(nFound) = j;     % j: the data entry contains the axis as jth dimension.
                indInValue(nFound) = ind(i);
            end
        end
    end 
end


end







