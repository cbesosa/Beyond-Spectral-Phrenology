
       
function obj = addDimension(obj, valueName, axisName, axisValue)

% add an extra dimension to a value, axisValue is a scalar.

    % create the new axis (need to be improved later):
    axisNamePre = obj.read(axisName);
    existFlag = 1;
    if class(axisNamePre) == "string"
        if axisNamePre == "nonexist"
            obj = obj.write(axisName, axisValue);
            existFlag = 0;
        end
    end
    if existFlag == 1
        axisValuePre = obj.read(axisName);
        if all(axisValuePre == axisValue) && length(axisValuePre) == 1
        else
            error("Yu: axisValue not compatible")
        end
    end
    
    % read old data:
    oldValue = obj.read(valueName);
    [oldNameShort, oldAxes, nameFull] = obj.nameID(valueName);
    
    % create new data:
    newValue = nan([1, size(oldValue)]);
    newValue(:) = oldValue(:);
    newName = [char(oldNameShort), '(', char(axisName)];
    for k = 1:length(oldAxes)
        newName = [newName, ',',char(oldAxes(k))];
    end
    newName = [newName, ')'];
    obj = obj.write(newName, newValue);


end
        

    
    


