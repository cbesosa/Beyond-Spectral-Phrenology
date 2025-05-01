
       
function obj = removeDimension(obj, valueName, axisName)

% add an extra dimension to a value, axisValue is a scalar.
    
    % read old axis:
    oldAxis = obj.read("axisName");
    if length(oldAxis) ~= 1
        error("Yu: only to remove singular dimension");
    end
    
    % read old data:
    oldValue = obj.read(valueName);
    [oldNameShort, oldAxes, nameFull] = obj.nameID(valueName);
    [~, Dim] = find(oldAxes == axisName);
    
    % find the size of new data:
    newAxes = oldAxes;
    newAxes(Dim) = [];
    [type, SD] = obj.dimension(valueName); 
    SD(Dim) = [];
    if length(SD) > 1
        if class(oldValue) == "cell"
            newValue = cell(SD);
        else
            newValue = nan(SD);
        end
    elseif length(SD) == 1
        if class(oldValue) == "cell"
            newValue = cell([SD,1]);
        else
            newValue = nan([SD,1]);
        end
    else
        error("Yu: please improve");
    end
    
    % create the new data:
    newValue(:) = oldValue(:); 
    
    % write back:
    newName = class_data.nameIDstr(oldNameShort, newAxes);
    obj = obj.write(newName, newValue);


end
        

    
    


