
        
function obj = mergeAxis(obj, axisGone, axisLeft)

% interp the data:
axisLeftValue = obj.read(axisLeft);
obj = obj.interpAxis(axisGone, axisLeftValue);

% update the accuCoeff to be max of axisGone and axisLeft
if class((obj.readAccuCoeff(axisGone))) ~= "string" && ...
   class((obj.readAccuCoeff(axisLeft))) ~= "string"   
        obj = obj.writeAccuCoeff(axisLeft, ...
              max(obj.readAccuCoeff(axisGone), obj.readAccuCoeff(axisLeft)));
else
    obj = obj.writeAccuCoeff(axisLeft, nan);        %  this is not the perfect way, change it later.
end

% find the affected data and replace names:
[nameShort, axisName] = obj.nameID;
Dim = 0*(1:length(nameShort));  % 0: the value entry does not contain the axisGone (exclude the axis itself).
for i = 1:length(nameShort)    
    newName(i) = string([char(nameShort(i))]); % construct the new name.
    if isempty(axisName{i})
        if strcmp(nameShort(i), axisGone)
            ind_axisGone = i;  
        end
    else
        newName(i) = string([char(newName(i)),'(']);
        for j = 1:length(axisName{i}) 
            if strcmp(axisName{i}(j), axisGone)
                Dim(i) = j;     % j: the data entry contains the axis as jth dimension.
                newName(i) = string([char(newName(i)), char(axisLeft),',']);
            else
                newName(i) = string([char(newName(i)), char(axisName{i}(j)),',']);
            end
        end
        temp = newName(i); temp = char(temp); 
        temp(end) = []; temp(end+1) = ')'; 
        newName(i) = string(temp);
    end 
end

%  write back into obj:
for i = 1:length(nameShort)    
    if Dim(i) > 0
        obj.name(i) = newName(i);
    end
end
obj = obj.remove(axisGone);


end
        

    
    


