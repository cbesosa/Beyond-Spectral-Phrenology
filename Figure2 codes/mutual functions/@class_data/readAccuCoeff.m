
function accuCoeff = readAccuCoeff(obj, name)

[nameShort, axis] = nameID(obj);

accuCoeffTable = obj.read("accuCoeff");
if class(accuCoeffTable) ~= "class_data" 
    if accuCoeffTable == "nonexist"
        accuCoeff = "nonexist";
        return;
    end
end

if strcmp(name, "all")
    name = nameShort;
end
if ischar(name)
    name = string(name);
end

for i = 1:length(name)
    foundName = find(nameShort == name(i));   % check if input name is correct.
    if isempty(foundName)
        error("Yu: name nonexist");
    end
    
%     ind = find(accuCoeffTable.name == name(i));
%     if length(ind) == 1
%         accuCoeff{i} = accuCoeffTable.coeff{ind};
%     elseif isempty(ind)
%         accuCoeff{i} = "nonexist";
%     else
%         error("Yu: unexpected");
%     end
    accuCoeff{i} = accuCoeffTable.read(name(i));   % change accuCoeff to be class_data obj.

end
        
if length(accuCoeff) == 1
    accuCoeff = accuCoeff{1};
end

end