
function obj = writeAccuCoeff(obj, name, accuCoeffValue)

found = 0;
[nameShort, axis] = nameID(obj);

accuCoeff = obj.read("accuCoeff");
if class(accuCoeff) ~= "class_data"
    if accuCoeff == "nonexist"
        clear accuCoeff;
    %     accuCoeff.name = strings(0);
    %     accuCoeff.coeff = {};
        accuCoeff = class_data;    % change accuCoeff to class_data obj.
    else
        error("Yu: unexpected");
    end
end

for i = 1:length(nameShort)
    if strcmp(name, nameShort(i)) || strcmp(name, "all")
%         ind = find(accuCoeff.name == nameShort(i));
%         if isempty(ind)
%             accuCoeff.name(end+1) = nameShort(i);
%             accuCoeff.coeff{end+1} = accuCoeffValue;
%             found = found+1;
%         elseif length(ind) == 1
%             accuCoeff.coeff{ind} = accuCoeffValue;
%             found = found+1;
%         else
%             error("Yu: unexpected");
%         end
        found = found+1;
        accuCoeff = accuCoeff.write(nameShort(i), accuCoeffValue, 'value');
        
    end
end

obj = obj.write("accuCoeff", accuCoeff,'table');

if found == 0
    error("Yu: name nonexist");
end
        

end

