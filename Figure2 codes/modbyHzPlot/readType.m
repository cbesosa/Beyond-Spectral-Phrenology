function [typeValue, typeValueStr] = readType(data, typeName_)
% if typeName

    stateNameAll = data.read('typeName');
    stateValueAll = data.read('cellType');

    if nargin < 2
        for i = 1:length(stateNameAll)
            typeValue.(stateNameAll(i)) = [stateValueAll{i, :}];
        end
        for j = 1:size(stateValueAll, 2)
            typeValueStr(j,1) = "";
            for i = 1:length(stateNameAll)
                if rem(i-1,4) == 0 && i>1   % max entrys per line.
                    typeValueStr(j,1) = string([char(typeValueStr(j,1)), 10]);
                end
                typeValueStr(j,1) = string([char(typeValueStr(j,1)), char(stateNameAll(i)), '= ', num2str(typeValue.(stateNameAll(i))(j), 3), '; ']);
            end
        end
        1;
    else
        ind = find(stateNameAll == string(typeName_));
        typeValue = [stateValueAll{ind, :}];
        for i = 1:length(typeValue)
            typeValueStr(i,1) = string([char(typeName_), '= ', num2str(typeValue(i), 3)]);
        end
        1;   
    end
    
end