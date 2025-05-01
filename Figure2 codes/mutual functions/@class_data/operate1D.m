
       
function obj2 = operate1D(obj1, valueNames, principalAxis, operator)

% valueNames = <string array>;
% axisName = <string>;
% operator = <function handle>;

% function [name, vOut] = operator(vIn, axis) 
% % name = <string array>, could include axis;  to distinguish value from axis, add () at end for singular value;
% % vOut = <cell array>
% % vIn = <cell array> of values, cause the results may depend on multiple values;   
% % axis = <array> same as axis values
% end

    if class(principalAxis) == "char"
        principalAxis = string(principalAxis);
    end
    if class(valueNames) == "char"
        valueNames = string(valueNames);
    end
    if class(operator) ~= "function_handle"
        error("Yu: ");
    end

    for i = 1:length(valueNames)
        [nameShort(i), axes{i}, nameFull(i)] = obj1.nameID(valueNames(i));
    end
    % check if axes are the same: (incapable of reorderred axis)
    for i = 2:length(valueNames)
        if length(axes{i}) ~= length(axes{i-1})
            error("Yu: axes dimension is not identical");
        end
        for j = 1:length(axes{i})
            if ~strcmp(axes{i}(j), axes{i-1}(j)) 
                error("Yu: axis names are not identical");
            end
        end
    end
    axes = axes{1};
    
    % construct a simplified class_data obj that is directly related to the necessary information: 
    data = class_data;
    for i = 1:length(axes)
        data = data.write(axes(i), obj1.read(axes(i)), 'value');
    end
    for i = 1:length(valueNames)
        data = data.write(nameFull(i), obj1.read(nameShort(i)), 'value');
    end

    
    % put principalAxis to be the first Axis:
    for i = 1:length(valueNames)
        data = data.permuteAxis(nameShort(i), axes(1), principalAxis); 
    end        
    [~, axesOri, ~] = data.nameID(nameShort(1));              % refresh the order of the axis.       
    % need here to make sure axes are in the same order.
    
    
    % find the input to operator:
    pricipalAxisValue = data.read(principalAxis);
    for i = 1:length(valueNames)
        vInAll{i} = data.read(nameShort(i));
    end   
    
    % find the size of the output:
    for i = 1:length(valueNames)
        vIn{i} = vInAll{i}(:, 1);
    end  
    [nOut, vOut] = operator(vIn, pricipalAxisValue);
    temp = class_data;
    for i = 1:length(nOut)
        temp = temp.write(nOut(i), vOut{i}, 'value');
    end
    [nameShort, axes, nameFull] = temp.nameID; 
    [type, SZ] = temp.dimension;     % SZ = size to distinguish.
    for i = 1:length(type)
        if type(i) == "value"
            SOD = size(vInAll{1});   % size original data, need to further ensure that vInAll{:} have compatible sizes.
            if SOD(2) == 1 && length(axes{i}) == 1
                SOD(2) = [];
            end
            if SOD(1) == 1 && length(axes{i}) == 1 && length(SOD) == 2
                SOD(1) = [];
            end
            if length(axes{i}) == 0 && length(SZ{i}) == 1
                if length(SOD) <= 2
                    SOD(1) = [];  SR{i} = [SZ{i}, SOD];
                else
                    SOD(1) = [];  SR{i} = [SOD];
                end
                if SZ{i} ~= 1
                    error("Yu: dimension and size are incompatible")
                end
            else
                SOD(1) = [];  SR{i} = [SZ{i}, SOD];    % size of result
            end
            
            if length(SR{i}) == 1
                vOutAll{i} = nan([SR{i},1]);
            else
                vOutAll{i} = nan(SR{i});
            end
        end
        1;
    end
        
    % perform the operator:
    for i = 1:length(type)
        if type(i) == "value"
            assignDim{i} = '';
            if length(axes{i}) == 0 && length(SZ{i}) == 1
                if SZ{i} ~= 1
                    error("Yu: dimension and size are incompatible")
                end
            else
                for k = 1:length(SZ{i})
                    assignDim{i} = [assignDim{i},':,'];
                end
            end
        end
    end
    
    % for j = 1:length(vInAll{1}(1, :))
    % for j = 1:size(vInAll{1}, 2)
    for j = 1:numel(vOutAll{i})
        for i = 1:length(valueNames)
            vIn{i} = vInAll{i}(:, j);
        end  
        [nOut, vOut] = operator(vIn, pricipalAxisValue);
%         if isnan(real(vOut{1}))
%             1;
%         end
        for i = 1:length(type)
            if type(i) == "value"
                eval(['vOutAll{i}(',assignDim{i},'j) = vOut{i};']);
            end
        end
    end
  
    
    % write into a data obj:
    % obj2 = class_data;
    obj2 = obj1;
    for i = 1:length(type)
        if type(i) == "value"
            name = char(nOut(i));  name(end) = ',';
            if length(axes{i}) == 0 && length(SZ{i}) == 1   % scalar value exception.
                name(end) = [];
            end
            for j = 2:length(axesOri)
                obj2 = obj2.write(axesOri(j), data.read(axesOri(j)), 'value');
                name = [name, char(axesOri(j)), ','];
            end
            if name(end) == ','
                name(end) = ')';
            else
                name(end+1) = ')';
            end
            obj2 = obj2.write(name, vOutAll{i}, 'value');
        elseif type(i) == "axis"
            obj2 = obj2.write(nOut(i), vOut{i}, 'value');
        end
    end
   
            
    return;
end
        

    
    


