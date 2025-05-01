
function obj = sortAxis(obj, axisName, direction)


% find the data:
[nameShort, axes, Dim] = axis2data(obj, axisName);

% construct the axis:
axis = obj.read(axisName);
[axis, ind] = sort(axis, direction);
obj = obj.write(axisName, axis, 'value');

% construct the data:
for i = 1:length(nameShort)
    obj = obj.permuteAxis(nameShort(i), axisName, axes{i}(1));
    data = obj.read(nameShort(i));
    if ~isempty(data)
        data(:,:) = data(ind,:);
        obj = obj.write(nameShort(i), data, 'value');
        obj = obj.permuteAxis(nameShort(i), axisName, axes{i}(1));
    end
end


end





