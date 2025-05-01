
function obj = concatenate(objArray, axisName)
% when facing repeated axis ticks, the code delete redundant axis-data paris without making sure the data associated with repeated axis ticks are identical.  
% objArray(1) is the template, with mergeTable function.

% find the data:
[nameShort, axes, Dim] = axis2data(objArray(1), axisName);


% find anchor index of each segment:
segSZ = nan(1, length(objArray));
indLeft = nan(1, length(objArray));  indRight = nan(1, length(objArray));

for i = 1:length(objArray)
    [type, segSZ(i)] = objArray(i).dimension(axisName);
    if i == 1
        indLeft(1) = 1;  indRight(1) = segSZ(i);
    else
        indLeft(i) = indRight(i-1) + 1;
        indRight(i) = indLeft(i) + segSZ(i) -1;
    end
end
tolSZ = sum(segSZ);


% construct the axis
if class(objArray(1).read(axisName)) == "double"
    axisValue = nan(1,tolSZ);
elseif class(objArray(1).read(axisName)) == "string"
    axisValue = repmat("", [1,tolSZ]);
else
    temp = objArray(1).read(axisName);
    axisValue = repmat(temp(1), [1,tolSZ]);
end
for i = 1:length(objArray)
    axisValue(indLeft(i):indRight(i)) = objArray(i).read(axisName);
end


% construct the data:
data = cell(1,length(nameShort));
for i = 1:length(nameShort)
    countf = fprintf(['class_data -> concatenate: ', num2str(i), '/', num2str(length(nameShort)), '(', char(nameShort(i)),')']);
    tempObj = objArray(1).permuteAxis(nameShort(i), axes{i}(1), axisName);
    sizeData = size(tempObj.read(nameShort(i)));
    sizeData(1) = indRight(end);
    if class(tempObj.read(nameShort(i))) == "double"
        if issparse(tempObj.read(nameShort(i)))
            data{i} = sparse(sizeData(1), sizeData(2));   
        else
            data{i} = nan(sizeData);    % nan default is slower than 0 default.
        end
    elseif class(tempObj.read(nameShort(i))) == "cell"
        data{i} = cell(sizeData);
    else
        data{i}(prod(sizeData)) = repmat(tempObj.read(nameShort(i)), sizeData); % is this right?
    end
    
    for j = 1:length(objArray)
        countf2 = fprintf(['; ', num2str(j), '/', num2str(length(objArray))]);
        tempObj = objArray(j).permuteAxis(nameShort(i), axes{i}(1), axisName);
        tempData = tempObj.read(nameShort(i));
        if ~isempty(tempData)
            data{i}(indLeft(j):indRight(j),:) = tempData((indLeft(j):indRight(j))-indLeft(j) +1, :);
        end
        fprintf(1, repmat('\b',1,countf2));
    end
    fprintf(1, repmat('\b',1,countf));
end


% write into obj:
obj = objArray(1);
table  = mergeTable(objArray);
obj.table = table;

obj = obj.write(axisName, axisValue, 'value');
for i = 1:length(nameShort)
    obj = obj.permuteAxis(nameShort(i), axes{i}(1), axisName);
    obj = obj.write(nameShort(i), data{i}, 'value');
end

% sortAxis into ascending order:
if class(obj.read(axisName)) == "double"
    obj = obj.sortAxis(axisName, 'ascend');
end

% delete repeated axis-data paris:
if class(obj.read(axisName)) == "double"
    axisValue = obj.read(axisName);
    axisValueNew = unique(axisValue);
    obj = obj.interpAxis(axisName, axisValueNew);
end


% permute back into original order:
for i = 1:length(nameShort)
    obj = obj.permuteAxis(nameShort(i), axes{i}(1), axisName);
end


end


function table = mergeTable(objArray)
% for table entries in double, concatenate will calculate the mean.
    table = objArray(1).table;

    for i = 1:size(table, 1)
        if isnumeric(table{i, 2})
            for j = 1:length(objArray)
                table{i, 2} = table{i, 2}*(j-1)/j + objArray(j).read(table{i, 1})/j;
            end
        else
            for j = 1:length(objArray)
                if ~isequal(table{i, 2}, objArray(j).read(table{i, 1}))
                    table{i, 2} = "";
                end
            end
        end
    end

end



















