% range: [left1, right1]
% segments are by default isometric in axis.
% marker can be N-D, & operation performs along non-principalAxis
       
function [obj, objArray] = segments(obj, principalAxis, minLength_, markerName_, range_)

if nargin == 2
    minLength_ = 0;
end

if nargin == 2 || nargin == 3
    [obj, objArray] = isoSegments(obj, principalAxis, minLength_);
    return;
end

if numel(range_) == 1
    range_ = [range_, range_];
end
    
if numel(range_) == 2
    range_ = range_(:)';
else
    error("Yu: unexpected")
end

oldAxis = obj.read(principalAxis);
marker = obj.read(markerName_);
[nameShort, axes, nameFull] = obj.nameID(markerName_);
ind = find(axes == principalAxis);
order = 1:length(axes);
order(1) = ind; order(ind) = 1;
if length(order) == 1
    marker = marker(:); 
    order = [1,2];
end
marker = permute(marker, order);
marker = marker(:,:);
valid = marker >= range_(1) & marker <= range_(2);
valid = all(valid, 2);
if numel(valid) ~= length(oldAxis)
    error("Yu: unexpected")
end

[ID, ~] = axisID(obj, principalAxis);
diffAxis = oldAxis(2:end) - oldAxis(1:end-1);
diffAxis = diffAxis(:);
diffBool = diffAxis > (ID.mode-0.001*abs(ID.mode)) & diffAxis < (ID.mode+0.001*abs(ID.mode));
% diffValid = valid(2:end) & valid(1:end-1) & diffBool;
% diffValid = [0,diffValid(:)',0];
% indLeft = find(diffValid(1:end-2) == 0 & diffValid(2:end-1) == 1);
% indRight = find(diffValid(1:end-1) == 1 & diffValid(2:end) == 0);
diffBool = [0,diffBool(:)',0];
valid = [0,valid(:)',0];
indLeft = find(valid(2:end-1) == 1 & (valid(1:end-2) == 0 | diffBool(1:end-1) == 0));
indRight = find(valid(2:end-1) == 1 & (valid(3:end) == 0 | diffBool(2:end) == 0));
if length(indLeft) ~= length(indRight)
    error("Yu: unexpected");
end

% eliminated the short segments:
indShort = [];
for i = 1:length(indLeft)
    LS = indRight(i) - indLeft(i); 
    LS = LS*abs(ID.mode);
    if LS < minLength_
        indShort = [indShort, i];
    end
end
indLeft(indShort) = []; indRight(indShort) = []; 


% construct new objs:
[obj, objArray] = constructNewObjs(obj, principalAxis, indLeft, indRight);

end

% get isolated segments:
function [obj, objArray] = isoSegments(obj, principalAxis, minLength_)
        
if nargin == 2
    minLength_ = 0;
end

% find isometric segments:
oldAxis = obj.read(principalAxis);
[AID, ~] = obj.axisID(principalAxis);

diffAxis = oldAxis(2:end) - oldAxis(1:end-1);
diffAxis = diffAxis(:);
diffBool = diffAxis > (AID.mode-0.001*abs(AID.mode)) & diffAxis < (AID.mode+0.001*abs(AID.mode));
diffBool = [0,diffBool(:)',0];
indLeft = find(diffBool(1:end-2) == 0 & diffBool(2:end-1) == 1);
indRight = find(diffBool(1:end-1) == 1 & diffBool(2:end) == 0);
if length(indLeft) ~= length(indRight)
    error("Yu: unexpected");
end

% eliminated the short segments:
indShort = [];
for i = 1:length(indLeft)
    LS = indRight(i) - indLeft(i); 
    LS = LS*abs(AID.mode);
    if LS < minLength_
        indShort = [indShort, i];
    end
end
indLeft(indShort) = []; indRight(indShort) = []; 


% construct objArray:
[obj, objArray] = constructNewObjs(obj, principalAxis, indLeft, indRight);
    
end


function [obj, objArray] = constructNewObjs(obj, principalAxis, indLeft, indRight)

oldAxis = obj.read(principalAxis);
if isempty(oldAxis)
    objArray = obj;
    return;
end

% construct objArray:
name = obj.name;
[nameShort, axes, Dim, indIV] = obj.axis2data(principalAxis); 
[~, ~, ~, indAxisIV] = obj.nameID(principalAxis); 
for i = 1:length(nameShort)   % permute the pricipalAxis to be the 1st axis.
    obj = obj.permuteAxis(nameShort(i), axes{i}(1), principalAxis);
    PermOrder{i} = 1:length(axes{i});
    PermOrder{i}(Dim(i)) = 1;
    PermOrder{i}(1) = Dim(i);
end

countPoints = 0;
indAllSegs = nan(1, sum(indRight-indLeft+1));
if isempty(indAllSegs)
    objArray = obj.cutAxis(principalAxis, []);
    obj = objArray;
    return;
end
[~, ind4sample] = min(indRight-indLeft);
objArraySample = obj.cutAxis(principalAxis, oldAxis([indLeft(ind4sample), indRight(ind4sample)]));    
objArraySample.name = name;
objArray(1:length(indLeft)) = objArraySample;

for i = 1:length(indLeft)
    count = fprintf(['class_data -> segments: ', num2str(i), '/', num2str(length(indLeft))]);

    indSegsArray = indLeft(i):indRight(i);
    indAllSegs(countPoints+1:countPoints+indRight(i)-indLeft(i)+1) = indSegsArray;
    countPoints = countPoints+indRight(i)-indLeft(i)+1;
    objArray(i).value{indAxisIV} = oldAxis(indSegsArray);
    for j = 1:length(nameShort)
        if isempty(obj.value{indIV(j)})
            objArray(i).value{indIV(j)} = [];
        else
            objArray(i).value{indIV(j)} = obj.value{indIV(j)}(indLeft(i):indRight(i),:);
        end
        % permute the pricipalAxis back to the original Dim:
        if length(axes{j}) > 1
            objArray(i).value{indIV(j)} = permute(objArray(i).value{indIV(j)}, PermOrder{j});
            1;
        end
    end
    fprintf(1, repmat('\b',1,count));
end


newAxis = oldAxis(indAllSegs);    
obj = obj.interpAxis(principalAxis, newAxis);
for i = 1:length(nameShort)   % permute back.
    obj = obj.permuteAxis(nameShort(i), axes{i}(1), principalAxis);
end


end










