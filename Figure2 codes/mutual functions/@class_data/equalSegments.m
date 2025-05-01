% range: [left1, right1]
% segments are by default isometric in axis.
% marker can be N-D, & operation performs along non-principalAxis
       
function [objArray, actualSegL] = equalSegments(obj, principalAxis, SegL, overlap)



% find isometric segments:
oldAxis = obj.read(principalAxis);
[AID, ~] = obj.axisID(principalAxis);
name = obj.name;

diffAxis = oldAxis(2:end) - oldAxis(1:end-1);
diffAxis = diffAxis(:);
diffBool = diffAxis > (AID.mode-0.0001*abs(AID.mode)) & diffAxis < (AID.mode+0.0001*abs(AID.mode));
diffBool = [0,diffBool(:)',0];
indLeft = find(diffBool(1:end-2) == 0 & diffBool(2:end-1) == 1);
indRight = find(diffBool(1:end-1) == 1 & diffBool(2:end) == 0);
if length(indLeft) ~= length(indRight)
    error("Yu: unexpected");
end

%  devide into SegL:
indLeftShort = [];
SegN = round(SegL/(abs(AID.mode)));
step = round(SegN *(1-overlap));
for i = 1:length(indLeft)
    for j = indLeft(i):step:indRight(i)-SegN+1
        indLeftShort = [indLeftShort, j];
    end
end
indLeft = indLeftShort; indRight = indLeftShort +SegN -1; 


% construct objArray:
[nameShort, axes, Dim, indIV] = obj.axis2data(principalAxis); 
[~, ~, ~, indAxisIV] = obj.nameID(principalAxis); 
for i = 1:length(nameShort)   % permute the pricipalAxis to be the 1st axis.
    obj = obj.permuteAxis(nameShort(i), axes{i}(1), principalAxis);
    PermOrder{i} = 1:length(axes{i});
    PermOrder{i}(Dim(i)) = 1;
    PermOrder{i}(1) = Dim(i);
end


objArraySample = obj.cutAxis(principalAxis, oldAxis([indLeft(1), indRight(1)]));    
objArraySample.name = name;
objArray(1:length(indLeft)) = objArraySample;
for i = 1:length(indLeft)
    count = fprintf(['class_data -> equalSegments: ', num2str(i), '/', num2str(length(indLeft))]);
    indSegsArray = indLeft(i):indRight(i);
    objArray(i).value{indAxisIV} = oldAxis(indSegsArray);
    for j = 1:length(nameShort)
        if ~isempty(objArray(i).value{indIV(j)})
            objArray(i).value{indIV(j)} = obj.value{indIV(j)}(indLeft(i):indRight(i),:);
        end
        % permute the pricipalAxis back to the original Dim:
        if length(axes{j}) > 1
            objArray(i).value{indIV(j)} = permute(objArray(i).value{indIV(j)}, PermOrder{j});
        end
    end
    fprintf(1, repmat('\b',1,count));
end

actualSegL = SegN *abs(AID.mode);

end
        

    
    


