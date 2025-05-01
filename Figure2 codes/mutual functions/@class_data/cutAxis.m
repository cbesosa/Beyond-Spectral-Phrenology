% range: [left1, right1; left2, right2; ....]
       
function [obj, objArray, segL] = cutAxis(obj, name, range)

if numel(range) == 2
    range = range(:)';
end
if numel(range) == 1
    range = [range,range];
end
if numel(range) == 0
    objArray = obj.interpAxis(name, []);
    obj = objArray;
end

oldAxis = obj.read(name);

ind = zeros(size(oldAxis));
segL = nan(1, size(range, 1));
for i = 1:size(range, 1)
   indArray{i} = oldAxis >= range(i,1) & oldAxis <= range(i,2);
   newAxisArray{i} = oldAxis(indArray{i});
   if isempty(newAxisArray{i})
       segL(i) = 0;
   else
       if class(newAxisArray{i}) == "string"
           segL(i) = nan;
       else
           segL(i) = newAxisArray{i}(end) - newAxisArray{i}(1);
       end
   end
   objArray(i) = obj.interpAxis(name, newAxisArray{i});
   
   ind = ind | indArray{i};
end
newAxis = oldAxis(ind);
obj = obj.interpAxis(name, newAxis);

end
        

    
    


