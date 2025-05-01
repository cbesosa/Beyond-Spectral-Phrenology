
       
function obj = permuteAxis(obj, dataName, axis1, axis2_)

if nargin == 3
    obj = moveAxisRightWard(obj, dataName, axis1);
    return;
end

if strcmp(axis1, axis2_)
    return;
end

[nameShort, axes, nameFull] = obj.nameID(dataName); 
data = obj.read(dataName);

ind1 = find(axes == string(axis1));
ind2 = find(axes == string(axis2_));

% resolve the same axis name problem:
ind1 = ind1(1);


order = 1:length(axes);
order(ind1) = ind2;
order(ind2) = ind1;
data = permute(data, order);

nameFullNew = [char(nameShort),'('];
for i = 1:length(axes)
    if i == ind1
        nameFullNew = [nameFullNew, char(axes(ind2)),','];
    elseif i == ind2
        nameFullNew = [nameFullNew, char(axes(ind1)),','];
    else
        nameFullNew = [nameFullNew, char(axes(i)),','];
    end
end
nameFullNew(end) = ')';

obj = obj.write(nameFullNew, data, 'value');
    
end
        

function obj = moveAxisRightWard(obj, dataName, axis1)

[nameShort, axes, nameFull] = obj.nameID(dataName); 
data = obj.read(dataName);

ind1 = find(axes == string(axis1));
if length(axes) == ind1
    return; 
end
ind2 = ind1 +1;

order = 1:length(axes);
order(ind1) = ind2;
order(ind2) = ind1;
data = permute(data, order);

nameFullNew = [char(nameShort),'('];
for i = 1:length(axes)
    if i == ind1
        nameFullNew = [nameFullNew, char(axes(ind2)),','];
    elseif i == ind2
        nameFullNew = [nameFullNew, char(axes(ind1)),','];
    else
        nameFullNew = [nameFullNew, char(axes(i)),','];
    end
end
nameFullNew(end) = ')';

obj = obj.write(nameFullNew, data, 'value');

end
    


