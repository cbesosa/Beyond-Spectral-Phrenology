% data identity
% show statistics and find the type of data.


function [ID, CF] = dataID(obj, dataName, axisName)

AID = obj.axisID(axisName);
if AID.type == "isometric" || AID.type == "isometricWithExceptions" 
else
    error("Yu: axisID non-suitable for dataID");
end
axisValue{0+1} = obj.read(axisName); 
data = obj.read(dataName); 
SD = size(data);
[nameShort, axes] = nameID(obj, dataName); 

% permute aimmed axis to be the 1st:
axesOrder = 1:length(axes);
ind = find(axes == axisName);
axesOrder(1) = ind;
axesOrder(ind) = 1;
if length(axesOrder) == 1; axesOrder = [1,2]; end
data = permute(data, axesOrder);
SD = SD(axesOrder);
axes = axes(axesOrder);


maxOrder = 3;
%  construct diff0(original), diff1 and diff2 (or maybe diff3) data:
for order = 0:maxOrder
    if order == 0
        diffData{order+1} = data;
        diff1AxisValidity{order+1} = nan;
    else
        if order == 1
            diff1Axis = axisValue{0+1}(2:end) - axisValue{0+1}(1:end-1);
            ind = find(diff1Axis <= AID.mode*1.0001 & diff1Axis >= AID.mode*0.9999);
            diffAxisValidity{order+1} = diff1Axis*0;
            diffAxisValidity{order+1}(ind) = 1;  % use 0, 1 to denote effectiveness of diffAxis. 
        else
            diffAxisValidity{order+1} = diffAxisValidity{order-1+1}(1:end-1) & diffAxisValidity{order-1+1}(2:end);
        end
        diffData{order+1} = diffData{order-1+1}(2:end, :) - diffData{order-1+1}(1:end-1, :); 
        axisValue{order+1} = (axisValue{order-1+1}(2:end) + axisValue{order-1+1}(1:end-1))/2;
        1;
    end
end


% construct scale analysis of diff1 and diff2 (or maybe diff3) data:
for order = 0:maxOrder
    N = 1:20;   % scales to be analyzed.
    scaleVar_Diff{order+1} = nan(length(N), SD(2:end));
    for i = 1:length(N)
        diff_diffData = diffData{order+1}(N(i)+1:end, :) - diffData{order+1}(1:end-N(i), :);
        if order+i >= 2
            diffAxisValidity{order+i+1} = diffAxisValidity{order+i-1+1}(1:end-1) & diffAxisValidity{order+i-1+1}(2:end);
        end
        diff_diffData = diff_diffData(diffAxisValidity{order+i+1}, :);
        scaleVar_Diff{order+1}(i, :) = var(diff_diffData, 'omitnan');
        scale(i) = AID.mode*N(i);
    end
    % scaleVar_Diff{order+1} = permute(scaleVar_Diff{order+1}, axesOrder);
end


% write the scale vs var result into a class_data obj:
scaleVar = class_data;
for order = 0:maxOrder
    valueName1{order+1} = ['diffData', num2str(order), '('];
    valueName2{order+1} = ['scaleVar_Diff', num2str(order), '('];
    for i = 1:length(axes)
        if axes(i) == axisName
            scaleVar = scaleVar.write(['originalAxis',num2str(order)], axisValue{order+1}, 'value');
            valueName1{order+1} = [valueName1{order+1}, 'originalAxis',num2str(order),', '];
            scaleVar = scaleVar.write('scale', scale, 'value');
            valueName2{order+1} = [valueName2{order+1}, 'scale, '];
        else
            scaleVar = scaleVar.write(axes(i), obj.read(axes(i)), 'value');
            valueName1{order+1} = [valueName1{order+1}, char(axes(i)), ', '];
            valueName2{order+1} = [valueName2{order+1}, char(axes(i)), ', '];
        end
    end
    valueName1{order+1}(end-1) = ')'; 
    valueName1{order+1}(end) = []; 
    valueName2{order+1}(end-1) = ')'; 
    valueName2{order+1}(end) = []; 
    scaleVar = scaleVar.write(valueName1{order+1}, diffData{order+1}, 'value');
    scaleVar = scaleVar.write(valueName2{order+1}, scaleVar_Diff{order+1}, 'value');
end

ID.scaleVar = scaleVar;

% visualization:
CF = dataID_plot(ID);

end



function CF = dataID_plot(ID)

maxOrder = 3;
Position = [79 790 1327 522]; 
for order = 0:maxOrder
    CF(order+1) = figure; CA = axes(CF(order+1));
    CF(order+1).Position = Position + [order*30, -order*170 , 0, 0]; 
    ID.scaleVar.plot(['originalAxis', num2str(order)],...
                     ['diffData', num2str(order)], CA);
    1; 
end

for order = 0:maxOrder
    CF(maxOrder+order+2) = figure; CA = axes(CF(maxOrder+order+2));
    CF(maxOrder+order+2).Position = Position + [400+order*30, -order*170 , 0, 0]; 
    ID.scaleVar.plot(['scale'],...
                     ['scaleVar_Diff', num2str(order)], CA);
    CA.YAxis.Limits(1) = 0;
    1; 
end


end









