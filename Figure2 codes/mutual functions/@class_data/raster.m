

function CP = raster(obj, axisName, DataName, CA, offSet_)

% offset_: offset value for plot2, isempty for default offset.


if nargin == 4
    offSet_ = [];
end

[nameShort, axes] = obj.nameID(DataName);
if length(axes) == 1
    CP = plot1(obj, DataName, CA);
elseif length(axes) == 2
    CP = plot2(obj, axisName, DataName, CA, offSet_);
else
    error("Yu: undeveloped");
end

end

function CP = plot1(obj, DataName, CA)

% read data and axis:
[nameShort, axes, nameFull] = obj.nameID(DataName);
data = obj.read(DataName);
AV = obj.read(axes);

% plot for 1D data:
data(data == 0) = nan;
data(data == 1) = 0;
CP = plot(CA, AV, data, '.'); 
hold(CA, 'on');

CA.XAxis.Limits = [min(AV), max(AV)];
CA.YAxis.Limits = [-0.5,0.5];

CA.Title.String = nameFull;
CA.XAxis.Label.String = axes;
CA.YAxis.Label.String = DataName;
legend(CP, DataName);

end

function CP = plot2(obj, axisName, DataName, CA, offSet_)

% general permutation for N-D data:
[nameShort, axes, nameFull] = obj.nameID(DataName);
data = obj.read(DataName);
for i = 1:length(axes)
    AV{i} = obj.read(axes(i));
end
ind = find(axes == axisName);
if length(ind) ~= 1; error("Yu: unexpected"); end
order = 1:length(axes);
order(1) = ind; 
order(ind) = 1;
AV = {AV{order}};
axes = axes(order);
data = permute(data, order);

% plot for 2D data:
data(data == 0) = nan;
data(data == 1) = 0;
colorMap = getNCLColmap('BlAqGrYeOrReVi200.rgb', length(AV{2}));
YLimits = [inf, -inf];
for i = 1:length(AV{2})
    if isempty(offSet_)
        offSet = (i-1)*1;
    else
        offSet = offSet_;
    end
    CP(i) = plot(CA, AV{1}, offSet + data(:,i), '.', 'color', colorMap(i,:)); 
    YLimits = [min(min(offSet + data(:,i))-0.5, YLimits(1)), max(max(offSet + data(:,i))+0.5, YLimits(2))];
    hold(CA, 'on');
    
    % plot(CA, AV{1}, offSet + data(:,i)*0, ':k', 'linewidth', 1.5);
end
legend(CP([1,end]), {[char(axes(2)),' = ',num2str(AV{2}(001))], ...
                     [char(axes(2)),' = ',num2str(AV{2}(end))]});
CA.XAxis.Limits = [min(AV{1}), max(AV{1})];
CA.YAxis.Limits = YLimits;

CA.Title.String = nameFull;
CA.XAxis.Label.String = axes(1);
CA.YAxis.Label.String = DataName;

end




