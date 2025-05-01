

function [CF, CA, CB, CP] = plot(obj, DataName, CA_, offSet_, plotArg_)

% offset_: offset value for plot2, isempty for default offset.
% plotArg_ = <cell>, arguments that go into plot function

if nargin == 2
    CA_ = [];
elseif nargin == 3
    offSet_ = [];
    plotArg_ = {};
elseif nargin == 4
    plotArg_ = {};
end

[nameShort, Axes] = obj.nameID(DataName);
if length(Axes) == 1
    if isempty(CA_); CF = figure('Position', [100, 100, 800, 700]); CA_ = axes(CF);  end
    if length(CA_) ~=1; error("Yu: n CA inconsistent"); end 
    CF = CA_.Parent;
    CP = plot1(obj, DataName, CA_, plotArg_);
elseif length(Axes) == 2
    if isempty(CA_); CF = figure('Position', [100, 100, 800, 700]); CA_ = axes(CF);  end
    if length(CA_) ~=1; error("Yu: n CA inconsistent"); end 
    CF = CA_.Parent;
    xAxisName = Axes(1);
    CP = plot2(obj, xAxisName, DataName, CA_, offSet_, plotArg_);
elseif length(Axes) > 2
    dataValue = obj.read(DataName);
    axis1 = obj.read(Axes(1));
    axis2 = obj.read(Axes(2));
    axis3 = obj.read(Axes(3));
    
    % check if dimension if CA is correct:
    nFig = length(dataValue(1,1,1,:));
    nAxes = length(dataValue(1,1,:,1));
    if ~isempty(CA_)
        if length(CA_) ~= length(axis3) || length(Axes) ~= 3
            error("Yu: n CA inconsistent");
        end
        CF = [];
        CA_ = {CA_};
    else
        % create CF and CA if not exist:
        figArr = class_figArr;
        figArr.X = ones(1,nAxes);  figArr.dx = [0.3 0.1*ones(1,nAxes-1) 0.3];
        figArr.Y = [1];     figArr.dy = [0.3 0.3]; 
        figArr.GM = 1:nAxes;
        figArr.AsRatio = 10/(nAxes*7.5+3);
        figArr.width = 0.8/figArr.AsRatio;
        figArr.screenOffset = [0.5, 0.8];
        figArr.figOut = 'name';
        for iFig = 1:nFig
            [CF(iFig), CA_{iFig}, CAnno{iFig}] = create(figArr);
        end
    end
    
    % plot each CA using the plot2 function:
    for iFig = 1:nFig
        for iCA = 1:nAxes
            dataValueSub = dataValue(:,:,iCA,iFig);
            dataSub = class_data;
            dataSub = dataSub.write(Axes(1), axis1);
            dataSub = dataSub.write(Axes(2), axis2);
            dataSub = dataSub.write(class_data.nameIDstr(nameShort, Axes(1:2)),dataValueSub);
            xAxisName = Axes(1);
            CP{iFig, iCA} = plot2(dataSub, xAxisName, DataName, CA_{iFig}(iCA), offSet_, plotArg_);
        end
    end
    if length(CA_) == 1; CA_ = CA_{1}; end
end

CA = CA_; 
CB = [];
return;
end

function CP = plot1(obj, DataName, CA, plotArg_)

% read data and axis:
[nameShort, axes, nameFull] = obj.nameID(DataName);
data = obj.read(DataName);
AV = obj.read(axes);

% plot for 1D data:
stdMax = max(std(data, 'omitnan'));
YLimits = [min(data), max(data)];
CP = plot(CA, AV, data, plotArg_{:}); 
hold(CA, 'on');

CA.XAxis.Limits = [min(AV), max(AV)];
if YLimits(2) > YLimits(1)
    CA.YAxis.Limits = YLimits;
end

CA.Title.String = nameFull;
CA.XAxis.Label.String = axes;
CA.YAxis.Label.String = DataName;
legend(CP, DataName, 'box', 'off');

end

function CP = plot2(obj, axisName, DataName, CA, offSet_, plotArg_)

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
stdMax = max(std(data, 'omitnan'));
colorMap = getNCLColmap('BlAqGrYeOrReVi200.rgb', length(AV{2}));
colorMap = getNCLColmap('NCV_rainbow2.rgb', length(AV{2}));
colorMap = getNCLColmap('helix1.rgb', length(AV{2}));
colorMap = getNCLColmap('thelix.rgb', length(AV{2}));
colorMap = getNCLColmap('GMT_split.rgb', length(AV{2}));
YLimits = [inf, -inf];
for i = 1:length(AV{2})
    if isempty(offSet_)
        offSet = (i-1)*stdMax*3.5;
    else
        offSet = offSet_;
    end
    CP(i) = plot(CA, AV{1}(:), offSet+data(:,i), 'color', colorMap(i,:), plotArg_{:}); 
    YLimits = [min(min(offSet + data(:,i)), YLimits(1)), max(max(offSet + data(:,i)), YLimits(2))];
    hold(CA, 'on');
    if offSet ~= 0
        plot(CA, AV{1}, offSet + data(:,i)*0, ':k', 'linewidth', 1.5);
    end
end

legendStr = {};
for i = 1:length(AV{2})
    legendStr{end+1} = [char(axes(2)),' = ',num2str(AV{2}(i))];
end
if length(AV{2}) >6
    CL = legend(CP([1,end]), {[char(axes(2)),' = ',num2str(AV{2}(001))], ...
                         [char(axes(2)),' = ',num2str(AV{2}(end))]},  'box', 'off');
else
    CL = legend(CP, legendStr, 'box', 'off');
end

CL.Box = 'off';
CL.EdgeColor = 'none';
CL.Color = 'none';
                 
CA.XAxis.Limits = [min(AV{1}), max(AV{1})];
CA.YAxis.Limits = YLimits;

CA.Title.String = nameFull;
CA.XAxis.Label.String = axes(1);
CA.YAxis.Label.String = DataName;

end




