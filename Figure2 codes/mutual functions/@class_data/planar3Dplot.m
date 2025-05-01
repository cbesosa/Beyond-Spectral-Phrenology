

function [CP, scaleFactor_] = planar3Dplot(obj, xAxisName, yAxisName, dataName, zDirection, CA, scaleFactor_, plotArg_)

% plot 3D data in 2D, squeeze the zAxis(Data) into either x or y axis, specified by zDirection.
% scaleFactor_: amplification of dx or dx to zDirection, scaleFactor_ = dz/dx or dz/dy.
% plotArg_ = <cell>, arguments that go into plot function


if nargin == 6
    scaleFactor_ = [];
    plotArg_ = {};
elseif nargin == 7
    plotArg_ = {};
end

[nameShort, axes] = obj.nameID(dataName);
if length(axes) ~= 2
    error("Yu: for 2D axes only");
else
    if xAxisName == axes(1) && yAxisName == axes(2)
    elseif xAxisName == axes(2) && yAxisName == axes(1)
        obj = obj.permuteAxis(dataName, axes(1), axes(2));
    else
        error("Yu: axes name miss match")
    end
        
    if zDirection == "y"
        [~, objSub] = obj.segments(xAxisName);
    elseif zDirection == "x"
        [~, objSub] = obj.segments(yAxisName);
    else
        error("Yu: unexpected.")
    end
    for i = 1:length(objSub)
        [CP, scaleFactor_] = plot2(objSub(i), dataName, zDirection, CA, scaleFactor_, plotArg_);
    end
    
end

end



function [CP, scaleFactor_] = plot2(obj, dataName, zDirection, CA, scaleFactor_, plotArg_)

[nameShort, axes] = obj.nameID(dataName);
dataValue = obj.read(dataName);
X = obj.read(axes(1));
Y = obj.read(axes(2));

% plot for 2D data:
scaleApprox = std(dataValue(:), 'omitnan');

if zDirection == "x"
    nLines = length(X);
elseif zDirection == "y"
    nLines = length(Y);
else
    error("Yu: zDirection wrong");
end
colorMap = getNCLColmap('BlAqGrYeOrReVi200.rgb', nLines);
colorMap = getNCLColmap('NCV_rainbow2.rgb', nLines);
colorMap = getNCLColmap('helix1.rgb', nLines);
colorMap = getNCLColmap('thelix.rgb', nLines);
colorMap = getNCLColmap('GMT_split.rgb', nLines);

if zDirection == "x"
    % XLimits = [inf, -inf];
    [AID, reAxis] = obj.axisID(axes(1));
    if isempty(scaleFactor_)
        scaleFactor_ = scaleApprox/AID.mode*3.0;
    end
    for i = 1:nLines
        subData = dataValue(i, :);
        subData = subData/scaleFactor_;
        if isempty(X)
            CP(i) = plot(CA, nan, nan, 'color', 'k', plotArg_{:}); 
        else
            CP(i) = plot(CA, X(i)+subData, Y, 'color', 'k', plotArg_{:}); 
        end
        % XLimits = [min(min(offSet + data(:,i)), YLimits(1)), max(max(offSet + data(:,i)), YLimits(2))];
        hold(CA, 'on');
        % plot(CA, AV{1}, offSet + data(:,i)*0, ':k', 'linewidth', 1.5);
    end
elseif zDirection == "y"
    % XLimits = [inf, -inf];
    [AID, reAxis] = obj.axisID(axes(2));
    if isempty(scaleFactor_)
        scaleFactor_ = scaleApprox/AID.mode*3.0;
    end
    for i = 1:nLines
        subData = dataValue(:, i);
        subData = subData/scaleFactor_;
        if isempty(X)
            CP(i) = plot(CA, nan, nan, 'color', 'k', plotArg_{:}); 
        else
            CP(i) = plot(CA, X, Y(i)+subData, 'color', 'k', plotArg_{:}); 
        end
        % YLimits = [min(min(offSet + data(:,i)), YLimits(1)), max(max(offSet + data(:,i)), YLimits(2))];
        hold(CA, 'on');
        % plot(CA, AV{1}, offSet + data(:,i)*0, ':k', 'linewidth', 1.5);
    end
end


% legend(CP([1,end]), {[char(axes(2)),' = ',num2str(AV{2}(001))], ...
%                      [char(axes(2)),' = ',num2str(AV{2}(end))]});
% CA.XAxis.Limits = [min(AV{1}), max(AV{1})];
% CA.YAxis.Limits = YLimits;

% CA.Title.String = nameFull;
CA.XAxis.Label.String = axes(1);
CA.YAxis.Label.String = axes(2);

end




