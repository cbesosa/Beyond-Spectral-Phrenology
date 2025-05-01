
function [CLines, CP, CB] = contourf(obj, axisName1, axisName2, DataName, CA, lines_)

if nargin == 5
    lines_ = [];
end
if isempty(lines_)
    CLines = [];
end


% read data:
X = obj.read(axisName1);  
XID = obj.axisID(axisName1);
dX = XID.mode;

Y = obj.read(axisName2);  
YID = obj.axisID(axisName2);
dY = YID.mode;

Data = obj.read(DataName);  

% permute Data so that the 1st index is axisName1
[nameShort, axis] = obj.nameID(DataName);
if strcmp(axis(1), axisName1) && strcmp(axis(2), axisName2)
elseif strcmp(axis(2), axisName1) && strcmp(axis(1), axisName2)
    Data = permute(Data, [2,1]);
else
    error("Yu: Uname is not consistent with axisName");
end


% pcolor:
CP(1) = pcolor(CA, X, Y, Data'); hold(CA, 'on');
shading(CA, 'flat');  % shading(CA, 'interp');
CA.Colormap = colMapPM(100, [0,0], [1,1]);
CA.YAxisLocation = 'left';
CA.XAxis.Label.String = axisName1;
CA.YAxis.Label.String = axisName2;

CAposRight = sum(CA.Position([1,3]));
CB = colorbar(CA);
CB.Position(1) = CAposRight + 0.01;
CB.Position(3) = 0.02;
CB.Label.String = DataName;

% contour:
if ~isempty(lines_)
    [CM, CP(2)] = contour(CA, X, Y, Data', lines_, 'color', 'm', 'ShowText','on');
    CLines = readContour(CM);
end

if contains(axisName1, 'time', 'IgnoreCase',true)
    datetickzoom(CA, 'x','HH:MM')
end


% rescale the axis so that grids are squares:
backUnits = get(CA, 'units');
set(CA, 'Units', 'pixels')
pixPos = get(CA, 'Position');
set(CA, 'Units', backUnits);
CA.XAxis.Limits = min(X(sum(~isnan(Data), 2)>0)) + [0, CA.XAxis.Limits(2) - CA.XAxis.Limits(1)];
CA.XAxis.Limits(2) = CA.XAxis.Limits(1) + ...
    pixPos(3)/pixPos(4)*(max(CA.YAxis.Limits)-min(CA.YAxis.Limits))*abs(dX/dY);


end









