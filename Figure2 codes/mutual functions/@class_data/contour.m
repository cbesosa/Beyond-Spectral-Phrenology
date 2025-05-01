

% the grid is rescaled to square, so that vectors provide intuitional
% sense if the grid is shown as square.

% scale_ = 0 means to normalize the amplitude.
% format of lines_ is consistent with Matlab contour.
function [CLines, CP] = contour(obj, axisName1, axisName2, DataName, CA, lines_)

if nargin == 5
    lines_ = [];
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


% plot:
[CM, CP] = contour(CA, X, Y, Data', lines_, 'color', 'm', 'ShowText','on'); hold(CA, 'on');
CLines = readContour(CM);

if contains(axisName1, 'time', 'IgnoreCase',true)
    datetickzoom(CA, 'x','HH:MM')
end


% rescale the axis so that grids are squares:
backUnits = get(CA, 'units');
set(CA, 'Units', 'pixels')
pixPos = get(CA, 'Position');
set(CA, 'Units', backUnits);
CA.XAxis.Limits(1) = min(X(sum(~isnan(Data), 2)>0));
CA.XAxis.Limits(2) = CA.XAxis.Limits(1) + ...
    pixPos(3)/pixPos(4)*(max(CA.YAxis.Limits)-min(CA.YAxis.Limits))*abs(dX/dY);


end









