

% the grid is rescaled to square, so that vectors provide intuitional
% sense if the grid is shown as square.

% scale_ = 0 means to normalize the amplitude.
function [CP, scaleY] = vector(obj, axisName1, axisName2, Uname, Vname, CA, scale_)

if nargin == 6
    scale_ = 1.7;
end

% read data:
X = obj.read(axisName1);  
XID = obj.axisID(axisName1);
dX = XID.mode;

Y = obj.read(axisName2);  
YID = obj.axisID(axisName2);
dY = YID.mode;

U = obj.read(Uname);  
V = obj.read(Vname);  
if scale_ == 0
    norm = sqrt(U.^2+V.^2);
    U = U./norm;  V = V./norm;
end

% permute U, V so that the 1st index is axisName1
[nameShort, axis] = obj.nameID(Uname);
if strcmp(axis(1), axisName1) && strcmp(axis(2), axisName2)
elseif strcmp(axis(2), axisName1) && strcmp(axis(1), axisName2)
    U = permute(U, [2,1]);
else
    error("Yu: Uname is not consistent with axisName");
end
[nameShort, axis] = obj.nameID(Vname);
if strcmp(axis(1), axisName1) && strcmp(axis(2), axisName2)
elseif strcmp(axis(2), axisName1) && strcmp(axis(1), axisName2)
    V = permute(V, [2,1]);
else
    error("Yu: Vname is not consistent with axisName");
end



% scale Y to make the grids to be squares:
scaleY = abs(dX/dY);
Y_scaled = Y.*scaleY;

% plot:
if scale_ == 0
    CP = quiver(CA, X, Y_scaled, U', V', 'k'); hold(CA, 'on');
else
    CP = quiver(CA, X, Y_scaled, U', V', scale_, 'k'); hold(CA, 'on');
end

if contains(axisName1, 'time', 'IgnoreCase',true)
    datetickzoom(CA, 'x','HH:MM');
end

CA.YAxis.TickValues = (min(Y):(max(Y)-min(Y))/4:max(Y))*scaleY;
CA.YAxis.TickLabels = min(Y):(max(Y)-min(Y))/4:max(Y);
CA.YAxis.Limits = [min(Y), max(Y)]*scaleY;

% rescale the axis so that grids are squares:
backUnits = get(CA, 'units');
set(CA, 'Units', 'pixels')
pixPos = get(CA, 'Position');
set(CA, 'Units', backUnits);
if sum(~isnan(U(:)+V(:))) ~= 0
    CA.XAxis.Limits(1) = min(X(sum(~isnan(U+V), 2)>0));
    CA.XAxis.Limits(2) = CA.XAxis.Limits(1) + ...
        pixPos(3)/pixPos(4)*(max(CA.YAxis.Limits)-min(CA.YAxis.Limits));
end

end









