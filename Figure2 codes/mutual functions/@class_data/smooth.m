% stdMeasure: sigma of measurement error. (assuming Gaussian)
% currently, designed for 1st order data. (random walk type)
%
% smooth a value along given axis.
        
function obj = smooth(obj, dataName, axisName, stdMeasure, stdWiener, N_)

if nargin == 5
    N_ = 41;     % number of point of the window, odd number.
end

% read necessary information:
[ID, reAxis] = obj.axisID(axisName);
if ID.type ~= "isometric"
    error("Yu: for uniform grid only.")
end
stdGrid = sqrt(stdWiener^2*ID.mode);
data = obj.read(dataName);

% permute the axis to be the 1st index:
[nameShort, axisNames] = obj.nameID(dataName); 
ind =  find(axisNames ==  axisName);
if length(ind) ~= 1
    error("Yu: check names error.")
end
order = 1:length(axisNames);
order(1) = ind; order(ind) = 1;
if length(order) > 1
    data = permute(data, order);
end

% the windowed average depending on 1st order data assumption: 
windowRight = 1./(stdMeasure^2+(1:(N_/2-1))*stdGrid^2);
windowLeft = fliplr(windowRight);  
mid = 1/(stdMeasure^2);
window = [windowLeft, mid, windowRight];
data = movmean_windowed(data, window);

% permute back and write:
if length(order) > 1
    data = permute(data, order);
end
[nameShort, axes, nameFull] = obj.nameID(dataName); 
nameShort = [char(nameShort), '_smooth'];
nameFull = class_data.nameIDstr(nameShort, axes);
obj = obj.write(nameFull, data);

end
        

function dataMean = movmean_windowed(data, window)
% For N-D arrays, movmean operates along the first array dimension.
% window is centered and odd in elements.

window = window(:);
if length(data) == numel(data)
    data = data(:);
end

% permute dimension and reshape:
SD = size(data);
order = 1:length(SD);
order(1) = order(end);
order(end) = 1;

data = permute(data, order);
SD_permuted = SD; 
SD_permuted(end) = SD(1);
SD_permuted(1) = SD(end);
data = reshape(data, prod(SD_permuted(1:end-1)), SD_permuted(end));

% calculate the result for range with full window:
dataMean = data*nan;
window = window/sum(window);
range1 = length(window)/2+0.5 : size(data,2)-length(window)/2+0.5;
for i = range1
    dataMean(:,i) = data(:,i-length(window)/2+0.5:i+length(window)/2-0.5)*...
                    window;
end

% calculate the result for range with partial window:
for i = setdiff(1:size(data,2), range1)
    nLeft = i -1;
    nLeft = min(nLeft, length(window)/2-0.5);
    nRight = size(data,2) -i;
    nRight = min(nRight, length(window)/2-0.5);
    partialWindow = window(length(window)/2+0.5 -nLeft : length(window)/2+0.5 +nRight);
    partialWindow = partialWindow/sum(partialWindow);
    dataMean(:,i) = data(:,i-nLeft:i+nRight)*partialWindow;
end


dataMean = reshape(dataMean, SD_permuted);
dataMean = permute(dataMean, order);

end





    

