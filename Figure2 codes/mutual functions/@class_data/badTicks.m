% 
% 
        
function [badTicks, reason, obj] = badTicks(obj, dataName, principalAxis)
% obj contains the badTicks of date with surfix of _badTick, type bool.

% read necessary information:
[ID, reAxis] = obj.axisID(principalAxis);
% if ID.type ~= "isometric"
%     error("Yu: for uniform grid only.")
% end
% stdGrid = sqrt(stdWiener^2*ID.mode);
axis = obj.read(principalAxis);
data = obj.read(dataName);

% permute the axis to be the 1st index:
[nameShort, axisNames] = obj.nameID(dataName); 
ind =  find(axisNames ==  principalAxis);
if length(ind) ~= 1
    error("Yu: check names error.")
end
order = 1:length(axisNames);
order(1) = ind; order(ind) = 1;
data = permute(data, order);

% the windowed average depending on 1st order data assumption: 
[badTickInd{1}, reasonSub{1}] = checkStd(data);
[badTickInd{2}, reasonSub{2}] = checkValueStdRatio(data);
[badTickInd{3}, reasonSub{3}] = checkCurvature(data);

badTickInds = [];
reason = [];
for i = 1:length(badTickInd)
    
    badTickInds = [badTickInds; badTickInd{i}(:)];
    reason = [reason; reasonSub{i}(:)];
end
badTickInds = unique(badTickInds);
badTicks = axis(badTickInds);

% construct obj with _badTick:
badTickBool = axis*0;
badTickBool(badTickInds) = 1;
[nameShort, axes, nameFull] = obj.nameID(dataName); 
nameFull = [char(nameShort), '_badTick(',char(principalAxis),')'];
obj = obj.write(nameFull, badTickBool);

1;
end
        

function [badTickInd, reason] = checkStd(data)
% For N-D arrays, movmean operates along the first array dimension.
% window is centered and odd in elements.

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
dataStd = nan(1,size(data,2));
for i = 1:size(data,2)
    dataStd(i) = std(data(:,i));
end
reciprocal = 1./dataStd;

criticalRange = mean(sortTrim(dataStd, [0.05, 0.95])) + [-1 1]*3*std(sortTrim(dataStd, [0.05, 0.95]));
[~, badTickInd1] = find(dataStd < criticalRange(1) | dataStd > criticalRange(2));
reason(1).str = "std out of 3 sigma";
reason(1).ticks = badTickInd1;

criticalRange = mean(sortTrim(reciprocal, [0.05, 0.95])) + [-1 1]*3*std(sortTrim(reciprocal, [0.05, 0.95]));
[~, badTickInd2] = find(reciprocal  < criticalRange(1) | reciprocal  > criticalRange(2));
reason(2).str = "std reciprocal out of 3 sigma";
reason(2).ticks = badTickInd2;

badTickInd = union(badTickInd1, badTickInd2);

end

function [badTickInd, reason] = checkValueStdRatio(data)
% For N-D arrays, movmean operates along the first array dimension.
% window is centered and odd in elements.

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
ValueStdRatio = nan(1,size(data,2));
for i = 1:size(data,2)
    ValueStdRatio(i) = mean(abs(data(:,i)))/std(data(:,i));
end
reciprocal = 1./ValueStdRatio;

criticalRange = mean(sortTrim(ValueStdRatio, [0.05, 0.95])) + [-1 1]*3*std(sortTrim(ValueStdRatio, [0.05, 0.95]));
[~, badTickInd1] = find(ValueStdRatio < criticalRange(1) | ValueStdRatio > criticalRange(2));
reason(1).str = "ValueStdRatio out of 3 sigma";
reason(1).ticks = badTickInd1;

criticalRange = mean(sortTrim(reciprocal, [0.05, 0.95])) + [-1 1]*3*std(sortTrim(reciprocal, [0.05, 0.95]));
[~, badTickInd2] = find(reciprocal  < criticalRange(1) | reciprocal  > criticalRange(2));
reason(2).str = "ValueStdRatio reciprocal out of 3 sigma";
reason(2).ticks = badTickInd2;

badTickInd = union(badTickInd1, badTickInd2);

end

function [badTickInd, reason] = checkCurvature(data)
% For N-D arrays, movmean operates along the first array dimension.
% window is centered and odd in elements.

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
featuredCurvature = nan(1,size(data,2));
for i = 3:size(data,2)-2  
    curvature1 = abs(2*data(:,i-1)-data(:,i-2)-data(:,i));
    curvature2 = abs(2*data(:,i)-data(:,i-1)-data(:,i+1));
    curvature3 = abs(2*data(:,i+1)-data(:,i)-data(:,i+2));  
    
    minCurvature = min( curvature1, curvature2 );  
    minCurvature = min( minCurvature, curvature3 );  
    featuredCurvature(i) = mean(minCurvature);                
end
for i = 2
    curvature2 = abs(2*data(:,i)-data(:,i-1)-data(:,i+1));
    curvature3 = abs(2*data(:,i+1)-data(:,i)-data(:,i+2));
    
    minCurvature = min( curvature2, curvature3 );  
    featuredCurvature(i) = mean(minCurvature);                
end
for i = size(data,2)-1  
    curvature1 = abs(2*data(:,i-1)-data(:,i-2)-data(:,i));
    curvature2 = abs(2*data(:,i)-data(:,i-1)-data(:,i+1));
    
    minCurvature = min( curvature1, curvature2 );  
    featuredCurvature(i) = mean(minCurvature);                
end
featuredCurvature(1) = mean(abs(data(:,1) - data(:,2)));
featuredCurvature(end) = mean(abs(data(:,end) - data(:,end-1)));


featuredCurvature = featuredCurvature/std(featuredCurvature);

criticalRange = mean(sortTrim(featuredCurvature, [0.05, 0.95])) + [-1 1]*6*std(sortTrim(featuredCurvature, [0.05, 0.95]));
[~, badTickInd1] = find(featuredCurvature > criticalRange(2));
reason(1).str = "featured Curvature exceed 6 sigma";
reason(1).ticks = badTickInd1;

badTickInd = union(badTickInd1, badTickInd1);

end




    

