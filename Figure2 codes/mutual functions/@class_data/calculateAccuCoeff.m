% for 1st order data:
% for higher dimensional AccuCoeff_old, perform on the first index. (similar to interp1)
function accuCoeff_new = calculateAccuCoeff(ticks_old, AccuCoeff_old, ticks_new)

if class(AccuCoeff_old) == "string" || class(AccuCoeff_old) == "char"
    if AccuCoeff_old == "nonexist"
        accuCoeff_new = "nonexist"; return;
    end
end

if isempty(setdiff(ticks_new, ticks_old))   % for fast operation for just cutting the axis.
    [~, iB] = ismember(ticks_new, ticks_old); 
    accuCoeff_new = AccuCoeff_old(iB); return;
end


if length(AccuCoeff_old(:))/length(AccuCoeff_old) == 1
    AccuCoeff_old = AccuCoeff_old(:);
end
[TO, ind] = sort(ticks_old, 'ascend');
ACO = AccuCoeff_old(ind,:);
[TN, ind] = sort(ticks_new, 'ascend');


% add a head and a tail so that ticks_new do not outbound:
if TO(1) > TN(1)
    TO = [TN(1)-9999*(TO(end)-TO(1)) , TO(:)'];
    ACO(2:size(ACO,1)+1, :) = ACO;
    ACO(1,:) = 9999*sqrt(TO(end)-TO(1));
end
if TO(end) <= TN(end)
    TO = [TO(:)', TN(end)+9999*(TO(end)-TO(1))];
    ACO(end+1,:) = 9999*sqrt(TO(end)-TO(1));
end

j1 = 1; j2 = 2;  i = 1;
sizeNew = size(ACO);  sizeNew(1) = length(TN);
ACN = nan(sizeNew);
while i <= length(TN)
    a = TN(i) - TO(j1);  b = TO(j2) - TN(i);
    if a >=0 && b >=0
        ACN(i,:) = sqrt( (ACO(j1,:).^2+a).*(ACO(j2,:).^2+b)./(ACO(j1,:).^2+a +ACO(j2,:).^2+b) );
        i = i+1;
    else
        j1 = j1+1; j2 = j2+1;
    end
end

ticks_new_reconstructed(ind) = TN;
accuCoeff_new = ACN * nan;
accuCoeff_new(ind,:) = ACN(:,:);

end





