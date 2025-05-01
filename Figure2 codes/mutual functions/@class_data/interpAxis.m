% name: name of the axis to be interpolated.
% value: value of the axis after interpolation.
% currently, designed for 1st order data. (random walk type)
%
% this function interpolate data on a replacement axis.
        
function obj = interpAxis(obj, name, value, interp1Arg_)

if nargin == 3
    interp1Arg_ = {};
end

% if the dimension after interpolation is singular, permute it to the first dimension.
if length(value) == 1 
    [nameShort, axes, ~, ind] = obj.axis2data(name); 
    for i = 1:length(nameShort)
        obj = obj.permuteAxis(nameShort(i), axes{i}(1), name);
    end
end

[nameShort, axisName] = obj.nameID;

% tic
% find the data:
for i = 1:length(nameShort)
    Dim{i} = [];             % []: the value entry does not contain the axis (exclude the axis itself).
    if isempty(axisName{i})
        if strcmp(nameShort(i), name)
            ind_axis = i;  
            axis_old = obj.read(nameShort(i));
%             accuCoeffAxis_old = obj.readAccuCoeff(nameShort(i));
%             if length(accuCoeffAxis_old) == 1 
%                 if  class(accuCoeffAxis_old) == "double"
%                     accuCoeffAxis_old = accuCoeffAxis_old + axis_old*0;
%                 end
%             end
        end
    else
        for j = 1:length(axisName{i}) 
            if strcmp(axisName{i}(j), name)
                Dim{i} = [Dim{i},j];     % j: the data entry contains the axis as jth dimension.
                value_old{i} = obj.read(nameShort(i));
%                 accuCoeffValue_old{i} = obj.readAccuCoeff(nameShort(i));
%                 if length(accuCoeffValue_old{i}) == 1 
%                     if class(accuCoeffValue_old{i}) == "double"
%                         accuCoeffValue_old{i} = accuCoeffValue_old{i} + value_old{i}*0;
%                     end
%                 end
            end
        end
    end 
end
% toc

%  interp the axis:
% tic
axis_new = value;
% accuCoeffAxis_new = class_data.calculateAccuCoeff(axis_old, accuCoeffAxis_old, axis_new);
% toc

% interp the corresponding values:
% tic
for i = 1:length(nameShort)
    if ~isempty(Dim{i})
        countf = fprintf(['class_data -> interpAxis: ', char(nameShort(i)), '.']);
        for iDim = 1:length(Dim{i})
            if iDim >=2; value_old{i} = value_new{i}; end
            % permute for higher dimentional data (put the principal axis to the first):
            if length(axisName{i}) >1
                order = 1:length(axisName{i});
                order(1) = Dim{i}(iDim); order(Dim{i}(iDim)) = 1;
                value_old{i} = permute(value_old{i}, order);
%                 accuCoeffValue_old{i} = permute(accuCoeffValue_old{i}, order);
            else
                value_old{i} = value_old{i}(:);
            end

            % tic
            if isempty(setdiff(axis_new, axis_old))   % for fast operation for just cutting data.
                [~, iB] = ismember(axis_new, axis_old); 
                sizeAfter = size(value_old{i});    sizeAfter(1) = length(iB);
                if class(value_old{i}) == "double"
                    if issparse(value_old{i})
                        value_new{i} = sparse(sizeAfter(1), sizeAfter(2));   
                    else
                        value_new{i} = nan(sizeAfter);    % nan default is slower than 0 default.
                    end
                elseif class(value_old{i}) == "cell"
                    value_new{i} = cell(sizeAfter);
                else
                    value_new{i}(prod(sizeAfter)) = eval([class(value_old{i}), '(42)']);
                    value_new{i} = reshape(value_new{i}, sizeAfter);
                end
                if prod(sizeAfter) ~= 0
                    value_new{i}(:,:) = value_old{i}(iB,:);
                else
                    value_new{i}(:,:) = [];
                end
            else
                value_new{i} = interp1(axis_old, value_old{i}, axis_new, interp1Arg_{:});
            end
            % toc
% 
%             accuCoeffValue_new{i} = class_data.calculateAccuCoeff(axis_old, accuCoeffValue_old{i}, axis_new);

            if length(axisName{i}) >1
                value_new{i} = permute(value_new{i}, order);
%                 accuCoeffValue_new{i} = permute(accuCoeffValue_new{i}, order);
            end
        end
        fprintf(1, repmat('\b',1,countf));
    end

end
% toc


% write into obj:
obj = obj.write(name, axis_new, 'value');
% if class(accuCoeffAxis_new) == "double"
%     obj = obj.writeAccuCoeff(name, accuCoeffAxis_new);
% end
for i = 1:length(nameShort)
    if ~isempty(Dim{i})
        obj = obj.write(nameShort(i), value_new{i}, 'value');
%         obj = obj.writeAccuCoeff(nameShort(i), accuCoeffValue_new{i});
    end
end

end
        

    
    


