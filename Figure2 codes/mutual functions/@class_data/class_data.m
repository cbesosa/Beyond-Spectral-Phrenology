
% % Example:
% data = class_data;
% data = data.write("readme", value, "table"); 
% data = data.write("name(axis1,axis2)", value, "value");   % for name of axis data, specify name only.
%                                                           % axis could be omitted for overwritting values.
%
% data = data.remove(["name1", "name2"]);  % (axis) could be ommitted to read/remove/keep data.
% data = data.keep(["name1", "name2"]); 
% 
% varargout = data.read(["name1", "name2"]); 
% data.info;
%

classdef class_data
% This is a universal class of how I store a dataset.

    properties
        value cell {class_data.valid_value} % data values, 1-D cell array
                                            % may contain any type of data
                                            
        name string         % names, are a 1-D string array;
                            % format of name: name_of_data(name_of_Dimension_1, name_of_Dimension_2, ...)
                            
        table cell          % table contains dataset related information   
                            % {entry1<string> , value1; 
                            %  entry2<string> , value2; ...}  nx2 cell array contains a table 
                            %  of any general information of the data storage.
    end
    
    methods
        function obj = class_data()
            value = [];
            name = [];
            table = [];
        end
    end
    
    methods   % data manipulation methods
        
        function [D, varargout] = read(obj, names) % read data from value or table
            % read every into a struct if names not specified.
            % if read multiple names at a same time then return a struct D.
            
            if nargin == 1
                if isempty(obj.table)
                    names = string([]);
                else
                    names = [obj.table{:,1}];
                end
                if ~isempty(obj.value)
                    names = [names, obj.name(:)'];
                end
            end
            
            if ischar(names)
                names = string(names);
            end
            
            D = struct();
            for iName = 1:length(names)
                namei = names(iName);
                nameExist = 0;

                for i = 1:size(obj.table, 1)
                    if strcmp(namei, obj.table{i,1})
                        varargout{iName} = obj.table{i,2};
                        D.(strrep(class_data.strNameID(namei),' ','')) = obj.table{i,2};
                        nameExist = 1;
                        break;
                    end
                end

                for i = 1:length(obj.name)
                    if nameCompare(namei, obj.name(i))
                        varargout{iName} = obj.value{i};  
                        D.(strrep(class_data.strNameID(namei),' ','')) = obj.value{i};
                        nameExist = 1;
                        break;
                    end
                end

                if nameExist == 0
                    varargout{1} = "nonexist";
                    D.(strrep(class_data.strNameID(namei),' ','')) = "nonexist";
                end
            end
            
            if length(names) == 1
                D = D.(strrep(class_data.strNameID(names),' ','')); 
            end
            
            1;
        end
        
        
        % for overwritting, axes can be omitted.
        function obj = write(obj, name, value, where_) % write data to value or table, specified by "where"
                                                       % where_<default> = "value".
            if isstring(name) || ischar(name)
            else
                error("Yu: name is not a tring");
            end
            
            % default to colume vectors for 1-N value (except with more than 1 axis):
            [nameShort, axes, nameFull] = obj.nameID(name);
            if nameShort ~= "nonexist"
                axes_dimension = length(axes);
            else
                [nameShort, axes] = strNameID(name);
                axes_dimension = length(axes);
            end
            
            if nargin == 3
                where_ = '';
                for i = 1:size(obj.table, 1)
                    if strcmp(obj.table{i,1}, name)
                        where_ = 'table';
                    end
                end
                for i = 1:length(obj.name)
                    [boolNameShort, boolAxes] = nameCompare(name, obj.name(i));
                    if boolNameShort
                        where_ = 'value';
                    end
                end
                if isempty(where_)
                    where_ = 'value';
                end
            end
            
            if numel(value) == max(size(value)) && axes_dimension <= 1 && ~(ischar(value) && strcmp(where_, 'table'))
                if numel(value) ~= 1
                    value = value(:);     % vectors are columes by default, except char array in table.
                end
            end
            
            
            if strcmp(where_, 'value') 
                for i = 1:length(obj.name)
                    [nameShort, axes] = strNameID(name);
                    [boolNameShort, boolAxes] = nameCompare(name, obj.name(i));
                    if boolNameShort && isempty(axes)
                        obj.value{i} = value;  
                        return;
                    elseif boolNameShort
                        obj.name{i} = char(name);
                        obj.value{i} = [];
                        obj.value{i} = value;  
                        return;
                    end
                end
                ind = length(obj.name);
                obj.name(ind+1) = name;
                obj.value{ind+1} = value;  
                return;
            elseif strcmp(where_, 'table') 
                for i = 1:size(obj.table, 1)
                    if strcmp(obj.table{i,1}, name)
                        obj.table{i,2} = value;
                        return;
                    end
                end
                ind = size(obj.table, 1);
                obj.table{ind+1,1} = name;
                obj.table{ind+1,2} = value;  
                return
            end

        end
        
        function obj = remove(obj, names)        % delete named values
            if ischar(names)
                names = string(names);
            end
            for i = 1:length(names)
                found = 1;  nFound = 0;
                while found == 1 && size(obj.table, 1) >= 1
                    for j = 1:size(obj.table, 1)
                        found = 0;
                        if strcmp(obj.table{j,1}, names(i))
                            obj.table(j,:) = [];
                            found = 1; nFound = nFound+1; break;
                        end
                    end
                end
                found = 1;
                while found == 1 && length(obj.name) >= 1
                    for j = 1:length(obj.name)
                        found = 0;
                        if nameCompare(names(i), obj.name(j))
                            obj.value(j) = [];
                            obj.name(j) = [];
                            found = 1; nFound = nFound+1; break;
                        end
                    end
                end
                if nFound > 1
                    warning('Yu: trying to remove duplicated entries')
                elseif nFound == 0
                    warning('Yu: trying to remove no existing entries')
                end
            end
        end
        
        function obj = keep(obj, names) % keep named values, remove all others
            if ~isstring(names); error("Yu: name has to be string array"); end
            indName = [];  indTable = [];
            for i = 1:length(names)
                found = 0;
                for j = 1:length(obj.name)
                    if nameCompare(names(i), obj.name(j))
                        indName = [indName, j];
                        found = found +1;
                    end
                end
                for j = 1:size(obj.table, 1)
                    if strcmp(obj.table{j,1}, names(i))
                        indTable = [indTable, j];
                        found = found +1;
                    end
                end
                if found == 0
                    warning('Yu: trying to keep no existing entries')
                elseif found > 1
                    warning('Yu: trying to keep duplicated entries')
                end
            end
            obj.name = obj.name(indName);
            obj.value = obj.value(indName);
            obj.table = obj.table(indTable, :);
        end
        
        function info(obj) % display information of the class_data obj  

            disp('table:');
            if isempty(obj.table)
                disp('empty');
            else
                for i = 1:size(obj.table, 1)
                    temp{i,1} = i;
                    temp{i,2} = obj.table{i,1};
                    if size(obj.table, 2) == 2
                        temp{i,3} = obj.table{i,2};
                    else
                        temp{i,3} = [];
                    end
                end
                disp(temp);
                clear('temp');
            end
            
            disp('name, value:');
            if isempty(obj.name)
                disp('empty');
            else
                for i = 1:length(obj.name)
                    temp{i,1} = i;
                    temp{i,2} = obj.name(i);
                    if i <= length(obj.value)
                        temp{i,3} = obj.value{i};
                    end
                end
                disp(temp);
                clear('temp');
            end
        end
    end
    
    
    
    methods   % administrative methods
        
        function [nameShort, axes, nameFull, indInValue] = nameID(obj, name_)
        % name identity
        % decompose the nameFull into nameShort and axes(axis in order).
        % return "nonexist" if not exist. 
            if length(obj.name) == 0
                axes = "nonexist";
                nameShort = "nonexist";
                nameFull = "nonexist";
                return;
            end
        
            for i = 1:length(obj.name)
                nameFull{i} = char(obj.name(i)); 
                nameFull{i}(nameFull{i} == ' ') = [];
                ind1 = find(nameFull{i} == '(');
                ind2 = find(nameFull{i} == ')');
                if isempty (ind1) && isempty (ind2)
                    nameShort(i) = string(nameFull{i});
                    axes{i} = [];
                elseif ind1 >1 && ind2 == length(nameFull{i}) && ind2>=ind1+2
                    nameShort(i) = string(nameFull{i}(1:ind1-1));
                    axisAll = nameFull{i}(ind1+1:ind2-1);
                    axisAll = split(axisAll,',');   % expect a cell array
                    for j = 1:length(axisAll)
                        axes{i}(j) = string(axisAll{j});
                    end
                elseif ind1 >1 && ind2 == length(nameFull{i}) && ind2 == ind1+1
                    nameShort(i) = string(nameFull{i}(1:ind1-1));
                    axes{i} = strings(0);
                end
                indInValue(i) = i;
            end
            
            for i = 1:length(obj.name)
                nameFullNew(i) = string(nameFull{i});
            end
            nameFull = nameFullNew;
            
            if nargin == 1
                return;
            elseif nargin == 2
                ind = find(nameShort == name_);
                if isempty(ind)
                    ind = find(nameFull == name_);
                end
                if isempty(ind)
                    axes = "nonexist";
                    nameShort = "nonexist";
                    nameFull = "nonexist";
                else
                    axes = axes{ind};
                    nameShort = nameShort(ind);
                    nameFull = nameFull(ind);
                    indInValue = indInValue(ind);
                end
            end
            
            return;
        end
          
    end
    
    methods   % children methods
        
        
    end
    
    
    methods (Static)
        function valid_value(value)
            
        end
        
        accuCoeff_new = calculateAccuCoeff(ticks_old, AccuCoeff_old, ticks_new)
        
        % merge objArray with same structure into a single obj:
        obj = concatenate(objArray, axisName)
        
        function [nameShort, axes] = strNameID(nameFull)
            [nameShort, axes] = strNameID(nameFull);
        end
        
        function nameFull = nameIDstr(nameShort, axes)
            nameFull = nameIDstr(nameShort, axes);
        end
    end

end


function [boolNameShort, boolAxes] = nameCompare(name1, name2)

% name1, name2 could be either name or name(axis)

    [nameShort1, axes1] = strNameID(name1);
    [nameShort2, axes2] = strNameID(name2);

    if strcmp(nameShort1, nameShort2)
        boolNameShort = 1;
    else
        boolNameShort = 0; 
    end
    
    boolAxes = 1;
    if length(axes1) == length(axes2)
        for i = 1:length(axes1)
            if strcmp(axes1(i), axes2(i))
            else
                boolAxes = 0;
            end
        end
    else
        boolAxes = 0;
    end

end


function [nameShort, axes] = strNameID(name)
% return empty if cannot decode the name.

    name = char(name);
    nameShort = '';
    axes = ""; axes(1) = [];
    for i = 1:length(name)
        if name(i) ~= '('
            nameShort = [nameShort, name(i)];
        else
            indAxes = 1;
            axes(indAxes) = "";
            for j = i+1:length(name)
                if name(j) == ','
                    indAxes = indAxes +1;
                    axes(indAxes) = "";
                elseif name(j) == ')' && j == length(name)
                    break;
                elseif name(j) ~= '(' && name(j) ~= ',' && name(j) ~= ')'
                    axes(indAxes) = string([char(axes(indAxes)), name(j)]);
                else
                    nameShort = [];
                    axes = ""; axes(1) = [];
                end
            end
            break;
        end
    end

end

function nameFull = nameIDstr(nameShort, axes)
% reconstruct nameFull from nameShort and axes.

    nameFull = char(nameShort);
    if isempty(axes)
        return;
    else
        nameFull = [nameFull, '('];
        for i = 1:length(axes)
            nameFull = [nameFull, char(axes(i)), ','];
        end
        nameFull(end) = ')';
    end
    
end


























