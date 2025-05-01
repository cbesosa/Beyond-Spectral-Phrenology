function [type, SZ] = dimension(obj, name_)
% SZ = <size>;

    
    [nameShort, axes, nameFull] = obj.nameID();
    if nargin == 1
        ind = 1:length(nameShort);
    else
        [~, ind] = ismember(name_, nameShort);
        [~, ind2] = ismember(name_, nameFull);
        ind(ind == 0) = ind2(ind == 0);
        if ismember(0, ind)
            error("Yu: woring name_ specified")
        end
    end

    for i = 1:length(obj.name)
        if ismember(i, ind)
            [nameShort, axes, nameFull] = obj.nameID(obj.name(i));

            % determine the type:
            if length(axes) == 0
                if contains(nameFull, '(') && contains(nameFull, ')')
                    type(i) = "value";
                else
                    type(i) = "axis";
                end
            else
                type(i) = "value";
            end

            % determine the size:
            SZ{i} = size(obj.value{i});
            if SZ{i}(2) == 1 && length(axes) <= 1
                SZ{i}(2) = [];
            end
            if SZ{i}(1) == 1 && length(axes) <= 1 && length(SZ{i}) == 2
                SZ{i}(1) = [];
            end
    %         if length(SZ{i}) == 2
    %             if SZ{i}(1) == 1 && SZ{i}(2) == 1 && length(axes) == 0
    %                 SZ{i} = 1;
    %             end
    %         end
        end
    end

    % add checking for dimension compatibility.
    
    
    
    
        
    %
    type = type(ind);
    SZ = SZ(ind);
    
    if nargin == 2 && length(ind) == 1
        SZ = SZ{:};
    end
    

    1;
end






