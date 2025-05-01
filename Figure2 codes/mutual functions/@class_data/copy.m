% range: [left1, right1; left2, right2; ....]
       
function obj1 = copy(obj1, obj2, name)

% copy value and associated axis in obj2 to obj1.

[nameShort, axes, nameFull, ind] = obj2.nameID(name); 

for i = 1:length(axes)
    axesValue = obj1.read(axes(i));
    if class(axesValue) == "string"
        if axesValue == "nonexist"
            obj1 = obj1.write(axes(i), obj2.read(axes(i)));
        end
    end 
end

obj1 = obj1.write(nameFull, obj2.read(nameShort));

end
        

    
    


