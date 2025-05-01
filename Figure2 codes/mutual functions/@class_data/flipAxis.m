
       
function obj = flipAxis(obj, name)

oldAxis = obj.read(name);
newAxis = oldAxis(end:-1:1);
obj = obj.interpAxis(name, newAxis);

end
        

    
    


