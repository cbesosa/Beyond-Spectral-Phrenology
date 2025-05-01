
function RV = subjectInfo(outPath, par)
% input is the animal session name.
% output is to create a dataobj, and record datapath and general meta info into table.

    data = class_data;
    
    % fillable entries without knowing a specific dataset:
    data = data.write("notes", '', "table");
    data = data.write("aux", struct(), "table");
    data = data.write("animalNameSession", par, "table");
    
    % dataset specific entries need to be filled:
    data = data.write("datasetName", "", "table");
    data = data.write("animalName", "", "table");
    data = data.write("animalSession", "", "table");
    data = data.write("dataPath", "", "table");
    data = data.write("dataCategory", "", "table");
    data = data.write("recordRate", nan, "table");
    data = data.write("LFP_std", nan, "table");

    % fiilling method for each dataset:
    if contains(par, 'IZ') && contains(par, 'sess')
        data = fill_subjectInfo_Zutshi(data);
    elseif contains(par, 'i01_maze')
        data = fill_subjectInfo_hc5(data);
    elseif contains(par, '_2') && length(char(par)) == 13
        data = fill_subjectInfo_LFPEEG(data);
    elseif all(isstrprop(char(par), 'digit')) && length(char(par)) == 4
        data = fill_subjectInfo_Burke_Nick(data);
    elseif contains(par, 'MD') && length(char(par)) == 11
        data = fill_subjectInfo_English(data);
    elseif contains(par, ["Con","Bon","Cor","Fiv","Fra","Mil","Ten","Dud", "Eig"])  % "Eig"  not sure where the probe.
        data = fill_subjectInfo_hc6(data);
    elseif contains(par, "e1") && (length(char(par)) == 14 || length(char(par)) == 15)
        data = fill_subjectInfo_000552(data);
    elseif contains(par, "539") && (length(char(par)) == 17)
        data = fill_subjectInfo_539(data);
    else
        error('Yu:')
    end
          
    
    save([char(outPath), 'data.mat'], 'data','-v7.3');
    RV = 1; return;
end


function data = fill_subjectInfo_hc6(data)

    animalNameSession = data.read("animalNameSession");

    newStr = split(animalNameSession, '_');
    animalName = newStr{1};
    animalSession = newStr{2};   % session is the day for this dataset.
    
    data = data.write("datasetName", "hc6", "table");
    data = data.write("animalName", string(animalName), "table");
    data = data.write("animalSession", string(animalSession), "table");
    
    dataPath = 'C:\1\Dropbox (UFL)\0. Data\(2024 02 08) hc6\data\';
    if animalName == "Con"
        dataPath = [dataPath, '_'];
    end
    dataPath = [dataPath, char(animalName), '\'];
    data = data.write("dataPath", string(dataPath), "table");

    if animalName == "Cor" || ...
       animalName == "Fiv" 
            dataCategory = "CA1/MEC";
    elseif animalName == "Bon" || ...
           animalName == "Fra" || ...
           animalName == "Con" || ...
           animalName == "Dud" || ...
           animalName == "Mil" || ...
           animalName == "Ten"
            dataCategory = "CA1/CA3";
    elseif animalName == "Eig" 
            dataCategory = "unknown";
    else
            error("animalName unrecognized");
    end
    data = data.write("dataCategory", dataCategory, "table");
    
end











