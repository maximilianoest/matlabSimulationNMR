function change_fileName_inConfigurationFileTo(newFileName ...
    ,configurationFilePath)

fileID = fopen(configurationFilePath);
txtLines = textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);
configurationFound = false;
for lineNr = 1:numel(txtLines{1,1})
    line = txtLines{1,1}{lineNr};
    if contains(line,'fileName=')
        txtLines{1,1}{lineNr} = sprintf('fileName=%s',newFileName);
        configurationFound = true;
    end
end

if ~configurationFound
    error('changeConfiguration:changeFileNameInConfigurationFile' ...
        ,'"fileName" was not part of the configuration file');
end

stringArray = string(txtLines{1,1});
fileID = fopen(configurationFilePath,'w');
fprintf(fileID,'%s\n',stringArray);
fclose(fileID);

end