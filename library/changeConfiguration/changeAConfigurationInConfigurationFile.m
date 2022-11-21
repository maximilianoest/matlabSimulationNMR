function configurationFound = changeAConfigurationInConfigurationFile( ...
    configurationKey,newValueAsString,configurationFilePath)

fileID = fopen(configurationFilePath);
txtLines = textscan(fileID,'%s','Delimiter','\n');
fclose(fileID);
configurationFound = false;
for lineNr = 1:numel(txtLines{1,1})
    line = txtLines{1,1}{lineNr};
    if contains(line,sprintf('%s=',configurationKey))
        txtLines{1,1}{lineNr} = sprintf('%s=%s',configurationKey ...
            ,newValueAsString);
        configurationFound = true;
    end
end

if ~configurationFound
    error('changeConfiguration:change%sInConfigurationFile' ...
        ,'"%s" was not part of the configuration file.' ...
        ,upper(configurationKey),upper(configurationKey));
end

stringArray = string(txtLines{1,1});
fileID = fopen(configurationFilePath,'w');
fprintf(fileID,'%s\n',stringArray);
fclose(fileID);

end