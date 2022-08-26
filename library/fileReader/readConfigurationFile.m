function [configuration] = readConfigurationFile(configurationFilePath)

configuration = {};
configFileId = fopen(configurationFilePath);
data = textscan(configFileId, '%s %s','Delimiter','=');
fclose(configFileId);

for element = 1:size(data{1},1)
    configurationName = data{1}{element};
    try
        configurationInformation = data{2}{element};
        if isnan(str2double(configurationInformation))
            configuration.(configurationName) = configurationInformation;
        else
            configuration.(configurationName) = ...
                str2double(configurationInformation);
        end 
    catch 
        warning('configuration:readConfigurationFile', ...
            ['The value for the variable name "%s" is non of the ' ...
            'specified data types the configuration file'], ...
            configurationName);
    end
end

end
