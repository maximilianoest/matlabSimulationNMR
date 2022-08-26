function checkIfFileExists(filePath,logFilePath)

if ~exist(filePath,'file')
    logMessage(sprintf('FILE CANNOT BE FOUND: %s',filePath) ...
        ,logFilePath);
    error('checkIfFileExists:DataFileNotFound' ...
        ,'File cannot be found %s. See log file.',filePath);
end

end
