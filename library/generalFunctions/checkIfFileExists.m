function checkIfFileExists(filePath,logFilePath)

if ~exist(filePath,'file')
    logMessage(sprintf('FILE: %s CANNOT BE FOUND!',filePath) ...
        ,path2LogFile, false);
    error('checkIfFileExists:DataFileNotFound' ...
        ,'File (%s) cannot be found. See log file.',filePath);
else
    logMessage(sprintf('File (%s) exists and will be loaded soon.' ...
        ,filePath),logFilePath,false);
end

end
