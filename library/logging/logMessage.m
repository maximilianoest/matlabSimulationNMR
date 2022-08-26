function logMessage(logInformation,path2LogFile,varargin)
configuration = readConfigurationFile('configMain.txt');
fileId = fopen(path2LogFile,'a');
if fileId == -1
    error('Cannot open log file.')
end

if numel(varargin) < 1
    printWithDate = true;
else
    printWithDate = varargin{1};
end

if printWithDate
    fprintf(fileId, '%s: %s\n', datestr(now,0),logInformation);
else
    fprintf(fileId, '%s\n',logInformation);
end
fclose(fileId);

if configuration.printToCommandWindow
    if printWithDate
        fprintf('%s: %s\n', datestr(now,0),logInformation);
    else
        fprintf('  %s\n',logInformation);
    end
end

end
