function [dataDirectory,resultsFileSavingPath,logFilePath] = ...
    setUpDependenciesBasedOn(configuration)

fileName = configuration.fileName;
whichLipid = getLipidNameFromFileName(fileName);
constituent = getConstituentFromFileName(fileName);
combinedName = [whichLipid constituent];

startingDate = datestr(date,'yyyymmdd');

if runningOnServer()
    dataDirectory = configuration.dataDirectoryOnServer;
else
    dataDirectory = configuration.dataDirectoryOnLocalMachine;
end
dataDirectory = sprintf('%s%s%s',dataDirectory,whichLipid,filesep);

resultsDirectory = sprintf('%s%s%s%s%s%s',pwd,filesep, ...
    'RESULTS',filesep,configuration.scriptFileToRun,filesep);
if ~isfolder(resultsDirectory)
    mkdir(resultsDirectory);
end

resultsFileSavingPath = sprintf('%s%s_Results_%s.mat',resultsDirectory ...
    ,startingDate,combinedName);
logFilePath = sprintf('%s%s_LogFile_%s.txt',resultsDirectory ...
    ,startingDate,combinedName); 

createNewLogFileAndRenameOldOne(logFilePath);

logMessage(sprintf(['The following directories are used: \n' ...
    '    Data is taken from directory: %s \n' ...
    '    File name is: %s \n' ...
    '    Results will be saved under: %s \n' ...
    '    Log file path: %s \n'],dataDirectory,fileName ...
    ,resultsFileSavingPath,logFilePath),logFilePath);

end