function [dataDirectory,resultsFileSavingPath,logFilePath] = ...
    setUpDependenciesBasedOn(configuration)

fileName = configuration.fileName;
lipidName = getLipidNameFromFileName(fileName);
constituent = getConstituentFromFileName(fileName);
combinedName = [lipidName constituent];

startingDate = datestr(date,'yyyymmdd');

if runningOnServer()
    dataDirectory = configuration.dataDirectoryOnServer;
else
    dataDirectory = configuration.dataDirectoryOnLocalMachine;
end

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

end