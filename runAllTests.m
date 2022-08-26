clc; clear all; close all;

allDirectories = genpath('testing\');
allDirectories = strsplit(allDirectories,';');
allDirectories = allDirectories(2:end-1);

for directoryNr = 1:length(allDirectories)
    directory = allDirectories{directoryNr};
    allFiles = dir(sprintf('%s%s%s',directory,filesep,'*m'));
    for fileNr = 1:length(allFiles)
        fileName = allFiles(fileNr).name;
        if ~isempty(fileName) 
            runtests(sprintf('%s%s%s',directory,filesep,fileName));
        end
    end
end
