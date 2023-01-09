clc; clear all; close all; fclose('all');

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile("constants.txt");
saving = 1;

logFileDirectory = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\angleDependency_lipidWater_bigSimulation\";

logFilesNames = [ ...
    "20221101_LogFile_DOPSwater_20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221025_LogFile_PLPCwater_20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221028_LogFile_PSMwater_20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"] ...
    + ".txt";

lipidNames = {};
estimatedR1 = [];
for logFileNr = 1:length(logFilesNames)
    logFilePath = logFileDirectory + logFilesNames(logFileNr);
    splittedFileName = strsplit(logFilesNames{logFileNr},"_");
    lipidNames{end+1} = splittedFileName{5}; %#ok<SAGROW>
    
    searchPattern = "=> nearestNeighbours8000";
    fileID = fopen(logFilePath);
    line = fgetl(fileID);
    estimatedR1ForALipid = [];
    while true
        if line == -1
            break;
        end
        
        if contains(line,searchPattern)
            line = fgetl(fileID);
            splittedLine = strsplit(line," ");
            estimatedR1ForALipid(end+1) = str2double(splittedLine{end}); %#ok<SAGROW>
        end
        line = fgetl(fileID);
    end
    estimatedR1(logFileNr,:) = estimatedR1ForALipid; %#ok<SAGROW>
    fclose(fileID);
end

initializeFigure();
atomAxis = 1:size(estimatedR1,2);
for lipidNr = 1:size(estimatedR1,1)
    plot(atomAxis,estimatedR1(lipidNr,:));
end

legend(lipidNames);
xlabel("Number of atoms");
ylabel("Realxation rate R1 $[Hz]$");

if saving
    saveFigureTo(logFileDirectory,"lipidWaterR1" ...
        ,"verificationEnoughAtoms","allLipids");
end





