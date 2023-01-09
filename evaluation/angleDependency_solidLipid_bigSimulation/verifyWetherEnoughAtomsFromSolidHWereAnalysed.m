clc; clear all; close all; fclose('all');

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile("constants.txt");
saving = 1;

logFileDirectory = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\angleDependency_solidLipid_bigSimulation\";

logFilesNames = [ ...
    "20221112_LogFile_DOPSlipid_20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns" ...
    "20221104_LogFile_PLPClipid_20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns" ...
    "20221108_LogFile_PSMlipid_20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns"] ...
    + ".txt";

lipidNames = {};
estimatedR1 = [];
for logFileNr = 1:length(logFilesNames)
    splittedFileName = strsplit(logFilesNames{logFileNr},"_");
    lipidNames{end+1} = splittedFileName{5}; %#ok<SAGROW>
    
    logFilePath = logFileDirectory + logFilesNames(logFileNr);
    searchPattern = "=> Estimated relaxation rate for calculated atoms:";
    
    fileID = fopen(logFilePath);
    line = fgetl(fileID);
    estimatedR1ForALipid = [];
    while true
        if line == -1
            break;
        end
        
        if contains(line,searchPattern)
            splittedLine = strsplit(line," ");
            estimatedR1ForALipid(end+1) = str2double(splittedLine(end)); %#ok<SAGROW>
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
    saveFigureTo(logFileDirectory,"solidLipidR1" ...
        ,"verificationEnoughAtoms","allLipids");
end





