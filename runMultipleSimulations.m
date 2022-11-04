clc; clear all; close all;
addpath(genpath('library'));
addpath(genpath('matFiles'));

scriptsToRun = ["angleDependency_solidLipid_bigSimulation"]; %#ok<NBRAK>
nearestNeighboursCases.lipids = "5500";
nearestNeighboursCases.lipidWaters = "8000";
configurationFilePath = sprintf('txtFiles%sconfigMain.txt',filesep);

% lipids.PLPC_lipidWater.fileNamesToAnalyse = [...
%     "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime10ns" ...
%     "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns"];
% lipids.PLPC_lipidWater.atomsCount = 10000;

% lipids.PLPC_lipidWater.fileNamesToAnalyse = [...
%     "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
% lipids.PLPC_lipidWater.atomsCount = 10000;
% 
% lipids.PSM_lipidWater.fileNamesToAnalyse = [ ...
%     "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
% lipids.PSM_lipidWater.atomsCount = 10000;
% 
% lipids.DOPS_lipidWater.fileNamesToAnalyse = [ ...
%     "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
% lipids.DOPS_lipidWater.atomsCount = 10000;

lipids.PLPC_lipid.fileNamesToAnalyse = [ ...
    "20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns.mat"];
lipids.PLPC_lipid.atomsCount = 8000;

lipids.PSM_lipid.fileNamesToAnalyse = [ ...
    "20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns.mat"];
lipids.PSM_lipid.atomsCount = 7900;

lipids.DOPS_lipid.fileNamesToAnalyse = [ ...
    "20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns.mat"];
lipids.DOPS_lipid.atomsCount = 7700;

checkIfFilesForMultiSimExist(lipids);

for lipidToAnalyse = string(fieldnames(lipids))'
    splittedFieldName = strsplit(lipidToAnalyse,'_');
    lipid = splittedFieldName{1};
    consituent = splittedFieldName{2};
    if strcmp(consituent,'lipid')
        change_nearestNeighbourCases_inConfigurationFileTo( ...
            nearestNeighboursCases.lipids,configurationFilePath);
    else
        change_nearestNeighbourCases_inConfigurationFileTo( ...
            nearestNeighboursCases.lipidWaters,configurationFilePath);
    end
    
    for newFileName = lipids.(lipidToAnalyse).fileNamesToAnalyse
        % change file name in configuration file
        change_fileName_inConfigurationFileTo(newFileName ...
            ,configurationFilePath);
        for scriptNr = 1:length(scriptsToRun)
            newScriptName = scriptsToRun(scriptNr);
            change_scriptName_inConfigurationFileTo(newScriptName ...
                ,configurationFilePath);
            % run main simulation with this configuration
            multiSimVariables = who;
            simulationMain;
            clearvars('-except',multiSimVariables{:});
        end
    end
end

function checkIfFilesForMultiSimExist(lipids)
dataDirectory = "/daten/a/Relaxation/Lipids/simulation_max/";
for lipidToAnalyse = string(fieldnames(lipids))'
    fileNamesCount = length(lipids.(lipidToAnalyse).fileNamesToAnalyse);
    whichLipid = strsplit(lipidToAnalyse,'_');
    whichLipid = whichLipid{1};
    for fileNR = 1:fileNamesCount
        fileName = lipids.(lipidToAnalyse).fileNamesToAnalyse(fileNR);
        filePath = sprintf('%s%s%s%s.mat',dataDirectory,whichLipid ...
            ,filesep,fileName);
        atomsCount = lipids.(lipidToAnalyse).atomsCount;
        fprintf('Number of hydrogen atoms in dataset: %i\n',atomsCount);
        if ~exist(filePath,'file')
            warning('File %s does not exist. Script will crash then.' ...
                ,filePath);
        end
    end
end
end
