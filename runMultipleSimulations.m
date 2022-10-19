clc; clear all; close all;
addpath(genpath('library'));
addpath(genpath('matFiles'));
% data was created but then the script crashed

scriptsToRun = ["validateMethodsForComplexCorrFunc"]; %#ok<NBRAK>
nearestNeighboursCases.lipids = "5500;2600;300";
nearestNeighboursCases.lipidWaters = "8000;5000;2000";
configurationFilePath = sprintf('txtFiles%sconfigMain.txt',filesep);


lipids.PLPC_lipidWater.fileNamesToAnalyse = [...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns" ...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
lipids.PLPC_lipidWater.atomsCount = 10000;

lipids.PSM_lipidWater.fileNamesToAnalyse = [ ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
lipids.PSM_lipidWater.atomsCount = 10000;

lipids.DOPS_lipidWater.fileNamesToAnalyse = [ ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"];
lipids.DOPS_lipidWater.atomsCount = 10000;

% lipids.PLPC_lipid.fileNamesToAnalyse = [ ...
%     "20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime500ns"]; % ...
% %     "20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime500ns" ...
% %     "20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt6ps_simTime500ns"];
% lipids.PLPC_lipid.atomsCount = 8000;
% 
% lipids.PSM_lipidWater.fileNamesToAnalyse = [ ...
%     "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns"];
% lipids.PSM_lipidWater.atomsCount = 10000;
% 
% lipids.PSM_lipid.fileNamesToAnalyse = [ ...
%     "20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime500ns"]; % ...
% %     "20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime500ns" ...
% %     "20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt6ps_simTime500ns"];
% lipids.PSM_lipid.atomsCount = 7900;
% 
% lipids.DOPS_lipidWater.fileNamesToAnalyse = [ ...
%     "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns"];
% lipids.DOPS_lipidWater.atomsCount = 10000;
% 
% lipids.DOPS_lipid.fileNamesToAnalyse = [ ...
%     "20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime500ns"]; % ...
% %     "20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime500ns" ...
% %     "20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt6ps_simTime500ns"];
% lipids.DOPS_lipid.atomsCount = 7700;

checkIfFilesForMultiSimExist(lipids);

for lipidToAnalyse = string(fieldnames(lipids))'
    splittedFieldName = strsplit(lipidToAnalyse,'_');
    lipid = splittedFieldName{1};
    consituent = splittedFieldName{2};
    % random sequence of atoms are used from the nearest neighbours
    % analysis of the lipids.
    load(sprintf('randomSequenceOfAtoms_%s%s',lipid,consituent));
    save(sprintf('matFiles%srandomSequenceOfAtoms.mat',filesep) ...
        ,'randomSequenceOfAtoms');
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

clear all;
createNewDataSetsFromOthers();

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
                ,fileName);
        end
    end
end
end
