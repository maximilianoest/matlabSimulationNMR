clc; clear all; close all;
addpath(genpath('library'));
configurationFilePath = sprintf('txtFiles%sconfigMain.txt',filesep);

lipids.PLPC.fileNamesToAnalyse = [...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime50ns" ...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns"];
lipids.PLPC.atomsCount = 10000;

lipids.PSM.fileNamesToAnalyse = [ ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime50ns" ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime50ns" ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns"];
lipids.PSM.atomsCount = 10000;

lipids.DOPS.fileNamesToAnalyse = [ ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime50ns" ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime50ns" ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt03ps_simTime25ns"];
lipids.DOPS.atomsCount = 10000;

for whichLipid = string(fieldnames(lipids))'
    % create random squence of atoms for one atom
    if ~strcmp(whichLipid,'PLPC')
        randomSequenceOfAtoms = randperm(lipids.(whichLipid).atomsCount);
        save(sprintf('matFiles%srandomSequenceOfAtoms.mat',filesep) ...
        ,'randomSequenceOfAtoms');
    end
    for newFileName = lipids.(whichLipid).fileNamesToAnalyse
        % change file name in configuration file
        change_fileName_inConfigurationFileTo(newFileName ...
            ,configurationFilePath);
        % run main simulation with this configuration
        simulationMain;
    end
end

