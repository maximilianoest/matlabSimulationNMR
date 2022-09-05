clc; clear all; close all;
addpath(genpath('library'));

configurationFilePath = 'txtFiles\configMain.txt';

lipids.PLPC.fileNamesToAnalyse = [...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime10ns" ...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns"];
lipids.PLPC.atomsCount = 10000;

lipids.PSM.fileNamesToAnalyse = [ ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime10ns" ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns"];
lipids.PSM.atomsCount = 10000;

lipids.DOPS.fileNamesToAnalyse = [ ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime10ns" ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns"];
lipids.DOPS.atomsCount = 10000;

for whichLipid = string(fieldnames(lipids))'
    % create random squence of atoms for one atom
    randomSequenceOfAtoms = randperm(lipids.(whichLipid).atomsCount);
    save('matFiles/randomSequenceOfAtoms.mat','randomSequenceOfAtoms');
    for newFileName = lipids.(whichLipid).fileNamesToAnalyse
        % change file name in configuration file
        change_fileName_inConfigurationFileTo(newFileName ...
            ,configurationFilePath);
        % run main simulation with this configuration
        simulationMain;
    end
end

