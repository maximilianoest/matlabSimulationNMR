clc; clear all; close all; fclose('all');

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions%s",createFilesepStringArray(5));

densDistribFileName = "20230405_myelinDistributionData";
densData = load(resultsDir + densDistribFileName + ".mat");

densData.avgDensity_locMolAtom = ...
    squeeze(mean(densData.density_timeLocMolAtom,1));
densData.avgFrequency_locMolAtom =  ...
    squeeze(mean(densData.frequency_timeLocMolAtom,1));
densData.avgDiscreteLocations = squeeze(mean( ...
    densData.discreteLocations,1));

densData = rmfield(densData,{'density_timeLocMolAtom' ...
    ,'frequency_timeLocMolAtom','discreteLocations'});

save(resultsDir + densDistribFileName + "_averaged",'-struct','densData');



