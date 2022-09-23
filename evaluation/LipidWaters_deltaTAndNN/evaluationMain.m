clc; clear all; close all;

addpath(genpath('../../library'));
resultsPath = '../../RESULTS/findDeltaTAndNN_lipidWater_longSimulationTime';
addpath(genpath(resultsPath));

allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsPath,filesep));

dataForEvaluation = struct();

for fileCounter = 1:size(allMatFilesInResultsDirectory,1)
fileName = allMatFilesInResultsDirectory(fileCounter,1).name;
if strcmp(fileName,"")
    error("Matfile name not found");
end
results = load(sprintf('%s%s%s',resultsPath,filesep,fileName));
simulationParameters = struct( ...
    'deltaT', results.deltaTInS ...
    ,'simulationDuration', results.simulationDurationInS ...
    ,'sumCorrFuncZerothOrder', results.sumCorrelationFunctionsSaverZerothOrder ...
    ,'sumCorrFuncFirstOrder', results.sumCorrelationFunctionSaverFirstOrder ...
    ,'sumCorrFuncSecondOrder', results.sumCorrelationFunctionSaverSecondOrder);

whichLipid = results.whichLipid;
if ~isfield(dataForEvaluation,whichLipid)
    dataForEvaluation.(whichLipid) = {};
end
dataForEvaluation.(whichLipid){end+1} = simulationParameters; 

    
end




