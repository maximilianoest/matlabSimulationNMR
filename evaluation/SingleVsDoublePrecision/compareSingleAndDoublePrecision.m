clc; clear all; close all;
fieldStrengthInT = 3; % main magnetic field strength
addpath(genpath('../../library'));
addpath(genpath('../../txtFiles'));
constants = readConstantsFile('constants.txt');
resultsDir = '../../RESULTS/validateMethodsForComplexCorrFunc';
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
savingDirectory = sprintf('%s%sPlots%s',resultsDir,filesep,filesep);
addpath(genpath(resultsDir));
dipoleDipoleConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
saving = 1;
theta = 1;
phi = 1;
nNCases.lipidWater = 'nearestNeighbours8000';
nNCases.lipid = 'nearestNeighbours5500';
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsDir,filesep));
figPosAndSize = [50 50 1400 700];

relevantResults = struct();

%% load data and put it into relevantResults:
for fileCounter = 1:length(allMatFilesInResultsDirectory)
    % create file path
    fileName = allMatFilesInResultsDirectory(fileCounter).name;
    filePath = sprintf('%s%s%s',resultsDir,filesep,fileName);
    
    % load data
    simulationResults = load(filePath);
    whichLipid = simulationResults.whichLipid;
    constituent = simulationResults.consituent;
    if strcmp(constituent,'lipidWater')
        whichCase = nNCases.lipidWater;
    else
        whichCase = nNCases.lipid;
    end
    corrFuncFirstOrderSingle = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverFirstOrderSingle.( ...
        whichCase)(theta,phi,:));
    corrFuncSecondOrderSingle = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverSecondOrderSingle.( ...
        whichCase)(theta,phi,:));
    corrFuncFirstOrderDouble = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverFirstOrder.( ...
        whichCase)(theta,phi,:));
    corrFuncSecondOrderDouble = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverSecondOrder.( ...
        whichCase)(theta,phi,:));
    
    % fill a struct
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid).constituent = constituent;
    else
        relevantResults.(whichLipid)(end+1).constituent = constituent;
    end
    
    relevantResults.(whichLipid)(end).deltaTInS ...
            = simulationResults.deltaTInS;
    relevantResults.(whichLipid)(end).corrFuncFirstOrderSingle = ...
        corrFuncFirstOrderSingle;
    relevantResults.(whichLipid)(end).corrFuncSecondOrderSingle = ...
        corrFuncSecondOrderSingle;
    relevantResults.(whichLipid)(end).corrFuncFirstOrderDouble = ...
        corrFuncFirstOrderDouble;
    relevantResults.(whichLipid)(end).corrFuncSecondOrderDouble = ...
        corrFuncSecondOrderDouble;
    relevantResults.(whichLipid)(end).simulationDurationInS = ...
        simulationResults.simulationDurationInS;
    relevantResults.(whichLipid)(end).fibreAnglesTheta = ...
        simulationResults.fibreAnglesTheta;
    relevantResults.(whichLipid)(end).fibreAnglesPhi = ...
        simulationResults.fibreAnglesPhi;
    relevantResults.(whichLipid)(end).matlabSimulationDate = ...
        simulationResults.matlabSimulationDate;
end

%% analyse data:
