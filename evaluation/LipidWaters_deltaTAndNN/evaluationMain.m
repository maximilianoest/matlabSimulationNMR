clc; clear all; close all;

addpath(genpath('functions'));
%% load results and determine R1
fieldStrengthInT = 3; % main magnetic field strength
addpath(genpath('../../library'));
addpath(genpath('../../txtFiles'));
constants = readConstantsFile('constants.txt');
resultsPath = '../../RESULTS/findDeltaTAndNN_lipidWater_longSimulationTime';
savingDirectory = sprintf('%s%sPlots%s',resultsPath,filesep,filesep);
addpath(genpath(resultsPath));
fieldstrength = 3; % fieldstrength in tesla
dipoleDipoleConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
% duration after which the correlation function is cutted
differentCorrelationFunctionDurations = [25e-9 10e-9 5e-9 1e-9];
histogramPlotRegion = [0.7 1];
binWidths = [1 3];
saving = 1;
theta = 1;
phi = 1;
nNCase = "nearestNeighbours8000";
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsPath,filesep));
figPosAndSize = [50 50 1400 700];

%% create relevant results
% 1. investigate whether the correlation function gets under zero?

relevantResults = struct();
for fileCounter = 1:size(allMatFilesInResultsDirectory,1)
    fileName = allMatFilesInResultsDirectory(fileCounter,1).name;
    if strcmp(fileName,"")
        error("Matfile not found.");
    end
    filePath = sprintf('%s%s%s',resultsPath,filesep,fileName);
    simulationResults = load(filePath);
    whichLipid = simulationResults.whichLipid;
    [corrFuncFirstOrder,corrFuncSecondOrder] =  ...
        loadFirstAndSecondOrderCorrelationFunctionsFromNNSimulation( ...
        filePath,nNCase);
    
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid).deltaTInS  ...
            = simulationResults.deltaTInS;
    else
        relevantResults.(whichLipid)(end+1).deltaTInS ...
            = simulationResults.deltaTInS;
    end
    relevantResults.(whichLipid)(end).corrFuncFirstOrder = ...
        corrFuncFirstOrder;
    relevantResults.(whichLipid)(end).corrFuncSecondOrder = ...
        corrFuncSecondOrder;
    relevantResults.(whichLipid)(end).simulationDurationInS = ...
        simulationResults.simulationDurationInS;
    relevantResults.(whichLipid)(end).fibreAnglesTheta = ...
        simulationResults.fibreAnglesTheta;
    relevantResults.(whichLipid)(end).fibreAnglesPhi = ...
        simulationResults.fibreAnglesPhi;
    relevantResults.(whichLipid)(end).matlabSimulationDate = ...
        simulationResults.matlabSimulationDate;
end

%% plot distributions to find whether correlation function is at 0
variationsOfCuttedCorrFunc = [40e-9 30e-9 25e-9 20e-9 15e-9 10e-9];
plottedFigures = plotCorrFuncDistributionsAtZero(relevantResults ...
    ,variationsOfCuttedCorrFunc,theta,phi,histogramPlotRegion ...
    ,binWidths,figPosAndSize);
if saving
    for figureNr = 1:length(plottedFigures)
        set(0, 'currentfigure', plottedFigures(figureNr).figure);
        saveFigureTo(savingDirectory ...
            ,plottedFigures(figureNr).whichLipid ...
            ,plottedFigures(figureNr).matlabSimulationDate ...
            ,'distributionsForDifferentCorrFuncDurations');
    end
end

%% smooth the curve an look at the correlation function
smoothingStep = 1500;
plottedFigures = plotSmoothedCorrelationFunctions(relevantResults ...
    ,smoothingStep,figPosAndSize,theta,phi);

if saving
    for figureNr = 1:length(plottedFigures)
        set(0, 'currentfigure', plottedFigures(figureNr).figure);
        saveFigureTo(savingDirectory ...
            ,plottedFigures(figureNr).whichLipid ...
            ,plottedFigures(figureNr).matlabSimulationDate ...
            ,'smoothedCorrelationFunctionForDifferentSimConfig');
    end
end

%% plot R1 over nearest neighbours
corrFuncDuration = 20e-9;
datasetsToPlot = [1:3];
nearestNeighbourCase = 'nearestNeighbours8000';
omega0 = fieldstrength*constants.gyromagneticRatioOfHydrogenAtom;
relevantResultsForAllNNCases = getRelevantResultsForAllNNCases( ...
    allMatFilesInResultsDirectory,resultsPath,omega0);

relevantResultsForAllNNCases ...
    = shortenCorrelationFunctionsForAllLipidsAndNNCases( ...
    relevantResultsForAllNNCases,corrFuncDuration);

relevantResultsForAllNNCases = calculateR1ForDifferentLipids( ...
    relevantResultsForAllNNCases);
plottedFigures = plotNNDependentR1RatesAndTheirCorrFunc( ...
    relevantResultsForAllNNCases,datasetsToPlot,nearestNeighbourCase ...
    ,figPosAndSize);

if saving
    for figureNr = 1:length(plottedFigures)
        set(0, 'currentfigure', plottedFigures(figureNr).figure);
        saveFigureTo(savingDirectory ...
            ,plottedFigures(figureNr).whichLipid ...
            ,plottedFigures(figureNr).matlabSimulationDate ...
            ,'R1AndCorrFuncDistributions');
    end
end
