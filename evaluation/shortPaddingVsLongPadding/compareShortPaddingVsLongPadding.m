clc; clear all; close all;
% While using the older code I have seen that the correlation function was
% calculated with a zeropadding of 2^nextpow(length+1). In this evaluation
% I want to find out whether it is possible to reduce this to
% 2^nextpow(length) to avoid a big oberhead in calculations.
% The data anaylsed in this script are collected from other simulations in
% order to avoid doubled simulations.

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath('functions'));
addpath(genpath(sprintf('..%s..%stxtFiles' ...
    ,createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf(['..%s..%sRESULTS' ...
    '%scompareShortVsLongZeroPaddingInCorrFunc'] ...
    ,createFilesepStringArray(3));
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
savingDirectory = sprintf('%s%sPlots%s',resultsDir ...
    ,createFilesepStringArray(2));
saving = 1;
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsDir,filesep));

fieldStrengthInT = 3; % main magnetic field strength
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrengthInT;
dipoleDipoleConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
theta = 2;
phi = 1;
nNCases.lipidWater = 'nearestNeighbours8000';
nNCases.lipid = 'nearestNeighbours5500';
averagingRegion = [0.98 1];

% [left bottom width height]
figPosAndSize = [50 50 1000 700];
widthAndHeight = [0.43 0.40];
subFigPosAndSize = [ ...
    0.07 0.55 widthAndHeight ; ...
    0.55 0.55 widthAndHeight ; ...
    0.07 0.06 widthAndHeight ; ...
    0.55 0.06 widthAndHeight];
clearvars widthAndHeight
relevantResults = struct();

configVariables = who;

%% load data and put it into relevantResults:
for fileCounter = 1:length(allMatFilesInResultsDirectory)
    % create file path
    fileName = allMatFilesInResultsDirectory(fileCounter).name;
    filePath = sprintf('%s%s%s',allMatFilesInResultsDirectory( ...
        fileCounter).folder,filesep,fileName);
    
    simulationResults = load(filePath);
    whichLipid = simulationResults.whichLipid;
    constituent = simulationResults.consituent;
    if strcmp(constituent,'lipidWater')
        whichNNCase = nNCases.lipidWater;
    else
        whichNNCase = nNCases.lipid;
    end
    
    try 
        corrFuncZerothOrder = squeeze( ...
            simulationResults.sumCorrelationFunctionsSaverZerothOrderSingle.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        corrFuncFirstOrder = squeeze( ...
            simulationResults.sumCorrelationFunctionSaverFirstOrderSingle.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        corrFuncSecondOrder = squeeze( ...
            simulationResults.sumCorrelationFunctionSaverSecondOrderSingle.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        zeroPaddingMethod = "long";
    catch
        corrFuncZerothOrder = squeeze( ...
            simulationResults.sumCorrFuncZerothOrderShortPadding.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        corrFuncFirstOrder = squeeze( ...
            simulationResults.sumCorrFuncFirstOrderShortPadding.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        corrFuncSecondOrder = squeeze( ...
            simulationResults.sumCorrFuncSecondOrderShortPadding.( ...
            whichNNCase)(theta,phi,:))./simulationResults.atomCounter;
        zeroPaddingMethod = "short";
    end
    
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid) = [];
    end
    relevantResults.(whichLipid)(end+1).constituent = string(constituent);
    relevantResults.(whichLipid)(end).deltaTInS ...
        = simulationResults.deltaTInS;
    relevantResults.(whichLipid)(end).zeroPaddingMethod ...
        = zeroPaddingMethod;
    relevantResults.(whichLipid)(end).corrFuncZerothOrder = ...
        corrFuncZerothOrder;
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
    relevantResults.(whichLipid)(end).atomCounter = ...
        simulationResults.atomCounter;
end

%% analyse data
corrFuncNames = ["corrFuncFirstOrder" ...
    "corrFuncSecondOrder"];

for whichLipid = string(fieldnames(relevantResults)')
    lipidData = relevantResults.(whichLipid);
    [~,indices] = sort([lipidData.zeroPaddingMethod]);
    lipidData = lipidData(indices);
    
    fig = initializeFigure();
    for subFigNr = 1:4
        subFigs(subFigNr) = initializeSubplot(fig,2,2,subFigNr); %#ok<SAGROW>
        subFigs(subFigNr).Position = subFigPosAndSize(subFigNr,:); %#ok<SAGROW>
    end
    legendEntries = {};
    specDensPlotTitleText = ["R1: " "R1: "];
    for datasetNr = 1:length(lipidData)
        dataset = lipidData(datasetNr);
        simulationDurationInS = dataset.simulationDurationInS;
        deltaTInS = dataset.deltaTInS;
        if dataset.constituent == "lipid"
            corrFuncIndex = 1;
            specDensIndex = 2;
            constituentIndex = 1;
            lgdVisibility = "On";
            corrFuncPlotTitleText = sprintf(['Lipid: %s, simDur: %.2d s' ...
                ', deltaT: %.2d s'],whichLipid,simulationDurationInS ...
                ,deltaTInS);
        elseif dataset.constituent == "lipidWater"
            corrFuncIndex = 3;
            specDensIndex = 4;
            constituentIndex = 2;
            lgdVisibility = "Off";
            corrFuncPlotTitleText = sprintf(['SimDur: %.2d s' ...
                ', deltaT: %.2d s'],simulationDurationInS,deltaTInS);
        else
            error("Unknown constituent");
        end
        
        if ~isfield(legendEntries,dataset.constituent)
            legendEntries(1).(dataset.constituent) = {};
            legendEntries(2).(dataset.constituent) = {};
        end
        
        % correlationFunctions
        timeAxis = 0:deltaTInS:simulationDurationInS;
        axes(subFigs(corrFuncIndex)); %#ok<LAXES>
        for corrFuncCaseNr = 1:length(corrFuncNames)
            corrFunc = dataset.(corrFuncNames(corrFuncCaseNr));
            absCorrFunc = abs(corrFunc);
            plot(timeAxis,absCorrFunc);
            legendEntries(1).(dataset.constituent){end+1} = ...
                sprintf("abs(%s), %s padding",corrFuncNames( ...
                corrFuncCaseNr),dataset.zeroPaddingMethod);
        end
        
        switch whichLipid
            case "PLPC"
                 axis([-inf inf -inf 0.5e4])
            case "PSM"
                 axis([-inf inf -inf 1e4])
            case "DOPS"
                 axis([-inf inf -inf 0.5e4])
        end
       
        % spectral densities
        axes(subFigs(specDensIndex)); %#ok<LAXES>
        for corrFuncCaseNr = 1:length(corrFuncNames)
            corrFuncFieldName = corrFuncNames(corrFuncCaseNr);
            corrFunc = dataset.(corrFuncFieldName);
            corrFuncLength = length(corrFunc);
            if corrFuncFieldName == "corrFuncFirstOrder"
                prefactor = -1;
                specDensOrder = "spec. dens. first order";
            elseif corrFuncFieldName == "corrFuncSecondOrder"
                prefactor = -2;
                specDensOrder = "spec. dens. second order";
            end
            specDens = 2*(deltaTInS*cumsum(corrFunc'.*exp(prefactor*1i ...
                *omega0*deltaTInS*(0:corrFuncLength-1))));
            plot(timeAxis,abs(real(specDens)));
            legendEntries(2).(dataset.constituent){end+1} = ...
                sprintf("abs(real(%s)), %s padding",specDensOrder ...
                ,dataset.zeroPaddingMethod);
        end
        
         % R1s in the title of spec dens plots:
        [specDensFirstOrder,specDensSecondOrder] = ...
            calculateSpectralDensities(dataset.corrFuncFirstOrder' ...
            ,dataset.corrFuncSecondOrder',omega0,deltaTInS ...
            ,averagingRegion);
        R1 = calculateR1WithSpectralDensity(specDensFirstOrder ...
            ,specDensSecondOrder,dipoleDipoleConstant);
        
        specDensPlotTitleText(constituentIndex) = sprintf( ...
            "%s %s: %.3f ",specDensPlotTitleText( ...
            constituentIndex),dataset.zeroPaddingMethod,R1);
        
        % titles, legends, axis labels
        lgd = legend(subFigs(corrFuncIndex),legendEntries(1).( ...
            dataset.constituent));
        lgd.Visible = lgdVisibility;
        
        lgd = legend(subFigs(specDensIndex),legendEntries(2).( ...
            dataset.constituent));
        lgd.Visible = lgdVisibility;
        title(subFigs(corrFuncIndex),corrFuncPlotTitleText);
        title(subFigs(specDensIndex),sprintf("%s" ...
            ,specDensPlotTitleText(constituentIndex)))
        ylabel(subFigs(corrFuncIndex),dataset.constituent);
        if dataset.constituent == "lipidWater"
            xlabel(subFigs(corrFuncIndex),"Lag time [s]");
            xlabel(subFigs(specDensIndex),"Upper limit expl. FT")
        end
    end
    
    if saving
        saveFigureTo(savingDirectory,whichLipid ...
            ,dataset.matlabSimulationDate ...
            ,'compareLongVsShortZeroPaddingLength');
    end
    
    
end
