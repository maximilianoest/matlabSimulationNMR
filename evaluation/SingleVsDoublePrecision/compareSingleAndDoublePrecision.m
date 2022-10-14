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
figPosAndSize = [50 50 1000 700];
configVariables = who;

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
   
    if strcmp(simulationResults.constituent,'lipid')
        arrayPos = 1;
    else
        arrayPos = 2;
    end
    relevantResults.(whichLipid)(arrayPos).constituent = constituent;
    relevantResults.(whichLipid)(arrayPos).deltaTInS ...
        = simulationResults.deltaTInS;
    relevantResults.(whichLipid)(arrayPos).corrFuncFirstOrderSingle = ...
        corrFuncFirstOrderSingle;
    relevantResults.(whichLipid)(arrayPos).corrFuncSecondOrderSingle = ...
        corrFuncSecondOrderSingle;
    relevantResults.(whichLipid)(arrayPos).corrFuncFirstOrderDouble = ...
        corrFuncFirstOrderDouble;
    relevantResults.(whichLipid)(arrayPos).corrFuncSecondOrderDouble = ...
        corrFuncSecondOrderDouble;
    relevantResults.(whichLipid)(arrayPos).simulationDurationInS = ...
        simulationResults.simulationDurationInS;
    relevantResults.(whichLipid)(arrayPos).fibreAnglesTheta = ...
        simulationResults.fibreAnglesTheta;
    relevantResults.(whichLipid)(arrayPos).fibreAnglesPhi = ...
        simulationResults.fibreAnglesPhi;
    relevantResults.(whichLipid)(arrayPos).matlabSimulationDate = ...
        simulationResults.matlabSimulationDate;
end
clearvars('-except',configVariables{:},'relevantResults');
%% analyse data:
corrFuncNames = ["corrFuncFirstOrder" "corrFuncSecondOrder"];

for whichLipid = string(fieldnames(relevantResults)')
    lipidData = relevantResults.(whichLipid);
    fig = initializeFigure('posAndSize',figPosAndSize);
    subplotCounter = 1;
    for datasetNr = 1:length(lipidData)
        constituentData = lipidData(datasetNr);
        for corrFuncCaseNr = 1:length(corrFuncNames)
            initializeSubplot(fig,2,2,subplotCounter);
            switch subplotCounter
                % [left bottom width height]
                case 1
                    left = 0.05;
                    bottom = 0.55;
                    
                case 2
                    left = 0.55;
                    bottom = 0.55;
                case 3
                    left = 0.05;
                    bottom = 0.055;
                case 4
                    left = 0.55;
                    bottom = 0.055;
            end
            set(gca,'Position',[left bottom 0.43 0.42]);
            subplotCounter = subplotCounter + 1;
            
            corrFuncSingle = constituentData.(corrFuncNames( ...
                corrFuncCaseNr)+"Single");
            corrFuncDouble = constituentData.(corrFuncNames( ...
                corrFuncCaseNr)+"Double");    
            difference = corrFuncDouble - corrFuncSingle;
            timeAxis = 0:constituentData.deltaTInS:constituentData ...
                .simulationDurationInS;
            lastPart = 80000;
            
            timeAxis = timeAxis(end-lastPart:end);
            plot(timeAxis,real(corrFuncDouble(end-lastPart:end)));
            plot(timeAxis,imag(corrFuncDouble(end-lastPart:end)));
            plot(timeAxis,real(corrFuncSingle(end-lastPart:end)),'--');
            plot(timeAxis,imag(corrFuncSingle(end-lastPart:end)),'--');
            plot(timeAxis,real(difference(end-lastPart:end)));
            plot(timeAxis,imag(difference(end-lastPart:end)),'--');

            axis([-inf inf -inf inf]);
            
            if corrFuncCaseNr == 1
                ylabel(constituentData.constituent);
            end
            
            if datasetNr == 2
                if corrFuncCaseNr == 1
                    xlabel('First order Corr. Func.');
                else
                    xlabel('Second order Corr. Func.');
                end
            end
                
            if datasetNr == 1 && corrFuncCaseNr == 1
                lgd = legend( ...
                    'Real part double','Imag. part double' ...
                    ,'Real part single','Imag. part single' ...
                    ,'Real part diff. dou. - sing.' ...
                    ,'Imag. part diff. dou. - sing.');
                lgd.Location = 'east';
                title(sprintf('Lipid: %s, Cumsum last part %i steps, double vs. single' ...
                    ,whichLipid,lastPart));
            else
                lgd = legend();
                lgd.Location = 'east';
                lgd.Visible = 'Off';
            end
        end
    end
    if saving
        saveFigureTo(savingDirectory,whichLipid ...
            ,constituentData.matlabSimulationDate ...
            ,'ComparisonCorrFuncSingleDouble');
    end
    
end




