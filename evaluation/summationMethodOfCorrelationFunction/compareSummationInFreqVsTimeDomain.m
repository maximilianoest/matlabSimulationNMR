clc; clear all; close all;
% I had the idea that the summation of the correlation function can already
% be performed in the frequency domain because the ifft is a linear
% transformation. This will be evaluated in this script with the lipid
% water data from the two different summation methods.
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath('functions'));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');
resultsDir = sprintf('..%s..%sRESULTS%scalculateCorrFuncWithShortPadding' ...
    ,createFilesepStringArray(3));
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
savingDirectory = sprintf('%s%sPlots%s',resultsDir ...
    ,createFilesepStringArray(2));
addpath(genpath(resultsDir));
dipoleDipoleConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
fieldStrengthInT = 3; % main magnetic field strength
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrengthInT;
saving = 1;
theta = 2;
phi = 1;
nNCases.lipidWater = 'nearestNeighbours8000';
nNCases.lipid = 'nearestNeighbours5500';
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsDir,filesep));
% [left bottom width height]
figPosAndSize = [50 50 1000 700];
widthAndHeight = [0.43 0.40];
subFigPosAndSize = [ ...
    0.05 0.55 widthAndHeight ; ...
    0.55 0.55 widthAndHeight ; ...
    0.05 0.055 widthAndHeight ; ...
    0.55 0.055 widthAndHeight];

timeDomainSum = "timeDomainSum";
freqDomainSum = "freqDomainSum";
configVariables = who;
datesToSearchFor = ["20221018" "20221019"];

relevantResults = struct();
%% load data and put it into relevantResults:
for fileCounter = 1:length(allMatFilesInResultsDirectory)
    % create file path
    fileName = allMatFilesInResultsDirectory(fileCounter).name;
    filePath = sprintf('%s%s%s',allMatFilesInResultsDirectory( ...
        fileCounter).folder,filesep,fileName);
    
    % load data
    simulationResults = load(filePath);
    whichLipid = simulationResults.whichLipid;
    constituent = simulationResults.consituent;
    if strcmp(constituent,'lipidWater')
        whichCase = nNCases.lipidWater;
    else
        whichCase = nNCases.lipid;
    end
    corrFuncZerothOrder = squeeze( ...
        simulationResults.sumCorrFuncZerothOrderShortPadding.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    corrFuncFirstOrder = squeeze( ...
        simulationResults.sumCorrFuncFirstOrderShortPadding.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    corrFuncSecondOrder = squeeze( ...
        simulationResults.sumCorrFuncSecondOrderShortPadding.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    
    additionalNote = 'No special information';
    if strcmp(constituent,'lipidWater')
        if strcmp(simulationResults.matlabSimulationDate ...
                ,datesToSearchFor(1))
            additionalNote = timeDomainSum;
        elseif strcmp(simulationResults.matlabSimulationDate ...
                ,datesToSearchFor(2))
            additionalNote = freqDomainSum;
        else
            error('Unknown date.');
        end
    end
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid) = [];
    end
    relevantResults.(whichLipid)(end+1).constituent = constituent;
    relevantResults.(whichLipid)(end).deltaTInS ...
        = simulationResults.deltaTInS;
    relevantResults.(whichLipid)(end).additionalNote = additionalNote;
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
clearvars('-except',configVariables{:},'relevantResults');


%% analyse data
clc
corrFuncNames = ["corrFuncZerothOrder" "corrFuncFirstOrder" ...
    "corrFuncSecondOrder"];
deltaTInSForSumMethodCompare = 3e-13;

for whichLipid = string(fieldnames(relevantResults)')
    lipidData = relevantResults.(whichLipid);
    lipidWaterData = getLipidWaterData(lipidData);
    lipidWaterDataSumMethodCompare = getLipidWaterDataWithSameDeltaT( ...
        lipidWaterData,deltaTInSForSumMethodCompare);
    fig = initializeFigure();
    for corrFuncCaseNr = 1:length(corrFuncNames)
        initializeSubplot(fig,2,2,corrFuncCaseNr);
        legendEntries = {};
        for datasetNr = 1:length(lipidWaterDataSumMethodCompare)
            dataset = lipidWaterDataSumMethodCompare(datasetNr);
            timeAxis = 0:dataset.deltaTInS:dataset.simulationDurationInS;
            corrFunc = dataset.(corrFuncNames(corrFuncCaseNr));
            plot(timeAxis,real(corrFunc));
            legendEntries{end+1} = sprintf('real: %s, atomCount: %i' ...
                ,dataset.additionalNote,dataset.atomCounter); %#ok<SAGROW>
            plot(timeAxis,imag(corrFunc));
            legendEntries{end+1} = sprintf('imag: %s, atomCount: %i' ...
                ,dataset.additionalNote,dataset.atomCounter); %#ok<SAGROW>
        end
        
        switch corrFuncNames(corrFuncCaseNr)
            case "corrFuncZerothOrder"
                axis([-inf inf -inf 2e3]);
            case "corrFuncFirstOrder"
                axis([-inf inf -inf 1e3]);
            case "corrFuncSecondOrder"
                axis([-inf inf -inf 2e3]);
        end
        ylabel(corrFuncNames(corrFuncCaseNr));
        
        lgd  = legend(legendEntries);
        if corrFuncCaseNr == 1
            title(sprintf( ...
                'Lipid: %s, dT in sec: %.2d, sim. Dur.: %.2d' ...
                ,whichLipid,deltaTInSForSumMethodCompare ...
                ,dataset.simulationDurationInS));
            
        end   
    end
    if saving
        saveFigureTo(savingDirectory,whichLipid ...
            ,dataset.matlabSimulationDate ...
            ,'ComparisonSummationInFreqVsInTimeDomain');
    end
end
    


function lipidWaterDataWithSameDeltaT = getLipidWaterDataWithSameDeltaT( ...
        lipidWaterData,givenDeltaTInS)
    fieldNames = fieldnames(lipidWaterData)';
    fieldNames{2,1} = {};
    lipidWaterDataWithSameDeltaT = struct(fieldNames{:});
    for datasetNr = 1:length(lipidWaterData)
        dataset = lipidWaterData(datasetNr);
        if dataset.deltaTInS == givenDeltaTInS
            lipidWaterDataWithSameDeltaT(end+1) = dataset; %#ok<AGROW>
        end
    end
end
   
function lipidWaterData = getLipidWaterData(lipidData)
fieldNames = fieldnames(lipidData)';
fieldNames{2,1} = {};
lipidWaterData = struct(fieldNames{:});
for datasetNr = 1:length(lipidData)
    dataset = lipidData(datasetNr);
    if string(dataset.constituent) == "lipidWater"
        lipidWaterData(end+1) = dataset; %#ok<AGROW>
    end
end

end
    
