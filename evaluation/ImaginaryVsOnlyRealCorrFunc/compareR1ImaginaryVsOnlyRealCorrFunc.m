clc; clear all; close all;
addpath(genpath(sprintf('..%s..%slibrary',createFilesepStringArray(2))));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');
resultsDir = sprintf('..%s..%sRESULTS%svalidateMethodsForComplexCorrFunc' ...
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
saving = 0;
theta = 1;
phi = 1;
nNCases.lipidWater = 'nearestNeighbours8000';
nNCases.lipid = 'nearestNeighbours5500';
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%ssingleDoublePrecisionComparison%s*.mat',resultsDir ...
    ,filesep,filesep));
% [left bottom width height]
figPosAndSize = [50 50 1000 700];
widthAndHeight = [0.43 0.40];
subFigPosAndSize = [ ...
    0.05 0.55 widthAndHeight ; ...
    0.55 0.55 widthAndHeight ; ...
    0.05 0.055 widthAndHeight ; ...
    0.55 0.055 widthAndHeight];
configVariables = who;

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
        simulationResults.sumCorrelationFunctionsSaverZerothOrderSingle.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    corrFuncFirstOrder = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverFirstOrderSingle.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    corrFuncSecondOrder = squeeze( ...
        simulationResults.sumCorrelationFunctionSaverSecondOrderSingle.( ...
        whichCase)(theta,phi,:))./simulationResults.atomCounter;
    
    if strcmp(simulationResults.constituent,'lipid')
        arrayPos = 1;
    else
        arrayPos = 2;
    end
    relevantResults.(whichLipid)(arrayPos).constituent = constituent;
    relevantResults.(whichLipid)(arrayPos).deltaTInS ...
        = simulationResults.deltaTInS;
    relevantResults.(whichLipid)(arrayPos).corrFuncZerothOrder = ...
        corrFuncZerothOrder;
    relevantResults.(whichLipid)(arrayPos).corrFuncFirstOrder = ...
        corrFuncFirstOrder;
    relevantResults.(whichLipid)(arrayPos).corrFuncSecondOrder = ...
        corrFuncSecondOrder;
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


%% analyse data
corrFuncNames = ["corrFuncFirstOrder" "corrFuncSecondOrder"];

for whichLipid = string(fieldnames(relevantResults)')
    lipidData = relevantResults.(whichLipid);
    fig = initializeFigure('posAndSize',figPosAndSize);
    subplotCounter = 1;
    for datasetNr = 1:length(lipidData)
        dataset = lipidData(datasetNr);
        initializeSubplot(fig,2,2,subplotCounter);
        set(gca,'Position',subFigPosAndSize(subplotCounter,:));
        
        deltaT = dataset.deltaTInS;
        simulationDuration = dataset.simulationDurationInS;
        for corrFuncCaseNr = 1:length(corrFuncNames)
            corrFunc = double(dataset.(corrFuncNames(corrFuncCaseNr)));
            if strcmp(dataset.constituent,'lipidWater')
                corrFunc = corrFunc(1:40000);
            end
            corrFuncLength = length(corrFunc);
            timeAxis = 0:deltaT:(corrFuncLength-1)*deltaT;
            quotientCorrFunc = imag(corrFunc)./real(corrFunc);
        
            plt = plot(timeAxis,quotientCorrFunc,'--');
            plot(timeAxis,real(corrFunc),'Color',plt.Color);
            
        end
        axis([-inf inf -1.5e4 1.5e4]);
        
        ylabel('Correlation function');
        if datasetNr == 1
            lgd = legend( ...
                'FO imag/real' ...
                ,'FO real' ...
                ,'SO imag/real' ...
                ,'SO real');
            lgd.Location = 'northeast';
            title(sprintf('Lipid: %s, influence of imaginary part' ...
                ,whichLipid));
        else
            lgd = legend();
            lgd.Visible = 'Off';
            title(sprintf('Lipid water: %s, influence of imaginary part' ...
                ,whichLipid));
        end
        
        subplotCounter = subplotCounter + 1;
        initializeSubplot(fig,2,2,subplotCounter);
        set(gca,'Position',subFigPosAndSize(subplotCounter,:));
        
        r1WholeCorrFunc = 0;
        r1OnlyReal = 0;
        for corrFuncCaseNr = 1:length(corrFuncNames)
            corrFunc = double(dataset.(corrFuncNames(corrFuncCaseNr)));
            if strcmp(dataset.constituent,'lipidWater')
                corrFunc = corrFunc(1:40000);
            end
            
            corrFuncLength = length(corrFunc);
            switch corrFuncNames(corrFuncCaseNr) 
                case 'corrFuncFirstOrder'
                    preFactor = 1;
                case 'corrFuncSecondOrder'
                    preFactor = 2;
                otherwise 
                    error('Not implemented.');
            end
            
            spectralDensityWholeCorrFunc = 2*(deltaT*cumsum(corrFunc' ... 
                .*exp(-preFactor*1i*omega0*deltaT*(0:corrFuncLength-1))));
            spectralDensityOnlyReal = 2*(deltaT*cumsum(real(corrFunc') ... 
                .*exp(-preFactor*1i*omega0*deltaT*(0:corrFuncLength-1))));
            r1WholeCorrFunc = r1WholeCorrFunc ...
                + dipoleDipoleConstant*3/2*abs(real(mean( ...
                spectralDensityWholeCorrFunc(end-5000:end))));
            r1OnlyReal = r1OnlyReal  ...
                + dipoleDipoleConstant*3/2*abs(real(mean( ...
                spectralDensityOnlyReal(end-5000:end))));
            
            plt = plot(timeAxis,real(spectralDensityWholeCorrFunc));
            plot(timeAxis,real(spectralDensityOnlyReal),':');
            
        end
        
        ylabel('Real(Spectral density)');
        xlabel('lag time [s]');
        title(sprintf(['Corr Func: R$_1$ = %.2f' ...
            ', Real(Corr. Func.): R$_1$ = %.2f'],r1WholeCorrFunc ...
            ,r1OnlyReal));
        if datasetNr == 1
            lgd = legend( ...
                'FO corrFunc' ...
                ,'FO real(corrFunc)' ...
                ,'SO corrFunc' ...
                ,'SO real(corrFunc)');
            lgd.Location = 'northeast';
        else
            lgd = legend();
            lgd.Visible = 'Off';
        end
        subplotCounter = subplotCounter + 1;
    end
    
    if saving
        saveFigureTo(savingDirectory,whichLipid ...
            ,dataset.matlabSimulationDate ...
            ,'realVsWholeCorrFunc');
    end
    
end
