clc; clear all; close all;
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
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
saving = 0;
theta = 2;
phi = 1;
nNCase = 'nearestNeighbours5500';
allMatFilesInResultsDirectory = dir(sprintf( ...
    '%s%s*.mat',resultsDir,filesep));
% [left bottom width height]
figPosAndSize = [50 50 1000 700];
widthAndHeight = [0.43 0.40];
subFigPosAndSize = [ ...
    0.07 0.55 widthAndHeight ; ...
    0.55 0.55 widthAndHeight ; ...
    0.07 0.06 widthAndHeight ; ...
    0.55 0.06 widthAndHeight];
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
    if strcmp(constituent,'lipidWater') ...
            || simulationResults.simulationDurationInS < 900e-9
        continue;
    end
    corrFuncZerothOrder = squeeze(...
        simulationResults.sumCorrFuncZerothOrderShortPadding.(nNCase) ...
        (theta,phi,:)./simulationResults.atomCounter)';
    corrFuncFirstOrder = squeeze(...
        simulationResults.sumCorrFuncFirstOrderShortPadding.(nNCase) ...
        (theta,phi,:)./simulationResults.atomCounter)';
    corrFuncSecondOrder = squeeze(...
        simulationResults.sumCorrFuncSecondOrderShortPadding.(nNCase) ...
        (theta,phi,:)./simulationResults.atomCounter)';
    
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid)(1).constituent = constituent;
    else
        relevantResults.(whichLipid)(end+1).constituent = constituent;
    end
    arrayPos = length(relevantResults.(whichLipid));
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
        rad2deg(simulationResults.fibreAnglesTheta);
    relevantResults.(whichLipid)(arrayPos).fibreAnglesPhi = ...
        rad2deg(simulationResults.fibreAnglesPhi);
    relevantResults.(whichLipid)(arrayPos).matlabSimulationDate = ...
        simulationResults.matlabSimulationDate;
end
clearvars('-except',configVariables{:},'relevantResults');


%% analyse data
corrFuncNames = ["corrFuncFirstOrder" "corrFuncSecondOrder"];

for whichLipid = string(fieldnames(relevantResults)')
    lipidData = relevantResults.(whichLipid);
    fig = initializeFigure('posAndSize',figPosAndSize);
    
    for subFigNr = 1:4
        subFigs(subFigNr) = initializeSubplot(fig,2,2,subFigNr); %#ok<SAGROW>
        subFigs(subFigNr).Position = subFigPosAndSize(subFigNr,:); %#ok<SAGROW>
    end
    lgd = legend();
    legendEntries = {};
    for datasetNr = 1:length(lipidData)
        dataset = lipidData(datasetNr);
        
        deltaT = dataset.deltaTInS;
        simulationDuration = dataset.simulationDurationInS;
        
        legendEntries{end+1} = sprintf("sim. Dur.: %.2d, dT: %.2d" ...
            ,deltaT,simulationDuration); %#ok<SAGROW>
        
        for corrFuncCaseNr = 1:length(corrFuncNames)
            if corrFuncNames(corrFuncCaseNr) == "corrFuncFirstOrder"
                corrFuncPlotIndex = 1;
                specDensPlotIndex = 3;
                preFactorSpecDens = -1;
            elseif corrFuncNames(corrFuncCaseNr) == "corrFuncSecondOrder"
                corrFuncPlotIndex = 2;
                specDensPlotIndex = 4;
                preFactorSpecDens = -2;
            else
                error("not implemented");
            end
            
            corrFunc = double(dataset.(corrFuncNames(corrFuncCaseNr( ...
                1:round(end*0.7))));
            corrFuncLength = length(corrFunc);
            timeAxis = 0:deltaT:(corrFuncLength-1)*deltaT;
        
            % correlation function
            axes(subFigs(corrFuncPlotIndex)); %#ok<LAXES>
            lgd = legend();
            plot(timeAxis,abs(corrFunc));
            axis([-inf inf -1.5e4 1.5e4]);
            
            if corrFuncPlotIndex == 1
                ylabel("abs(corr. Func. FO)");
                title(sprintf("Lipid: %s, $\\theta$: %.2f, $\\phi$: %.2f" ...
                    ,whichLipid,dataset.fibreAnglesTheta(theta) ...
                    ,dataset.fibreAnglesPhi(phi)));
                lgd.String = legendEntries;
                lgd.Visible = 'On';
            elseif corrFuncPlotIndex == 2
                ylabel("abs(corr. Func. SO)");
                lgd.Visible = 'Off';
            else
                error("not implemented");
            end
            xlabel("lag time [s]");
            
            % spectral density
            specDens = 2*(deltaT*cumsum(corrFunc.*exp( ...
                preFactorSpecDens*1i*omega0*deltaT*(0:corrFuncLength-1))));
            axes(subFigs(specDensPlotIndex)); %#ok<LAXES>
            lgd = legend();
            plot(timeAxis,abs(real(specDens)));
            
            if specDensPlotIndex == 3
                ylabel("abs(real(Spec. Dens. FO))");
                lgd.Visible = 'Off';
            elseif specDensPlotIndex == 4
                ylabel("abs(reall(Spec. Dens. SO))");
                lgd.Visible = 'Off';
            else
                error("not implemented");
            end
            xlabel("upper limit of FT [s]")
            
        end
        
        
%         ylabel('Correlation function');
%         if datasetNr == 1
%             lgd = legend( ...
%                 'FO imag/real' ...
%                 ,'FO real' ...
%                 ,'SO imag/real' ...
%                 ,'SO real');
%             lgd.Location = 'northeast';
%             title(sprintf('Lipid: %s, influence of imaginary part' ...
%                 ,whichLipid));
%         else
%             lgd = legend();
%             lgd.Visible = 'Off';
%             title(sprintf('Lipid water: %s, influence of imaginary part' ...
%                 ,whichLipid));
%         end
%         
%         subplotCounter = subplotCounter + 1;
%         initializeSubplot(fig,2,2,subplotCounter);
%         set(gca,'Position',subFigPosAndSize(subplotCounter,:));
%         
%         r1WholeCorrFunc = 0;
%         r1OnlyReal = 0;
%         for corrFuncCaseNr = 1:length(corrFuncNames)
%             corrFunc = double(dataset.(corrFuncNames(corrFuncCaseNr)));
%             if strcmp(dataset.constituent,'lipidWater')
%                 corrFunc = corrFunc(1:40000);
%             end
%             
%             corrFuncLength = length(corrFunc);
%             switch corrFuncNames(corrFuncCaseNr) 
%                 case 'corrFuncFirstOrder'
%                     preFactor = 1;
%                 case 'corrFuncSecondOrder'
%                     preFactor = 2;
%                 otherwise 
%                     error('Not implemented.');
%             end
%             
%             spectralDensityWholeCorrFunc = 2*(deltaT*cumsum(corrFunc' ... 
%                 .*exp(-preFactor*1i*omega0*deltaT*(0:corrFuncLength-1))));
%             spectralDensityOnlyReal = 2*(deltaT*cumsum(real(corrFunc') ... 
%                 .*exp(-preFactor*1i*omega0*deltaT*(0:corrFuncLength-1))));
%             r1WholeCorrFunc = r1WholeCorrFunc ...
%                 + dipoleDipoleConstant*3/2*abs(real(mean( ...
%                 spectralDensityWholeCorrFunc(end-5000:end))));
%             r1OnlyReal = r1OnlyReal  ...
%                 + dipoleDipoleConstant*3/2*abs(real(mean( ...
%                 spectralDensityOnlyReal(end-5000:end))));
%             
%             plt = plot(timeAxis,real(spectralDensityWholeCorrFunc));
%             plot(timeAxis,real(spectralDensityOnlyReal),':');
%             
%         end
%         
%         ylabel('Real(Spectral density)');
%         xlabel('lag time [s]');
%         title(sprintf(['Corr Func: R$_1$ = %.2f' ...
%             ', Real(Corr. Func.): R$_1$ = %.2f'],r1WholeCorrFunc ...
%             ,r1OnlyReal));
%         if datasetNr == 1
%             lgd = legend( ...
%                 'FO corrFunc' ...
%                 ,'FO real(corrFunc)' ...
%                 ,'SO corrFunc' ...
%                 ,'SO real(corrFunc)');
%             lgd.Location = 'northeast';
%         else
%             lgd = legend();
%             lgd.Visible = 'Off';
%         end
%         subplotCounter = subplotCounter + 1;
%     end
%     
%     if saving
%         saveFigureTo(savingDirectory,whichLipid ...
%             ,dataset.matlabSimulationDate ...
%             ,'realVsWholeCorrFunc');
%     end
    
    end
end
