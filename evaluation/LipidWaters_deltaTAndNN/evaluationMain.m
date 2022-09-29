clc; clear all; close all;

%% load results and determine R1
fieldStrengthInT = 3; % main magnetic field strength
addpath(genpath('../../library'));
addpath(genpath('../../txtFiles'));
constants = readConstantsFile('constants.txt');
resultsPath = '../../RESULTS/findDeltaTAndNN_lipidWater_longSimulationTime';
savingPath = sprintf('%s%sPlots%s',resultsPath,filesep,filesep);
addpath(genpath(resultsPath));
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrengthInT;
dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
differentCorrelationFunctionDurations = [25e-9 5e-9 2e-9 1e-9]; % duration after which the correlation function is cutted
offsetReductionRegion = [0.9 1];
saving = 1;
allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsPath,filesep));

for corrFuncDurNr = 1:size(differentCorrelationFunctionDurations,2)
    correlationFunctionDuration = differentCorrelationFunctionDurations( ...
        corrFuncDurNr);
    
    relevantResults = struct();
    
    for fileCounter = 1:size(allMatFilesInResultsDirectory,1)
        fileName = allMatFilesInResultsDirectory(fileCounter,1).name;
        if strcmp(fileName,"")
            error("Matfile name not found.");
        end
        filePath = sprintf('%s%s%s',resultsPath,filesep,fileName);
        load(filePath,'whichLipid','deltaTInS','simulationDurationInS');
        if simulationDurationInS < correlationFunctionDuration
            error(['The simulation duration was shorter than the' ...
                ' given correlation function duration.']);
        end
        lags = round((correlationFunctionDuration-deltaTInS)/deltaTInS);
        datasetResults = loadAndPostprocessResults(filePath ...
            ,lags,offsetReductionRegion);
        datasetResults.r1_NN_theta_phi = ...
            calculateR1ForNearestNeighbourCases(...
            datasetResults.corrFuncFirstOrder ...
            ,datasetResults.corrFuncSecondOrder,omega0,deltaTInS ...
            ,dipolDipolConstant);
        datasetResults.correlationFunctionDuration = ...
            correlationFunctionDuration;
        
        if ~isfield(relevantResults,whichLipid)
            relevantResults.(whichLipid) = datasetResults;
        else
            relevantResults.(whichLipid)(end+1) = datasetResults;
        end
    end
    disp('postprocessing finished');
%     clearvars -except relevantResults savingPath saving
    
    %% Subplot R1 at the top corr func at the bottom for every lipid for
    theta = 1;
    phi = 1;
    nNCase = "nearestNeighbours8000";
    
    for whichLipid = fieldnames(relevantResults)'
        lipidData = relevantResults.(whichLipid{1});
        fig = initializeFigure();
        initializeSubplot(fig,2,1,1);
        legendEntries = {};
        title(sprintf(['Lipid: %s, $\\theta$: %.2f, $\\phi$: %.2f' ...
            ', corrFuncLength: %.2d sec'],whichLipid{1} ...
            ,lipidData(1).fibreAnglesTheta(theta) ...
            ,lipidData(1).fibreAnglesPhi(phi) ...
            ,lipidData(1).correlationFunctionDuration));
        ylabel('R$_1$ [Hz]');
        xlabel('Nearest Neighbours');
        for datasetNr = 1:size(lipidData,2)
            dataset = lipidData(datasetNr);
            plot(dataset.nearestNeighbourCases ...
                ,dataset.r1_NN_theta_phi(:,theta,phi)');
            legendEntries{end+1} = sprintf(['simDur: %.1d' ...
                ', deltaT: %.1d'],dataset.simulationDuration ...
                ,dataset.deltaT); %#ok<SAGROW>
        end
        legend(legendEntries,'location','northwest');
        initializeSubplot(fig,2,1,2);
        legendEntries = {};
        title(sprintf(['Normalized Correlation functions %s'],nNCase));
        plotHandles = [];
        xlabel('Correlation time [s]');
        for datasetNr = 1:size(lipidData,2)
            dataset = lipidData(datasetNr);
            corrFuncFirstOrder = squeeze(dataset.corrFuncFirstOrder.( ...
                nNCase)(theta,phi,:))';
            timeAxis = [0:length(corrFuncFirstOrder)-1]*dataset.deltaT;
            plot(timeAxis,corrFuncFirstOrder,':');
            legendEntries{end+1} = sprintf(['first order, simDur: %.2d, deltaT:' ...
                ' %.2d, corrFuncLength: %.2d'],dataset.simulationDuration ...
                ,dataset.deltaT,dataset.simulationDuration); %#ok<SAGROW>
            corrFuncSecondOrder = squeeze(dataset.corrFuncSecondOrder.( ...
                nNCase)(theta,phi,:))';
            plot(timeAxis,corrFuncSecondOrder,'--');
            legendEntries{end+1} = sprintf(['second order, simDur: %.2d, deltaT:' ...
                ' %.2d, corrFuncLength: %.2d'],dataset.simulationDuration ...
                ,dataset.deltaT,dataset.simulationDuration); %#ok<SAGROW>
        end
        legend(legendEntries);
        if saving
            savingName = sprintf('%s%s',savingPath,filesep);
            fileSavingName = sprintf( ...
                '%s%s_R1sAndCorrFunc_%.0d' ...
                ,savingName,whichLipid{1} ...
                ,lipidData(1).correlationFunctionDuration);
            if ~exist(savingName,'dir')
                mkdir(savingName)
            end
            print(gcf,fileSavingName,'-dpng','-r300');
        end
    end
    
    
end




% %% Present/Plot results
% 
% % 1. plot r1 rates
% for whichLipid = fieldnames(relevantResults)'
%     lipidData = relevantResults.(whichLipid{1});
%     initializeFigure();
%     legendEntries = {};
%     title(sprintf('Lipid: %s', whichLipid{1}));
%     for datasetNr = 1:size(lipidData,2)
%         dataset = lipidData(datasetNr);
%         r1_NN_theta_phi = dataset.r1_NN_theta_phi;
%         for thetaNr = 1:size(dataset.fibreAnglesTheta,2)
%             for phiNr = 1:size(dataset.fibreAnglesPhi,2)
%                 plot(dataset.nearestNeighbourCases ...
%                     ,r1_NN_theta_phi(:,thetaNr,phiNr)');
%                 legendEntries{end+1} = sprintf(['simDuration: %.2d, dt:' ...
%                     ' %.2d, $\\theta$: %.2f, $\\phi$: %.2f'] ...
%                     , dataset.simulationDuration, dataset.deltaT ...
%                     , dataset.fibreAnglesTheta(thetaNr) ...
%                     , dataset.fibreAnglesPhi(phiNr)); %#ok<SAGROW>
%             end
%         end
%     end
%     legend(legendEntries);
%     if saving
%         savingName = sprintf('%srelaxationRates%s',savingPath,filesep);
%         fileSavingName = sprintf('%s%s_R1ForDifferentSimConfigs',savingName,whichLipid{1});
%         if ~exist(savingName,'dir')
%             mkdir(savingName)
%         end
%         print(gcf,fileSavingName,'-dpng','-r300');
%     end
% end
% 
% %% 2. plot the correlation functions
% % plot correlation functions for everey NN case but different sim configs.
% whichCorrFunc = 'corrFuncSecondOrder';
% thetaNr = 1;
% phiNr = 1;
% 
% for whichLipid = fieldnames(relevantResults)'
%     lipidData = relevantResults.(whichLipid{1});
%     for whichNNCase = fieldnames(lipidData(1).(whichCorrFunc))'
%         if  getNNFromNNCase(whichNNCase) < 8000
%             continue
%         end
%         initializeFigure();
%         title(sprintf('Lipid: %s, %s',whichLipid{1}, whichNNCase{1}));
%         legendEntries = {};
%         for datasetNr = 1:size(lipidData,2)
%             dataset = lipidData(datasetNr);
%             corrFunc = squeeze(dataset.(whichCorrFunc).( ...
%                 whichNNCase{1})(thetaNr,phiNr,:))';
%             timeAxis = 0:dataset.deltaT:(size(corrFunc,2)-1)*dataset.deltaT;
%             plot(timeAxis,corrFunc/corrFunc(1));
%             legendEntries{end+1} = sprintf('simDuration: %.2d, deltaT: %.2d' ...
%                 ,dataset.simulationDuration,dataset.deltaT); %#ok<SAGROW>
%         end
%         legend(legendEntries);
%         if saving
%             savingName = sprintf('%scorrelationFunctions%s',savingPath,filesep);
%             fileSavingName = sprintf( ...
%                 '%s%s_CorrelationFunctionsForDifferentSimConfigs' ...
%                 ,savingName,whichLipid{1});
%             if ~exist(savingName,'dir')
%                 mkdir(savingName)
%             end
%             print(gcf,fileSavingName,'-dpng','-r300');
%         end
%     end
% end

