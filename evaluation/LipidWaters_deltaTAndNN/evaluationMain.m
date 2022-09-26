clc; clear all; close all;

%% load results and determine R1
addpath(genpath('../../library'));
addpath(genpath('../../txtFiles'));
constants = readConstantsFile('constants.txt');
resultsPath = '../../RESULTS/findDeltaTAndNN_lipidWater_longSimulationTime';
savingPath = sprintf('%s%sPlots%s',resultsPath,filesep,filesep);
addpath(genpath(resultsPath));
fieldStrengthInT = 3; % main magnetic field strength
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrengthInT;
dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);

allMatFilesInResultsDirectory = dir( ...
    sprintf('%s%s*.mat',resultsPath,filesep));

relevantResults = struct();

for fileCounter = 1:size(allMatFilesInResultsDirectory,1)
    fileName = allMatFilesInResultsDirectory(fileCounter,1).name;
    
    if strcmp(fileName,"")
        error("Matfile name not found. ");
    end
    
    filePath = sprintf('%s%s%s',resultsPath,filesep,fileName);
    load(filePath,'whichLipid');
    if ~isfield(relevantResults,whichLipid)
        relevantResults.(whichLipid) = {};
    end
    relevantResults.(whichLipid){end+1} = loadResultsAndDetermineR1( ...
        filePath,omega0,dipolDipolConstant);
    relevantResults.(whichLipid){end}.('fieldStrengthInT') = fieldStrengthInT;
    
end

clearvars -except relevantResults savingPath
%% Present/Plot results

% 1. plot r1 rates
colors = getColorsFromColorConfigurationFile();
for whichLipid = fieldnames(relevantResults)'
    lipidData = relevantResults.(whichLipid{1});
    initializeFigure();
    legendEntries = {};
    for datasetNr = 1:size(lipidData,2)
        dataset = lipidData{datasetNr};
        r1_theta_phi = dataset.r1_theta_phi;
        title(sprintf('Lipid: %s', whichLipid{1}));
        for thetaNr = 1:size(dataset.fibreAnglesTheta,2)
            for phiNr = 1:size(dataset.fibreAnglesPhi,2)
                r1_NN = [];
                for nearestNeighbourCase = fieldnames(dataset.r1_theta_phi)'
                    r1_NN(end+1) = r1_theta_phi.(nearestNeighbourCase{1})(thetaNr,phiNr); %#ok<SAGROW>
                end
                plot(dataset.nearestNeighbourCases,r1_NN);
                legendEntries{end+1} = sprintf('simDuration: %.2d, dt: %.2d, $\\theta$: %.2f, $\\phi$: %.2f', dataset.simulationDuration, dataset.deltaT, dataset.fibreAnglesTheta(thetaNr), dataset.fibreAnglesPhi(phiNr)); %#ok<SAGROW>
            end
        end
    end
    legend(legendEntries);
    savingName = sprintf('%srelaxationRates%s',savingPath,filesep); 
    fileSavingName = sprintf('%s%s_R1ForDifferentSimConfigs',savingName,whichLipid{1});
    if ~exist(savingName,'dir')
        mkdir(savingName)
    end
    print(gcf,fileSavingName,'-dpng','-r300');
end


% 2. plot the correlation functions 



