clc; clear all; close all;

% ATTENTION: with mean positions, it is possible that more than one atom is
% at one position. Thus the distribution and density is not right.
% Therefore, a snapshot of the production run is used.

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf( ...
    '..%s..%sRESULTS%sangleDependency_lipidWater_bigSimulation' ...
    ,createFilesepStringArray(3));
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
addpath(genpath(resultsDir));

savingDirectory = sprintf('%s%sPlots%s',resultsDir ...
    ,createFilesepStringArray(2));
saving = 1;

%% get relevant results from files
directory = "C:\Users\maxoe\Documents\Gromacs\testDatasets\";
lipidNames = ["DOPS" "PLPC" "PSM"];
filesToLoad = [ ...
    "20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns" ...
    "20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns" ...
    "20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt1ps_simTime20ns"];
    
relevantInformation = struct();
fig = initializeFigure();
legendEntries = {};
avogadro = 6.02214076e23;
molMassWaterKG = 18.015268/1000;
massWaterKG = molMassWaterKG/avogadro;

for datasetNr = 1:length(lipidNames)
    relevantInformation(datasetNr).lipidName = lipidNames(datasetNr);
    getSimulationConfigurationFromFileName(filesToLoad(datasetNr));
    disp("Loading data.");
    simData = load(sprintf("%s%s%s%s.mat",directory,lipidNames(datasetNr) ...
        ,filesep,filesToLoad(datasetNr)));
    
    positionsInNM = simData.trajectories;
    positionsAroundZero = positionsInNM - mean(positionsInNM,1);
    positionsInM = positionsAroundZero * constants.nanoMeter;
    
    disp("Plotting histogram.");
    initializeSubplot(fig,2,1,1);
    title("Distributions of water molecules in model");
    histogram(positionsInM(:,1,10000),'BinWidth',1e-10,'FaceAlpha',0.3);
    xlabel('Mean Position $[nm]$');
    ylabel('Frequency $[\#]$');
    legendEntries{end+1} = sprintf("%s",whichLipid); %#ok<SAGROW>
    legend(legendEntries);
    
    disp("Plotting density.");
    initializeSubplot(fig,2,1,2);
    averagingTimePoints = 1:100:size(positionsInM,3);
    densityVector = zeros(1,100);
    densityDistribution = zeros(length(averagingTimePoints) ...
        ,length(densityVector));
    densityLocations = zeros(length(averagingTimePoints) ...
        ,length(densityVector));
    for timePointNr = 1:length(averagingTimePoints)
        densityVector = zeros(1,length(densityVector));
        timePoint = averagingTimePoints(timePointNr);
        positionsXInM = positionsInM(:,1,timePoint);
        positionsYInM = positionsInM(:,2,timePoint);
        positionsZInM = positionsInM(:,3,timePoint);
        
        densityLocations(timePointNr,:) = linspace(min(positionsXInM) ...
            ,max(positionsXInM),length(densityVector));
        for atomNr = 1:length(positionsXInM)
            xPosition = positionsXInM(atomNr);
            for locationNr = 1:size(densityLocations,2)-1
                if locationNr == 1
                    lowerBoundDecrease = - constants.nanoMeter;
                    upperBoundIncrease = 0;
                elseif locationNr == length(densityLocations( ...
                        timePointNr,:))-1
                    lowerBoundDecrease = 0;
                    upperBoundIncrease = constants.nanoMeter;
                else
                    lowerBoundDecrease = 0;
                    upperBoundIncrease = 0;
                end
                lowerBound = densityLocations(timePointNr,locationNr)  ...
                    + lowerBoundDecrease;
                upperBound = densityLocations(timePointNr,locationNr ...
                    + 1) + upperBoundIncrease;
                if xPosition >= lowerBound && xPosition < upperBound
                    densityVector(locationNr) ...
                        = densityVector(locationNr) + 1;
                    continue;
                end
            end
        end
        dX = mean(diff(densityLocations(timePointNr,:)));
        dY = max(positionsYInM) - min(positionsYInM);
        dZ = max(positionsZInM) - min(positionsZInM);
        
        waterMoleculeCount = densityVector/2;
        waterMoleculesMass = waterMoleculeCount * massWaterKG;
        dVolume = dX*dY*dZ;
        
        densityDistribution(timePointNr,:) = waterMoleculesMass / dVolume;
    end
    avgDensityDistribution = mean(densityDistribution,1);
    avgDensityLocations = mean(densityLocations,1);
    plot(avgDensityLocations,avgDensityDistribution);
    xlabel('Position $[nm]$');
    ylabel('Density $[kg/m^3]$');
    lgd = legend();
    lgd.Visible = "Off";
    
end

%%
if saving
    saveFigureTo(savingDirectory,"LocationDistributions",""...
        ,"lipidWater");
end

