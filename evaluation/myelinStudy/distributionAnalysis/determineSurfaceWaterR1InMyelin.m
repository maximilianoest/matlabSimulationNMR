clc; clear all; close all; fclose('all');


addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions%s",createFilesepStringArray(5));
if ~exist(resultsDir,'dir')
    warning('The results directory does not exist but will be created.');
    mkdir(resultsDir)
end

%% set up
windowSize = 21;
surfWaterDensBorder = 0.93;


%% density stuff
densDistribPathFolder = "C:\Users\maxoe\Google Drive\Promotion" ...
    + "\Simulation\RESULTS\wholeMyelin_densityDistributions\";
densDistribFileName = "20230404_myelinDistributionData";
densData = load(densDistribPathFolder + densDistribFileName + ".mat");

avgLocations = mean(densData.discreteLocations(:,2:end-1),1);
avgDeltaVol = mean(densData.dVolumeForTimeStep);

waterIndex = densData.moleculeNameCollection == "WATER";
hydrogenIndex = densData.atomNameCollection == "H";

densDistribWater = squeeze(sum(mean( ...
    densData.density_timeLocMolAtom(:,:,waterIndex,:),1),3));
densDistribMembr = squeeze(sum(mean( ...
    densData.density_timeLocMolAtom(:,:,~waterIndex,:),1),3));

modelWaterMass = sum(sum(densDistribWater)) * avgDeltaVol;
modelMembrMass = sum(sum(densDistribMembr)) * avgDeltaVol;
waterMassRatio = modelWaterMass/modelMembrMass;

fprintf("<strong> In this model: </strong> \n");
fprintf("  Water mass: %.4d kg.\n  Membrane mass: %.4d kg.\n" ...
    ,modelWaterMass,modelMembrMass);
fprintf("  => water mass ratio: %.f\n",waterMassRatio);


densDistribWaterH = meanFilter(squeeze(sum(mean( ...
    densData.density_timeLocMolAtom(:,2:end-1,waterIndex,hydrogenIndex) ...
    ,1),3)),windowSize);
densDistribMembrH = meanFilter(squeeze(sum(mean( ...
    densData.density_timeLocMolAtom(:,2:end-1,~waterIndex,hydrogenIndex) ...
    ,1),3)),windowSize);

surfWaterIndices = densDistribWaterH < ...
    surfWaterDensBorder * mean(densDistribWaterH(end/2-100:end/2+100));
surfWaterBorders = find(diff(surfWaterIndices) ~= 0);
if length(surfWaterBorders) ~= 2
    error("The borders for surface water are not chosen correctly.");
end

densSurfWater = nan(1,length(densDistribWaterH));
densSurfWater(surfWaterIndices) = densDistribWaterH(surfWaterIndices);

initializeFigure();
plot(avgLocations,densDistribWaterH);
plot(avgLocations,densSurfWater);
plot(avgLocations,densDistribMembrH);
plotVerticalLine(avgLocations(surfWaterBorders(1)),max(densDistribWaterH));
plotVerticalLine(avgLocations(surfWaterBorders(2)),max(densDistribWaterH));
legend("Free H$_{2}$O","Surf. H$_{2}$O","Membrane",'Location','north');
title("Hydrogen atoms density");
xlabel('x position $[m]$');
ylabel('density $[kg/m^3]$');


%% relaxation stuff

r1DataFolder = "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\wholeMyelin_relaxationRates\";
r1DataFileName = "20230405_compartmentAndCrossR1FromSimulations";
r1Data = load(r1DataFolder + r1DataFileName + ".mat");










%% functions in use
function plt = plotVerticalLine(xPos,yValue)
plt = plot([xPos xPos],[0 yValue],'--k','LineWidth',0.7);
end









