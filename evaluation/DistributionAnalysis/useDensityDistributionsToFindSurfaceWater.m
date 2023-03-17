clc; clear; close all; fclose('all');
% Description:
% In the orientation-dependency paper the authors said that in the
% cross-region is where both densities of solid lipid and lipid water H are
% larger than 0. They found that this region are 32% of the membran
% associated hydrogen nuceli. Here I am interested in the part of the lipid
% water pool that is interacting with the solid lipid part. Therefore, I
% divide the density of the interacting part by the whole water pool.

%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

%% paths to load
densityDistributionsFolder = "C:\Users\maxoe\Google Drive\Promotion" ...
    + "\Simulation\RESULTS\densityDistributions\";
distributionDateID = "20230309";
lipidNames = ["DOPS" "PLPC" "PSM"];
fileNames = distributionDateID + "densityDistribution_" ...
    + lipidNames + ".mat";

matlabSimResultsFolder = "C:\Users\maxoe\Google Drive\Promotion\" ...
    + "Simulation\RESULTS\angleDependency_lipidWater_bigSimulation\";
matFileNames = [ ...
    "20221101_Results_DOPSwater_20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221025_Results_PLPCwater_20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221028_Results_PSMwater_20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ] ...
    + ".mat";
corrFuncFilePaths = matlabSimResultsFolder + matFileNames;

%% set up
% fieldstrength
fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 18; % °C
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);

windowSize = 9;
averageDensityRegion = 50;
bulkWaterFractionForBorder = 0.95;
fprintf("Bulk water density fraction for border: %.4f\n" ...
    ,bulkWaterFractionForBorder);
fourthSubplotText  = sprintf("Bulk water density fraction: %.4f \n" ...
    ,bulkWaterFractionForBorder);

saving = 1;

%% load data

fig = initializeFigure();
for fileNr = 1:length(fileNames)
    lipidName = lipidNames(fileNr);
    fprintf(" ---- %s \n",lipidName);
    distribData = load(sprintf("%s%s",densityDistributionsFolder ...
        ,fileNames(fileNr)));
    if lipidName ~= distribData.lipidName
        error("Not the same lipid!");
    end
    
    hDensityLipid = squeeze(mean( ...
        distribData.densityDistributionOnlyHdydrogen(: ...
        ,distribData.moleculesToCompare == lipidName,2:end-1),1))';
    hDensityWater = squeeze(mean( ...
        distribData.densityDistributionOnlyHdydrogen(: ...
        ,distribData.moleculesToCompare == "SOL",2:end-1),1))';
    locations = distribData.avgLocations(2:end-1);
    
    meanFilteredHDensityLipid = meanFilter(hDensityLipid,windowSize);
    meanFilteredHDensityWater = meanFilter(hDensityWater,windowSize);
    
    hydrogenAtomWeight = distribData.atomWeightCollection( ...
        distribData.atomNameCollection == "H");
    hFrequencyLipid = meanFilteredHDensityLipid*mean( ...
        distribData.dVolumeForTimeStep)/hydrogenAtomWeight;
    hFrequencyWater = meanFilteredHDensityWater*mean( ...
        distribData.dVolumeForTimeStep)/hydrogenAtomWeight;
    
    allHsInLipidPool = sum(hFrequencyLipid);
    allHsHsInWaterPool = sum(hFrequencyWater);
    
    fprintf("All Hs in lipid: %.2f \nAll Hs in water: %.2f \n" ...
        ,allHsInLipidPool,allHsHsInWaterPool);
    
    bulkWaterHDensity = mean(hDensityWater(end/2-averageDensityRegion ...
        : end/2+averageDensityRegion));
    fprintf("Bulk water hydrogen atom density: %.4f. kg/m^3. \n" ...
        + "Density border: %.4f \n" ...
        ,bulkWaterHDensity,bulkWaterHDensity*bulkWaterFractionForBorder);
    
    surfWater = meanFilteredHDensityWater  ...
        < bulkWaterHDensity*bulkWaterFractionForBorder;
    borders = find(abs(diff(surfWater)) == 1);
    
    surfaceWaterFraction = sum(meanFilteredHDensityWater(surfWater)) ...
        /sum(meanFilteredHDensityWater);
    freeWaterFraction = sum(meanFilteredHDensityWater(~surfWater)) ...
        /sum(meanFilteredHDensityWater);
    
    fprintf("Interacting water fraction: %.4f \n" ...
        +"Free water fraction: %.4f\n",surfaceWaterFraction ...
        ,freeWaterFraction);
    
    
    initializeSubplot(fig,2,2,fileNr);
    title(lipidName);
    
    plot(locations(1:end-1),meanFilteredHDensityLipid(1:end-1));
    plot(locations(1:end-1),meanFilteredHDensityWater(1:end-1));
    
    plotVerticalLine(locations(borders(1)),120);
    plotVerticalLine(locations(borders(end)),120);
    
    ylabel("Density [$\frac{kg}{m^3}$]");
    xlabel("Location in [m]");
    xticks(-5e-9 : 1e-9 : 5e-9);
    lgd = legend("Lipid","Water");
    lgd.Location = 'best';
    
    if fileNr ~= 1
        lgd.Visible = 'off';
    end
    
    
    fourthSubplotText = sprintf("%sLipid: %s\n  surf water: " ...
        +" %.4f \n  free water: %.4f \n" ...
        ,fourthSubplotText,lipidName,surfaceWaterFraction ...
        ,freeWaterFraction);
    
    
    corrFuncData = load(corrFuncFilePaths(fileNr));
    if corrFuncData.whichLipid ~= lipidName
        error("Not the same lipids are compared with each other.");
    end
    deltaTInS = corrFuncData.deltaTInS;
    corrFuncFirstOrder = corrFuncData.sumCorrFuncFirstOrder ...
        .nearestNeighbours8000 / corrFuncData.atomCounter;
    corrFuncSecondOrder = corrFuncData.sumCorrFuncSecondOrder ...
        .nearestNeighbours8000 / corrFuncData.atomCounter;
    
    [specDensFirstOrder, specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder ...
        ,corrFuncSecondOrder,omega0,deltaTInS,[0.9 1]);
    
    dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
        *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
        /(constants.nanoMeter^6);
    totalWaterPoolR1 = calculateR1WithSpectralDensity( ...
        specDensFirstOrder,specDensSecondOrder,dipolDipolConstant);
    fprintf("  Total water pool R1 = %.4f \n",totalWaterPoolR1);
    
   
    fprintf("  Surface/Free water fraction: %.4f / %.4f\n" ...
        ,surfaceWaterFraction,freeWaterFraction);
    
    surfaceWaterR1 = totalWaterPoolR1 / surfaceWaterFraction ...
        - freeWaterFraction / surfaceWaterFraction ...
        * freeWaterR1;
    
    fprintf("  Surface water R1: %.4f\n",surfaceWaterR1);
    
    fourthSubplotText = sprintf("%sR1surf: %.4f\n \n",fourthSubplotText ...
        ,surfaceWaterR1);
    
end
%%
ax = initializeSubplot(fig,2,2,4);
text(0.05,0.5,fourthSubplotText,'interpreter','latex');
set(ax,'visible','off');
lgd = legend();
lgd.Visible = 'Off';

if saving
    saveFigureTo(densityDistributionsFolder ...
        ,datestr(now,"yyyymmdd") + "_allLipids" ...
        ,"DensityDistributionsAndSurfWater","onlyHydrogen_bulkBorder" ...
        + num2str(round(bulkWaterFractionForBorder*100)),true);
end

%% functions in use
function plt = plotVerticalLine(xPos,yValue)
plt = plot([xPos xPos],[0 yValue],'--k','LineWidth',0.7);
end


