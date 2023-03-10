clc; clear all; close all;

%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

%% data paths
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

%% set up parameter
hCountBorder = 0.05;
fprintf("H atoms count border: %2.f \n",hCountBorder);

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 18; % °C
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);

% fieldstrength
fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;

% equation to determine lipid water R1 based on pool fraction factor:
% R1_tot = f_surf * R1_surf + f_free * R1_free
% 
% f_free(f) = #_H,free / #_H,LW = 1 - #_H,surf / #_H,SL * f
% f_surf(f) = #_H,surf / #_H,LW = #_H,surf / #_H,SL * f = 1 - f_free(f)

largestPossiblePoolFractionFactor = 8;
saving = 1;

%% evaluation
for lipidNr = 1:length(fileNames)
    lipid = lipidNames(lipidNr);
    fprintf("Lipid %s\n",lipid);
    distribData = load(densityDistributionsFolder + fileNames(lipidNr));
    avgDensDistribOfH = squeeze(mean( ...
        distribData.densityDistributionOnlyHdydrogen,1));
    hFrequencyDistributionInWaterPool = avgDensDistribOfH( ...
        distribData.moleculesToCompare == "SOL",:) ...
        /distribData.atomWeightCollection( ...
        distribData.atomNameCollection == "H") ...
        * mean(distribData.dVolumeForTimeStep);
    hFrequencyDistributionInLipidPool = avgDensDistribOfH( ...
        distribData.moleculesToCompare == lipid,:) ...
        /distribData.atomWeightCollection( ...
        distribData.atomNameCollection == "H") ...
        * mean(distribData.dVolumeForTimeStep);
    
    overallHsInWaterPool = sum(hFrequencyDistributionInWaterPool);
    overallHsInLipidPool = sum(hFrequencyDistributionInLipidPool);
    poolFractionFactorInSim = overallHsInLipidPool/overallHsInWaterPool;
    fprintf("  Pool fraction factor in MD sim: %.4f\n" ...
        ,poolFractionFactorInSim);
    
    surfaceWaterHDistrib = hFrequencyDistributionInWaterPool ...
        .* (hFrequencyDistributionInLipidPool > hCountBorder ...
        & hFrequencyDistributionInWaterPool > hCountBorder);
    whereIsSurfWater = surfaceWaterHDistrib > 0;
    borders = find(diff(whereIsSurfWater) ~= 0);
    overallHsInSurfacePool = sum(surfaceWaterHDistrib);
    
    freeWaterHDistrib = zeros(1,length(hFrequencyDistributionInWaterPool));
    freeWaterHDistrib(borders(2)+1:borders(3)) = ...
        hFrequencyDistributionInWaterPool(borders(2)+1:borders(3));
    whereIsFreeWater = freeWaterHDistrib > 0;
    
    avgLocations = distribData.avgLocations;
    
    % 1. determine R1^free (Tsukiashi2018) - check
    
    % 2. determine R1^surf (with equation) - check
    
    corrFuncsData = load(corrFuncFilePaths(lipidNr));
    if corrFuncsData.whichLipid ~= lipid
        error("wrong dataset!");
    end
    deltaTInS = corrFuncsData.deltaTInS;
    corrFuncFirstOrder = corrFuncsData.sumCorrFuncFirstOrder ...
        .nearestNeighbours8000 / corrFuncsData.atomCounter;
    corrFuncSecondOrder = corrFuncsData.sumCorrFuncSecondOrder ...
        .nearestNeighbours8000 / corrFuncsData.atomCounter;
    
    [specDensFirstOrder, specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder ...
        ,corrFuncSecondOrder,omega0,deltaTInS,[0.9 1]);
    
    dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
        *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
        /(constants.nanoMeter^6);
    totalWaterPoolR1 = calculateR1WithSpectralDensity( ...
        specDensFirstOrder,specDensSecondOrder,dipolDipolConstant);
    fprintf("  Total water pool R1 = %.4f \n",totalWaterPoolR1);
    
    
    surfaceInteractionFraction = sum(surfaceWaterHDistrib) ...
        / overallHsInWaterPool;
    freeWaterInteractionFraction = sum(freeWaterHDistrib) ...
        / overallHsInWaterPool;
    fprintf("  Surface/Free water fraction: %.4f / %.4f\n" ...
        ,surfaceInteractionFraction,freeWaterInteractionFraction);
    
    surfaceWaterR1 = totalWaterPoolR1 / surfaceInteractionFraction ...
        - freeWaterInteractionFraction / surfaceInteractionFraction ...
        * freeWaterR1;
    
    fprintf("  Surface water R1: %.4f\n",surfaceWaterR1);
    
    % 3. determine f_max and R1_LW for different f till f_max
    
    % right wing:
    rightWingLipidHDistrib = hFrequencyDistributionInLipidPool( ...
        length(hFrequencyDistributionInLipidPool)/2+1:end);
    rightWingWaterHDistrib = hFrequencyDistributionInWaterPool( ...
        length(hFrequencyDistributionInWaterPool)/2+1:end);
    rightWingLocations = avgLocations(length(avgLocations)/2+1:end);
    
    fig1 = initializeFigure();
    initializeSubplot(fig1,2,1,1);
    title(sprintf("Wings %s",lipid));
    ylabel("Frequency H atoms");
    xlabel("Location from center [m]");
    [rightR1_LwBelowF_max,rightFBelowF_max,rightR1SurfForLargeF ...
        ,rightLargerF] = determineR1SurfForWing( ...
        "right",rightWingWaterHDistrib,rightWingLipidHDistrib ...
        ,rightWingLocations,hCountBorder,largestPossiblePoolFractionFactor ...
        ,surfaceWaterR1,freeWaterR1);
    
    leftWingLipidHDistrib = fliplr(hFrequencyDistributionInLipidPool( ...
        1:length(hFrequencyDistributionInLipidPool)/2));
    leftWingWaterHDistrib = fliplr(hFrequencyDistributionInWaterPool( ...
        1:length(hFrequencyDistributionInWaterPool)/2));
    leftWingLocations = abs(fliplr(avgLocations(1:length(avgLocations)/2)));
    
    [leftR1_LwBelowF_max,leftFBelowF_max,leftR1SurfForLargeF ...
        ,leftLargerF] = determineR1SurfForWing( ...
        "left",leftWingWaterHDistrib,leftWingLipidHDistrib ...
        ,leftWingLocations,hCountBorder,largestPossiblePoolFractionFactor ...
        ,surfaceWaterR1,freeWaterR1);
    legend("right surf","right free","right lipid","left surf" ...
        ,"left free","left lipid","Location","east");
    
    initializeSubplot(fig1,2,1,2);
    plot(rightFBelowF_max,rightR1_LwBelowF_max);
    plot(rightLargerF,rightR1SurfForLargeF,':');
    plot(leftFBelowF_max,leftR1_LwBelowF_max);
    plot(leftLargerF,leftR1SurfForLargeF,'--');
    
    % averaged over both wings
    fForAverage = round([0:0.01:largestPossiblePoolFractionFactor]*100);
    rightWingF = round([rightFBelowF_max rightLargerF]*100);
    rightWingR1 = [rightR1_LwBelowF_max rightR1SurfForLargeF];
    leftWingF = round([leftFBelowF_max leftLargerF]*100);
    leftWingR1 = [leftR1_LwBelowF_max leftR1SurfForLargeF];
    
    averagedF = [];
    averagedR1 = [];
    for fCounter = 1:length(fForAverage)
        f = fForAverage(fCounter);
        foundCounter = 0;
        summedR1 = 0;
        if any(rightWingF == f)
            summedR1 = rightWingR1(rightWingF == f);
            foundCounter = foundCounter + 1;
        end
        if any(leftWingF == f)
            summedR1 = summedR1 + leftWingR1(leftWingF == f);
            foundCounter = foundCounter + 1;
        end
        
        if foundCounter > 0
            averagedF(end+1) = f/100; %#ok<SAGROW>
            averagedR1(end+1) = summedR1/foundCounter; %#ok<SAGROW>
        end
        
    end
    
    plot(averagedF,averagedR1);
    
    legend("right surf","right free","left surf","left free" ...
        ,"averaged","Location","east");
    ylabel("R$_{1,LW}$ [Hz]");
    xlabel("Pool fraction factor f");
    
    if saving
        saveFigureTo(matlabSimResultsFolder + "plots" + filesep ...
            ,lipid,corrFuncsData.matlabSimulationDate + "_" ...
            + corrFuncsData.gromacsSimulationDate ...
            ,"poolFractionDependentR1",true);
        save(matlabSimResultsFolder  ...
            + corrFuncsData.matlabSimulationDate + "_" ...
            + corrFuncsData.gromacsSimulationDate + "_" + lipid ...
            + "_poolFractionDependentLipidWaterR1.mat" ...
            ,'averagedF','averagedR1');
            
    end
    
end

function [R1_LwBelowF_max,fBelowF_max,r1SurfForLargeF ...
    ,largerF] = determineR1SurfForWing(whichWing ...
    ,waterHDistrib,lipidHDistrib,locations,hCountBorder ...
    ,largestPossiblePoolFractionFactor,surfaceWaterR1,freeWaterR1)
fprintf(" - Analysing the %s wing.\n",whichWing);

overallLipidHAtomsInWing = sum(lipidHDistrib);
overallWaterHAtomsInWing = sum(waterHDistrib);

f_sim = overallLipidHAtomsInWing/overallWaterHAtomsInWing;
fprintf("  f_sim: %.4f \n",f_sim);

surfaceWaterHDistrib = waterHDistrib .* ( ...
    lipidHDistrib > hCountBorder & waterHDistrib > hCountBorder);
whereIsSurfWater = surfaceWaterHDistrib > 0;
borders = find(diff(whereIsSurfWater) ~= 0);
borders = [borders(1) borders(end)];
overallHsInSurfacePool = sum(surfaceWaterHDistrib);
surfaceHCumulated = cumsum(surfaceWaterHDistrib,'reverse');

freeWaterHDistrib = zeros(1,length(waterHDistrib));
freeWaterHDistrib(1:borders(1)) = ...
    waterHDistrib(1:borders(1));

plot(locations(surfaceWaterHDistrib ~= 0) ...
    ,surfaceWaterHDistrib(surfaceWaterHDistrib ~= 0));
plot(locations(freeWaterHDistrib ~= 0) ...
    ,freeWaterHDistrib(freeWaterHDistrib ~= 0));
plot(locations(1:end-2),lipidHDistrib(1:end-2));

f_max = round((overallLipidHAtomsInWing/overallHsInSurfacePool)*100)/100;
fprintf("  f_max: %.4f \n",f_max);

fBelowF_max = 0:0.01:f_max;
f_surf = zeros(1,length(fBelowF_max));
f_free = zeros(1,length(fBelowF_max));

for fCounter = 1:length(fBelowF_max) 
    f = fBelowF_max(fCounter);
    f_surf(fCounter) = overallHsInSurfacePool / overallLipidHAtomsInWing * f;
    f_free(fCounter) = 1-f_surf(fCounter);
end

R1_LwBelowF_max = f_surf * surfaceWaterR1 + f_free * freeWaterR1;

fAboveF_max = f_max + 0.01 : 0.01 : largestPossiblePoolFractionFactor;
deltaZ = locations(borders(2)) - locations(borders(1)+1);

largerF = [];
r1SurfForLargeF = [];
for fCounter = 1:length(fAboveF_max)
    f = fAboveF_max(fCounter);
    howManyWaterHAtoms = overallLipidHAtomsInWing / f;
    
    % right wing
    lastHigher = find(surfaceHCumulated > howManyWaterHAtoms,1,'last');
    if isempty(lastHigher)
        continue;
    end
    z_null = locations(lastHigher)-locations(borders(1)+1);
    largerF(end+1) = f;   %#ok<AGROW>
    r1SurfForLargeF(end+1) = (freeWaterR1*(deltaZ - z_null) ...
        - (surfaceWaterR1 - freeWaterR1)*(deltaZ - z_null^3/deltaZ^2) ...
        +(3/2*surfaceWaterR1 - 2*freeWaterR1)*(deltaZ - z_null^2/deltaZ)) ...
        / (1/2 * deltaZ - z_null + 1/2 * z_null^2/deltaZ);  %#ok<AGROW>
    
end

end

