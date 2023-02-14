clc; clear all; close all;

% Description:
% Based on the number of hydrogen atoms in the water pool and in the solid
% lipid pool and the distributions of H atoms, the surface water fraction
% and free water fraction is determined. With this information the pool
% fraction dependent R1 of the lipid water pool is determined.
% Basically, the number of surface water H atoms stays constant if the
% number of H atoms in the lipids stays constant. Thus, the free water
% fraction must be somehow inversely proportional to the pool fraction
% factor.
%
% 1. get distributions (cacluated by createDensityDistributions.m)
% 2. 

%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 18; % °C
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);

% equation for determining surface water R1
% R1_tot = f_surf * R1_surf + f_free * R1_free
% => R1_surf = R1_tot / f_surf - f_free / f_surf * R1_free    (1)
% 
% equation to determine lipid water R1 based on pool fraction factor:
% R1_tot = f_surf * R1_surf + f_free * R1_free
% 
% f_free(f) = #_H,free / #_H,LW = 1 - #_H,surf / #_H,SL * f
% f_surf(f) = #_H,surf / #_H,LW = #_H,surf / #_H,SL * f = 1 - f_free(f)
% 
% #_H,i: number of atoms in i pool or part of pool i
% f = M_0,SL / M_0,LW = #_H,SL / #_H,LW : pool fraction factor
fprintf("Fromula: R1_LW(f) = (1 - #_H,surf / #_H,SL * f) * R1_free" ...
    + " + #_H,surf / #_H,SL * f * R1_surf \n");

% fieldstrength
fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;
%% data
resultsDirectory = sprintf( ...
    "..%s..%sRESULTS%sdensityDistributions%sPlots%s" ...
    ,createFilesepStringArray(5));
corrFuncDirectories = sprintf("..%s..%sRESULTS" ...
    + "%sangleDependency_lipidWater_bigSimulation%s" ...
    ,createFilesepStringArray(4));
corrFuncDirectories = corrFuncDirectories ...
    + ["20221101_Results_DOPSwater_20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221025_Results_PLPCwater_20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221028_Results_PSMwater_20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns"] ...
    + ".mat";

lipidNames = ["DOPS","PLPC","PSM"];
gromacsSimulationDate = "20221122";
fileNames = gromacsSimulationDate + "densityDistribution_" ...
    + lipidNames + ".mat";

filePaths = resultsDirectory + fileNames;
checkIfAllFilesExist(filePaths);

fig = initializeFigure();
for lipidNr = 1:length(lipidNames) 
    distributionData = load(filePaths(lipidNr));
    lipid = lipidNames(lipidNr);
    avgDensityDistributionHAtoms = squeeze(mean( ...
        distributionData.densityDistributionOnlyHdydrogen,1));
    numberOfHsInWaterPool = avgDensityDistributionHAtoms( ...
        distributionData.moleculesToCompare == "SOL",:) ...
        /distributionData.atomWeightCollection( ...
        distributionData.atomNameCollection == "H") ...
        * mean(distributionData.dVolumeForTimeStep);
    numberOfHsInLipidPool = avgDensityDistributionHAtoms( ...
        distributionData.moleculesToCompare == lipid,:) ...
        /distributionData.atomWeightCollection( ...
        distributionData.atomNameCollection == "H") ...
        * mean(distributionData.dVolumeForTimeStep);
    
    numberOfWaterHsInSurfaceRegion = numberOfHsInWaterPool ...
        .* (numberOfHsInLipidPool > 0 & numberOfHsInWaterPool > 0);
    numberOfWaterHsInFreeRegion = numberOfHsInWaterPool ...
        .* ~(numberOfHsInLipidPool > 0 & numberOfHsInWaterPool > 0);
    
    
    ax = initializeSubplot(fig,2,2,lipidNr);
    plot(distributionData.avgLocations,numberOfHsInLipidPool);
    plot(distributionData.avgLocations,numberOfHsInWaterPool);
    plot(distributionData.avgLocations,numberOfWaterHsInSurfaceRegion,"--");
    plot(distributionData.avgLocations,numberOfWaterHsInFreeRegion,"--");
    lgd = legend("lipid pool", "water pool", "surface water", "free water");
    lgd.Location = "south";
    xlabel("Location [m]");
    ylabel("Number of atoms");
    title(sprintf("%s",lipid));
    
    
    overallHsInWaterPool = sum(numberOfHsInWaterPool);
    overallHsInLipidPool = sum(numberOfHsInLipidPool);
    
    surfaceInteractionFraction = sum(numberOfWaterHsInSurfaceRegion) ...
        / overallHsInWaterPool;
    freeWaterInteractionFraction = sum(numberOfWaterHsInFreeRegion) ...
        / overallHsInWaterPool;
    
    poolFractionFactor = overallHsInLipidPool/overallHsInWaterPool;
    
    fprintf("Lipid %s\n",lipidNames(lipidNr));
    fprintf("  Surface/Free water fraction: %.4f / %.4f\n" ...
        ,surfaceInteractionFraction,freeWaterInteractionFraction);
    fprintf("  Pool fraction factor for simulation: %.4f\n" ...
        ,poolFractionFactor);
    
    corrFuncsData = load(corrFuncDirectories(lipidNr));
    deltaTInS = corrFuncsData.deltaTInS;
    corrFuncFirstOrder = corrFuncsData.sumCorrFuncFirstOrder ...
        .nearestNeighbours8000 / corrFuncsData.atomCounter;
    corrFuncSecondOrder = corrFuncsData.sumCorrFuncSecondOrder ...
        .nearestNeighbours8000 / corrFuncsData.atomCounter;
    
    [specDensFirstOrder, specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder,corrFuncSecondOrder ...
        ,omega0,deltaTInS,[0.9 1]);
    
    dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
        *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
        /(constants.nanoMeter^6);
    overallPoolR1 = calculateR1WithSpectralDensity( ...
        specDensFirstOrder,specDensSecondOrder,dipolDipolConstant);
    
    fprintf("  Overall pool R1 = %.4f \n",overallPoolR1);
    
    surfaceWaterR1 = overallPoolR1 / surfaceInteractionFraction ...
        - freeWaterInteractionFraction / surfaceInteractionFraction ...
        * freeWaterR1;
    
    fprintf("  Surface water R1 = %.4f \n",surfaceWaterR1);
    
    % ---- control whether function (1) is correct
    f_free = (1 - sum(numberOfWaterHsInSurfaceRegion) ...
        / overallHsInLipidPool * poolFractionFactor);
    f_surf = sum(numberOfWaterHsInSurfaceRegion) ...
        / overallHsInLipidPool * poolFractionFactor; 
    overallWaterR1_control = f_surf * surfaceWaterR1 ...
        +  f_free * freeWaterR1;
    if abs(overallWaterR1_control - overallPoolR1) > 0.00001
       error("Results are not equal");
    end
    fprintf("  Formula: R1_LW(f) = (1 - %.4f * f) * %.4f + %.4f * f * %.4f \n" ...
        ,sum(numberOfWaterHsInSurfaceRegion)/overallHsInLipidPool ...
        ,freeWaterR1,sum(numberOfWaterHsInSurfaceRegion) ...
        /overallHsInLipidPool,surfaceWaterR1);
    fprintf("  Largest possible pool fraction factor: f_max = %.4f \n" ...
        ,overallHsInLipidPool/sum(numberOfWaterHsInSurfaceRegion));
    
    
    fprintf("\n\n");
    
    
    
    

end




%% functions

function checkIfAllFilesExist(filePaths)

if isempty(filePaths)
    error("No file paths were given.");
end

for filePath = filePaths
    if ~exist(filePath,'file')
        error("File: %s does not exist.",filePath);
    end
end
end









