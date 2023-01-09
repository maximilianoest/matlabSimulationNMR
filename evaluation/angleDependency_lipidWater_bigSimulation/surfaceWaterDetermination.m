clc; clear all; close all;
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile("constants.txt");
saving = 1;

fieldStrength = 3; % Tesla
gyromagneticRatio = constants.gyromagneticRatioOfHydrogenAtom;
omega0 = fieldStrength * gyromagneticRatio;

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 20; % °C
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);

densityDirectory = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\densityDistributions\Plots\";
lipids = ["DOPS" "PLPC" "PSM"];

correlationDirectory = "C:\Users\maxoe\Google Drive\Promotion" ...
    + "\Simulation\RESULTS\angleDependency_lipidWater_bigSimulation\";
resultsFiles = [ ...
    "20221101_Results_DOPSwater_20220110_DOPS_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ...
    "20221025_Results_PLPCwater_20220804_PLPC_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ... 
    "20221028_Results_PSMwater_20220804_PSM_TIP4_Monolayer_50water_water_H_whole_dt02ps_simTime25ns" ] ...
    + ".mat";
fileNames = "20221122densityDistribution_" + lipids ...
    + "_withInteractionFractions" + ".mat";

fourthSubplotText = "";

for lipidNr = 1:length(lipids)
    load(densityDirectory + fileNames(lipidNr));
    hydrogenAtomsInSurfaceRegion = numberOfHsInWaterPool ...
        .* (numberOfHsInLipidPool > 0.1 & numberOfHsInWaterPool > 0.1);
    hydrogenAtomsInFreeWaterRegion = numberOfHsInWaterPool ...
        .* ~(numberOfHsInLipidPool > 0.1 & numberOfHsInWaterPool > 0.1);
    
    fig = initializeFigure();
    ax = initializeSubplot(fig,2,2,1);
    plot(avgLocations,hydrogenAtomsInSurfaceRegion);
    plot(avgLocations,hydrogenAtomsInFreeWaterRegion);
    xticks([-5e-9 : 1e-9 : 5e-9]); %#ok<NBRAK>
    xticklabels([-5 : 1 : 5]);
    xlabel("Postions [nm]")
    legend("Surface region", "Free water region");
    title("Number of H atoms");
    
    numberOfHAtomsInSurfaceRegion = sum(hydrogenAtomsInSurfaceRegion);
    numberOfHAtomsInFreeWaterRegion = sum(hydrogenAtomsInFreeWaterRegion);
    
    surfaceWaterFraction = numberOfHAtomsInSurfaceRegion ...
        / sum(numberOfHsInWaterPool);
    freeWaterFraction = numberOfHAtomsInFreeWaterRegion ...
        / sum(numberOfHsInWaterPool);
    
    fprintf("Lipid %s: \n   f_SW = %.4f ; f_FW = %.4f \n" ...
        ,lipids(lipidNr),surfaceWaterFraction,freeWaterFraction);
    
    load(correlationDirectory + resultsFiles(lipidNr));
    corrFuncFirstOrder = sumCorrFuncFirstOrder.nearestNeighbours8000 ...
        / atomCounter;
    corrFuncSecondOrder = sumCorrFuncSecondOrder.nearestNeighbours8000 ...
        / atomCounter;
    ax = initializeSubplot(fig,2,2,2);
    timeAxis = [0:deltaTInS:(length(corrFuncFirstOrder)-1) * deltaTInS];
    plot(timeAxis,abs(corrFuncFirstOrder));
    plot(timeAxis,abs(corrFuncSecondOrder));
    legend("1st order", "2nd order");
    title(sprintf("abs(Correlation functions)"));
    
    [specDensFirstOrder, specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder,corrFuncSecondOrder ...
        ,omega0,deltaTInS,[0 1]);
    ax = initializeSubplot(fig,2,2,3);
    plot(timeAxis,real(specDensFirstOrder));
    plot(timeAxis,real(specDensSecondOrder));
    title("real(Spectral density)");
    legend("1st order", "2nd order");
    xlabel("Upper integration limit");
    
    dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
        *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
        /(constants.nanoMeter^6);
    overallPoolR1 = calculateR1WithSpectralDensity( ...
        specDensFirstOrder(round(0.9*end):end),specDensSecondOrder(round(0.9*end):end) ...
        ,dipolDipolConstant);
    
    fprintf("   overall pool R1 = %.4f \n",overallPoolR1);
    
    surfaceWaterR1 = overallPoolR1 / surfaceWaterFraction ...
        - freeWaterFraction / surfaceWaterFraction * freeWaterR1;
    
    fprintf("   surface water R1 = %.4f \n",surfaceWaterR1);
    lipids = ["DOPS" "PLPC" "PSM"];
    
    fourthSubplotText = sprintf("Lipid: %s\n  surface water fraction: " ...
       +"%.4f \n  free water fraction: %.4f \n" ...
       +"  overall R1: %.4f\n  surface R1: %.4f\n  free water R1: %.4f\n\n" ...
       ,lipids(lipidNr),surfaceWaterFraction,freeWaterFraction ...
       ,overallPoolR1,surfaceWaterR1,freeWaterR1);
    
    ax = initializeSubplot(fig,2,2,4);
    text(0.05,0.5,fourthSubplotText,'interpreter','latex');
    set(ax,'visible','off');
    lgd = legend();
    lgd.Visible = 'Off';
    
    if saving
        saveFigureTo(densityDirectory,"surfaceR1",lipids(lipidNr) ...
            ,matlabSimulationDate);
    end
    
end





