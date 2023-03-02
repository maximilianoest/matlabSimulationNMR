clc; clear all; close all;

% Description:
% Based on the number of hydrogen atoms in the water pool and in the solid
% lipid pool and the distributions of H atoms, the surface water fraction
% and free water fraction is determined. With this information the pool
% fraction dependent R1 of the lipid water pool is determined.
% For PLPC, the maximum pool fraction factor in
% lipidWater_poolFractionFactorDependentR1.m is exceeded. Thus, the number
% of water molecules at the surface need to be reduced. Therefore, the
% maximum pool fraction factor can be increased and R1_surf have to be
% adapted. Thus, it is assumed that the lipid head of PLPC come very close
% together and the water between them is decreased. Furthermore, a table is
% created where R1_surf is given for a f_max


%% add directories
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

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

savingPath = sprintf("..%s..%sRESULTS%sdensityDistributions%s" ...
    ,createFilesepStringArray(4));
saving = 0;

%% configuration

% set border where lipid part is cutted
lipidHCountBorders = [0:0.1:90]; %#ok<NBRAK>

% based on Tsukiashi2018
% T1 vs. temperature: y = 0.058071*x + 1.7623 with x = Temperature in °C
temperature = 18; % °C
freeWaterR1 =  1/(0.058071 * temperature + 1.7623);
fprintf("Free water R1: %.4f \n", freeWaterR1);

% % equation for determining surface water R1
% % R1_tot = f_surf * R1_surf + f_free * R1_free
% % => R1_surf = R1_tot / f_surf - f_free / f_surf * R1_free    (1)
% % 
% % equation to determine lipid water R1 based on pool fraction factor:
% % R1_tot = f_surf * R1_surf + f_free * R1_free
% % 
% % f_free(f) = #_H,free / #_H,LW = 1 - #_H,surf / #_H,SL * f
% % f_surf(f) = #_H,surf / #_H,LW = #_H,surf / #_H,SL * f = 1 - f_free(f)
% % 
% % #_H,i: number of atoms in i pool or part of pool i
% % f = M_0,SL / M_0,LW = #_H,SL / #_H,LW : pool fraction factor
% formula = "R1_LW(f) = (1 - #_H,surf / #_H,SL * f) * R1_free" ...
%     + " + #_H,surf / #_H,SL * f * R1_surf \n";
% fprintf("Fromula: " + formula);

% fieldstrength
fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;

% % information in subplot
% fourthSubplotText = sprintf("Field strength: %.2f Tesla\n" ...
%     + "Lipid H count border: %.2f \n", fieldStrength,lipidHCountBorder);

fig = initializeFigure();
legendEntries = {};
for lipidNr = 1:length(lipidNames)
    distributionData = load(filePaths(lipidNr));
    lipid = lipidNames(lipidNr);
    fprintf("Lipid %s\n",lipidNames(lipidNr));
    
    avgDensityDistributionHAtoms = squeeze(mean( ...
        distributionData.densityDistributionOnlyHdydrogen,1));
    
    hFrequencyDistributionInWaterPool = avgDensityDistributionHAtoms( ...
        distributionData.moleculesToCompare == "SOL",:) ...
        /distributionData.atomWeightCollection( ...
        distributionData.atomNameCollection == "H") ...
        * mean(distributionData.dVolumeForTimeStep);
    hFrequencyDistributionInLipidPool = avgDensityDistributionHAtoms( ...
        distributionData.moleculesToCompare == lipid,:) ...
        /distributionData.atomWeightCollection( ...
        distributionData.atomNameCollection == "H") ...
        * mean(distributionData.dVolumeForTimeStep);
    
    overallHsInWaterPool = sum(hFrequencyDistributionInWaterPool);
    overallHsInLipidPool = sum(hFrequencyDistributionInLipidPool);
    poolFractionFactorInSim = overallHsInLipidPool/overallHsInWaterPool;
    fprintf("  Pool fraction factor in MD sim: %.4f\n" ...
        ,poolFractionFactorInSim);
    
    corrFuncsData = load(corrFuncDirectories(lipidNr));
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
    
    % lipid H border count is set to 0
    hFrequencyInSurfaceWater = hFrequencyDistributionInWaterPool ...
        .* (hFrequencyDistributionInLipidPool > 0 ...
        & hFrequencyDistributionInWaterPool > 0);
    hFrequencyInFreeWater = hFrequencyDistributionInWaterPool ...
        .* ~(hFrequencyDistributionInLipidPool > 0 ...
        & hFrequencyDistributionInWaterPool > 0);
    
    surfaceInteractionFraction = sum(hFrequencyInSurfaceWater) ...
        / overallHsInWaterPool;
    freeWaterInteractionFraction = sum(hFrequencyInFreeWater) ...
        / overallHsInWaterPool;
     fprintf("  Surface/Free water fraction: %.4f / %.4f\n" ...
        ,surfaceInteractionFraction,freeWaterInteractionFraction);
    
    surfaceWaterR1 = totalWaterPoolR1 / surfaceInteractionFraction ...
        - freeWaterInteractionFraction / surfaceInteractionFraction ...
        * freeWaterR1;
    
    fprintf("  Surface water R1 = %.4f \n",surfaceWaterR1);
    
    % ---- control whether function (1) is correct
    f_free = (1 - sum(hFrequencyInSurfaceWater) ...
        / overallHsInLipidPool * poolFractionFactorInSim);
    f_surf = sum(hFrequencyInSurfaceWater) ...
        / overallHsInLipidPool * poolFractionFactorInSim;
    overallWaterR1_control = f_surf * surfaceWaterR1 ...
        +  f_free * freeWaterR1;
    if abs(overallWaterR1_control - totalWaterPoolR1) > 0.00001
        error("Results are not equal");
    end
    filledFormula = sprintf("R1_LW(f) = (1 - %.4f * f) * %.4f + %.4f " ...
        + "* f * %.4f" ...
        ,sum(hFrequencyInSurfaceWater)/overallHsInLipidPool ...
        ,freeWaterR1,sum(hFrequencyInFreeWater) ...
        /overallHsInLipidPool,surfaceWaterR1);
    fprintf("  Formula: %s \n",filledFormula);
    
    maximumPoolFractionFactorWithFreeWater = overallHsInLipidPool ...
        /sum(hFrequencyInSurfaceWater);
    fprintf("  Largest possible pool fraction factor: f_max = %.4f" ...
        + " for zero H border count\n" ...
        ,maximumPoolFractionFactorWithFreeWater);
    fprintf("  There are %.0f H atoms in the surface region and %.0f H" ...
        + " atoms in the solid lipid region. \n" ...
        ,sum(hFrequencyInSurfaceWater),overallHsInLipidPool);
    
    
    % reducing surface water step by step and determine according pool
    % fraction factor and R1_surf. The equation of surface water and free
    % water still assumes that there is only free water and surface water,
    % although some of the free water interacts with the lipid.
    poolFractionFactorsForReducedSurfaceWater = ...
        zeros(1,length(lipidHCountBorders));
    surfaceWaterR1ForReducedSurfaceWater = ...
        zeros(1,length(lipidHCountBorders));
    for borderIndex = 1:length(lipidHCountBorders)
        lipidHCount = lipidHCountBorders(borderIndex);
        hFrequencyInSurfaceWater = hFrequencyDistributionInWaterPool ...
            .* (hFrequencyDistributionInLipidPool > lipidHCount ...
            & hFrequencyDistributionInWaterPool > 0);
        hFrequencyInFreeWater = hFrequencyDistributionInWaterPool ...
            .* ~(hFrequencyDistributionInLipidPool > lipidHCount ...
            & hFrequencyDistributionInWaterPool > 0);
        
        surfaceWaterInteractionFraction ...
            = sum(hFrequencyInSurfaceWater)/ overallHsInWaterPool;
        freeWaterInteractionFraction = sum(hFrequencyInFreeWater) ...
            / overallHsInWaterPool;
        
        surfaceWaterR1ForReducedSurfaceWater(borderIndex) =  ...
            totalWaterPoolR1 / surfaceWaterInteractionFraction ...
            - freeWaterInteractionFraction ...
            / surfaceWaterInteractionFraction * freeWaterR1;
        
        poolFractionFactorsForReducedSurfaceWater(borderIndex) ...
            = overallHsInLipidPool/sum(hFrequencyInSurfaceWater);
        
    end
    
    initializeSubplot(fig,2,2,lipidNr);
    plot(lipidHCountBorders,poolFractionFactorsForReducedSurfaceWater);
    plot(lipidHCountBorders,surfaceWaterR1ForReducedSurfaceWater)
    title(sprintf("Lipid: %s",lipid));
    legend("pool fraction factor", "R1$_{surf}$")
    xlabel("Number of H atoms in lipid");
    ylabel("Value of parameter");
    lgd = legend();
    lgd.Location = "northwest";
    
    initializeSubplot(fig,2,2,4);
    plot(poolFractionFactorsForReducedSurfaceWater ...
        ,surfaceWaterR1ForReducedSurfaceWater);
    legendEntries{end+1} = lipid; %#ok<SAGROW>
    
       
    save(sprintf("..%s..%sResults%slipidWater_poolFractionFactors%s%s" ...
        + "_%s_R1surfForGivenPoolFractionFactor.mat" ...
        ,createFilesepStringArray(4),datestr(now,"yyyymmdd"),lipid) ...
        ,'lipid','poolFractionFactorsForReducedSurfaceWater' ...
        ,'surfaceWaterR1ForReducedSurfaceWater','filledFormula');

    
end
lgd = legend(legendEntries);
xlabel("Pool fraction factor")
ylabel("R1$_{surf}$");
lgd.Location = "northwest";

saveFigureTo(sprintf("..%s..%sResults%slipidWater_" ...
        + "poolFractionFactors%s",createFilesepStringArray(4)) ...
        ,"allLipidWaters",datestr(now,"yyyymmdd") ...
        ,"poolFractionFactorDependentSurfaceR1",false);


% ax = initializeSubplot(fig,2,2,4);
% text(0.05,0.5,fourthSubplotText,'interpreter','latex');
% set(ax,'visible','off');
% lgd = legend();
% lgd.Visible = 'Off';
% 
% if saving
%     saveFigureTo(savingPath,"allLipids",datestr(now,"yyyymmdd") ...
%         ,sprintf("poolFractionDependentR1_lipidBorderHCount%d" ...
%         ,lipidHCountBorder));
% end










