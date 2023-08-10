clc; clear all; close all; fclose('all');

%% paths


addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath("../molecularDynamicSimulationAnalysis/"));
addpath(genpath("../distributionAnalysis/"));
addpath(genpath("../brainConstitutionAnalysis/"));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));

runAllScriptsAgain = 0;
if runAllScriptsAgain
    run("brainConstitution"); %#ok<UNRCH>
    
    run("determineR1OfMyelinBasedOnCompartments");
    run("determineSurfaceWaterAndHistologR1InMyelin");
end


resultsPath = ...
    "C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\wholeBrain\";
constants = readConstantsFile('constants.txt');
saving = 1;

axisFontSize = 10;
legendFontSize = 10;

dataPointsCountLowerLimit = 1;
fieldStrengthsToIgnoreInWM = 0;
fieldStrengthsToIgnoreInGM = 0;

minimumDataPointsForFitting = 2;
fieldStrengthAxis = 0.01:0.01:7.1;
highestR1InPlot = 20;
lowestFieldStrengthForMdSimFit = 0.349;

lowFieldLimit = 0.75;

% ==== relaxation rates ====
csvR1DataFilePath = resultsPath + "ObservedT1Values.CSV";
observedRelaxationTimesTable = readtable(csvR1DataFilePath);

r1DataFolderName = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\wholeMyelin_relaxationRates\";
r1DataFileName = "solidMyelinAndMyelinWater_histCompartmentAndCrossR1";
r1Data = load(r1DataFolderName + r1DataFileName);
lowerBoundR1_SPFraction = 0.01;
upperBoundR1_SPFraction = 1.5;
r1_SPSamplingPoints = 500;

% ==== exchange rates ====
csvExchangeRateDataFilePath = resultsPath + "ObservedExchangeRates.CSV";
exchangeRatesTable = readtable(csvExchangeRateDataFilePath);

exchangeRateIncreaseFactor = 1.10;
exchangeRateDecreaseFactor = 0.90;

% ==== constitution data ====
wmWaterFraction = r1Data.constitutionData.wmWaterContent;
wmNonWaterFraction = r1Data.constitutionData.wmNonWaterContent;
wmLipidFraction = wmNonWaterFraction ...
    * r1Data.constitutionData.wmLipidContent;
wmProteinFraction = wmNonWaterFraction ...
    * r1Data.constitutionData.wmProteinContent;

gmWaterFraction = r1Data.constitutionData.gmWaterContent;
gmNonWaterFraction = r1Data.constitutionData.gmNonWaterContent;
gmLipidFraction = gmNonWaterFraction ...
    * r1Data.constitutionData.gmLipidContent;
gmProteinFraction = gmNonWaterFraction ...
    * r1Data.constitutionData.gmProteinContent;

% ==== concluded data ====
ESIVR = r1Data.ESIVR ...
    * r1Data.fractionOfSurfaceWaterRelativeToHistMyelinWaterAmount;
wmSurfaceWaterFraction = wmNonWaterFraction * ESIVR;
wmLipidWaterFraction = wmLipidFraction * ESIVR;
wmProteinWaterFraction = wmProteinFraction * ESIVR;
wmFreeWaterFraction = wmWaterFraction - wmLipidWaterFraction ...
    - wmProteinWaterFraction;

gmSurfaceWaterFraction = gmNonWaterFraction * ESIVR;
gmLipidWaterFraction = gmLipidFraction * ESIVR;
gmProteinWaterFraction = gmProteinFraction * ESIVR;
gmFreeWaterFraction = gmWaterFraction - gmLipidWaterFraction ...
    - gmProteinWaterFraction;

%% investigate relaxation rates in literature
[wmFieldStrengthArray,r1WMCollection,~,~] ...
   = getCollectionsFromTable(observedRelaxationTimesTable ...
    ,dataPointsCountLowerLimit,fieldStrengthsToIgnoreInWM);
    
[~,~,gmFieldStrengthArray,r1GMCollection] ...
    = getCollectionsFromTable(observedRelaxationTimesTable ...
    ,dataPointsCountLowerLimit,fieldStrengthsToIgnoreInGM);
showDistributionOfValues(wmFieldStrengthArray,r1WMCollection ...
    ,gmFieldStrengthArray,r1GMCollection,axisFontSize,legendFontSize);

if saving
    saveFigureTo(r1DataFolderName,datestr(now,'yyyymmdd') ...
        ,"distributionOfR1Values","Literature");
end

r1WMCollection = ignoreOutliers(r1WMCollection,0.25,0.75);
r1GMCollection = ignoreOutliers(r1GMCollection,0.25,0.75);

r1Data.wmFieldStrengthArray = wmFieldStrengthArray;
r1Data.r1WMColletion = r1WMCollection;
r1Data.gmFieldStrengthArray = gmFieldStrengthArray;
r1Data.r1GMCollection = r1GMCollection;

[r1WMAvg,r1WMSTD,r1WMCounts] = getAveragSTDandCountsOfR1Collection( ...
    r1WMCollection,observedRelaxationTimesTable,wmFieldStrengthArray);
r1Data.r1WMAvg = r1WMAvg;
r1Data.r1WMSTD = r1WMSTD;
r1Data.r1WMCounts = r1WMCounts;

[r1GMAvg,r1GMSTD,r1GMCounts] = getAveragSTDandCountsOfR1Collection( ...
    r1GMCollection,observedRelaxationTimesTable,gmFieldStrengthArray);
r1Data.r1GMAVg = r1GMAvg;
r1Data.r1GMSTD = r1GMSTD;
r1Data.r1GMCounts = r1GMCounts;

% ==== parameter set up ====
% -> time set up
lowFieldDeltaT = 1e-4;
highFieldDeltaT = 1e-3;
lowFieldDuration = 6;
highFieldDuration = 20;
deltaTsForWM = ones(1,length(r1WMCollection));
deltaTsForWM(wmFieldStrengthArray > 1.6) = highFieldDeltaT;
deltaTsForWM(~(wmFieldStrengthArray > 1.6)) = lowFieldDeltaT;
durationsForWM = ones(1,length(r1WMCollection));
durationsForWM(wmFieldStrengthArray > 1.6) = highFieldDuration;
durationsForWM(~(wmFieldStrengthArray > 1.6)) = lowFieldDuration;
timeStepCountForWM = durationsForWM./deltaTsForWM;

deltaTsForGM = ones(1,length(r1GMCollection));
deltaTsForGM(gmFieldStrengthArray > 1.6) = highFieldDeltaT;
deltaTsForGM(~(gmFieldStrengthArray > 1.6)) = lowFieldDeltaT;
durationsForGM = ones(1,length(r1GMCollection));
durationsForGM(gmFieldStrengthArray > 1.6) = highFieldDuration;
durationsForGM(~(gmFieldStrengthArray > 1.6)) = lowFieldDuration;
timeStepCountForGM = durationsForGM./deltaTsForGM;

%% fitting of literature values (only >= 0.35T)
wmFieldStrengthArrayForFitting = wmFieldStrengthArray(r1WMCounts ...
    >= minimumDataPointsForFitting);
r1WMCollectionForFitting = r1WMCollection(r1WMCounts  ...
    >= minimumDataPointsForFitting);

r1Data.wmFieldStrengthArrayForFitting = wmFieldStrengthArrayForFitting;
r1Data.r1WMCollectionForFitting = r1WMCollectionForFitting;

[r1WMAvgForFitting,r1WMStdForFitting,~] ...
    = getAveragSTDandCountsOfR1Collection(r1WMCollectionForFitting);
    
[wmFactor,wmExponent,wmCovariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    wmFieldStrengthArrayForFitting,r1WMAvgForFitting,[0.5 0.5] ...
    ,r1WMStdForFitting.^2);
r1Data.wmFactor = wmFactor;
r1Data.wmExponent = wmExponent;
r1Data.wmCovariance = wmCovariance;

gmFieldStrengthArrayForFitting = gmFieldStrengthArray(r1GMCounts ...
    >= minimumDataPointsForFitting);
r1GMCollectionForFitting = r1GMCollection(r1GMCounts ...
    >= minimumDataPointsForFitting);

r1Data.gmFieldStrengthArrayForFitting = gmFieldStrengthArrayForFitting;
r1Data.r1GMCollectionForFitting = r1GMCollectionForFitting;

[r1GMAvgForFitting,r1GMStdForFitting,~] ...
    = getAveragSTDandCountsOfR1Collection(r1GMCollectionForFitting);

[gmFactor,gmExponent,gmCovariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    gmFieldStrengthArrayForFitting,r1GMAvgForFitting,[0.5 0.5] ...
    ,r1GMStdForFitting.^2);
r1Data.gmFactor = gmFactor;
r1Data.gmExponent = gmExponent;
r1Data.gmCovariance = gmCovariance;

fprintf("<strong>FITTING RESULTS OF TISSUE R1: </strong> \n");
fprintf("WM: %.3f B_0 ^(-%.3f) \n variance/covariance matrix:\n" ...
    ,wmFactor,wmExponent) 
fprintf(" %.7f %.7f \n",wmCovariance);
fprintf("GM: %.3f B_0 ^(-%.3f) \n variance/covariance matrix:\n" ...
    ,gmFactor,gmExponent) 
fprintf(" %.7f %.7f \n",gmCovariance);

%% plotting literature values (Observed R1 in WM and GM)
literatureFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
literatureLegendEntries = {};

wmSubplot = initializeSubplot(literatureFigure,2,1,1);
ylabel("R$_{1,WM}$ [Hz]");
xticks(0:0.5:7);

literaturePlt(1) = errorbar(wmFieldStrengthArray,r1WMAvg,r1WMSTD ...
    ,'LineStyle','none','LineWidth',1.5,'Color','g');
literatureLegendEntries{end+1} = "R$_{1,WM}$ $\pm$ STD";
literaturePlt(2) = errorbar(wmFieldStrengthArrayForFitting ...
    ,r1WMAvgForFitting,r1WMStdForFitting,'LineStyle','none' ...
    ,'LineWidth',1.5);
literatureLegendEntries{end+1} = "R$_{1,WM}$ $\pm$ STD for fit";
literaturePlt(3) = plot(fieldStrengthAxis,wmFactor*fieldStrengthAxis ...
    .^(-wmExponent));
literatureLegendEntries{end+1} = "fitted R$_{1,WM}$ $\pm$ STD";
drawSTDRegionAroundPowerFunction(wmFactor,wmExponent...
    ,wmCovariance,fieldStrengthAxis,literaturePlt(3).Color);

legend(literaturePlt(1:3),literatureLegendEntries);
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);

gmSubplot = initializeSubplot(literatureFigure,2,1,2);
literatureLegendEntries = {};
xlabel("Field strengths [T]");
xticks(0:0.5:7);
ylabel("R$_{1,GM}$ [Hz]");

literaturePlt(4) = errorbar(gmFieldStrengthArray,r1GMAvg,r1GMSTD ...
    ,'LineStyle','none','LineWidth',1.5,'Color','g');
literatureLegendEntries{end+1} = "R$_{1,GM}$ $\pm$ STD";
literaturePlt(5) = errorbar(gmFieldStrengthArrayForFitting ...
    ,r1GMAvgForFitting,r1GMStdForFitting,'LineStyle','none' ...
    ,'LineWidth',1.5);
literatureLegendEntries{end+1} = "R$_{1,WM}$ $\pm$ STD for fit";
literaturePlt(6) = plot(fieldStrengthAxis,gmFactor*fieldStrengthAxis ...
    .^(-gmExponent));
literatureLegendEntries{end+1} = "fitted R$_{1,GM}$ $\pm$ STD";
drawSTDRegionAroundPowerFunction(gmFactor,gmExponent...
    ,gmCovariance,fieldStrengthAxis,literaturePlt(6).Color);

legend(literaturePlt(4:6),literatureLegendEntries);
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);

%% fitting R1_SL and R1_LW to power function
[factorR1_SL,exponentR1_SL,covarianceR1_SL] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    r1Data.fieldStrengths(r1Data.fieldStrengths ...
    > lowestFieldStrengthForMdSimFit) ...
    ,r1Data.r1Eff_SM(r1Data.fieldStrengths ...
    > lowestFieldStrengthForMdSimFit),[0.5 0.5]);
r1Data.factorR1_SL = factorR1_SL;
r1Data.exponentR1_SL = exponentR1_SL;
r1Data.covarianceR1_SL = covarianceR1_SL;

[factorR1_LW,exponentR1_LW,covarianceR1_LW] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    r1Data.fieldStrengths(r1Data.fieldStrengths ...
    > lowestFieldStrengthForMdSimFit),r1Data.r1Hist_MW( ...
    r1Data.fieldStrengths > lowestFieldStrengthForMdSimFit),[0.5 0.5]);
r1Data.factorR1_LW = factorR1_LW;
r1Data.exponentR1_LW = exponentR1_LW;
r1Data.covarianceR1_LW = covarianceR1_LW;

%% plotting R1_SL and R1_LW to figure
mdSimFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
indicesOfInterest = (r1Data.fieldStrengths ...
    >= fieldStrengthAxis(1)) & (r1Data.fieldStrengths ...
    <= fieldStrengthAxis(end));

simResultsSubplot = initializeSubplot(mdSimFigure,2,1,1);
legendEntriesSimResultsSubplot = {};

mdSimPlt(1) = plot(r1Data.fieldStrengths(indicesOfInterest) ...
    ,r1Data.r1Eff_SM(indicesOfInterest));
legendEntriesSimResultsSubplot{end+1} ...
    = "R$^{eff}_{1,SL}$ MD sim.";
mdSimPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SL*fieldStrengthAxis.^(-exponentR1_SL));
drawSTDRegionAroundPowerFunction(factorR1_SL,exponentR1_SL...
    ,covarianceR1_SL,fieldStrengthAxis,mdSimPlt(2).Color);
legendEntriesSimResultsSubplot{end+1} ...
    = "R$^{eff}_{1,SL}$ $\pm$ STD fit";

mdSimPlt(3) = plot(r1Data.fieldStrengths(indicesOfInterest) ...
    ,r1Data.r1Hist_MW(indicesOfInterest));
legendEntriesSimResultsSubplot{end+1} ...
    = "R$^{hist}_{1,LW}$ MD Sim.";
mdSimPlt(4) = plot(fieldStrengthAxis ...
    ,factorR1_LW*fieldStrengthAxis.^(-exponentR1_LW));
drawSTDRegionAroundPowerFunction(factorR1_LW,exponentR1_LW...
    ,covarianceR1_LW,fieldStrengthAxis,mdSimPlt(4).Color);
legendEntriesSimResultsSubplot{end+1} ...
    = "R$^{hist}_{1,LW}$ $\pm$ STD fit ";

ylabel("R$_{1}$ [Hz]");

legend(mdSimPlt(1:4),legendEntriesSimResultsSubplot);
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
drawnow;
%% determine R1_SP with IE model
fprintf("<strong> === INTERMEDIATE EXCHANGE MODEL ===</strong>\n");

exchangeRate = mean(exchangeRatesTable.exchangeRateWholePoolNorm(:));
macroMolecularPoolFraction = mean(...
    exchangeRatesTable.macromolecularPoolfraction(:));

if wmFieldStrengthArrayForFitting == gmFieldStrengthArrayForFitting 
    spFieldStrengthForWMFitting = wmFieldStrengthArrayForFitting;
    spFieldStrengthForGMFitting = gmFieldStrengthArrayForFitting;
    spIndicesForWMFitting = false(1,length(wmFieldStrengthArray));
    spIndicesForGMFitting = false(1,length(gmFieldStrengthArray));
    for fieldStrengthNr = 1:length(spFieldStrengthForWMFitting)
        spIndicesForWMFitting((wmFieldStrengthArray ...
            >= spFieldStrengthForWMFitting(fieldStrengthNr) - 0.0001) ...
            & (wmFieldStrengthArray <= spFieldStrengthForWMFitting( ...
            fieldStrengthNr) + 0.0001)) = true;
        spIndicesForGMFitting((gmFieldStrengthArray ...
            >= spFieldStrengthForGMFitting(fieldStrengthNr) - 0.0001) ...
            & (gmFieldStrengthArray <= spFieldStrengthForGMFitting( ...
            fieldStrengthNr) + 0.0001)) = true;
        
    end
else
    error("Not implemented");
end

% ==== white matter ====
fprintf("<strong>WHITE MATTER </strong>\n");
exchange_waterToSolid_WM = exchangeRate/(wmWaterFraction);
exchange_lipidToWater_WM = exchangeRate ...
    /(wmNonWaterFraction * r1Data.constitutionData.wmLipidContent);
exchange_proteinToWater_WM = exchangeRate ...
    /(wmNonWaterFraction * r1Data.constitutionData.wmProteinContent);
r1Data.exchange_waterToSolid_Wm = exchange_waterToSolid_WM;
r1Data.exchange_lipidToWater_Wm = exchange_lipidToWater_WM;
r1Data.exchange_proteinToWater_Wm = exchange_proteinToWater_WM;

[bestGuessTissueR1_WM_IE,bestGuessR1_SP_WM_IE] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    wmFieldStrengthArray,r1WMCollection,r1Data,wmSurfaceWaterFraction ...
    ,wmWaterFraction,wmLipidFraction,wmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_WM,exchange_lipidToWater_WM ...
    ,exchange_proteinToWater_WM,timeStepCountForWM,deltaTsForWM);
r1Data.r1_SPWM_IEModel = bestGuessR1_SP_WM_IE;
r1Data.bestGuessTissueR1_WM_IEModel = bestGuessTissueR1_WM_IE;

[r1_SP_WMAvg_IE,r1_SP_WMStd_IE,~] = getAveragSTDandCountsOfR1Collection( ...
    bestGuessR1_SP_WM_IE);
r1Data.r1_SP_WMAvg_IEModel = r1_SP_WMAvg_IE;
r1Data.r1_SP_WMStd_IEModel = r1_SP_WMStd_IE;

[factorR1_SP_WM_IE,exponentR1_SP_WM_IE,covarianceR1_SP_WM_IE] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForWMFitting,r1_SP_WMAvg_IE(spIndicesForWMFitting) ...
    ,[0.5 0.5],r1_SP_WMStd_IE(spIndicesForWMFitting));
r1Data.factorR1_SPWM_IEModel = factorR1_SP_WM_IE;
r1Data.exponentR1_SPWM_IEModel = exponentR1_SP_WM_IE;
r1Data.covarianceR1_SPWM_IEModel = covarianceR1_SP_WM_IE;

% ==== gray matter ====
fprintf("<strong>GRAY MATTER </strong>\n");
exchange_waterToSolid_GM = exchangeRate ...
    /(gmWaterFraction);
exchange_lipidToWater_GM = exchangeRate ...
    /(gmNonWaterFraction * r1Data.constitutionData.gmLipidContent);
exchange_proteinToWater_GM = exchangeRate ...
    /(gmNonWaterFraction * r1Data.constitutionData.gmProteinContent);
r1Data.exchange_waterToSolid_GM = exchange_waterToSolid_GM;
r1Data.exchange_lipidToWater_GM = exchange_lipidToWater_GM;
r1Data.exchange_proteinToWater_GM = exchange_proteinToWater_GM;

[bestGuessTissueR1_GM_IE,bestGuessR1_SP_GM_IE] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    gmFieldStrengthArray,r1GMCollection,r1Data,gmSurfaceWaterFraction ...
    ,gmWaterFraction,gmLipidFraction,gmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_GM,exchange_lipidToWater_GM ...
    ,exchange_proteinToWater_GM,timeStepCountForGM,deltaTsForGM);
r1Data.r1_SPGM_IEModel = bestGuessR1_SP_GM_IE;
r1Data.bestGuessTissueR1_GM_IEModel = bestGuessTissueR1_GM_IE;

[r1_SP_GMAvg_IE,r1_SP_GMStd_IE,~] = getAveragSTDandCountsOfR1Collection( ...
    bestGuessR1_SP_GM_IE);
r1Data.r1_SP_GMAvg_IEModel = r1_SP_GMAvg_IE;
r1Data.r1_SP_GMStd_IEModel = r1_SP_GMStd_IE;

[factorR1_SP_GM_IE,exponentR1_SP_GM_IE,covarianceR1_SP_GM_IE] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForGMFitting,r1_SP_GMAvg_IE(spIndicesForGMFitting) ...
    ,[0.5 0.5],r1_SP_GMStd_IE(spIndicesForGMFitting));
r1Data.factorR1_SPGM_IEModel = factorR1_SP_GM_IE;
r1Data.exponentR1_SPGM_IEModel = exponentR1_SP_GM_IE;
r1Data.covarianceR1_SPGM_IEModel = covarianceR1_SP_GM_IE;

%% plotting R1_SP dependent on WM and GM tissue
set(0,'CurrentFigure', mdSimFigure); 
ieModelResultsSubPlt = initializeSubplot(mdSimFigure,2,1,2);
cla(ieModelResultsSubPlt);

ieModelResultsPlt(1) = plot(fieldStrengthAxis,factorR1_SP_WM_IE ...
    *fieldStrengthAxis.^(-exponentR1_SP_WM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE,exponentR1_SP_WM_IE ...
    ,covarianceR1_SP_WM_IE,fieldStrengthAxis,ieModelResultsPlt(1).Color);
ieModelResultsPlt(2) = errorbar(wmFieldStrengthArray,r1_SP_WMAvg_IE ...
    ,r1_SP_WMStd_IE,'LineStyle','none','LineWidth',1.5);
ieModelResultsPlt(3) = errorbar(wmFieldStrengthArray( ...
    spIndicesForWMFitting),r1_SP_WMAvg_IE(spIndicesForWMFitting) ...
    ,r1_SP_WMStd_IE(spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5);

ieModelResultsPlt(4) = plot(fieldStrengthAxis,factorR1_SP_GM_IE ...
    *fieldStrengthAxis.^(-exponentR1_SP_GM_IE),'r');
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE,exponentR1_SP_GM_IE ...
    ,covarianceR1_SP_GM_IE,fieldStrengthAxis,ieModelResultsPlt(2).Color);
ieModelResultsPlt(5) = errorbar(gmFieldStrengthArray,r1_SP_GMAvg_IE ...
    ,r1_SP_GMStd_IE,'LineStyle','none','LineWidth',1.5);
ieModelResultsPlt(6) = errorbar(gmFieldStrengthArray( ...
    spIndicesForGMFitting),r1_SP_GMAvg_IE(spIndicesForGMFitting) ...
    ,r1_SP_GMStd_IE(spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Color','g');

axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
xlabel("Field strength [T]");
ylabel("R$_{1,SP}$");
legend(ieModelResultsPlt(1:6),"Fit for WM","For all B$_0$ in WM" ...
    ,"B$_0$ for fitting in WM","Fit for GM" ...
    ,"For all B$_0$ in GM","B$_0$ for fitting in GM");
drawnow;

%% Averaging

% First average R1_SP and then perform the fit
onlyWMDataPoints = zeros(1,length(wmFieldStrengthArray));
r1_SP_AvgTissue_IE = zeros(1,length(wmFieldStrengthArray));
r1_SP_StatisticalError = zeros(1,length(wmFieldStrengthArray));
for fieldStrengthNr = 1:length(wmFieldStrengthArray)
    gmFieldStrengthIndex = find((wmFieldStrengthArray(fieldStrengthNr) ...
        - 0.0001 <= gmFieldStrengthArray) .* (wmFieldStrengthArray( ...
        fieldStrengthNr) + 0.0001 >= gmFieldStrengthArray));
    if isempty(gmFieldStrengthIndex)
        r1_SP_AvgTissue_IE(fieldStrengthNr) = r1_SP_GMAvg_IE( ...
            fieldStrengthNr);
        r1_SP_StatisticalError(fieldStrengthNr) = ...
            r1_SP_WMStd_IE(fieldStrengthNr);
        onlyWMDataPoints(fieldStrengthNr) = wmFieldStrengthArray( ...
            fieldStrengthNr);
    else
        r1_SP_AvgTissue_IE(fieldStrengthNr) = (r1_SP_GMAvg_IE( ...
            gmFieldStrengthIndex) + r1_SP_WMAvg_IE(fieldStrengthNr))/2;
        r1_SP_StatisticalError(fieldStrengthNr) = sqrt(1/( ...
            1/(r1_SP_WMStd_IE(fieldStrengthNr)^2) ...
            + 1/(r1_SP_GMStd_IE(gmFieldStrengthIndex)^2)));
    end
end

r1Data.r1_SP_AvgTissue_IEModel = r1_SP_AvgTissue_IE;
r1Data.r1_SP_StatisticalError = r1_SP_StatisticalError;

if wmFieldStrengthArray(~logical(onlyWMDataPoints)) ~= gmFieldStrengthArray
    error("There is GM data point that isn't found in WM. " ...
        + "This case is not implemented yet");
end

if spFieldStrengthForWMFitting == spFieldStrengthForGMFitting
    [factorR1_SPAvg_IE,exponentR1_SPAvg_IE,covarianceR1_SPAvg_IE] ...
        = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
        spFieldStrengthForWMFitting,r1_SP_AvgTissue_IE( ...
        spIndicesForWMFitting),[0.5 0.5],r1_SP_StatisticalError( ...
        spIndicesForWMFitting).^2);
    
    r1Data.factorR1_SPAvg_IEModel = factorR1_SPAvg_IE;
    r1Data.exponentR1_SPAvg_IEModel = exponentR1_SPAvg_IE;
    r1Data.covarianceR1_SPAvg_IEModel = covarianceR1_SPAvg_IE;
else
    error("Case not implemented yet.");
end

averageSpPowerFunction = (factorR1_SP_WM_IE*fieldStrengthAxis.^( ...
    -exponentR1_SP_WM_IE) + factorR1_SP_GM_IE*fieldStrengthAxis.^( ...
    -exponentR1_SP_GM_IE))/2;

r1_SP_systematicTissueBasedError = ones(1,length(r1_SP_AvgTissue_IE));
r1_SP_systematicTissueBasedError(logical(onlyWMDataPoints)) = abs( ...
    r1_SP_WMAvg_IE(logical(onlyWMDataPoints)) ...
    - r1_SP_AvgTissue_IE(logical(onlyWMDataPoints)));
r1Data.r1_SP_systematicTissueBasedError = r1_SP_systematicTissueBasedError;

r1_SP_systematicTissueBasedError(~logical(onlyWMDataPoints)) = (abs( ...
    r1_SP_WMAvg_IE(~logical(onlyWMDataPoints)) ...
    - r1_SP_AvgTissue_IE(~logical(onlyWMDataPoints))) ...
    + abs(r1_SP_GMAvg_IE - r1_SP_AvgTissue_IE(~logical( ...
    onlyWMDataPoints))))/2;
r1Data.r1_SP_systematicTissueBasedError = r1_SP_systematicTissueBasedError;

%% plotting comparison of fit avg. and avg. fit
initializeFigure('axisFontSize',14,'legendFontSize',12);
title("Comparison of avg. power functions of SP");
plot(fieldStrengthAxis,averageSpPowerFunction);
plot(fieldStrengthAxis,factorR1_SPAvg_IE*fieldStrengthAxis.^( ...
    -exponentR1_SPAvg_IE));
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
xlabel("Field strength [T]");
ylabel("R$_{1,SP}^{Avg.}$ [Hz]");
legend("Avg. power function","Fit to avg. data points $>$ 0.35T");
drawnow;
if saving
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"AvgPowerFuncVsFitToAvgDataPointsForR1" ...
        ,"SP");
end

%% plotting average results
set(0, 'CurrentFigure', mdSimFigure); 
mdSimFigure.CurrentAxes = simResultsSubplot;

mdSimPlt(5) = plot(fieldStrengthAxis,factorR1_SPAvg_IE ...
    *fieldStrengthAxis.^(-exponentR1_SPAvg_IE));
drawSTDRegionAroundPowerFunction(factorR1_SPAvg_IE,exponentR1_SPAvg_IE ...
    ,covarianceR1_SPAvg_IE,fieldStrengthAxis,mdSimPlt(5).Color);
legendEntriesSimResultsSubplot{end+1} = "R$_{1,SP}^{Avg} \pm$ STD fit";
legend(mdSimPlt(1:5),legendEntriesSimResultsSubplot);
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
drawnow;

%% plotting STD at WM and GM to show overlap
additionalFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
r1SPStdSubplot = initializeSubplot(additionalFigure,2,2,1);
cla(r1SPStdSubplot);
r1SPStdLegendEntries = {};
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE,exponentR1_SP_WM_IE ...
    ,covarianceR1_SP_WM_IE,fieldStrengthAxis,mdSimPlt(1).Color,0.5);
r1SPStdLegendEntries{end+1} = "STD WM";

drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE,exponentR1_SP_GM_IE ...
    ,covarianceR1_SP_GM_IE,fieldStrengthAxis,mdSimPlt(2).Color,0.5);
r1SPStdLegendEntries{end+1} = "STD GM";
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);

legend(r1SPStdLegendEntries);
xlabel("Field Strength [T]");
ylabel("R$_{1,SP}$ [Hz]");
drawnow;

%% plotting ratio of R1_SP and R1_SL

set(0,'CurrentFigure',additionalFigure);
r1SPSLRatioSubPlot = initializeSubplot(additionalFigure,2,2,3);
cla(r1SPSLRatioSubPlot);
lgd = legend();
lgd.Visible = 'off';

plot(fieldStrengthAxis,(factorR1_SL*fieldStrengthAxis.^(-exponentR1_SL)) ...
    ./(factorR1_SPAvg_IE*fieldStrengthAxis.^(-exponentR1_SPAvg_IE)));

xlabel("Field strength [T]");
ylabel("${}^{R_{1,SL}}{\mskip -5mu/\mskip -3mu}_{R_{1,SP}}$ [a.u.]");
drawnow;


%% determine R1_SP with larger exchange rates
fprintf("<strong>WHITE MATTER WITH INCREASED EXCHANGE RATES</strong>\n");
exchangeRate_incr = exchangeRate * exchangeRateIncreaseFactor;
[exchange_waterToSolid_WMincr,exchange_lipidToWater_WMincr ...
    ,exchange_proteinToWater_WMincr] ...
    = determinePoolSizeDependentExchangeRates(exchangeRate_incr ...
    ,wmWaterFraction,wmLipidFraction,wmProteinFraction);
r1Data.exchange_waterToSolid_WMincr = exchange_waterToSolid_WMincr;
r1Data.exchange_lipidToWater_WMincr = exchange_lipidToWater_WMincr;
r1Data.exchange_proteinToWater_WMincr = exchange_proteinToWater_WMincr;

[bestGuessTissueR1_WM_IE_incr,bestGuessR1_SP_WM_IE_incr] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    wmFieldStrengthArray,r1WMCollection,r1Data,wmSurfaceWaterFraction ...
    ,wmWaterFraction,wmLipidFraction,wmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_WMincr,exchange_lipidToWater_WMincr ...
    ,exchange_proteinToWater_WMincr,timeStepCountForWM,deltaTsForWM);
r1Data.r1_SPWM_IEModel_incr = bestGuessR1_SP_WM_IE_incr;
r1Data.bestGuessTissueR1_WM_IEModel_incr = bestGuessTissueR1_WM_IE_incr;

[r1_SP_WMAvg_IE_incr,r1_SP_WMStd_IE_incr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_WM_IE_incr);
r1Data.r1_SP_WMAvg_IEModel_incr = r1_SP_WMAvg_IE_incr;
r1Data.r1_SP_WMStd_IEModel_incr = r1_SP_WMStd_IE_incr;

[factorR1_SP_WM_IE_incr,exponentR1_SP_WM_IE_incr ...
    ,covarianceR1_SP_WM_IE_incr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForWMFitting,r1_SP_WMAvg_IE_incr( ...
    spIndicesForWMFitting),[0.5 0.5],r1_SP_WMStd_IE_incr( ...
    spIndicesForWMFitting));
r1Data.factorR1_SPWM_IEModel_incr = factorR1_SP_WM_IE_incr;
r1Data.exponentR1_SPWM_IEModel_incr = exponentR1_SP_WM_IE_incr;
r1Data.covarianceR1_SPWM_IEModel_incr = covarianceR1_SP_WM_IE_incr;

fprintf("<strong>WHITE MATTER WITH DECREASED EXCHANGE RATES</strong>\n");
exchangeRate_decr = exchangeRate * exchangeRateDecreaseFactor;
[exchange_waterToSolid_WMdecr,exchange_lipidToWater_WMdecr ...
    ,exchange_proteinToWater_WMdecr] ...
    = determinePoolSizeDependentExchangeRates(exchangeRate_decr ...
    ,wmWaterFraction,wmLipidFraction,wmProteinFraction);
r1Data.exchange_waterToSolid_WMdecr = exchange_waterToSolid_WMdecr;
r1Data.exchange_lipidToWater_WMdecr = exchange_lipidToWater_WMdecr;
r1Data.exchange_proteinToWater_WMdecr = exchange_proteinToWater_WMdecr;

[bestGuessTissueR1_WM_IE_decr,bestGuessR1_SP_WM_IE_decr] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    wmFieldStrengthArray,r1WMCollection,r1Data,wmSurfaceWaterFraction ...
    ,wmWaterFraction,wmLipidFraction,wmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_WMdecr,exchange_lipidToWater_WMdecr ...
    ,exchange_proteinToWater_WMdecr,timeStepCountForWM,deltaTsForWM);
r1Data.r1_SPWM_IEModel_decr = bestGuessR1_SP_WM_IE_decr;
r1Data.bestGuessTissueR1_WM_IEModel_decr = bestGuessTissueR1_WM_IE_decr;

[r1_SP_WMAvg_IE_decr,r1_SP_WMStd_IE_decr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_WM_IE_decr);
r1Data.r1_SP_WMAvg_IEModel_decr = r1_SP_WMAvg_IE_decr;
r1Data.r1_SP_WMStd_IEModel_decr = r1_SP_WMStd_IE_decr;

[factorR1_SP_WM_IE_decr,exponentR1_SP_WM_IE_decr ...
    ,covarianceR1_SP_WM_IE_decr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForWMFitting,r1_SP_WMAvg_IE_decr( ...
    spIndicesForWMFitting),[0.5 0.5],r1_SP_WMStd_IE_decr( ...
    spIndicesForWMFitting));
r1Data.factorR1_SPWM_IEModel_decr = factorR1_SP_WM_IE_decr;
r1Data.exponentR1_SPWM_IEModel_decr = exponentR1_SP_WM_IE_decr;
r1Data.covarianceR1_SPWM_IEModel_decr = covarianceR1_SP_WM_IE_decr;

fprintf("<strong>GRAY MATTER WITH INCREASED EXCHANGE RATES</strong>\n");
[exchange_waterToSolid_GMincr,exchange_lipidToWater_GMincr ...
    ,exchange_proteinToWater_GMincr] ...
    = determinePoolSizeDependentExchangeRates(exchangeRate_incr ...
    ,gmWaterFraction,gmLipidFraction,gmProteinFraction);
r1Data.exchange_waterToSolid_GMincr = exchange_waterToSolid_GMincr;
r1Data.exchange_lipidToWater_GMincr = exchange_lipidToWater_GMincr;
r1Data.exchange_proteinToWater_GMincr = exchange_proteinToWater_GMincr;

[bestGuessTissueR1_GM_IE_incr,bestGuessR1_SP_GM_IE_incr] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    gmFieldStrengthArray,r1GMCollection,r1Data,gmSurfaceWaterFraction ...
    ,gmWaterFraction,gmLipidFraction,gmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_GMincr,exchange_lipidToWater_GMincr ...
    ,exchange_proteinToWater_GMincr,timeStepCountForGM,deltaTsForGM);
r1Data.r1_SPGM_IEModel_incr = bestGuessR1_SP_GM_IE_incr;
r1Data.bestGuessTissueR1_GM_IEModel_incr = bestGuessTissueR1_GM_IE_incr;

[r1_SP_GMAvg_IE_incr,r1_SP_GMStd_IE_incr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_GM_IE_incr);
r1Data.r1_SP_GMAvg_IEModel_incr = r1_SP_GMAvg_IE_incr;
r1Data.r1_SP_GMStd_IEModel_incr = r1_SP_GMStd_IE_incr;

[factorR1_SP_GM_IE_incr,exponentR1_SP_GM_IE_incr ...
    ,covarianceR1_SP_GM_IE_incr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForGMFitting,r1_SP_GMAvg_IE_incr( ...
    spIndicesForGMFitting),[0.5 0.5],r1_SP_GMStd_IE_incr( ...
    spIndicesForGMFitting));
r1Data.factorR1_SPGM_IEModel_incr = factorR1_SP_GM_IE_incr;
r1Data.exponentR1_SPGM_IEModel_incr = exponentR1_SP_GM_IE_incr;
r1Data.covarianceR1_SPGM_IEModel_incr = covarianceR1_SP_GM_IE_incr;

fprintf("<strong>GRAY MATTER WITH DECREASED EXCHANGE RATES</strong>\n");
[exchange_waterToSolid_GMdecr,exchange_lipidToWater_GMdecr ...
    ,exchange_proteinToWater_GMdecr] ...
    = determinePoolSizeDependentExchangeRates(exchangeRate_decr ...
    ,gmWaterFraction,gmLipidFraction,gmProteinFraction);
r1Data.exchange_waterToSolid_GMdecr = exchange_waterToSolid_GMdecr;
r1Data.exchange_lipidToWater_GMdecr = exchange_lipidToWater_GMdecr;
r1Data.exchange_proteinToWater_GMdecr = exchange_proteinToWater_GMdecr;

[bestGuessTissueR1_GM_IE_decr,bestGuessR1_SP_GM_IE_decr] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel( ...
    gmFieldStrengthArray,r1GMCollection,r1Data,gmSurfaceWaterFraction ...
    ,gmWaterFraction,gmLipidFraction,gmProteinFraction ...
    ,lowerBoundR1_SPFraction,upperBoundR1_SPFraction,r1_SPSamplingPoints ...
    ,exchange_waterToSolid_GMdecr,exchange_lipidToWater_GMdecr ...
    ,exchange_proteinToWater_GMdecr,timeStepCountForGM,deltaTsForGM);
r1Data.r1_SPGM_IEModel_decr = bestGuessR1_SP_GM_IE_decr;
r1Data.bestGuessTissueR1_GM_IEModel_decr = bestGuessTissueR1_GM_IE_decr;

[r1_SP_GMAvg_IE_decr,r1_SP_GMStd_IE_decr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_GM_IE_decr);
r1Data.r1_SP_GMAvg_IEModel_decr = r1_SP_GMAvg_IE_decr;
r1Data.r1_SP_GMStd_IEModel_decr = r1_SP_GMStd_IE_decr;

[factorR1_SP_GM_IE_decr,exponentR1_SP_GM_IE_decr ...
    ,covarianceR1_SP_GM_IE_decr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForGMFitting,r1_SP_GMAvg_IE_decr( ...
    spIndicesForGMFitting),[0.5 0.5],r1_SP_GMStd_IE_decr( ...
    spIndicesForGMFitting));
r1Data.factorR1_SPGM_IEModel_decr = factorR1_SP_GM_IE_decr;
r1Data.exponentR1_SPGM_IEModel_decr = exponentR1_SP_GM_IE_decr;
r1Data.covarianceR1_SPGM_IEModel_decr = covarianceR1_SP_GM_IE_decr;

%% Find systematic exchange rate based error

r1_SP_systematicPositiveExchangeRateBasedError = abs( ...
    r1_SP_AvgTissue_IE - r1_SP_WMAvg_IE_decr);
r1Data.r1_SP_systematicPositiveExchangeRateBasedError ...
    = r1_SP_systematicPositiveExchangeRateBasedError;

lowerR1SPValues = ones(1,length(r1_SP_AvgTissue_IE));
lowerR1SPValues(~logical(onlyWMDataPoints)) = r1_SP_GMAvg_IE_incr;
lowerR1SPValues(logical(onlyWMDataPoints)) = r1_SP_WMAvg_IE_incr( ...
    logical(onlyWMDataPoints));
r1_SP_systematicNegativeExchangeRateBasedError ...
    = ones(1,length(r1_SP_AvgTissue_IE));
r1_SP_systematicNegativeExchangeRateBasedError( ...
    logical(onlyWMDataPoints)) = abs( ...
    r1_SP_AvgTissue_IE(logical(onlyWMDataPoints)) ...
    - r1_SP_AvgTissue_IE(logical(onlyWMDataPoints)));
r1_SP_systematicNegativeExchangeRateBasedError( ...
    ~logical(onlyWMDataPoints)) = (abs( ...
    r1_SP_AvgTissue_IE(~logical(onlyWMDataPoints)) ...
    - r1_SP_AvgTissue_IE(~logical(onlyWMDataPoints))) ...
    + abs(r1_SP_AvgTissue_IE(~logical(onlyWMDataPoints)) ...
    - r1_SP_GMAvg_IE_incr))/2;
r1Data.r1_SP_systematicNegativeExchangeRateBasedError ... 
    = r1_SP_systematicNegativeExchangeRateBasedError;


fprintf("<strong>Error Analyis</strong>\n");
for index = find(spIndicesForWMFitting)
   fprintf("%.4fT: %.4f +/- %.4f + %.4f - %.4f Hz \n" ...
       ,wmFieldStrengthArray(index) ...
       ,r1_SP_AvgTissue_IE(index) ...
       ,r1_SP_StatisticalError(index) ...
       ,r1_SP_systematicTissueBasedError(index) ...
       + r1_SP_systematicPositiveExchangeRateBasedError(index) ...
       ,r1_SP_systematicTissueBasedError(index) ...
       + r1_SP_systematicNegativeExchangeRateBasedError(index));
   
    
end

%% Plotting results of different exchange rates
set(0,'CurrentFigure',additionalFigure);
r1WmChangedExchangeRaterSubPlt = initializeSubplot(additionalFigure,2,2,2);
cla(r1WmChangedExchangeRaterSubPlt);
ylabel("R$_{1,SP}$ based on WM");
xlabel("Field strength [T]");
r1WmChangedExchangeRaterSubPlt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE*fieldStrengthAxis.^(-exponentR1_SP_WM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE,exponentR1_SP_WM_IE ...
    ,covarianceR1_SP_WM_IE,fieldStrengthAxis ...
    ,r1WmChangedExchangeRaterSubPlt(1).Color);

r1WmChangedExchangeRaterSubPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE_incr*fieldStrengthAxis.^(-exponentR1_SP_WM_IE_incr));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE_incr ...
    ,exponentR1_SP_WM_IE_incr,covarianceR1_SP_WM_IE_incr ...
    ,fieldStrengthAxis,r1WmChangedExchangeRaterSubPlt(2).Color);

r1WmChangedExchangeRaterSubPlt(3) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE_decr*fieldStrengthAxis.^(-exponentR1_SP_WM_IE_decr));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE_decr ...
    ,exponentR1_SP_WM_IE_decr,covarianceR1_SP_WM_IE_decr ...
    ,fieldStrengthAxis,r1WmChangedExchangeRaterSubPlt(3).Color);
axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);

legend(r1WmChangedExchangeRaterSubPlt(1:3) ...
    ,"unchanged","increased","decreased");
drawnow;

r1GmChangedExchangeRaterSubPlt = initializeSubplot(additionalFigure,2,2,4);
cla(r1GmChangedExchangeRaterSubPlt);
ylabel("R$_{1,SP}$ based on GM");
xlabel("Field strength [T]");
r1GmChangedExchangeRaterSubPlt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE*fieldStrengthAxis.^(-exponentR1_SP_GM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE,exponentR1_SP_GM_IE ...
    ,covarianceR1_SP_GM_IE,fieldStrengthAxis ...
    ,r1GmChangedExchangeRaterSubPlt(1).Color);

r1GmChangedExchangeRaterSubPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE_incr*fieldStrengthAxis.^(-exponentR1_SP_GM_IE_incr));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE_incr ...
    ,exponentR1_SP_GM_IE_incr,covarianceR1_SP_GM_IE_incr ...
    ,fieldStrengthAxis,r1GmChangedExchangeRaterSubPlt(2).Color);

r1GmChangedExchangeRaterSubPlt(3) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE_decr*fieldStrengthAxis.^(-exponentR1_SP_GM_IE_decr));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE_decr ...
    ,exponentR1_SP_GM_IE_decr,covarianceR1_SP_GM_IE_decr ...
    ,fieldStrengthAxis,r1GmChangedExchangeRaterSubPlt(3).Color);

axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
lgd = legend();
lgd.Visible = 'off';
drawnow;
%% determine R1_SP with FE model
fprintf("<strong> === FAST EXCHANGE MODEL ===</strong>\n");

% ==== White matter ====
bestGuessR1_SP_WM_FE = determineR1_SPWithFastExchangeMode( ...
    wmFieldStrengthArray,r1WMCollection,r1Data,wmProteinFraction ...
    ,wmLipidFraction,wmProteinWaterFraction,wmLipidWaterFraction ...
    ,wmFreeWaterFraction,lowFieldLimit);
r1Data.r1_SPWM_FEModel = bestGuessR1_SP_WM_FE;

[r1_SP_WMAvg_FE,r1_SP_WMStd_FE,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_WM_FE);
r1Data.r1_SP_WMAvg_FEModel = r1_SP_WMAvg_FE;
r1Data.r1_SP_WMStd_FEModel = r1_SP_WMStd_FE;

[factorR1_SP_WM_FE,exponentR1_SP_WM_FE,covarianceR1_SP_WM_FE] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForWMFitting,r1_SP_WMAvg_FE(spIndicesForWMFitting) ...
    ,[0.5 0.5],r1_SP_WMStd_FE(spIndicesForWMFitting));
r1Data.factorR1_SPWM_FEModel = factorR1_SP_WM_FE;
r1Data.exponentR1_SPWM_FEModel = exponentR1_SP_WM_FE;
r1Data.covarianceR1_SPWM_FEModel = covarianceR1_SP_WM_FE;

% ==== Gray matter ====
bestGuessR1_SP_GM_FE = determineR1_SPWithFastExchangeMode( ...
    gmFieldStrengthArray,r1GMCollection,r1Data,gmProteinFraction ...
    ,gmLipidFraction,gmProteinWaterFraction,gmLipidWaterFraction ...
    ,gmFreeWaterFraction,lowFieldLimit);
r1Data.r1_SPGM_FEModel = bestGuessR1_SP_GM_FE;

[r1_SP_GMAvg_FE,r1_SP_GMStd_FE,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_GM_FE);
r1Data.r1_SP_GMAvg_FEModel = r1_SP_GMAvg_FE;
r1Data.r1_SP_GMStd_FEModel = r1_SP_GMStd_FE;

[factorR1_SP_GM_FE,exponentR1_SP_GM_FE,covarianceR1_SP_GM_FE] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    spFieldStrengthForGMFitting,r1_SP_GMAvg_FE(spIndicesForGMFitting) ...
    ,[0.5 0.5],r1_SP_GMStd_FE(spIndicesForGMFitting));
r1Data.factorR1_SPGM_FEModel = factorR1_SP_GM_FE;
r1Data.exponentR1_SPGM_FEModel = exponentR1_SP_GM_FE;
r1Data.covarianceR1_SPGM_FEModel = covarianceR1_SP_GM_FE;

%% plotting FE model results
FEModelResultsPlt = initializeFigure( ...
    'axisFontSize',axisFontSize,'legendFontSize',legendFontSize);
title("FE model results");
ylabel("R$_{1,SP}$ [Hz]");
xlabel("Field strength [T]");
feModelPlt(1) = errorbar(wmFieldStrengthArray,r1_SP_WMAvg_FE ...
    ,r1_SP_WMStd_FE,'LineStyle','none','LineWidth',1.5);
feModelPlt(2) = plot(fieldStrengthAxis,factorR1_SP_WM_FE* ...
    fieldStrengthAxis.^(- exponentR1_SP_WM_FE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_FE,exponentR1_SP_WM_FE...
    ,covarianceR1_SP_WM_FE,fieldStrengthAxis,feModelPlt(2).Color);

feModelPlt(3) = errorbar(gmFieldStrengthArray,r1_SP_GMAvg_FE ...
    ,r1_SP_GMStd_FE,'LineStyle','none','LineWidth',1.5,'Color','g');
feModelPlt(4) = plot(fieldStrengthAxis,factorR1_SP_GM_FE ...
    *fieldStrengthAxis.^(- exponentR1_SP_GM_FE));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_FE,exponentR1_SP_GM_FE...
    ,covarianceR1_SP_GM_FE,fieldStrengthAxis,feModelPlt(4).Color);

axis([0 fieldStrengthAxis(end) 0 highestR1InPlot])
lgd = legend();
lgd.Visible = 'off';

legend(feModelPlt(1:4),"FE WM $\pm$ STD","FE WM fit $\pm$ STD" ...
    ,"FE GM $\pm$ STD","FE GM fit $\pm$ STD");
drawnow;

%% Writing results to console
fprintf("\n \n \n");
fprintf("<strong>***** RESULTS *****</strong> \n");

fprintf("1. Literature \n");
fprintf("  R1_WM(B_0) = %.4f * B_0 ^(-%.4f) \n",wmFactor,wmExponent);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",wmCovariance);

fprintf("  R1_GM(B_0) = %.4f * B_0 ^(-%.4f) \n",gmFactor,gmExponent);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",gmCovariance);
fprintf("\n\n");

fprintf("2. MD sim. \n");
fprintf("  R1_SL(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SL,exponentR1_SL);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SL);

fprintf("  R1_LW(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_LW,exponentR1_LW);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_LW);
fprintf("\n\n");

fprintf("3. Model predictions \n");
fprintf(" Based on WM: \n");
fprintf("  R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_WM_IE ...
    ,exponentR1_SP_WM_IE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_WM_IE);
fprintf("\n");

fprintf(" Based on GM: \n");
fprintf("  IE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_GM_IE ...
    ,exponentR1_SP_GM_IE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_GM_IE);
fprintf("\n");

fprintf(" On average \n");
fprintf("  IE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SPAvg_IE ...
    ,exponentR1_SPAvg_IE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_GM_IE);
fprintf("\n");


%% saving

if saving
    save(r1DataFolderName + "SM_MW_SP_histCompartmentAndCrossR1" ...
        ,'-struct','r1Data');
    
    set(0, 'CurrentFigure', literatureFigure);
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"ObservedR1ValuesFit","WMandGM");
    
    set(0, 'CurrentFigure', FEModelResultsPlt); 
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"FEModelResultsForR1","SP");
    
    set(0, 'CurrentFigure', mdSimFigure); 
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"SL_LW_SP" ...
        ,"PredictedValuesAndFieldStrengthBehvior");
    
    set(0,'CurrentFigure',additionalFigure);
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"Overlap_Ration_IncrAndDecrExchRats" ...
        ,"ofR1_SP");
end

%% functions

function [exchange_waterToSolid,exchange_lipidToWater...
    ,exchange_proteinToWater] ...
    = determinePoolSizeDependentExchangeRates( ...
    exchangeRateNormalizedToOverallPoolSize,waterFraction ...
    ,lipidFraction,proteinFraction)

exchange_waterToSolid = exchangeRateNormalizedToOverallPoolSize ...
    /waterFraction;
exchange_lipidToWater = exchangeRateNormalizedToOverallPoolSize ...
    /lipidFraction;
exchange_proteinToWater = exchangeRateNormalizedToOverallPoolSize ...
    /proteinFraction;

end

function [r1_SPCollection] = determineR1_SPWithFastExchangeMode( ...
    fieldStrengthArray,r1Collection,r1Data,solidProteinContent ...
    ,solidLipidContent,proteinWaterContent,lipidWaterContent ...
    ,freeWaterContent,lowFieldLimit)

r1_SPCollection = cell(1,length(fieldStrengthArray));
notLipidContent = solidProteinContent + lipidWaterContent ...
    + proteinWaterContent + freeWaterContent;

for fieldStrengthNr = 1:length(fieldStrengthArray)
    fieldStrength = fieldStrengthArray(fieldStrengthNr);
    fprintf(" == %.2f Tesla ==  \n",fieldStrength);
    
    fieldStrengthIndex = ...
        r1Data.fieldStrengths > fieldStrength - 0.0001 ...
        & r1Data.fieldStrengths < fieldStrength + 0.0001;
    if isempty(fieldStrengthIndex)
        error("No field strength found-");
    end
    
    r1_SL = r1Data.r1Eff_SM(fieldStrengthIndex);
    r1_LW = r1Data.r1Hist_MW(fieldStrengthIndex);
    
    for dataPointNr = 1:length(r1Collection{fieldStrengthNr})
        if fieldStrength < lowFieldLimit
            
            r1_SPCollection{fieldStrengthNr}(dataPointNr) ...
                = notLipidContent/solidProteinContent ...
                *(r1Collection{fieldStrengthNr}(dataPointNr) ...
                - lipidWaterContent/notLipidContent*r1_LW ...
                - proteinWaterContent/notLipidContent*r1_LW ...
                - freeWaterContent/notLipidContent*r1Data.r1_free);
        else
            r1_SPCollection{fieldStrengthNr}(dataPointNr) ...
                = 1/solidProteinContent ...
                *(r1Collection{fieldStrengthNr}(dataPointNr) ...
                - solidLipidContent*r1_SL...
                - lipidWaterContent*r1_LW ...
                - proteinWaterContent*r1_LW ...
                - freeWaterContent*r1Data.r1_free);
        end
    end
    
    fprintf("  avg. R1_SP = %.4f\n",mean( ...
        r1_SPCollection{fieldStrengthNr}));
   
end

end

function [bestGuessTissueR1,bestGuessR1_SP] ...
    = findBestGuessR1sForDifferentDataPointsWithIEModel(fieldStrengthArray ...
    ,r1Collection,r1Data,surfaceWaterFraction,waterFraction,lipidFraction ...
    ,proteinFraction,lowerBoundR1_SPFraction,upperBoundR1_SPFraction ...
    ,r1_SPSamplingPoints,exchange_waterToSolid,exchange_lipidToWater ...
    ,exchange_proteinToWater,timeStepCounts,deltaTArray)

bestGuessTissueR1 = cell(1,length(fieldStrengthArray));
bestGuessR1_SP = cell(1,length(fieldStrengthArray));

for fieldStrengthNr = 1:length(fieldStrengthArray)
    fieldStrength = fieldStrengthArray(fieldStrengthNr);
    fprintf(" == %.2f Tesla ==  \n",fieldStrength);
    deltaT = deltaTArray(fieldStrengthNr);
    fprintf("     dT = %.4d \n",deltaT);
    timeStepCount = timeStepCounts(fieldStrengthNr);
    fprintf("     time step count = %.2d \n",timeStepCount);
    
    fieldStrengthIndex = ...
        r1Data.fieldStrengths > fieldStrength - 0.0001 ...
        & r1Data.fieldStrengths < fieldStrength + 0.0001;
    if isempty(fieldStrengthIndex)
        error("No field strength found");
    end
    
    r1_SL = r1Data.r1Eff_SM(fieldStrengthIndex);
    fprintf(" R1_SL = %.4f \n",r1_SL);
    r1_LW = r1Data.r1Surf_MW(fieldStrengthIndex);
    r1_water = surfaceWaterFraction/waterFraction * r1_LW ...
        + (1-surfaceWaterFraction)/waterFraction*r1Data.r1_free;
    
    r1_SPArray = linspace(lowerBoundR1_SPFraction*r1_SL ...
        ,upperBoundR1_SPFraction*r1_SL,r1_SPSamplingPoints);
    
    for dataPointNr = 1:length(r1Collection{fieldStrengthNr})
        r1 = r1Collection{fieldStrengthNr}(dataPointNr);
        waterR1_wm = findWaterR1ForGivenR1_SP(r1_SPArray,r1_SL ...
            ,r1_water,waterFraction,lipidFraction,proteinFraction ...
            ,exchange_waterToSolid,exchange_lipidToWater ...
            ,exchange_proteinToWater,timeStepCount,deltaT);
        
        [~,bestGuessIndex] = min(abs(waterR1_wm - r1));
        bestGuessTissueR1{fieldStrengthNr}(dataPointNr) ...
            = waterR1_wm(bestGuessIndex);
        bestGuessR1_SP{fieldStrengthNr}(dataPointNr) ...
            = r1_SPArray(bestGuessIndex);
    end
    
    fprintf("  avg. Obs./Calc. R1 = %.4f / %.4f\n" ...
        ,mean(r1Collection{fieldStrengthNr}) ...
        ,mean(bestGuessTissueR1{fieldStrengthNr}));
    fprintf("  avg. R1_SP = %.4f\n",mean(bestGuessR1_SP{fieldStrengthNr}));
    
end
end


function [r1Avg,r1STD,r1Counts] ...
    = getAveragSTDandCountsOfR1Collection(r1Collection,varargin)

if ~isempty(varargin) && length(varargin) ~= 2
    error("Too less or too many input arguments");
elseif ~isempty(varargin) && length(varargin) == 2
    relaxationTimesTable = varargin{1};
    fieldStrengthArray = varargin{2};
end

variableName = inputname(1);
r1Avg = [];
r1STD = [];
r1Counts = [];

for dataPointNr = 1:length(r1Collection) 
    r1 = r1Collection{dataPointNr};
    r1Counts(end+1) = length(r1); %#ok<AGROW>
    if r1Counts(end) == 1 && isempty(varargin)
        r1Avg(end+1) = mean(r1); %#ok<AGROW>
        r1STD(end+1) = nan; %#ok<AGROW>
        continue;
    elseif r1Counts(end) == 1 && istable(varargin{1})
        fieldStrength = fieldStrengthArray(dataPointNr);
        indexInTable = ...
            relaxationTimesTable.fieldStrength > (fieldStrength - 0.001) ...
            & relaxationTimesTable.fieldStrength < (fieldStrength + 0.001);
        if variableName == "r1WMCollection"
            standardDev = relaxationTimesTable.T1_WMSTD(indexInTable);
        elseif variableName == "r1GMCollection"
            standardDev = relaxationTimesTable.T1_GMSTD(indexInTable);
        else
            error("Not known variablen name");
        end
        r1STD(end+1) = standardDev; %#ok<AGROW>
    elseif ~isempty(varargin) && ~istable(varargin{1})
        error("Not implemented!");
    else
        r1STD(end+1) = std(r1); %#ok<AGROW>
    end
    r1Avg(end+1) = mean(r1); %#ok<AGROW>
end

r1STD(isnan(r1STD)) = max(r1STD)*1.5;

end

function [factor,exponent,covariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction(fieldStrengths ...
    ,relaxationRates,startValues,variance,varargin)

powerFunction = @(parameters,fieldStrengths)(parameters(1) ...
        *fieldStrengths.^(-parameters(2)));

if isempty(varargin)

    [parameters,~,~,covariance] = nlinfit(fieldStrengths ...
        ,relaxationRates,powerFunction,startValues,'Weights'...
        ,1./variance);

else
    error("Not implemented yet.");
    
end

factor = parameters(1);
exponent = parameters(2);

end

function [factor,exponent,covariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengths,relaxationRates,startValues,varargin)

powerFunction = @(parameters,fieldStrengths)(parameters(1) ...
        *fieldStrengths.^(-parameters(2)));

if isempty(varargin)

    [parameters,~,~,covariance] = nlinfit(fieldStrengths ...
        ,relaxationRates,powerFunction,startValues);

else
    error("Not implemented yet.");
    
end

factor = parameters(1);
exponent = parameters(2);

end

function showDistributionOfValues(wmFieldStrengthArray,r1WMCollection ...
    ,gmFieldStrengthArray,r1GMCollection,axisFontSize,legendFontSize)

r1WMData = createMatrixOfDataPoints(r1WMCollection);
r1GMData = createMatrixOfDataPoints(r1GMCollection);

[r1AvgWMData,~,wmCounts] ...
    = getAveragSTDandCountsOfR1Collection(r1WMCollection);
[r1AvgGMData,~,gmCounts] ...
    = getAveragSTDandCountsOfR1Collection(r1GMCollection);

fig = initializeFigure('legend',false,'axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
initializeSubplot(fig,2,1,1);
labels = strsplit(sprintf("%.3f (#%i) \n" ...
    ,[wmFieldStrengthArray;wmCounts]),'\n');
boxplot(r1WMData','Labels',labels(1:end-1));
plot(r1AvgWMData,'g*');
xtickangle(45);
ylabel("R$_{1,WM}$ [Hz]");
lgd = legend();
lgd.Visible = 'off';

initializeSubplot(fig,2,1,2);
labels = strsplit(sprintf("%.3f (#%i) \n" ...
    ,[gmFieldStrengthArray;gmCounts]),'\n');
boxplot(r1GMData','Labels',labels(1:end-1));
plot(r1AvgGMData,'g*');
xlabel("Field strength [T]");
xtickangle(45);
ylabel("R$_{1,GM}$ [Hz]");
lgd = legend();
lgd.Visible = 'off';

end

function allDataMatrix = createMatrixOfDataPoints(r1Collection)
allDataMatrix = nan(length(r1Collection),200);

for fieldStrengthNr = 1:length(r1Collection)
   allDataMatrix(fieldStrengthNr,1:length(r1Collection{ ...
       fieldStrengthNr}))= r1Collection{fieldStrengthNr};
end

end


function r1CollectionWithoutOutliers = ignoreOutliers(r1Collection ...
    ,lowerQuantile,upperQuantile)

r1CollectionWithoutOutliers = cell(1,length(r1Collection));

for fieldStrengthNr = 1:length(r1Collection)
    Q1 = quantile(r1Collection{fieldStrengthNr},lowerQuantile);
    Q3 = quantile(r1Collection{fieldStrengthNr},upperQuantile);
    Spread = 1.5*(Q3-Q1);
    MaxValue = Q3 + Spread;
    MinValue = Q1 - Spread;
    
    r1CollectionWithoutOutliers{fieldStrengthNr} ...
        = r1Collection{fieldStrengthNr}( ...
        r1Collection{fieldStrengthNr} <= MaxValue ...
        & r1Collection{fieldStrengthNr} >= MinValue);
    
end

end

