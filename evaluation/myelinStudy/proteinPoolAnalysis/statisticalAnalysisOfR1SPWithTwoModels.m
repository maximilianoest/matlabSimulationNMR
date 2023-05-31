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
dataPointsCountLowerLimit = 1;
fieldStrengthsToIgnoreInWM = [0.5 0.7 1.0];
fieldStrengthsToIgnoreInGM = [0.05 0.15 0.2 0.28 0.5 0.7 1.0];
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
r1_SPSamplingPoints = 100;

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
ESIVR = r1Data.ESIVR;
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
    ,gmFieldStrengthArray,r1GMCollection);
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


%% fitting of literature values
[wmFactor,wmExponent,wmCovariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    wmFieldStrengthArray,r1WMAvg,[0.5 0.5],r1WMSTD.^2);
r1Data.wmFactor = wmFactor;
r1Data.wmExponent = wmExponent;
r1Data.wmCovariance = wmCovariance;
[gmFactor,gmExponent,gmCovariance] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    gmFieldStrengthArray,r1GMAvg,[0.5 0.5],r1GMSTD.^2);
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


%% determine R1_SP with IE model
fprintf("<strong> === INTERMEDIATE EXCHANGE MODEL ===</strong>\n");
% ==== parameter set up ====
% -> time set up
durationForEachFieldStrength = 15; % seconds
deltaTForEachFieldStrength = 1e-3; % seconds
deltaT = deltaTForEachFieldStrength(1);
duration = durationForEachFieldStrength(1);
timeStepCount = duration/deltaT;
timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;

exchangeRate = mean(exchangeRatesTable.exchangeRateWholePoolNorm(:));
macroMolecularPoolFraction = mean(...
    exchangeRatesTable.macromolecularPoolfraction(:));

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
    ,exchange_proteinToWater_WM,timeStepCount,deltaT);
r1Data.r1_SPWM_IEModel = bestGuessR1_SP_WM_IE;
r1Data.bestGuessTissueR1_WM_IEModel = bestGuessTissueR1_WM_IE;

[r1_SP_WMAvg_IE,r1_SP_WMStd_IE,~] = getAveragSTDandCountsOfR1Collection( ...
    bestGuessR1_SP_WM_IE);
r1Data.r1_SP_WMAvg_IEModel = r1_SP_WMAvg_IE;
r1Data.r1_SP_WMStd_IEModel = r1_SP_WMStd_IE;

[factorR1_SP_WM_IE,exponentR1_SP_WM_IE,covarianceR1_SP_WM_IE] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction(wmFieldStrengthArray ...
    ,r1_SP_WMAvg_IE,[0.5 0.5],r1_SP_WMStd_IE);
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
    ,exchange_proteinToWater_GM,timeStepCount,deltaT);
r1Data.r1_SPGM_IEModel = bestGuessR1_SP_GM_IE;
r1Data.bestGuessTissueR1_GM_IEModel = bestGuessTissueR1_GM_IE;

[r1_SP_GMAvg_IE,r1_SP_GMStd_IE,~] = getAveragSTDandCountsOfR1Collection( ...
    bestGuessR1_SP_GM_IE);
r1Data.r1_SP_GMAvg_IEModel = r1_SP_GMAvg_IE;
r1Data.r1_SP_GMStd_IEModel = r1_SP_GMStd_IE;

[factorR1_SP_GM_IE,exponentR1_SP_GM_IE,covarianceR1_SP_GM_IE] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction(gmFieldStrengthArray ...
    ,r1_SP_GMAvg_IE,[0.5 0.5],r1_SP_GMStd_IE);
r1Data.factorR1_SPGM_IEModel = factorR1_SP_GM_IE;
r1Data.exponentR1_SPGM_IEModel = exponentR1_SP_GM_IE;
r1Data.covarianceR1_SPGM_IEModel = covarianceR1_SP_GM_IE;

%% determine R1SP with larger exchange rates
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
    ,exchange_proteinToWater_WMincr,timeStepCount,deltaT);
r1Data.r1_SPWM_IEModel_incr = bestGuessR1_SP_WM_IE_incr;
r1Data.bestGuessTissueR1_WM_IEModel_incr = bestGuessTissueR1_WM_IE_incr;

[r1_SP_WMAvg_IE_incr,r1_SP_WMStd_IE_incr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_WM_IE_incr);
r1Data.r1_SP_WMAvg_IEModel_incr = r1_SP_WMAvg_IE_incr;
r1Data.r1_SP_WMStd_IEModel_incr = r1_SP_WMStd_IE_incr;

[factorR1_SP_WM_IE_incr,exponentR1_SP_WM_IE_incr ...
    ,covarianceR1_SP_WM_IE_incr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    wmFieldStrengthArray,r1_SP_WMAvg_IE_incr,[0.5 0.5] ...
    ,r1_SP_WMStd_IE_incr);
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
    ,exchange_proteinToWater_WMdecr,timeStepCount,deltaT);
r1Data.r1_SPWM_IEModel_decr = bestGuessR1_SP_WM_IE_decr;
r1Data.bestGuessTissueR1_WM_IEModel_decr = bestGuessTissueR1_WM_IE_decr;

[r1_SP_WMAvg_IE_decr,r1_SP_WMStd_IE_decr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_WM_IE_decr);
r1Data.r1_SP_WMAvg_IEModel_decr = r1_SP_WMAvg_IE_decr;
r1Data.r1_SP_WMStd_IEModel_decr = r1_SP_WMStd_IE_decr;

[factorR1_SP_WM_IE_decr,exponentR1_SP_WM_IE_decr ...
    ,covarianceR1_SP_WM_IE_decr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    wmFieldStrengthArray,r1_SP_WMAvg_IE_decr,[0.5 0.5] ...
    ,r1_SP_WMStd_IE_decr);
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
    ,exchange_proteinToWater_GMincr,timeStepCount,deltaT);
r1Data.r1_SPGM_IEModel_incr = bestGuessR1_SP_GM_IE_incr;
r1Data.bestGuessTissueR1_GM_IEModel_incr = bestGuessTissueR1_GM_IE_incr;

[r1_SP_GMAvg_IE_incr,r1_SP_GMStd_IE_incr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_GM_IE_incr);
r1Data.r1_SP_GMAvg_IEModel_incr = r1_SP_GMAvg_IE_incr;
r1Data.r1_SP_GMStd_IEModel_incr = r1_SP_GMStd_IE_incr;

[factorR1_SP_GM_IE_incr,exponentR1_SP_GM_IE_incr ...
    ,covarianceR1_SP_GM_IE_incr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    gmFieldStrengthArray,r1_SP_GMAvg_IE_incr,[0.5 0.5] ...
    ,r1_SP_GMStd_IE_incr);
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
    ,exchange_proteinToWater_GMdecr,timeStepCount,deltaT);
r1Data.r1_SPGM_IEModel_decr = bestGuessR1_SP_GM_IE_decr;
r1Data.bestGuessTissueR1_GM_IEModel_decr = bestGuessTissueR1_GM_IE_decr;

[r1_SP_GMAvg_IE_decr,r1_SP_GMStd_IE_decr,~] ...
    = getAveragSTDandCountsOfR1Collection(bestGuessR1_SP_GM_IE_decr);
r1Data.r1_SP_GMAvg_IEModel_decr = r1_SP_GMAvg_IE_decr;
r1Data.r1_SP_GMStd_IEModel_decr = r1_SP_GMStd_IE_decr;

[factorR1_SP_GM_IE_decr,exponentR1_SP_GM_IE_decr ...
    ,covarianceR1_SP_GM_IE_decr] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    gmFieldStrengthArray,r1_SP_GMAvg_IE_decr,[0.5 0.5] ...
    ,r1_SP_GMStd_IE_decr);
r1Data.factorR1_SPGM_IEModel_decr = factorR1_SP_GM_IE_decr;
r1Data.exponentR1_SPGM_IEModel_decr = exponentR1_SP_GM_IE_decr;
r1Data.covarianceR1_SPGM_IEModel_decr = covarianceR1_SP_GM_IE_decr;


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
    fitFieldStrengthBehaviorOfR1ToPowerFunction(wmFieldStrengthArray(1:end) ...
    ,r1_SP_WMAvg_FE(1:end),[0.5 0.5],r1_SP_WMStd_FE(1:end));
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
    = fitFieldStrengthBehaviorOfR1ToPowerFunction(gmFieldStrengthArray(1:end) ...
    ,r1_SP_GMAvg_FE(1:end),[0.5 0.5],r1_SP_GMStd_FE(1:end));
r1Data.factorR1_SPGM_FEModel = factorR1_SP_GM_FE;
r1Data.exponentR1_SPGM_FEModel = exponentR1_SP_GM_FE;
r1Data.covarianceR1_SPGM_FEModel = covarianceR1_SP_GM_FE;

%% fitting R1_SL and R1_LW to power function
[factorR1_SL,exponentR1_SL,covarianceR1_SL] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    r1Data.fieldStrengths,r1Data.r1Eff_SM,[0.5 0.5]);
r1Data.factorR1_SL = factorR1_SL;
r1Data.exponentR1_LW = exponentR1_SL;
r1Data.convarianceR1_SL = covarianceR1_SL;

[factorR1_LW,exponentR1_LW,covarianceR1_LW] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    r1Data.fieldStrengths,r1Data.r1Hist_MW,[0.5 0.5]);
r1Data.factorR1_LW = factorR1_LW;
r1Data.exponentR1_LW = exponentR1_LW;
r1Data.covarianceR1_LW = covarianceR1_LW;

%% determine an average R1_SP based on WM and GM in IE model




%% plotting results of differnt models

% ==== plotting of literature values
legendEntries = {};
fieldStrengthAxis = 0.01:0.01:7.1;

fig = initializeFigure();
initializeSubplot(fig,2,2,1);
xlabel("Field strengths [T]");
ylabel("R$_{1,literature}$ [Hz]");
title("$\textbf{a}$ Literature values");

plt(1) = errorbar(wmFieldStrengthArray,r1WMAvg,r1WMSTD ...
    ,'LineStyle','none','LineWidth',1.5);
plt(2) = errorbar(gmFieldStrengthArray,r1GMAvg,r1GMSTD ...
    ,'LineStyle','none','LineWidth',1.5);
legendEntries{end+1} = "avg WM $\pm$ STD";
legendEntries{end+1} = "avg GM $\pm$ STD";

plt(3) = plot(fieldStrengthAxis,wmFactor*fieldStrengthAxis.^( ...
    -wmExponent));
drawSTDRegionAroundPowerFunction(wmFactor,wmExponent...
    ,wmCovariance,fieldStrengthAxis,plt(3).Color);

plt(4) = plot(fieldStrengthAxis,gmFactor*fieldStrengthAxis.^(-gmExponent));
drawSTDRegionAroundPowerFunction(gmFactor,gmExponent...
    ,gmCovariance,fieldStrengthAxis,plt(4).Color);
legendEntries{end+1} = "fit WM";
legendEntries{end+1} = "fit GM";

legend(plt(1:4),legendEntries);

% ==== molecular dynamics simulation results ====

initializeSubplot(fig,2,2,2);
skip = 10;
plt(1) = plot(r1Data.fieldStrengths(1:skip:end) ...
    ,r1Data.r1Eff_SM(1:skip:end),'*');
plt(2) = plot(r1Data.fieldStrengths, factorR1_SL  ...
    * r1Data.fieldStrengths .^(-exponentR1_SL),'Color', plt(1).Color);
plt(3) = plot(r1Data.fieldStrengths(1:skip:end) ...
    ,r1Data.r1Hist_MW(1:skip:end),'*');
plt(4) = plot(r1Data.fieldStrengths, factorR1_LW ...
    * r1Data.fieldStrengths .^(-exponentR1_LW),'Color',plt(3).Color);
axis([0 inf 0 20]);

legend(plt(1:4),"SL MD Sim","SL fit","LW MD sim","LW fit");
xlabel("Field strength [T]");
ylabel("R$_{1,simulation}$ [Hz]");
title("$\textbf{b}$ MD simulation results")


% ==== model predictions
initializeSubplot(fig,2,2,3);
title("$\textbf{c}$ IE model results");
ylabel("R$_{1,SP}$ [Hz]");
xlabel("Field strength [T]");

% ==== IE model ====
plt(1) = errorbar(wmFieldStrengthArray,r1_SP_WMAvg_IE,r1_SP_WMStd_IE ...
    ,'LineStyle','none','LineWidth',1.5);
plt(2) = plot(fieldStrengthAxis,factorR1_SP_WM_IE*fieldStrengthAxis.^( ...
    - exponentR1_SP_WM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE,exponentR1_SP_WM_IE...
    ,covarianceR1_SP_WM_IE,fieldStrengthAxis,plt(2).Color);

plt(3) = errorbar(gmFieldStrengthArray,r1_SP_GMAvg_IE,r1_SP_GMStd_IE ...
    ,'LineStyle','none','LineWidth',1.5);
plt(4) = plot(fieldStrengthAxis,factorR1_SP_GM_IE*fieldStrengthAxis.^( ...
    - exponentR1_SP_GM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE,exponentR1_SP_GM_IE...
    ,covarianceR1_SP_GM_IE,fieldStrengthAxis,plt(4).Color);

axis([0 inf 0 25])
legend(plt(1:4),"WM $\pm$ STD","WM fit","GM $\pm$ STD","GM fit");

% ==== FE model ====
initializeSubplot(fig,2,2,4);
title("$\textbf{d}$ FE model results");
ylabel("R$_{1,SP}$ [Hz]");
xlabel("Field strength [T]");
plt(5) = errorbar(wmFieldStrengthArray,r1_SP_WMAvg_FE,r1_SP_WMStd_FE ...
    ,'LineStyle','none','LineWidth',1.5);
plt(6) = plot(fieldStrengthAxis,factorR1_SP_WM_FE*fieldStrengthAxis.^( ...
    - exponentR1_SP_WM_FE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_FE,exponentR1_SP_WM_FE...
    ,covarianceR1_SP_WM_FE,fieldStrengthAxis,plt(6).Color);

plt(7) = errorbar(gmFieldStrengthArray,r1_SP_GMAvg_FE,r1_SP_GMStd_FE ...
    ,'LineStyle','none','LineWidth',1.5);
plt(8) = plot(fieldStrengthAxis,factorR1_SP_GM_FE*fieldStrengthAxis.^( ...
    - exponentR1_SP_GM_FE));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_FE,exponentR1_SP_GM_FE...
    ,covarianceR1_SP_GM_FE,fieldStrengthAxis,plt(8).Color);

axis([0 inf 0 25])
lgd = legend();
lgd.Visible = 'off';

% legend(plt(1:4),"IE WM $\pm$ STD","IE WM fit","IE GM $\pm$ STD" ...
%     ,"IE GM fit","FE WM $\pm$ STD","FE WM fit","FE GM $\pm$ STD" ...
%     ,"FE GM fit");

if saving
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"TissueR1_ComparmentR1_R1","SP");
end

%% Plotting results of different exchange rates

fig = initializeFigure();
initializeSubplot(fig,2,1,1);
title("WM");
ylabel("R$_{1,SP}$");
xlabel("Field strength [T]");
plt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE*fieldStrengthAxis.^(-exponentR1_SP_WM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE,exponentR1_SP_WM_IE ...
    ,covarianceR1_SP_WM_IE,fieldStrengthAxis,plt(1).Color);

plt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE_incr*fieldStrengthAxis.^(-exponentR1_SP_WM_IE_incr));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE_incr ...
    ,exponentR1_SP_WM_IE_incr,covarianceR1_SP_WM_IE_incr ...
    ,fieldStrengthAxis,plt(2).Color);

plt(3) = plot(fieldStrengthAxis ...
    ,factorR1_SP_WM_IE_decr*fieldStrengthAxis.^(-exponentR1_SP_WM_IE_decr));
drawSTDRegionAroundPowerFunction(factorR1_SP_WM_IE_decr ...
    ,exponentR1_SP_WM_IE_decr,covarianceR1_SP_WM_IE_decr ...
    ,fieldStrengthAxis,plt(3).Color);
axis([0 inf 0 25]);

legend(plt(1:3),"unchanged","increased","decreased");

initializeSubplot(fig,2,1,2);
title("GM");
ylabel("R$_{1,SP}$");
xlabel("Field strength [T]");
plt(4) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE*fieldStrengthAxis.^(-exponentR1_SP_GM_IE));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE,exponentR1_SP_GM_IE ...
    ,covarianceR1_SP_GM_IE,fieldStrengthAxis,plt(4).Color);

plt(5) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE_incr*fieldStrengthAxis.^(-exponentR1_SP_GM_IE_incr));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE_incr ...
    ,exponentR1_SP_GM_IE_incr,covarianceR1_SP_GM_IE_incr ...
    ,fieldStrengthAxis,plt(5).Color);

plt(6) = plot(fieldStrengthAxis ...
    ,factorR1_SP_GM_IE_decr*fieldStrengthAxis.^(-exponentR1_SP_GM_IE_decr));
drawSTDRegionAroundPowerFunction(factorR1_SP_GM_IE_decr ...
    ,exponentR1_SP_GM_IE_decr,covarianceR1_SP_GM_IE_decr ...
    ,fieldStrengthAxis,plt(6).Color);


axis([0 inf 0 25]);
lgd = legend();
lgd.Visible = 'off';

if saving
    saveFigureTo(r1DataFolderName ...
        ,datestr(now,'yyyymmdd'),"ComparingIncreasedAndDecreasedExchOn" ...
        ,"R1_SP");
end


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
fprintf("  IE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_WM_IE ...
    ,exponentR1_SP_WM_IE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_WM_IE);

fprintf("  FE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_WM_FE ...
    ,exponentR1_SP_WM_FE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_WM_FE);
fprintf("\n");

fprintf(" Based on GM: \n");
fprintf("  IE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_GM_IE ...
    ,exponentR1_SP_GM_IE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_GM_IE);

fprintf("  FE: R1_SP(B_0) = %.4f * B_0 ^(-%.4f) \n",factorR1_SP_GM_FE ...
    ,exponentR1_SP_GM_FE);
fprintf("  variance/covariance matrix:\n")
fprintf("  %.7f %.7f \n",covarianceR1_SP_GM_FE);




%% saving

if saving
    save(r1DataFolderName + "SM_MW_SP_histCompartmentAndCrossR1" ...
        ,'-struct','r1Data'); %#ok<UNRCH>
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
    ,exchange_proteinToWater,timeStepCount,deltaT)

bestGuessTissueR1 = cell(1,length(fieldStrengthArray));
bestGuessR1_SP = cell(1,length(fieldStrengthArray));

for fieldStrengthNr = 1:length(fieldStrengthArray)
    fieldStrength = fieldStrengthArray(fieldStrengthNr);
    fprintf(" == %.2f Tesla ==  \n",fieldStrength);
    
    fieldStrengthIndex = ...
        r1Data.fieldStrengths > fieldStrength - 0.0001 ...
        & r1Data.fieldStrengths < fieldStrength + 0.0001;
    if isempty(fieldStrengthIndex)
        error("No field strength found");
    end
    
    r1_SL = r1Data.r1Eff_SM(fieldStrengthIndex);
    fprintf(" R1_SL = %.4f \n",r1_SL);
    r1_LW = r1Data.r1Hist_MW(fieldStrengthIndex);
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
    ,gmFieldStrengthArray,r1GMCollection)

r1WMData = createMatrixOfDataPoints(r1WMCollection);
r1GMData = createMatrixOfDataPoints(r1GMCollection);

[r1AvgWMData,~,wmCounts] ...
    = getAveragSTDandCountsOfR1Collection(r1WMCollection);
[r1AvgGMData,~,gmCounts] ...
    = getAveragSTDandCountsOfR1Collection(r1GMCollection);

fig = initializeFigure('legend',false);
initializeSubplot(fig,2,1,1);
labels = strsplit(sprintf("%.3f (#%i) \n" ...
    ,[wmFieldStrengthArray;wmCounts]),'\n');
title("WM");
boxplot(r1WMData','Labels',labels(1:end-1));
plot(r1AvgWMData,'g*');
xlabel("Field strength [T]");
ylabel("R$_{1,WM}$ [Hz]");
lgd = legend();
lgd.Visible = 'off';

initializeSubplot(fig,2,1,2);
labels = strsplit(sprintf("%.3f (#%i) \n" ...
    ,[gmFieldStrengthArray;gmCounts]),'\n');
title("GM");
boxplot(r1GMData','Labels',labels(1:end-1));
plot(r1AvgGMData,'g*');
xlabel("Field strength [T]");
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




