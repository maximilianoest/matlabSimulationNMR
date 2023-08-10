clc; clear all; close all;


addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath("../molecularDynamicSimulationAnalysis/"));
addpath(genpath("../distributionAnalysis/"));
addpath(genpath("../brainConstitutionAnalysis/"));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));

saving = true;
scriptImagesFolderName = "C:\Users\maxoe\Google Drive\Promotion\Literatur\Lipid and Protein R1 Relaxation Dispersion in the Human Brain\images\";
dataForPaperFolderName = "C:\Users\maxoe\Google Drive\Promotion\Literatur\Lipid and Protein R1 Relaxation Dispersion in the Human Brain\data\";



colorArray = cell(1,7);
fig = initializeFigure();
for i = 1:7
    plt = plot(1,1);
    colorArray{i} = plt.Color;
    
end
close(fig);

densDistribDirectory = sprintf("..%s..%s..%sRESULTS%swholeMyelin" ...
    + "_densityDistributions%s",createFilesepStringArray(5));
densDistribFileName = "20230405_myelinDistributionData_averaged";
densData = load(densDistribDirectory + densDistribFileName + ".mat");


load("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\wholeMyelin_relaxationRates\SM_MW_SP_histCompartmentAndCrossR1.mat");

% === plot configuration ===
axisFontSize = 14;
legendFontSize = 12;
fieldStrengthAxis = 0.001:0.01:7.01;
highestR1InPlot = 100;
highestR1InLitPlot = 8;
lowestFieldStrength = min(wmFieldStrengthArray) ...
    - min(wmFieldStrengthArray)*0.1;

[r1WMAvgForFitting,r1WMStdForFitting,~] ...
    = getAveragSTDandCountsOfR1Collection(r1WMCollectionForFitting);
[r1GMAvgForFitting,r1GMStdForFitting,~] ...
    = getAveragSTDandCountsOfR1Collection(r1GMCollectionForFitting);

spIndicesForWMFitting = getIndicesInFieldStrengthArry( ...
    wmFieldStrengthArray,wmFieldStrengthArrayForFitting);
spIndicesForGMFitting = getIndicesInFieldStrengthArry( ...
    gmFieldStrengthArray,gmFieldStrengthArrayForFitting);

% ==== denstiy data stuff
windowSize = 21;
surfWaterDensBorder = 0.95;
borderShiftOffset = 0;

%% plotting literature values (Observed R1 in WM and GM)
literatureFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
literatureLegendEntries = {};

wmSubplot = initializeSubplot(literatureFigure,2,1,1);
title("$\textbf{a)}$");
ylabel("R$_{1,WM}$ [Hz]");
xlabel("Field strength [T]");
xticks(0:0.5:7);

literaturePlt(1) = errorbar(wmFieldStrengthArray,r1WMAvg,r1WMSTD ...
    ,'LineStyle','none','LineWidth',1.5,'Color','g');
literatureLegendEntries{end+1} = "R$_{1,WM}$";
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
title("$\textbf{b)}$");

literaturePlt(4) = errorbar(gmFieldStrengthArray,r1GMAVg,r1GMSTD ...
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

%% plotting literature values (Observed R1 in WM and GM) 2.0 -> loglog
literatureFigure2_0 = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);

ylabel("R$_{1}$, logarithmic [Hz]");
set(gca,'YScale','log','XScale','log');
xlabel("Field strength, logarithmic [T]");
xticks([0:0.2:1 2:fieldStrengthAxis(end)]);
yticks([0:0.2:1 2:10])
axis([lowestFieldStrength ...
    fieldStrengthAxis(end)+1 0.4 highestR1InLitPlot]);

% == WM ==
literaturePlt2_0(1) = plot(fieldStrengthAxis,wmFactor*fieldStrengthAxis ...
    .^(-wmExponent),'Color',colorArray{1});
drawSTDRegionAroundPowerFunction(wmFactor,wmExponent...
    ,wmCovariance,fieldStrengthAxis,literaturePlt2_0(1).Color);
literaturePlt2_0(2) = errorbar(wmFieldStrengthArray( ...
    spIndicesForWMFitting),r1WMAvg(spIndicesForWMFitting) ...
    ,r1WMSTD(spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_0(1).Color);
literaturePlt2_0(3) = plot(wmFieldStrengthArray(~spIndicesForWMFitting) ...
    ,r1WMAvg(~spIndicesForWMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',colorArray{1},'Marker','o');

% == GM ==
literaturePlt2_0(4) = plot(fieldStrengthAxis,gmFactor*fieldStrengthAxis ...
    .^(-gmExponent),'Color',colorArray{2});
drawSTDRegionAroundPowerFunction(gmFactor,gmExponent...
    ,gmCovariance,fieldStrengthAxis,literaturePlt2_0(4).Color);
literaturePlt2_0(5) = errorbar(gmFieldStrengthArray( ...
    spIndicesForGMFitting),r1GMAVg(spIndicesForGMFitting) ...
    ,r1GMSTD(spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_0(4).Color);
literaturePlt2_0(6) = plot(gmFieldStrengthArray(~spIndicesForGMFitting) ...
    ,r1GMAVg(~spIndicesForGMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',colorArray{2},'Marker','o');

legend(literaturePlt2_0([1 4]),"WM","GM");


%% plotting literature values (Observed R1 in WM and GM) 2.1 -> ylog
literatureFigure2_1 = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
literatureLegendEntries = {};

ylabel("R$_{1}$ [Hz]");
set(gca,'YScale','log');
xlabel("Field strength [T]");
xticks([0:fieldStrengthAxis(end)]);
yticks([0:0.2:1 2:10])
axis([min(wmFieldStrengthArray) - min(wmFieldStrengthArray)*0.1 ...
    fieldStrengthAxis(end) 0.4 highestR1InLitPlot]);

% == WM ==
literaturePlt2_1(1) = plot(fieldStrengthAxis,wmFactor*fieldStrengthAxis ...
    .^(-wmExponent),'Color',colorArray{1});
drawSTDRegionAroundPowerFunction(wmFactor,wmExponent...
    ,wmCovariance,fieldStrengthAxis,literaturePlt2_1(1).Color);
literaturePlt2_1(2) = errorbar(wmFieldStrengthArray( ...
    spIndicesForWMFitting),r1WMAvg(spIndicesForWMFitting) ...
    ,r1WMSTD(spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_1(1).Color);
literaturePlt2_1(3) = plot(wmFieldStrengthArray(~spIndicesForWMFitting) ...
    ,r1WMAvg(~spIndicesForWMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',literaturePlt2_1(1).Color,'Marker','o');

% == GM ==
literaturePlt2_1(4) = plot(fieldStrengthAxis,gmFactor*fieldStrengthAxis ...
    .^(-gmExponent),'Color',colorArray{2});
drawSTDRegionAroundPowerFunction(gmFactor,gmExponent...
    ,gmCovariance,fieldStrengthAxis,literaturePlt2_1(4).Color);
literaturePlt2_1(5) = errorbar(gmFieldStrengthArray( ...
    spIndicesForGMFitting),r1GMAVg(spIndicesForGMFitting) ...
    ,r1GMSTD(spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_1(4).Color);
literaturePlt2_1(6) = plot(gmFieldStrengthArray(~spIndicesForGMFitting) ...
    ,r1GMAVg(~spIndicesForGMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',literaturePlt2_1(4).Color,'Marker','o');

legend(literaturePlt2_1([1 4]),"WM","GM");

%% plotting literature values (Observed R1 in WM and GM) 2.2 -> xlog
literatureFigure2_2 = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
literatureLegendEntries = {};

ylabel("R$_{1}$ [Hz]");
set(gca,'XScale','log');
xlabel("Field strength[T]");
% xticks([0:fieldStrengthAxis(end)]);
xticks([0.02:0.04:0.1 0.2:0.2:1 2:2:10])
axis([min(wmFieldStrengthArray) - min(wmFieldStrengthArray)*0.1 ...
    fieldStrengthAxis(end)+1 0.4 highestR1InLitPlot]);

% == WM ==
literaturePlt2_2(1) = plot(fieldStrengthAxis,wmFactor*fieldStrengthAxis ...
    .^(-wmExponent),'Color',colorArray{1});
drawSTDRegionAroundPowerFunction(wmFactor,wmExponent...
    ,wmCovariance,fieldStrengthAxis,literaturePlt2_2(1).Color);
literaturePlt2_2(2) = errorbar(wmFieldStrengthArray( ...
    spIndicesForWMFitting),r1WMAvg(spIndicesForWMFitting) ...
    ,r1WMSTD(spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_2(1).Color);
literaturePlt2_2(3) = plot(wmFieldStrengthArray(~spIndicesForWMFitting) ...
    ,r1WMAvg(~spIndicesForWMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',literaturePlt2_1(1).Color,'Marker','o');

% == GM ==
literaturePlt2_2(4) = plot(fieldStrengthAxis,gmFactor*fieldStrengthAxis ...
    .^(-gmExponent),'Color',colorArray{2});
drawSTDRegionAroundPowerFunction(gmFactor,gmExponent...
    ,gmCovariance,fieldStrengthAxis,literaturePlt2_2(4).Color);
literaturePlt2_2(5) = errorbar(gmFieldStrengthArray( ...
    spIndicesForGMFitting),r1GMAVg(spIndicesForGMFitting) ...
    ,r1GMSTD(spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',literaturePlt2_2(4).Color);
literaturePlt2_2(6) = plot(gmFieldStrengthArray(~spIndicesForGMFitting) ...
    ,r1GMAVg(~spIndicesForGMFitting),'LineStyle','none','LineWidth' ...
    ,1.5,'Color',literaturePlt2_2(4).Color,'Marker','o');

legend(literaturePlt2_2([1 4]),"WM","GM");

%% plotting R1_SL and average R1_SP
mdSimFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
indicesOfInterest = (fieldStrengths ...
    >= fieldStrengthAxis(1)) & (fieldStrengths ...
    <= fieldStrengthAxis(end));
skip = 20;
fieldStrengthsOfInterest = fieldStrengths(indicesOfInterest);
% skipArray = [1:5 6:10:95 96:100:length(fieldStrengthsOfInterest)];
skipArray = 1:skip:length(fieldStrengthsOfInterest);
r1_SLOfInterest = r1Eff_SM(indicesOfInterest);

slSimResultsSubplot = initializeSubplot(mdSimFigure,2,1,1);

% = R1_SL
slMdSimPlt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SL*fieldStrengthAxis.^(-exponentR1_SL));
drawSTDRegionAroundPowerFunction(factorR1_SL,exponentR1_SL...
    ,covarianceR1_SL,fieldStrengthAxis,colorArray{1});

slMdSimPlt(2) = plot(fieldStrengthsOfInterest(skipArray) ...
    ,r1_SLOfInterest(skipArray),'Marker','*','LineStyle','none' ...
    ,'Color',slMdSimPlt(1).Color);

% = R1_SP_avg
slMdSimPlt(3) = plot(fieldStrengthAxis,factorR1_SPAvg_IEModel ...
    *fieldStrengthAxis.^(-exponentR1_SPAvg_IEModel),'Color',colorArray{2});
drawSTDRegionAroundPowerFunction(factorR1_SPAvg_IEModel ...
    ,exponentR1_SPAvg_IEModel,covarianceR1_SPAvg_IEModel ...
    ,fieldStrengthAxis,slMdSimPlt(3).Color);
slMdSimPlt(4) = errorbar(wmFieldStrengthArray(spIndicesForWMFitting) ...
    , r1_SP_AvgTissue_IEModel(spIndicesForWMFitting) ...
    , r1_SP_StatisticalError(spIndicesForWMFitting) ...
    + r1_SP_systematicTissueBasedError(spIndicesForWMFitting)  ...
    + r1_SP_systematicNegativeExchangeRateBasedError( ...
    spIndicesForWMFitting), r1_SP_StatisticalError( ...
    spIndicesForWMFitting) + r1_SP_systematicTissueBasedError( ...
    spIndicesForWMFitting) ...
    + r1_SP_systematicPositiveExchangeRateBasedError( ...
    spIndicesForWMFitting),'Color',slMdSimPlt(3).Color,'Marker','*' ...
    ,'LineStyle','none','LineWidth',get(gcf,'DefaultLineLineWidth'));

legend(slMdSimPlt([1 3]),"R$_{1,SL}$","R$_{1,SP}$");
axis([min(wmFieldStrengthArray) - min(wmFieldStrengthArray)*0.1 ...
    fieldStrengthAxis(end) 0 highestR1InPlot]);

ylabel("R$_{1}$ [Hz]");
set(gca,'YScale','log');
% xticks([0.02:0.04:0.1 0.2:0.2:1 2:2:10])
xlabel("Field strength [T]");
% yticks([0.2:0.4:1 2:4:10 20:40:80])
title("$\textbf{a)}$");

% tissue based R1_SP
set(0,'CurrentFigure', mdSimFigure);
ieModelResultsSubPlt = initializeSubplot(mdSimFigure,2,1,2);
cla(ieModelResultsSubPlt);

ieModelResultsPlt(1) = plot(fieldStrengthAxis,factorR1_SPWM_IEModel ...
    *fieldStrengthAxis.^(-exponentR1_SPWM_IEModel),'Color',colorArray{4});
drawSTDRegionAroundPowerFunction(factorR1_SPWM_IEModel ...
    ,exponentR1_SPWM_IEModel,covarianceR1_SPWM_IEModel ...
    ,fieldStrengthAxis,ieModelResultsPlt(1).Color);
ieModelResultsPlt(2) = errorbar(wmFieldStrengthArray( ...
    spIndicesForWMFitting),r1_SP_WMAvg_IEModel(spIndicesForWMFitting) ...
    ,r1_SP_WMStd_IEModel(spIndicesForWMFitting),'*','LineWidth' ...
    ,get(gcf,'DefaultLineLineWidth'),'Color',ieModelResultsPlt(1).Color);
ieModelResultsPlt(3) = plot(wmFieldStrengthArray( ...
    ~spIndicesForWMFitting),r1_SP_WMAvg_IEModel(~spIndicesForWMFitting) ...
    ,'LineStyle','none','Marker','o' ...
    ,'Color',ieModelResultsPlt(1).Color);

ieModelResultsPlt(4) = plot(fieldStrengthAxis,factorR1_SPGM_IEModel ...
    *fieldStrengthAxis.^(-exponentR1_SPGM_IEModel),'Color',colorArray{5});
drawSTDRegionAroundPowerFunction(factorR1_SPGM_IEModel ...
    ,exponentR1_SPGM_IEModel,covarianceR1_SPGM_IEModel ...
    ,fieldStrengthAxis,ieModelResultsPlt(4).Color);
ieModelResultsPlt(5) = errorbar(gmFieldStrengthArray( ...
    spIndicesForGMFitting),r1_SP_GMAvg_IEModel(spIndicesForGMFitting) ...
    ,r1_SP_GMStd_IEModel(spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',get(gcf,'DefaultLineLineWidth'),'Marker','*' ...
    ,'Color',ieModelResultsPlt(4).Color);
ieModelResultsPlt(6) = plot(gmFieldStrengthArray( ...
    ~spIndicesForGMFitting),r1_SP_GMAvg_IEModel(~spIndicesForGMFitting) ...
    ,'Marker','o','LineStyle','none','Color',ieModelResultsPlt(4).Color);


axis([min(wmFieldStrengthArray) - min(wmFieldStrengthArray)*0.1 ...
    fieldStrengthAxis(end) 0 40]);
xlabel("Field strength [T]");
set(gca,'YScale','log');
% xticks([0.02:0.04:0.1 0.2:0.2:1 2:2:10])
ylabel("R$_{1,SP}$ [Hz]");
title("$\textbf{b)}$");
legend(ieModelResultsPlt([1 4]),"WM","GM");
drawnow;

%% plotting STD at WM and GM to show overlap
% additionalFigure = initializeFigure('axisFontSize',axisFontSize ...
%     ,'legendFontSize',legendFontSize);
% r1SPStdSubplot = initializeSubplot(additionalFigure,2,2,1);
% cla(r1SPStdSubplot);
% r1SPStdLegendEntries = {};
% drawSTDRegionAroundPowerFunction(factorR1_SPWM_IEModel ...
%     ,exponentR1_SPWM_IEModel,covarianceR1_SPWM_IEModel ...
%     ,fieldStrengthAxis,slMdSimPlt(1).Color,0.5);
% r1SPStdLegendEntries{end+1} = "STD WM";
%
% drawSTDRegionAroundPowerFunction(factorR1_SPGM_IEModel ...
%     ,exponentR1_SPGM_IEModel,covarianceR1_SPGM_IEModel ...
%     ,fieldStrengthAxis,slMdSimPlt(2).Color,0.5);
% r1SPStdLegendEntries{end+1} = "STD GM";
% axis([0 fieldStrengthAxis(end) 0 highestR1InPlot]);
%
% legend(r1SPStdLegendEntries);
% xlabel("Field Strength [T]");
% ylabel("R$_{1,SP}$ [Hz]");
% title("$\textbf{a)}$");
% drawnow;

%% plotting R1_LW and corresponding fit

appendixFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
lwSimResultsSubplot = initializeSubplot(mdSimFigure,2,2,1:2);

lwFieldStrengths = fieldStrengths(indicesOfInterest);
lwRelaxationRates = r1Hist_MW(indicesOfInterest);

lwMdSimPlt(1) = plot(lwFieldStrengths(1:skip:end) ...
    ,lwRelaxationRates(1:skip:end),'Color',colorArray{1},'Marker','*' ...
    ,'LineStyle','none');
lwMdSimPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_LW*fieldStrengthAxis.^(-exponentR1_LW) ...
    ,'Color',lwMdSimPlt(1).Color);
drawSTDRegionAroundPowerFunction(factorR1_LW,exponentR1_LW...
    ,covarianceR1_LW,fieldStrengthAxis,lwMdSimPlt(1).Color);

ylabel("R$_{1,LW}$ [Hz]");
xlabel("Field strength [T]");
title("$\textbf{a)}$");

lgd = legend();
lgd.Visible = 'off';
axis([0 fieldStrengthAxis(end) 0 3]);
drawnow;

% Plotting results of different exchange rates
set(0,'CurrentFigure',appendixFigure);
r1WmChangedExchangeRaterSubPlt = initializeSubplot(appendixFigure,2,2,3);
cla(r1WmChangedExchangeRaterSubPlt);
ylabel("R$_{1,SP}$ based on WM");
xlabel("Field strength [T]");
% set(gca,'YScale','log');
title("$\textbf{b)}$");
r1WmChangedExchangeRaterSubPlt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SPWM_IEModel*fieldStrengthAxis.^(-exponentR1_SPWM_IEModel));
drawSTDRegionAroundPowerFunction(factorR1_SPWM_IEModel ...
    ,exponentR1_SPWM_IEModel,covarianceR1_SPWM_IEModel,fieldStrengthAxis ...
    ,r1WmChangedExchangeRaterSubPlt(1).Color);

r1WmChangedExchangeRaterSubPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SPWM_IEModel_incr*fieldStrengthAxis.^( ...
    -exponentR1_SPWM_IEModel_incr));
drawSTDRegionAroundPowerFunction(factorR1_SPWM_IEModel_incr ...
    ,exponentR1_SPWM_IEModel_incr,covarianceR1_SPWM_IEModel_incr ...
    ,fieldStrengthAxis,r1WmChangedExchangeRaterSubPlt(2).Color);

r1WmChangedExchangeRaterSubPlt(3) = plot(fieldStrengthAxis ...
    ,factorR1_SPWM_IEModel_decr*fieldStrengthAxis.^( ...
    -exponentR1_SPWM_IEModel_decr));
drawSTDRegionAroundPowerFunction(factorR1_SPWM_IEModel_decr ...
    ,exponentR1_SPWM_IEModel_decr,covarianceR1_SPWM_IEModel_decr ...
    ,fieldStrengthAxis,r1WmChangedExchangeRaterSubPlt(3).Color);
axis([0 fieldStrengthAxis(end) 0 15]);

legend(r1WmChangedExchangeRaterSubPlt(1:3) ...
    ,"k = 1.44Hz (default)","k = 1.58Hz (+ 10 \%)","k = 1.30Hz (- 10 \%)");
drawnow;

r1GmChangedExchangeRaterSubPlt = initializeSubplot(appendixFigure,2,2,4);
cla(r1GmChangedExchangeRaterSubPlt);
ylabel("R$_{1,SP}$ based on GM");
xlabel("Field strength [T]");
title("$\textbf{c)}$");
r1GmChangedExchangeRaterSubPlt(1) = plot(fieldStrengthAxis ...
    ,factorR1_SPGM_IEModel*fieldStrengthAxis.^(-exponentR1_SPGM_IEModel));
drawSTDRegionAroundPowerFunction(factorR1_SPGM_IEModel ...
    ,exponentR1_SPGM_IEModel,covarianceR1_SPGM_IEModel ...
    ,fieldStrengthAxis,r1GmChangedExchangeRaterSubPlt(1).Color);

r1GmChangedExchangeRaterSubPlt(2) = plot(fieldStrengthAxis ...
    ,factorR1_SPGM_IEModel_incr*fieldStrengthAxis.^( ...
    -exponentR1_SPGM_IEModel_incr));
drawSTDRegionAroundPowerFunction(factorR1_SPGM_IEModel_incr ...
    ,exponentR1_SPGM_IEModel_incr,covarianceR1_SPGM_IEModel_incr ...
    ,fieldStrengthAxis,r1GmChangedExchangeRaterSubPlt(2).Color);

r1GmChangedExchangeRaterSubPlt(3) = plot(fieldStrengthAxis ...
    ,factorR1_SPGM_IEModel_decr*fieldStrengthAxis.^( ...
    -exponentR1_SPGM_IEModel_decr));
drawSTDRegionAroundPowerFunction(factorR1_SPGM_IEModel_decr ...
    ,exponentR1_SPGM_IEModel_decr,covarianceR1_SPGM_IEModel_decr ...
    ,fieldStrengthAxis,r1GmChangedExchangeRaterSubPlt(3).Color);

legend(r1GmChangedExchangeRaterSubPlt(1:3) ...
    ,"k = 1.44Hz (default)","k = 1.58Hz (+ 10 \%)","k = 1.30Hz (- 10 \%)");
axis([0 fieldStrengthAxis(end) 0 15]);
% lgd = legend();
% lgd.Visible = 'off';
drawnow;

%% plotting FE model results
FEModelResultsPlt = initializeFigure( ...
    'axisFontSize',axisFontSize,'legendFontSize',legendFontSize);

ylabel("R$_{1,SP}^{FE}$ [Hz]");
set(gca,'YScale','log');
xlabel("Field strength [T]");
feModelPlt(1) = errorbar(wmFieldStrengthArray(spIndicesForWMFitting) ...
    ,r1_SP_WMAvg_FEModel(spIndicesForWMFitting) ...
    ,r1_SP_WMStd_FEModel(spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','*','Color',colorArray{1});
feModelPlt(2) = plot(wmFieldStrengthArray(~spIndicesForWMFitting) ...
    ,r1_SP_WMAvg_FEModel(~spIndicesForWMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','o','Color',feModelPlt(1).Color);
feModelPlt(3) = plot(fieldStrengthAxis,factorR1_SPWM_FEModel* ...
    fieldStrengthAxis.^(- exponentR1_SPWM_FEModel),'Color' ...
    ,feModelPlt(1).Color);
drawSTDRegionAroundPowerFunction(factorR1_SPWM_FEModel ...
    ,exponentR1_SPWM_FEModel,covarianceR1_SPWM_FEModel ...
    ,fieldStrengthAxis,feModelPlt(1).Color);

feModelPlt(4) = errorbar(gmFieldStrengthArray(spIndicesForGMFitting) ...
    ,r1_SP_GMAvg_FEModel(spIndicesForGMFitting) ...
    ,r1_SP_GMStd_FEModel(spIndicesForGMFitting) ...
    ,'LineStyle','none','Marker','*','LineWidth',1.5,'Color' ...
    ,colorArray{2});
feModelPlt(5) = plot(gmFieldStrengthArray(~spIndicesForGMFitting) ...
    ,r1_SP_GMAvg_FEModel(~spIndicesForGMFitting),'LineStyle','none' ...
    ,'LineWidth',1.5,'Marker','o','Color',feModelPlt(4).Color);
feModelPlt(6) = plot(fieldStrengthAxis,factorR1_SPGM_FEModel ...
    *fieldStrengthAxis.^(- exponentR1_SPGM_FEModel),'Color' ...
    ,feModelPlt(4).Color);
drawSTDRegionAroundPowerFunction(factorR1_SPGM_FEModel ...
    ,exponentR1_SPGM_FEModel,covarianceR1_SPGM_FEModel ...
    ,fieldStrengthAxis,feModelPlt(4).Color);

axis([min(wmFieldStrengthArray) - min(wmFieldStrengthArray)*0.1 ...
    fieldStrengthAxis(end) 0 40]);
lgd = legend();
lgd.Visible = 'off';

legend(feModelPlt([3 6]),"WM","GM");
drawnow;

%% plotting density distribution and MD model
fprintf("<strong> Density distribution: </strong> \n");

densFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
bilayerImageSubplot = initializeSubplot(densFigure,2,2,1);
set(bilayerImageSubplot,'Position',get(bilayerImageSubplot,'Position')  ...
    .* [1, 1, 1.1, 1.1]);
bilayerImage = imread('C:\Users\maxoe\Google Drive\Promotion\Literatur\Lipid and Protein R1 Relaxation Dispersion in the Human Brain\images\BilayerMDModel.jpg');
imshow(bilayerImage,'border','tight');
lgd = legend();
lgd.Visible = 'off';
title("$\textbf{a)}$");

monolayerImageSubplot = initializeSubplot(densFigure,2,2,2);
set(monolayerImageSubplot,'Position',get(monolayerImageSubplot,'Position')  ...
    .* [0.9, 0.85, 1.42, 1.42]);
monolayerImage = imread('C:\Users\maxoe\Google Drive\Promotion\Literatur\Lipid and Protein R1 Relaxation Dispersion in the Human Brain\images\MonolayerMDModel.jpg');
imshow(monolayerImage,'border','tight');
lgd = legend();
lgd.Visible = 'off';
title("$\textbf{b)}$");

locations = densData.avgDiscreteLocations(:,2:end-1);
dVol = mean(densData.dVolumeForTimeStep);

waterIndex = (densData.moleculeNameCollection == "WATER");
membraneIndices = (densData.moleculeNameCollection ~= "WATER") ...
    .* (densData.moleculeNameCollection ~= "POT") ...
    .* (densData.moleculeNameCollection ~= "CLA");
hydrogenIndex = densData.atomNameCollection == "H";

densDistribWater = squeeze(sum( ...
    densData.avgDensity_locMolAtom(:,waterIndex,:),2));
densDistribMembr = squeeze(sum( ...
    densData.avgDensity_locMolAtom(:,~waterIndex,:),2));

modelWaterMass = sum(sum(densDistribWater)) * dVol;
modelMembrMass = sum(sum(densDistribMembr)) * dVol;

densDistribWaterH = meanFilter(sum( ...
    densData.avgDensity_locMolAtom(2:end-1,waterIndex,hydrogenIndex) ...
    ,2),windowSize);
densDistribMembrH = meanFilter(sum( ...
    densData.avgDensity_locMolAtom(2:end-1,~waterIndex,hydrogenIndex) ...
    ,2),windowSize);
modelWaterHMass = sum(sum(densDistribWaterH))*dVol;
modelMembrHMass = sum(sum(densDistribMembrH))*dVol;


waterMassRatio = modelWaterMass/(modelWaterMass+modelMembrMass);
fprintf("  Water mass: %.4d kg.\n  Membrane mass: %.4d kg.\n" ...
    ,modelWaterMass,modelMembrMass);
fprintf("  => water mass ratio: %.4f\n",waterMassRatio);

waterMassRatioH = modelWaterHMass/(modelWaterHMass+modelMembrHMass);
fprintf("  Water H mass: %.4d kg.\n  Membrane H mass: %.4d kg.\n" ...
    ,modelWaterHMass,modelMembrHMass);
fprintf("  => water mass ratio: %.4f\n",waterMassRatioH);

surfWaterBorderIndices = densDistribWaterH < ...
    surfWaterDensBorder * mean(densDistribWaterH(end/2-100:end/2+100));
surfWaterBorders = find(diff(surfWaterBorderIndices) ~= 0) ...
    + [borderShiftOffset, -borderShiftOffset];
surfWaterBorderIndices(1:surfWaterBorders(1)) = 1;
surfWaterBorderIndices(surfWaterBorders(2):end) = 1;
if length(surfWaterBorders) ~= 2
    error("The borders for surface water are not chosen correctly.");
end

densSubplot = initializeSubplot(densFigure,2,2,3:4);
densSurfWaterH = nan(1,length(densDistribWaterH));
densSurfWaterH(surfWaterBorderIndices) = densDistribWaterH( ...
    surfWaterBorderIndices);

densPlt(1) = plot(locations,densDistribWaterH,'Color',colorArray{1});
densPlt(2) = plot(locations,densSurfWaterH,'Color',colorArray{2});
densPlt(3) = plot(locations,densDistribMembrH,'Color',colorArray{3});
densPlt(4) = plotVerticalLine(locations(surfWaterBorders(1)) ...
    ,max(densDistribWaterH));
densPlt(5) = plotVerticalLine(locations(surfWaterBorders(2)) ...
    ,max(densDistribWaterH));
legend("Free H$_{2}$O","Surf. H$_{2}$O","Membrane",'Location','south');
xlabel('x position $[m]$');
ylabel('$\rho$ $[kg/m^3]$');
title("$\textbf{c)}$");

%% plot density distribtuion
densityFigure = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
densSurfWaterH = nan(1,length(densDistribWaterH));
densSurfWaterH(surfWaterBorderIndices) = densDistribWaterH( ...
    surfWaterBorderIndices);

densPlt(1) = plot(locations,densDistribWaterH,'Color',colorArray{1});
densPlt(2) = plot(locations,densSurfWaterH,'Color',colorArray{2});
densPlt(3) = plot(locations,densDistribMembrH,'Color',colorArray{3});
densPlt(4) = plotVerticalLine(locations(surfWaterBorders(1)) ...
    ,max(densDistribWaterH));
densPlt(5) = plotVerticalLine(locations(surfWaterBorders(2)) ...
    ,max(densDistribWaterH));
legend("Free H$_{2}$O","Surf. H$_{2}$O","Membrane",'Location','south');
xlabel('x position $[m]$');
ylabel('$\rho$ $[kg/m^3]$');

%% plot cross and auto relaxation rates

crossFig = initializeFigure('axisFontSize',axisFontSize ...
    ,'legendFontSize',legendFontSize);
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

r1Auto_SL_ofInterest = r1Auto_SM(indicesOfInterest);
r1Auto_LW_ofInterest = r1Auto_MW(indicesOfInterest);

r1Cross_SL_ofInterest = r1Cross_SM(indicesOfInterest);
r1Cross_LW_ofInterest = r1Cross_MW(indicesOfInterest);

crossPlt(1) = plot(fieldStrengthsOfInterest,r1Auto_SL_ofInterest);
crossPlt(2) = plot(fieldStrengthsOfInterest,r1Auto_LW_ofInterest);

crossPlt(3) = plot(fieldStrengthsOfInterest,r1Cross_SL_ofInterest);
crossPlt(4) = plot(fieldStrengthsOfInterest,r1Cross_LW_ofInterest);

legend("$R_{1,SL}^{auto}$", "$R_{1,LW}^{auto}$", "$R_{1,SL}^{cross}$" ...
    , "$R_{1,LW}^{cross}$");
axis([ 0 inf 0 0.13])


%% saving

if saving
    
    set(0, 'CurrentFigure', literatureFigure);
    saveFigureTo(scriptImagesFolderName ...
        ,"ObservedR1Values","WithFit","WMandGM");
    
    set(0, 'CurrentFigure', literatureFigure2_0);
    saveFigureTo(scriptImagesFolderName ...
        ,"ObservedR1Values","WithFit","WMandGM_loglog");
    
    set(0, 'CurrentFigure', literatureFigure2_1);
    saveFigureTo(scriptImagesFolderName ...
        ,"ObservedR1Values","WithFit","WMandGM_ylog");
    
    set(0, 'CurrentFigure', literatureFigure2_2);
    saveFigureTo(scriptImagesFolderName ...
        ,"ObservedR1Values","WithFit","WMandGM_xlog");
    
    set(0, 'CurrentFigure', FEModelResultsPlt);
    saveFigureTo(scriptImagesFolderName ...
        ,"FEModelResultsForR1","SP","withFit");
    
    set(0, 'CurrentFigure', mdSimFigure);
    saveFigureTo(scriptImagesFolderName ...
        ,"SL","LW_SP","PredictedValuesAndFieldStrengthBehvior");
    
    set(0,'CurrentFigure',appendixFigure);
    saveFigureTo(scriptImagesFolderName ...
        ,"AppendixS1","LWR1","IncreasedAndDecreasedExchRates");
    
    set(0,'CurrentFigure',densFigure);
    saveFigureTo(scriptImagesFolderName ...
        ,"DensityDistribution","HAtoms","Monolayer");
    
    set(0, 'CurrentFigure', densityFigure);
    saveFigureTo(scriptImagesFolderName ...
        ,"DensityDistribution","H","Atoms");
    
    set(0, 'CurrentFigure', crossFig);
    saveFigureTo(scriptImagesFolderName ...
        ,"CrossAndAutoRelaxation","SL","LW");
    
end

%% save data of interest
if saving
    fileID = fopen(sprintf("%s%s_informationForPaper.txt" ...
        ,dataForPaperFolderName,datestr(now,'yyyymmdd')),'w');
    
    formatSpecification = "%s: %.4f * B_0 ^(-%.4f) sigma_a = %.4e " ...
        + "sigma_b = %.4e \n\n";
    fprintf(fileID,formatSpecification ...
        ,"R1_WM",wmFactor,wmExponent,sqrt(wmCovariance(1,1)) ...
        ,sqrt(wmCovariance(2,2)));
    
    fprintf(fileID,formatSpecification ...
        ,"R1_GM",gmFactor,gmExponent,sqrt(gmCovariance(1,1)) ...
        ,sqrt(gmCovariance(2,2)));
    
    writeR1ToTxtFile(fileID,"SL");
    writeR1ToTxtFile(fileID,"LW");
    writeR1ToTxtFile(fileID,"SPAvg_IEModel");
    writeR1ToTxtFile(fileID,"SPWM_IEModel");
    writeR1ToTxtFile(fileID,"SPGM_IEModel");
    
    fprintf(fileID,"%20s | %4s +/- stat.  + syst. - syst. error: \n","field strengths", "R1_SP");
    fprintf(fileID,"%20.3f | %4.3f +/- %.4f +%.4f -%.4f \n", [ ...
        wmFieldStrengthArray(spIndicesForWMFitting); ...
        r1_SP_AvgTissue_IEModel(spIndicesForWMFitting); ...
        r1_SP_StatisticalError(spIndicesForWMFitting); ...
        r1_SP_systematicTissueBasedError(spIndicesForWMFitting) ...
        + r1_SP_systematicPositiveExchangeRateBasedError( ...
        spIndicesForWMFitting); ...
        r1_SP_systematicTissueBasedError(spIndicesForWMFitting) ...
        + r1_SP_systematicNegativeExchangeRateBasedError( ...
        spIndicesForWMFitting)]);
    
    fclose(fileID);
    
end

%% functions
function plt = plotVerticalLine(xPos,yValue)
plt = plot([xPos xPos],[0 yValue],'--k','LineWidth',0.7);
end

function writeR1ToTxtFile(fileID,appendix)
formatSpecification = "%s: %.4f * B_0 ^(-%.4f) sigma_a = %.4e " ...
    + "sigma_b = %.4e \n\n";
fprintf(fileID,formatSpecification ...
    ,"R1_" + appendix,evalin('base',"factorR1_" + appendix) ...
    ,evalin('base',"exponentR1_" + appendix) ...
    ,evalin('base',"sqrt(covarianceR1_" + appendix + "(1,1))") ...
    ,evalin('base',"sqrt(covarianceR1_" + appendix + "(2,2))"));
end

function logicalIndices = getIndicesInFieldStrengthArry(fieldStrengthArray ...
    ,fieldStrengthsOfInterest)

for fieldStrengthNr = 1:length(fieldStrengthsOfInterest)
    indices(fieldStrengthNr) = find((fieldStrengthArray ...
        > fieldStrengthsOfInterest(fieldStrengthNr) - 0.000001) ...
        .* (fieldStrengthArray < fieldStrengthsOfInterest(fieldStrengthNr) ...
        + 0.000001)); %#ok<AGROW>
    if isempty(indices(fieldStrengthNr))
        error("Field strength %.4f not found." ...
            ,fieldStrengthsOfInterest(fieldStrengthNr));
    end
    
end

logicalIndices = false(1,length(fieldStrengthArray));
logicalIndices(indices) = true;

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


%%

% %% plotting ratio of R1_SP and R1_SL
%
% set(0,'CurrentFigure',additionalFigure);
% r1SPSLRatioSubPlot = initializeSubplot(additionalFigure,2,2,2);
% cla(r1SPSLRatioSubPlot);
% lgd = legend();
% lgd.Visible = 'off';
%
% plot(fieldStrengthAxis,(factorR1_SL*fieldStrengthAxis.^(-exponentR1_SL)) ...
%     ./(factorR1_SPAvg_IEModel*fieldStrengthAxis.^( ...
%     -exponentR1_SPAvg_IEModel)));
%
% xlabel("Field strength [T]");
% ylabel("${}^{R_{1,SL}}{\mskip -5mu/\mskip -3mu}_{R_{1,SP}}$ [a.u.]");
% title("$\textbf{b)}$");
% drawnow;


