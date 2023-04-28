clc; clear all; close all; fclose('all');

%% set dependencies
addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');


%% dependencies
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
constitutionFolder = "/wholeMyelin_brainConstitution/";

%% constitution data
constitutionFileName = "20230426_wmAndGMCompositionBasedOnLiterature";
constData = load(resultsDir + constitutionFolder + constitutionFileName);

wmWaterFraction = constData.wmWaterContent;
wmNonWaterFraction = constData.wmNonWaterContent;
wmLipidFraction = wmNonWaterFraction * constData.wmLipidContent;
wmProteinFraction = wmNonWaterFraction * constData.wmProteinContent;

%% set up
timeStepCount = 10000;
deltaT = 0.001;

initialWaterMagnetization = -1; % IR experiment
equilibriumWaterMagnetization = 1;

initialSLMagnetization = 0;
equilibriumSLMagnetization = 1;

initialSPMagnetization = 0;
equilibriumSPMagnetization = 1;

%% model assumptions
% according to Wang2020, the exchange rate is:
% "k_w and k_m can be equivalently described by f and k, where f is the MP
% pool fraction and k is the magnetization exchange rate as a fraction of
% the entire proton pool (both MP and WP)
wang2020FieldStrengths = [3 7]; % T
wang2020ExchangeRates_k = [1.5 1.38]; % normalized on whol MP and WP pool
wang2020PoolFractionFactors_f = [0.281 0.289]; % M_0,m / ( M_0,m + M_0,w )

exchange_macroToWater_normPoolSize = mean(wang2020ExchangeRates_k)  ...
    /mean(wang2020PoolFractionFactors_f);
exchange_waterToMacro_normPoolSize = mean(wang2020ExchangeRates_k)  ...
    /mean(1-wang2020PoolFractionFactors_f);

exchange_lipidToWater = exchange_macroToWater_normPoolSize*wmLipidFraction;
exchange_waterToLipid = exchange_waterToMacro_normPoolSize*wmWaterFraction;

exchange_proteinToWater = exchange_macroToWater_normPoolSize*wmProteinFraction;
exchange_waterToProtein = exchange_waterToMacro_normPoolSize*wmWaterFraction;


%% r1 rates
r1RatesFileName = "20230427_SM_WM_SP_histCompartmentAndCrossR1_fitted";
r1Data = load(resultsDir + r1RatesFolder + r1RatesFileName);
ESIVR = r1Data.ESIVR;

r1_SL = r1Data.r1Eff_SM(r1Data.fieldStrengths > 2.999 ...
    & r1Data.fieldStrengths < 3.0001);

r1_MW = r1Data.r1Hist_MW(r1Data.fieldStrengths > 2.999 ...
    & r1Data.fieldStrengths < 3.0001);
r1_free = r1Data.r1_free;

r1_wMWater = wmNonWaterFraction * ESIVR * r1_MW ...
    + wmWaterFraction*r1_free;

r1ObservedWm = r1Data.r1Observed_WM(r1Data.litFieldStrengths > 2.9999 ...
    & r1Data.litFieldStrengths < 3.0001);


%% finding the right R1_SP
r1_SPArray =  0.1*r1_SL : 0.05 : 5 * r1_SL;
r1Wm_waterExchanging = zeros(1,length(r1_SPArray));

for r1_SPCounter = 1:length(r1_SPArray)
    r1_SP = r1_SPArray(r1_SPCounter);
    r1Wm_waterExchanging(r1_SPCounter) = determineWaterR1(r1_SP,r1_SL ...
        ,r1_wMWater,timeStepCount,deltaT ...
        ,initialSLMagnetization,initialSPMagnetization ...
        ,initialWaterMagnetization,exchange_waterToLipid ...
        ,exchange_waterToProtein,exchange_lipidToWater ...
        ,exchange_proteinToWater);
    
    if mod(r1_SPCounter,100) == 0
        fprintf("r1_SP = %.4f \n",r1_SP);
    end
    
end

%% find best estimate for r1_waterExchanging
difference = abs(r1Wm_waterExchanging - r1ObservedWm);
[~,bestGuessIndex] = min(difference);



%% plotting
fig = initializeFigure();
initializeSubplot(fig,2,1,1);
plot(r1_SPArray, r1Wm_waterExchanging)
xlabel("R$_{1,SP}$");
ylabel("Predicted R$_{1,w}$");
lgd = legend();
lgd.Visible = 'off';

initializeSubplot(fig,2,1,2);
writeIntoPlotWindow(fig,[2,1,2],sprintf("Water R$_1$ = %.4f Hz \n" ...
    + "Water T$_1$ = %.4f s \n\nLiterature: %.4f Hz / %.4f s" ...
    ,r1Wm_waterExchanging(bestGuessIndex) ...
    ,1/r1Wm_waterExchanging(bestGuessIndex),r1ObservedWm ...
    ,1/r1ObservedWm));
%% functions

function r1_exchangingWater = determineWaterR1(r1_SP,r1_SL ...
        ,r1_water,timeStepCount,deltaT ...
        ,initialSLMagnetization,initialSPMagnetization ...
        ,initialWaterMagnetization,exchange_waterToLipid ...
        ,exchange_waterToProtein,exchange_lipidToWater ...
        ,exchange_proteinToWater)

% preallocation
slMagnetization = zeros(1,timeStepCount);
slMagnetization(1) = initialSLMagnetization;
equilibriumSLMagnetization = 1;

spMagnetization = zeros(1,timeStepCount);
spMagnetization(1) = initialSPMagnetization(1);
equilibriumSPMagnetization = 1;

waterMagnetization = zeros(1,timeStepCount);
waterMagnetization(1) = initialWaterMagnetization;
equilibriumWaterMagnetization = 1;

% simulation
oldStep = 1;
timeAxis = zeros(1,timeStepCount);
for newStep = 2:timeStepCount
    slMagnetization(newStep) = slMagnetization(oldStep) + deltaT * ( ...
        (equilibriumSLMagnetization - slMagnetization(oldStep))*r1_SL ...
        - exchange_lipidToWater*slMagnetization(oldStep) ...
        + exchange_waterToLipid*waterMagnetization(oldStep));
    
    spMagnetization(newStep) = spMagnetization(oldStep) + deltaT * ( ...
        (equilibriumSPMagnetization - spMagnetization(oldStep))*r1_SP ...
        - exchange_proteinToWater*spMagnetization(oldStep) ...
        + exchange_waterToProtein*waterMagnetization(oldStep));
    
    waterMagnetization(newStep) = waterMagnetization(oldStep) + deltaT *( ...
        (equilibriumWaterMagnetization - waterMagnetization(oldStep))*r1_water ...
        - exchange_waterToLipid*waterMagnetization(oldStep) ...
        - exchange_waterToProtein*waterMagnetization(oldStep) ...
        + exchange_lipidToWater*slMagnetization(oldStep) ...
        + exchange_proteinToWater*spMagnetization(oldStep));
    
    timeAxis(newStep) = timeAxis(oldStep) + deltaT;
    
    oldStep = newStep;
end

opts = optimset('Display','off');

relaxationFunction = @(relaxationRate,timeAxis)( ...
    1 - 2*exp(-relaxationRate*timeAxis));

r1_exchangingWater = lsqcurvefit(relaxationFunction,1,timeAxis ...
    ,waterMagnetization,[],[],opts);
        
end








