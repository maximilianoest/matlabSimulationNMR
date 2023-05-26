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

%% r1 rates
r1RatesFileName = "20230427_SM_WM_SP_histCompartmentAndCrossR1_fitted";
r1Data = load(resultsDir + r1RatesFolder + r1RatesFileName);
ESIVR = r1Data.ESIVR;
fieldStrengthOfInterest = 0.35; % Tesla
simFieldStrengthIndex = ...
    r1Data.fieldStrengths > fieldStrengthOfInterest - 0.001 ...
    & r1Data.fieldStrengths < fieldStrengthOfInterest + 0.001;
litFieldStrengthIndex =  ...
    r1Data.litFieldStrengths > fieldStrengthOfInterest - 0.001 ...
    & r1Data.litFieldStrengths < fieldStrengthOfInterest + 0.001;

r1_SL = r1Data.r1Eff_SM(simFieldStrengthIndex);

r1_MW = r1Data.r1Hist_MW(simFieldStrengthIndex);
r1_free = r1Data.r1_free;

r1_wMWater = wmNonWaterFraction * ESIVR * r1_MW ...
    + wmWaterFraction*r1_free;

r1ObservedWm = r1Data.r1Observed_WM(litFieldStrengthIndex);

%% model assumptions
% according to Wang2020, the exchange rate is:
% "k_w and k_m can be equivalently described by f and k, where f is the MP
% pool fraction and k is the magnetization exchange rate as a fraction of
% the entire proton pool (both MP and WP)
wang2020FieldStrengths = [3 7]; % T
wang2020ExchangeRates_k = [1.5 1.38]; % normalized on whol MP and WP pool
wang2020PoolFractionFactors_f = mean([0.281 0.289]); % M_0,m / ( M_0,m + M_0,w )

% k_Slope = (wang2020ExchangeRates_k(2) - wang2020ExchangeRates_k(1)) ...
%     /(wang2020FieldStrengths(2) - wang2020FieldStrengths(1));
% k_Offset = wang2020ExchangeRates_k(1) - (k_Slope * wang2020FieldStrengths(1));
% fieldDependentExchangeRate =  k_Slope * fieldStrengthOfInterest + k_Offset;


[factor,exponent] = fitFieldStrengthBehaviorToPowerFunction( ...
    wang2020FieldStrengths,wang2020ExchangeRates_k,[1 1]);
fieldDependentExchangeRate = factor * fieldStrengthOfInterest^(exponent);

fieldDependentExchangeRate = mean(wang2020ExchangeRates_k);

exchange_macroToWater_normPoolSize = fieldDependentExchangeRate  ...
    /mean(wang2020PoolFractionFactors_f);
exchange_waterToMacro_normPoolSize = fieldDependentExchangeRate ...
    /mean(1-wang2020PoolFractionFactors_f);

exchange_lipidToWater = exchange_macroToWater_normPoolSize*wmLipidFraction;
exchange_waterToLipid = exchange_waterToMacro_normPoolSize ...
    *wmWaterFraction*constData.wmLipidContent;

exchange_proteinToWater = exchange_macroToWater_normPoolSize ...
    *wmProteinFraction;
exchange_waterToProtein = exchange_waterToMacro_normPoolSize ...
    *wmWaterFraction*constData.wmLipidContent;

%% set up
duration = 3; % seconds
deltaT = 0.001; % seconds
timeStepCount = duration/deltaT;

initialWaterMagnetization = -1; % IR experiment
equilibriumWaterMagnetization = 1;

initialSLMagnetization = 0;
equilibriumSLMagnetization = 1;

initialSPMagnetization = 0;
equilibriumSPMagnetization = 1;

r1_SPArray =  linspace(0,2*r1_SL,200);

exchangeRates = linspace(0.5*mean(wang2020ExchangeRates_k) ...
    ,3*mean(wang2020ExchangeRates_k),50);

%% finding the right R1_SP

[r1Wm_waterExchanging,~,~,~] = determineWaterR1(r1_SPArray,r1_SL ...
    ,r1_wMWater,timeStepCount,deltaT ...
    ,initialSLMagnetization,initialSPMagnetization ...
    ,initialWaterMagnetization,exchange_waterToLipid ...
    ,exchange_waterToProtein,exchange_lipidToWater ...
    ,exchange_proteinToWater);

timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;


%% find best estimate for r1_waterExchanging
difference = abs(r1Wm_waterExchanging - r1ObservedWm);
[~,bestGuessIndex] = min(difference);

[~,bestGuessWaterMagnetization,bestGuessSLMagnetization ...
    ,bestGuessSPMagnetization] = determineWaterR1( ...
    r1_SPArray(bestGuessIndex),r1_SL,r1_wMWater,timeStepCount,deltaT ...
    ,initialSLMagnetization,initialSPMagnetization ...
    ,initialWaterMagnetization,exchange_waterToLipid ...
    ,exchange_waterToProtein,exchange_lipidToWater ...
    ,exchange_proteinToWater);

%% find best exchange rate

r1_SP = r1Data.r1_SP_GM(simFieldStrengthIndex);
[r1Wm_waterExchangingBasedOnK,~,~,~] ...
    = determineWaterR1BasedOnFieldDependentExchangeRate( ...
    r1_SL,initialSLMagnetization,r1_SP,initialSPMagnetization ...
    ,r1_wMWater,initialWaterMagnetization,timeStepCount,deltaT ...
    ,exchangeRates,wang2020PoolFractionFactors_f,constData);

difference = abs(r1Wm_waterExchangingBasedOnK - r1ObservedWm);
[~,bestGuessIndex_k] = min(difference);

bestWorkingExchangeRate = exchangeRates(bestGuessIndex_k);
[~,bestGuessWaterMagnetization_k,bestGuessSLMagnetization_k ...
    ,bestGuessSPMagnetization_k] ...
    = determineWaterR1BasedOnFieldDependentExchangeRate( ...
    r1_SL,initialSLMagnetization,r1_SP,initialSPMagnetization ...
    ,r1_wMWater,initialWaterMagnetization,timeStepCount,deltaT ...
    ,bestWorkingExchangeRate,wang2020PoolFractionFactors_f,constData);

%% plotting to command window
fprintf("\n\n<strong>Solid protein R1 fit </strong> \n");
fprintf(" Literature tissue R1: %.4f \n",r1ObservedWm);
fprintf(" Determined water R1: %.4f \n",r1_SPArray(bestGuessIndex));

fprintf("<strong>Exchange rate fit </strong> \n");
fprintf(" Determined exchange rate: %.4f \n",bestWorkingExchangeRate);
fprintf(" Used R1_SP: %.4f \n",r1_SP);

%% plotting
fig = initializeFigure();
initializeSubplot(fig,2,2,3);
plot(r1_SPArray, r1Wm_waterExchanging)
xlabel("R$_{1,SP}$");
ylabel("Predicted R$_{1,w}$");
lgd = legend();
lgd.Visible = 'off';

initializeSubplot(fig,2,2,1);
plot(timeAxis,1-2*exp(-timeAxis*r1Wm_waterExchanging(bestGuessIndex)));
plot(timeAxis,1-2*exp(-timeAxis*r1ObservedWm));
plot(timeAxis,bestGuessWaterMagnetization);
plot(timeAxis,bestGuessSLMagnetization);
plot(timeAxis,bestGuessSPMagnetization);
plot(timeAxis,1-exp(-timeAxis*r1_SPArray(bestGuessIndex)));
plot(timeAxis,1-exp(-timeAxis*r1_SL));
legend("R$_1^{fit}$","R$_1^{obs.}$","water magnetization" ...
    ,"lipid magnetization","protein magnetization" ...
    ,"Best guess R$_{1,SP}$","R$_{1,SL}$",'location','southeast');
xlabel("Time [sec]");
ylabel("Magnetization [a.u.]");

initializeSubplot(fig,2,2,2);
plot(timeAxis,1-2*exp(-timeAxis*r1Wm_waterExchangingBasedOnK( ...
    bestGuessIndex_k)));
plot(timeAxis,1-2*exp(-timeAxis*r1ObservedWm));
plot(timeAxis,bestGuessWaterMagnetization_k);
plot(timeAxis,bestGuessSLMagnetization_k);
plot(timeAxis,bestGuessSPMagnetization_k);
plot(timeAxis,1-exp(-timeAxis*r1_SP));
plot(timeAxis,1-exp(-timeAxis*r1_SL));
legend("R$_1^{fit}$","R$_1^{obs.}$","water magnetization" ...
    ,"lipid magnetization","protein magnetization" ...
    ,"R$_{1,SP}$","R$_{1,SL}$",'location','southeast');
xlabel("Time [sec]");
ylabel("Magnetization [a.u.]");

initializeSubplot(fig,2,2,4);
writeIntoPlotWindow(fig,[2,2,4],sprintf("Water R$_1$ = %.4f Hz \n" ...
    + "Water T$_1$ = %.4f s \n\nLiterature: %.4f Hz / %.4f s \n\n" ...
    + "R$_{1,SP}$ = %.4f \nR$_{1,SL}$ = %.4f \n\n" ...
    + "Exchange rate: %.4f \n\n" ...
    ,r1Wm_waterExchanging(bestGuessIndex) ...
    ,1/r1Wm_waterExchanging(bestGuessIndex),r1ObservedWm ...
    ,1/r1ObservedWm,r1_SPArray(bestGuessIndex),r1_SL ...
    ,bestWorkingExchangeRate));
%% functions

function [r1_exchangingWater,waterMagnetization,slMagnetization ...
    ,spMagnetization] = determineWaterR1(r1_SPArray,r1_SL ...
    ,r1_water,timeStepCount,deltaT ...
    ,initialSLMagnetization,initialSPMagnetization ...
    ,initialWaterMagnetization,exchange_waterToLipid ...
    ,exchange_waterToProtein,exchange_lipidToWater ...
    ,exchange_proteinToWater)

% preallocation
equilibriumSLMagnetization = 1;
slMagnetization = zeros(1,timeStepCount);

equilibriumSPMagnetization = 1;
spMagnetization = zeros(1,timeStepCount);

equilibriumWaterMagnetization = 1;
waterMagnetization = zeros(1,timeStepCount);

r1_exchangingWater = zeros(1,length(r1_SPArray));

timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;

for r1_SPCounter = 1:length(r1_SPArray)
    r1_SP = r1_SPArray(r1_SPCounter);
    slMagnetization(1) = initialSLMagnetization;
    spMagnetization(1) = initialSPMagnetization;
    waterMagnetization(1) = initialWaterMagnetization;
    
    oldStep = 1;
    % solving Bloch equations numerically
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
        oldStep = newStep;
    end
    
    opts = optimset('Display','off');
    
    relaxationFunction = @(relaxationRate,timeAxis)( ...
        1 - 2*exp(-relaxationRate*timeAxis));
    
    r1_exchangingWater(r1_SPCounter) = lsqcurvefit(relaxationFunction ...
        ,1,timeAxis,waterMagnetization,[],[],opts);
    
    if mod(r1_SPCounter,100) == 0
        fprintf("r1_SP = %.4f \n",r1_SP);
    end
end


        
end

function [r1_exchangingWater,waterMagnetization,slMagnetization ...
    ,spMagnetization] ...
    = determineWaterR1BasedOnFieldDependentExchangeRate(...
    r1_SL,initialSLMagnetization,r1_SP,initialSPMagnetization ...
    ,r1_water,initialWaterMagnetization,timeStepCount,deltaT ...
    ,exchangeRates,macroMolecularPoolFraction,constData)


wmWaterFraction = constData.wmWaterContent;
wmNonWaterFraction = constData.wmNonWaterContent;
wmLipidFraction = wmNonWaterFraction * constData.wmLipidContent;
wmProteinFraction = wmNonWaterFraction * constData.wmProteinContent;


equilibriumSLMagnetization = 1;
slMagnetization = zeros(1,timeStepCount);

equilibriumSPMagnetization = 1;
spMagnetization = zeros(1,timeStepCount);

equilibriumWaterMagnetization = 1;
waterMagnetization = zeros(1,timeStepCount);

r1_exchangingWater = zeros(1,length(exchangeRates));

timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;

for exchangeRateCounter = 1:length(exchangeRates)
    rate_macroToWater = exchangeRates(exchangeRateCounter) ...
        /macroMolecularPoolFraction;
    rate_waterToMacro = exchangeRates( ...
        exchangeRateCounter)/(1 - macroMolecularPoolFraction);
    
    exchange_lipidToWater = rate_macroToWater*wmLipidFraction;
    exchange_waterToLipid = rate_waterToMacro ...
        *wmWaterFraction*constData.wmLipidContent;
    
    exchange_proteinToWater = rate_macroToWater*wmProteinFraction;
    exchange_waterToProtein = rate_waterToMacro ...
        *wmWaterFraction*constData.wmLipidContent;
    
    slMagnetization(1) = initialSLMagnetization;
    spMagnetization(1) = initialSPMagnetization;
    waterMagnetization(1) = initialWaterMagnetization;
    
    oldStep = 1;
    % solving Bloch equations numerically
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
        oldStep = newStep;
    end
    
    opts = optimset('Display','off');
    
    relaxationFunction = @(relaxationRate,timeAxis)( ...
        1 - 2*exp(-relaxationRate*timeAxis));
    
    r1_exchangingWater(exchangeRateCounter) = lsqcurvefit( ...
        relaxationFunction,1,timeAxis,waterMagnetization,[],[],opts);
    
    if mod(exchangeRateCounter,100) == 0
        fprintf("Exchange Rate = %.4f \n",r1_SP);
    end
end




end


function [factor,exponent] = fitFieldStrengthBehaviorToPowerFunction( ...
   fieldStrengths,exchangeRates,startValues,varargin)

opts = optimset('Display','off');

powerFunction = @(parameters,fieldStrengths)(parameters(1) ...
        *fieldStrengths.^(-parameters(2)));

if isempty(varargin)

    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,exchangeRates,[],[],opts);
    
elseif length(varargin{:}) == 2
    
    lowerLimitFactor = 0;
    upperLimitFactor = inf;
    
    varargin = sort(varargin{:});
    lowerLimitExponent = varargin(1);
    upperLimitExponent = varargin(2);
    
    lowerLimits = [lowerLimitFactor lowerLimitExponent];
    upperLimits = [upperLimitFactor upperLimitExponent];
    
    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,exchangeRates,lowerLimits,upperLimits,opts);
    
else
    error("Not implemented yet.");
    
end

factor = parameters(1);
    exponent = parameters(2);

end









