clc; clear all; close all; fclose('all');

addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

%% set up system
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
constitutionFolder = "/wholeMyelin_brainConstitution/";
saving = 1;

% ==== constitution data ====
constitutionFileName = "20230426_wmAndGMCompositionBasedOnLiterature";
constData = load(resultsDir + constitutionFolder + constitutionFileName);

wmWaterFraction = constData.wmWaterContent;
wmNonWaterFraction = constData.wmNonWaterContent;
wmLipidFraction = wmNonWaterFraction * constData.wmLipidContent;
wmProteinFraction = wmNonWaterFraction * constData.wmProteinContent;

gmWaterFraction = constData.gmWaterContent;
gmNonWaterFraction = constData.gmNonWaterContent;
gmLipidFraction = gmNonWaterFraction * constData.gmLipidContent;
gmProteinFraction = gmNonWaterFraction * constData.gmProteinContent;

% ==== exchange parameters based on Wang 2020 ====
% according to Wang2020, the exchange rate is:
% "k_w and k_m can be equivalently described by f and k, where f is the MP
% pool fraction and k is the magnetization exchange rate as a fraction of
% the entire proton pool (both MP and WP)
wang2020FieldStrengths = [3 7]; % T
wang2020ExchangeRates_k = [1.5 1.38]; % normalized on whol MP and WP pool
wang2020PoolFractionFactors_f = mean([0.281 0.289]); % M_0,m / ( M_0,m + M_0,w )

% ==== relaxation rates ====
r1RatesFileName = "20230508_SM_WM_SP_histCompartmentAndCrossR1_fitted";
r1Data = load(resultsDir + r1RatesFolder + r1RatesFileName);
ESIVR = r1Data.ESIVR;
wmSurfaceWaterFraction = wmNonWaterFraction * ESIVR;
gmSurfaceWaterFraction = gmNonWaterFraction * ESIVR;
% -> field strengths to fit
fieldStrengthsToInvestigate = r1Data.litFieldStrengths;
r1Wm = r1Data.r1Observed_WM;
r1Gm = r1Data.r1Observed_GM;
r1_free = r1Data.r1_free;

% ==== parameter set up ====
% -> time set up
durationForEachFieldStrength = ones(1,5)*15; % seconds
deltaTForEachFieldStrength = ones(1,5)*5e-3; % seconds


% -> fitting configuration
numberOfPointsToFit = 1000;
highestR1_SPFactor = 1.5;
lowestR1_SPFactor = 0.5;

highestFieldStrengthForPlotting = 7.1;


% -> experiment set up
initialWaterMagnetization = -1; % IR experiment
equilibriumWaterMagnetization = 1;

initialSLMagnetization = 0;
equilibriumSLMagnetization = 1;

initialSPMagnetization = 0;
equilibriumSPMagnetization = 1;

%% testing exchange rates

deltaT = deltaTForEachFieldStrength(1);
duration = durationForEachFieldStrength(1);
timeStepCount = duration/deltaT;
timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;

exchange_waterToSolid_Wm = mean(wang2020ExchangeRates_k) ...
    /(wmWaterFraction);
exchange_lipidToWater_Wm = mean(wang2020ExchangeRates_k) ...
    /(wmNonWaterFraction * constData.wmLipidContent);
exchange_proteinToWater_Wm = mean(wang2020ExchangeRates_k) ...
    /(wmNonWaterFraction * constData.wmProteinContent);

[waterMagnetization_wm,slMagnetization_wm,spMagnetization_wm] ...
    = determineMagnetizationsWhenOnlyExchanging(exchange_waterToSolid_Wm ...
    ,exchange_lipidToWater_Wm,exchange_proteinToWater_Wm,timeStepCount ...
    ,deltaT);

exchange_waterToSolid_Gm = mean(wang2020ExchangeRates_k) ...
    /(gmWaterFraction);
exchange_lipidToWater_Gm = mean(wang2020ExchangeRates_k) ...
    /(gmNonWaterFraction * constData.gmLipidContent);
exchange_proteinToWater_Gm = mean(wang2020ExchangeRates_k) ...
    /(gmNonWaterFraction * constData.gmProteinContent);

[waterMagnetization_gm,slMagnetization_gm,spMagnetization_gm] ...
    = determineMagnetizationsWhenOnlyExchanging(exchange_waterToSolid_Gm ...
    ,exchange_lipidToWater_Gm,exchange_proteinToWater_Gm,timeStepCount ...
    ,deltaT);

initializeFigure();
plot(timeAxis,waterMagnetization_wm);
plot(timeAxis,slMagnetization_wm);
plot(timeAxis,spMagnetization_wm);

plot(timeAxis,waterMagnetization_gm,'--');
plot(timeAxis,slMagnetization_gm,'--');
plot(timeAxis,spMagnetization_gm,'--');

legend("Water WM","SL WM","SP WM","Water GM","SL GM","SP GM",'Location' ...
    ,'east');
xlabel("Time [sec]");
ylabel("Magnetization [a.u.]");

%% find best working R1_SP

bestGuessTissueR1_wm = zeros(1,length(fieldStrengthsToInvestigate));
bestGuessR1_SP_wm = zeros(1,length(fieldStrengthsToInvestigate));

bestGuessTissueR1_gm = zeros(1,length(fieldStrengthsToInvestigate));
bestGuessR1_SP_gm = zeros(1,length(fieldStrengthsToInvestigate));

for fieldStrengthCounter = 1:length(fieldStrengthsToInvestigate)
    fieldStrength = fieldStrengthsToInvestigate(fieldStrengthCounter);
    fprintf("<strong>Field strength %.2f T </strong> \n",fieldStrength);
    fieldStrengthIndex_sim = find(r1Data.fieldStrengths ...
        > (fieldStrength - 0.001) & r1Data.fieldStrengths ...
        < (fieldStrength + 0.001));
    r1_SL = r1Data.r1Eff_SM(fieldStrengthIndex_sim);
    fprintf(" R1_SL: %.4f\n", r1_SL);
    r1_SPArray = linspace(lowestR1_SPFactor*r1_SL ...
        ,highestR1_SPFactor*r1_SL,numberOfPointsToFit);
    r1_MW = r1Data.r1Hist_MW(fieldStrengthIndex_sim);
    
    deltaT = deltaTForEachFieldStrength(fieldStrengthCounter);
    duration = durationForEachFieldStrength(fieldStrengthCounter);
    timeStepCount = round(duration/deltaT);
    timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;
    
    % --- white matter
    r1_wmWater = wmSurfaceWaterFraction/wmWaterFraction * r1_MW ...
        + (1-wmSurfaceWaterFraction)/wmWaterFraction*r1_free;
    
    [waterR1_wm,~] = findWaterR1ForGivenR1_SP(r1_SPArray,r1_SL,r1_wmWater ...
        ,wmWaterFraction,wmLipidFraction,wmProteinFraction ...
        ,exchange_waterToSolid_Wm,exchange_lipidToWater_Wm ...
        ,exchange_proteinToWater_Wm,timeStepCount,deltaT);
    
    [~,bestGuessIndex] = min(abs(waterR1_wm - r1Wm(fieldStrengthCounter)));
    bestGuessTissueR1_wm(fieldStrengthCounter) = waterR1_wm( ...
        bestGuessIndex);
    bestGuessR1_SP_wm(fieldStrengthCounter) = r1_SPArray(bestGuessIndex);
    [~,bestGuessWaterMagnetization_wm(fieldStrengthCounter,:)] ...
        = findWaterR1ForGivenR1_SP( ...
        bestGuessR1_SP_wm(fieldStrengthCounter),r1_SL,r1_wmWater ...
        ,wmWaterFraction,wmLipidFraction,wmProteinFraction ...
        ,exchange_waterToSolid_Wm,exchange_lipidToWater_Wm ...
        ,exchange_proteinToWater_Wm,timeStepCount,deltaT); %#ok<SAGROW>
    
    % --- gray matter 
    r1_gmWater = gmSurfaceWaterFraction/gmWaterFraction * r1_MW ...
        + (1-gmSurfaceWaterFraction)/gmWaterFraction*r1_free;
    
    [waterR1_gm,~] = findWaterR1ForGivenR1_SP(r1_SPArray,r1_SL,r1_gmWater ...
        ,gmWaterFraction,gmLipidFraction,gmProteinFraction ...
        ,exchange_waterToSolid_Gm,exchange_lipidToWater_Gm ...
        ,exchange_proteinToWater_Gm,timeStepCount,deltaT);
    
    [~,bestGuessIndex] = min(abs(waterR1_gm - r1Gm(fieldStrengthCounter)));
    bestGuessTissueR1_gm(fieldStrengthCounter) = waterR1_gm( ...
        bestGuessIndex);
    bestGuessR1_SP_gm(fieldStrengthCounter) = r1_SPArray(bestGuessIndex);
    [~,bestGuessWaterMagnetization_gm(fieldStrengthCounter,:)] ....
        = findWaterR1ForGivenR1_SP( ...
        bestGuessR1_SP_gm(fieldStrengthCounter),r1_SL,r1_gmWater ...
        ,gmWaterFraction,gmLipidFraction,gmProteinFraction ...
        ,exchange_waterToSolid_Gm,exchange_lipidToWater_Gm ...
        ,exchange_proteinToWater_Gm,timeStepCount,deltaT); %#ok<SAGROW>
    
    fprintf(" Tissue R1 in literature (WM/GM): %.4f / %.4f \n" ...
        + " Tissue R1 determined here (WM/GM): %.4f / %.4f \n" ...
        + " R1_SP based on WM/GM: %.4f / %.4f\n" ...
        ,r1Wm(fieldStrengthCounter),r1Gm(fieldStrengthCounter) ...
        ,bestGuessTissueR1_wm(fieldStrengthCounter) ...
        ,bestGuessTissueR1_gm(fieldStrengthCounter) ...
        ,bestGuessR1_SP_wm(fieldStrengthCounter) ...
        ,bestGuessR1_SP_gm(fieldStrengthCounter));
    
end

%% fit R1_SP to power function

[factorWm,exponentWm] = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    fieldStrengthsToInvestigate,bestGuessR1_SP_wm,[6 1]);
[factorGm,exponentGm] = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    fieldStrengthsToInvestigate,bestGuessR1_SP_gm,[6 1]);

fieldStrengths = r1Data.fieldStrengths(r1Data.fieldStrengths ...
    < highestFieldStrengthForPlotting);

fig = initializeFigure(); %#ok<NASGU>
plot(fieldStrengthsToInvestigate,bestGuessR1_SP_wm,'*');
plot(fieldStrengths,factorWm*fieldStrengths.^(-exponentWm),'LineWidth',1);

plot(fieldStrengthsToInvestigate,bestGuessR1_SP_gm,'*');
plot(fieldStrengths,factorGm*fieldStrengths.^(-exponentGm),'LineWidth',1);

plot(fieldStrengthsToInvestigate,r1Data.r1Prot_WM,'o');
plot(fieldStrengths,r1Data.spFitParamsWM_AllFieldStrengths(1)*fieldStrengths.^( ...
    -r1Data.spFitParamsWM_AllFieldStrengths(2)),'--','LineWidth',1);

plot(fieldStrengthsToInvestigate,r1Data.r1Prot_GM,'o');
plot(fieldStrengths,r1Data.spFitParamsGM_AllFieldStrengths(1)*fieldStrengths.^( ...
    -r1Data.spFitParamsGM_AllFieldStrengths(2)),'--','LineWidth',1);

axis([0 inf 0 25])
legend("WM$_{IE}$", "WM$_{IE}$ fitted", "GM$_{IE}$", "GM$_{IE}$ fitted" ...
    ,"WM$_{FE}$","WM$_{FE} fitted$","GM$_{FE}$","GM$_{FE}$ fitted");
xlabel("Field strength [T]");
ylabel("R$_{1,SP}$ [Hz]");

if saving
    saveFigureTo(resultsDir + r1RatesFolder,'solidProteinR1' ...
        ,datestr(now,'yyyymmdd'),'FEvsNotFE');
end

%% magnetization part for checking
fig = initializeFigure();
legendEntries = {};
for fieldStrengthCounter = 1:length(fieldStrengthsToInvestigate)
   fieldStrength = fieldStrengthsToInvestigate(fieldStrengthCounter);
   
   deltaT = deltaTForEachFieldStrength(fieldStrengthCounter);
   duration = durationForEachFieldStrength(fieldStrengthCounter);
   timeStepCount = round(duration/deltaT);
   timeAxis = 0:deltaT:(timeStepCount-1)*deltaT;
   
   
   initializeSubplot(fig,2,1,1);
   wmBlochPlot(fieldStrengthCounter) = plot(timeAxis ...
       ,bestGuessWaterMagnetization_wm(fieldStrengthCounter,:)); %#ok<SAGROW>
   wmFitPlot = plot(timeAxis,relaxationFunction(bestGuessTissueR1_wm( ...
       fieldStrengthCounter),{timeAxis,wmWaterFraction}));
   wmFitPlot.Color = wmBlochPlot(fieldStrengthCounter).Color;
   wmFitPlot.LineStyle = '--';
   
   
   initializeSubplot(fig,2,1,2);
   gmBlochPlot = plot(timeAxis,bestGuessWaterMagnetization_gm( ...
       fieldStrengthCounter,:));
   gmFitPlot = plot(timeAxis,relaxationFunction(bestGuessTissueR1_gm( ...
       fieldStrengthCounter),{timeAxis,gmWaterFraction}));
   gmFitPlot.Color = gmBlochPlot.Color;
   gmFitPlot.LineStyle = '--';

   legendEntries{end+1} = sprintf("%.2f T",fieldStrength); %#ok<SAGROW>

end



initializeSubplot(fig,2,1,1);
title("White matter water ( -- num. Bloch equation, - - Fitted R$_1$)");
ylabel("Magnetization [a.u.]");
lgd1 = legend(wmBlochPlot,legendEntries,'Location','southeast');
axis([0 10 0 1]);

initializeSubplot(fig,2,1,2);
title("Gray matter water");
ylabel("Magnetization [a.u.]");
xlabel("Time [sec]");
lgd2 = legend();
lgd2.Visible = 'off';
axis([0 10 0 1]);

fprintf("<strong>Fitting results </strong>\n");
fprintf(" IE Model: \n -> WM: %.4f * B_0 ^( - %.4f) \n" ...
    + " -> GM: %.4f * B_0 ^( - %.4f) \n",factorWm,exponentWm,factorGm ...
    ,exponentGm);
fprintf(" FE Model: \n -> WM: %.4f * B_0 ^( - %.4f) \n" ...
    + " -> GM: %.4f * B_0 ^( - %.4f) \n" ...
    ,r1Data.spFitParamsWM_AllFieldStrengths(1) ...
    ,r1Data.spFitParamsWM_AllFieldStrengths(2) ...
    ,r1Data.spFitParamsGM_AllFieldStrengths(1) ...
    ,r1Data.spFitParamsGM_AllFieldStrengths(2));




%% functions

function [waterMagnetization,slMagnetization,spMagnetization] ...
    = determineMagnetizationsWhenOnlyExchanging(exchangeRate_waterToMacro ...
    ,exchangeRate_lipidToWater,exchangeRate_proteinToWater ...
    ,timeStepCount,deltaT)

slMagnetization = zeros(1,timeStepCount);
spMagnetization = zeros(1,timeStepCount);
waterMagnetization = zeros(1,timeStepCount);

slMagnetization(1) = 0;
spMagnetization(1) = 0;
waterMagnetization(1) = 1;

oldStep = 1;
% solving Bloch equations numerically
for newStep = 2:timeStepCount
    slMagnetization(newStep) = slMagnetization(oldStep) + deltaT * ( ...
        - exchangeRate_lipidToWater*slMagnetization(oldStep) ...
        + exchangeRate_waterToMacro*waterMagnetization(oldStep));
    
    spMagnetization(newStep) = spMagnetization(oldStep) + deltaT * ( ...
        - exchangeRate_proteinToWater*spMagnetization(oldStep) ...
        + exchangeRate_waterToMacro*waterMagnetization(oldStep));
    
    waterMagnetization(newStep) = waterMagnetization(oldStep) + deltaT *( ...
        - exchangeRate_waterToMacro*waterMagnetization(oldStep) ...
        - exchangeRate_waterToMacro*waterMagnetization(oldStep) ...
        + exchangeRate_lipidToWater*slMagnetization(oldStep) ...
        + exchangeRate_proteinToWater*spMagnetization(oldStep));
    oldStep = newStep;
end


end

function magnetization =  relaxationFunction(relaxationRate,inputArguments)
magnetization = ...
    inputArguments{2}*(1 - 1*exp(-relaxationRate*inputArguments{1}));
end

function [waterR1,waterMagnetization] ...
    = findWaterR1ForGivenR1_SP(r1_SPArray,r1_SL,r1_water ...
    ,waterFraction,solidLipidFraction,solidProteinFraction ...
    ,exchange_waterToSolid,exchange_lipidToWater,exchange_proteinToWater ...
    ,timeStepCount,deltaT)

initialWaterMagnetization = 0;
initialSLMagnetization = 0;
initialSPMagnetization = 0;

waterR1 = zeros(1,length(r1_SPArray));
timeAxis = 0:deltaT:(timeStepCount - 1)*deltaT;

equilibriumSLMagnetization = solidLipidFraction;
slMagnetization = zeros(1,timeStepCount);

equilibriumSPMagnetization = solidProteinFraction;
spMagnetization = zeros(1,timeStepCount);

equilibriumWaterMagnetization = waterFraction;
waterMagnetization = zeros(1,timeStepCount);

opts = optimset('Display','off');

for r1_SPCounter = 1:length(r1_SPArray)
    r1_SP = r1_SPArray(r1_SPCounter);
    
    slMagnetization(1) = initialSLMagnetization;
    spMagnetization(1) = initialSPMagnetization;
    waterMagnetization(1) = initialWaterMagnetization;
    
    oldStep = 1;
    % solving Bloch equations numerically
    for newStep = 2:timeStepCount
        slMagnetization(newStep) = slMagnetization(oldStep) + deltaT  ...
            * ((equilibriumSLMagnetization - slMagnetization(oldStep)) ...
            *r1_SL - exchange_lipidToWater*slMagnetization(oldStep) ...
            + exchange_waterToSolid*waterMagnetization(oldStep));
        
        spMagnetization(newStep) = spMagnetization(oldStep) + deltaT  ...
            * ((equilibriumSPMagnetization - spMagnetization(oldStep)) ...
            *r1_SP - exchange_proteinToWater*spMagnetization(oldStep) ...
            + exchange_waterToSolid*waterMagnetization(oldStep));
        
        waterMagnetization(newStep) = waterMagnetization(oldStep)  ...
            + deltaT *((equilibriumWaterMagnetization  ...
            - waterMagnetization(oldStep))*r1_water ...
            - exchange_waterToSolid*waterMagnetization(oldStep) ...
            - exchange_waterToSolid*waterMagnetization(oldStep) ...
            + exchange_lipidToWater*slMagnetization(oldStep) ...
            + exchange_proteinToWater*spMagnetization(oldStep));
        oldStep = newStep;
    end
    
    waterR1(r1_SPCounter) = lsqcurvefit(@relaxationFunction ...
        ,1,{timeAxis equilibriumWaterMagnetization} ...
        ,waterMagnetization,[],[],opts);
    
end
end


function [factor,exponent] = fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
   fieldStrengths,relaxationRates,startValues,varargin)

opts = optimset('Display','off');

powerFunction = @(parameters,fieldStrengths)(parameters(1) ...
        *fieldStrengths.^(-parameters(2)));

if isempty(varargin)

    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,relaxationRates,[],[],opts);
    
elseif length(varargin{:}) == 2
    
    lowerLimitFactor = 0;
    upperLimitFactor = inf;
    
    varargin = sort(varargin{:});
    lowerLimitExponent = varargin(1);
    upperLimitExponent = varargin(2);
    
    lowerLimits = [lowerLimitFactor lowerLimitExponent];
    upperLimits = [upperLimitFactor upperLimitExponent];
    
    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,relaxationRates,lowerLimits,upperLimits,opts);
    
else
    error("Not implemented yet.");
    
end

factor = parameters(1);
    exponent = parameters(2);

end









