clc; clear all; close all; fclose('all');

%% set dependencies
addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

%% set up
saving = 1;

%% directories
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));

%% relaxation rates 
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
r1RatesFileName  ...
    = "20230418_histCompartmentAndCrossR1WithProteinR1";
r1Data =  load(resultsDir + r1RatesFolder + r1RatesFileName);

%% fitting field strength behavior of comparments to power function

[factorSM,exponentSM] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Eff_SM,[5.98 0.45]);

[factorMW,exponentMW] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Hist_MW,[1.53 0.35]);

fprintf("<strong>My fitting results: a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSM,exponentSM);
fprintf("  MW: a = %.4f, b = %.4f\n",factorMW,exponentMW);

%% Fitting on GM and WM separately at different field strengths and on average
% Due to fast relaxation of lipids and high lipid content in WM a
% calculation of R1_SP produces negative relaxation rates in WM. This is,
% due to small lipid content, not the case for GM. That's why smaller field
% strengths for fitting of GM R1_SP can be used and for GM a fast-exchang
% model is more realistic than for WM. Therefore, R1 calculated in WM with
% field strengths >= 1.5T and in GM for field strengths >= 0.55T are used
% for fitting to the power function.

protFieldStrengths = r1Data.litFieldStrengths;
r1_ProteinWM = r1Data.r1Prot_WM;
r1_ProteinGM = r1Data.r1Prot_GM;

smallestFieldStrengthForFittingWM = 0.3;
smallestFieldStrengthForFittingGM = 0.3;

fieldStrengthsForFittingWM = protFieldStrengths(protFieldStrengths ...
    > smallestFieldStrengthForFittingWM);
r1_ProteinWMForFitting = r1_ProteinWM(protFieldStrengths ...
    > smallestFieldStrengthForFittingWM);

fieldStrengthsForFittingGM = protFieldStrengths(protFieldStrengths ...
    > smallestFieldStrengthForFittingGM);
r1_ProteinGMForFitting = r1_ProteinGM(protFieldStrengths ...
    > smallestFieldStrengthForFittingGM);

% Fitting WM for field strengths > 1.4T
[factorProtWM,exponentProtWM] = ...
    fitFieldStrengthBehaviorToPowerFunction(fieldStrengthsForFittingWM ...
    ,r1_ProteinWMForFitting,[factorSM exponentSM]);

% Fitting GM for field strengths > 0.5T
[factorProtGM,exponentProtGM] = ...
    fitFieldStrengthBehaviorToPowerFunction(fieldStrengthsForFittingGM ...
    ,r1_ProteinGMForFitting,[factorSM exponentSM]);

% fitting AVG for field strengths > 0.5T
% r1_ProteinGMForFitting(1:2)
fieldStrengthsForAvgFitting = fieldStrengthsForFittingGM;
r1_AvgProtein = (r1_ProteinWMForFitting ...
    + r1_ProteinGMForFitting)/2;

[avgFactorProt,avgExponentProt] = fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForAvgFitting,r1_AvgProtein,[factorSM exponentSM]);

% Why can smaller field strengths be used for fitting R1_SP in GM? This
% assumes that a fast-exchange model is valid for GM. This can be explained
% by a high protein content (67%) and a low lipid content (37%) in GM.
% Thus, the highly field-dependent relaxation rate of lipids does not
% influence the measured signal much and can be ignored. The results of the
% fitting seem to give valid results when compared to Leuze2017. Their
% investigations revealed that no contrast between WM and GM can be seen
% after clearing. Thus, although protein content in WM and GM is different,
% no contrast is observed in T1-relaxation. Therefore, proteins have
% clearly not that much influence on water than lipids. 

%% Fit R1_SP with high and low field strengths
fieldStrengthsForHighAvgFitting = fieldStrengthsForAvgFitting( ...
    fieldStrengthsForAvgFitting > 0.4);
r1_AvgHighField = r1_AvgProtein(fieldStrengthsForAvgFitting > 0.4);

[highFAvgFactor,highFAvgExponent] = fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForHighAvgFitting,r1_AvgHighField,[factorSM,factorMW]);

fig1 = initializeFigure();
initializeSubplot(fig1,2,1,1);
plot(fieldStrengthsForAvgFitting,r1_AvgProtein,'*');
plot(r1Data.fieldStrengths,avgFactorProt ...
    *r1Data.fieldStrengths.^(-avgExponentProt));
plot(r1Data.fieldStrengths,highFAvgFactor ...
    *r1Data.fieldStrengths.^(-highFAvgExponent));
legend("Data points","All data points fitting","High field fitting only");
xlabel("Field strength [T]");
ylabel("R$_{1,SP}$");

initializeSubplot(fig1,2,1,2);
writeIntoPlotWindow(fig1,[2,1,2],sprintf( ...
    "Fitted power functions:\n\n" ...
    + "R$_{1,SP,AVG}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP,AVG}^{high}$ = %.3f*B$_{0}^{-%.3f}$" ...
    ,avgFactorProt,avgExponentProt,highFAvgFactor,highFAvgExponent));


if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"solidProteinR1FittedToPowerFunction","",false); 
end
%% Reproduce fitting results from Felix by only fitting higher field strengths

fieldStrengthIndices = ... 
    r1Data.fieldStrengths > 1.4999 & r1Data.fieldStrengths < 1.5001 ...
    | r1Data.fieldStrengths > 2.999 & r1Data.fieldStrengths < 3.0001 ...
    | r1Data.fieldStrengths > 6.999 & r1Data.fieldStrengths < 7.0001;

[factorSM_red,exponentSM_red] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths(fieldStrengthIndices ) ...
    ,r1Data.r1Eff_SM(fieldStrengthIndices),[factorSM,exponentSM]);

[factorMW_red,exponentMW_red] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths(fieldStrengthIndices) ...
    ,r1Data.r1Hist_MW(fieldStrengthIndices),[factorMW,exponentMW]);


fprintf("<strong>Reproducing Schyboll2019: a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSM_red,exponentSM_red);
fprintf("  MW: a = %.4f, b = %.4f\n",factorMW_red,exponentMW_red);

factorSchybollSM = 5.98;
factorSchybollMW = 1.53;
exponentSchybollSM = 0.45;
exponentSchybollMW = 0.31;

lowFieldStrengthIndices = ...
    r1Data.fieldStrengths > 0.34999 & r1Data.fieldStrengths < 0.35001 ...
    | r1Data.fieldStrengths > 0.54999 & r1Data.fieldStrengths < 0.550001;

fieldsForSchyboll = [r1Data.fieldStrengths(lowFieldStrengthIndices) ...
    r1Data.fieldStrengths(fieldStrengthIndices)];
r1ForSchyboll = [r1Data.r1Eff_SM(lowFieldStrengthIndices) ...
    factorSchybollSM * [1.5 3 7].^(-exponentSchybollSM)];
[factorSM_SchybollExpanded,exponentSM_SchybollExpanded] ...
    = fitFieldStrengthBehaviorToPowerFunction(fieldsForSchyboll ...
    ,r1ForSchyboll,[factorSM,exponentSM]);
    
fprintf("<strong>Schyboll2019: a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSchybollSM,exponentSchybollSM);
fprintf("  MW: a = %.4f, b = %.4f\n",factorSchybollMW,exponentSchybollMW);

fig2 = initializeFigure();
initializeSubplot(fig2,2,1,1);
title("Reproduced behavior");
plot(r1Data.fieldStrengths(fieldStrengthIndices) ...
    ,r1Data.r1Eff_SM(fieldStrengthIndices),'*');
plot(r1Data.fieldStrengths ...
    ,factorSM_red*r1Data.fieldStrengths.^(-exponentSM_red));
plot(r1Data.fieldStrengths ...
    ,factorSchybollSM*r1Data.fieldStrengths.^(-exponentSchybollSM));
plot(fieldsForSchyboll,r1ForSchyboll,'*');
plot(r1Data.fieldStrengths ...
    ,factorSM_SchybollExpanded ...
    *r1Data.fieldStrengths.^(-exponentSM_SchybollExpanded));
legend("Data points","Fitted","Schyboll2019" ...
    ,"Data points Schyboll + my data","Fitted Schyboll expanded");
xlabel("Field strength [T]");
ylabel("R$_{1,SM}$");

initializeSubplot(fig2,2,1,2);
writeIntoPlotWindow(fig1,[2,1,2],sprintf( ...
    "Fitted power functions:\n\n" ...
    + "R$_{1,SM}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SM}^{high}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP}^{Schyboll}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP}^{Schyboll,expanded}$ = %.3f*B$_{0}^{-%.3f}$" ...
    ,factorSM,exponentSM,factorSM_red,exponentSM_red,factorSchybollSM ...
    ,exponentSchybollSM,factorSM_SchybollExpanded ...
    ,exponentSM_SchybollExpanded));


if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"","compareSchyboll2019ToMyResults",false); 
end


%% transfer behavior 
% transfer field strength behavior of solid myelin to field strength
% behavior of solide protein for high field strengths because for low field
% strengths a fast exchange model cannot be assumed andymore

% smallestFieldStrengthForFitting = 1.4; %T
% 
% protFieldStrengths = r1Data.litFieldStrengths;
% r1_ProteinWM = r1Data.r1Prot_WM;
% r1_ProteinGM = r1Data.r1Prot_GM;
% 
% fieldStrengthsForFitting = protFieldStrengths(protFieldStrengths ...
%     > smallestFieldStrengthForFitting);
% r1_ProteinWMForFitting = r1_ProteinWM(protFieldStrengths ...
%     > smallestFieldStrengthForFitting);
% r1_ProteinGMForFitting = r1_ProteinGM(protFieldStrengths ...
%     > smallestFieldStrengthForFitting);
% 
% % Fitting calulcated R1s of Protein pool to power function and using only
% % higher fields
% 
% [factorProtWM,exponentProtWM] = fitFieldStrengthBehaviorToPowerFunction(...
%     fieldStrengthsForFitting,r1_ProteinWMForFitting,[factorSM exponentSM]);
% 
% [factorProtGM,exponentProtGM] = fitFieldStrengthBehaviorToPowerFunction(...
%     fieldStrengthsForFitting,r1_ProteinGMForFitting,[factorSM exponentSM]);
% 
% % -> Fit works, but not so well. There is a large difference between GM and
% % WM. The values for the fitted parameters work not too good. They are more
% % like the values from Schyboll  2019 and less like the ones from solid
% % myelin. The reason for that is that the low field is not in the fitting
% % data which leads to results that are not so well.
% 
% % Can fitting the factor and not the exponent be a helpfull approach?
% % Assume that change in exponent in fit from caluclated data is due to
% % missing low field strengths (compare to Schyboll2019)
% 
% factorProtWMBasedOnSM = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
%     fieldStrengthsForFitting,r1_ProteinWMForFitting,exponentSM,factorSM);
% 
% factorProtGMBasedOnSM = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
%     fieldStrengthsForFitting,r1_ProteinGMForFitting,exponentSM,factorSM);
% 
% % The curves from SM and SP look very similiar now. Maybe try to fit the
% % a curve to the datapoints where the exponent can be fitted with limited
% % range.
% 
% exponentBound = exponentSM * [0.8 1.2];
% 
% [factorProtWM_boundedExp,exponentProtWM_boundedExp]  ...
%     = fitFieldStrengthBehaviorToPowerFunction(...
%     fieldStrengthsForFitting,r1_ProteinWMForFitting ...
%     ,[factorSM exponentSM],exponentBound);
% 
% [factorProtGM_boundedExp,exponentProtGM_boundedExp] ...
%     = fitFieldStrengthBehaviorToPowerFunction( ...
%     fieldStrengthsForFitting,r1_ProteinWMForFitting ...
%     ,[factorSM exponentSM],exponentBound);

% The problem with this approach is that the exponential fitting parameter
% is running into the exponentBound. Therefore this approach is not so
% suitable either.


%% The best approach to determine R1_prot
% Calculate low-field values is not possible due to a not met condition of
% fast exchange between the pools because the solid myelin pool is
% relaxing too fast. Therefore, a fitting with data from field strengths > 1.5T
% should be the best approach. Although we know that this approach
% underestimates the two parameters in the power function, approaches where
% only one parameter is fitted would predict protein R1 values that are
% very close to what we observe in SM. This is not valid if we look at what
% Leuze2017 observed, where no contrast between WM and GM was seen for
% cleared tissue. 

%% plotting
fig = initializeFigure();

% solid myelin (SM) fit
initializeSubplot(fig,2,2,1);
plot(r1Data.fieldStrengths,r1Data.r1Eff_SM,'*');
plot(r1Data.fieldStrengths,factorSM*r1Data.fieldStrengths.^(-exponentSM));

title(sprintf("R$_{1}$: R$_{1,SM}$ = %.3f * B$_{0}^{-%.3f}$" ...
    ,factorSM,exponentSM));
title("Solid myelin (SM)");
legend("Caclulated","Fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

% myelin water (MW) fit
initializeSubplot(fig,2,2,2);
plot(r1Data.fieldStrengths,r1Data.r1Hist_MW,'*');
plot(r1Data.fieldStrengths,factorMW*r1Data.fieldStrengths.^(-exponentMW));

title("Myelin water (MW)");
legend("Caclulated","Fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

% WM & GM solid protein (SP) fit
initializeSubplot(fig,2,2,3);
% -> WM
plot(fieldStrengthsForFittingWM,r1_ProteinWMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtWM*r1Data.fieldStrengths ...
    .^(-exponentProtWM));
% -> GM
plot(fieldStrengthsForFittingGM,r1_ProteinGMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtGM*r1Data.fieldStrengths ...
    .^(-exponentProtGM));
% -> average
plot(fieldStrengthsForAvgFitting,r1_AvgProtein,'*');
plot(r1Data.fieldStrengths,avgFactorProt*r1Data.fieldStrengths ...
    .^(-avgExponentProt));

title("Solid protein (SP)");
legend("Caclulated SP in WM","Fitted SP in WM","Caclulated SP in GM" ...
    ,"Fitted SP in GM","Avg. calculated","Avg fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

writeIntoPlotWindow(fig,[2,2,4],sprintf( ...
    "R$_{1,SM}$ = %.3f*B$_{0}^{-%.3f}$\n" + ...
    "R$_{1,MW}$ = %.3f*B$_{0}^{-%.3f}$\n" + ...
    "R$_{1,SP,WM}$ = %.3f*B$_{0}^{-%.3f}$\n" + ...
    "R$_{1,SP,GM}$ = %.3f*B$_{0}^{-%.3f}$\n" + ...
    "R$_{1,SP,AVG}$ = %.3f*B$_{0}^{-%.3f}$"...
    ,factorSM,exponentSM,factorMW,exponentMW,factorProtWM ...
    ,exponentProtWM,factorProtGM,exponentProtGM ...
    ,avgFactorProt,avgExponentProt));


if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"solidMyelin_MyelinWater_SolidProtein","R1FittedToPowerFunction",false); 
end

%% functions

function [factor,exponent] = fitFieldStrengthBehaviorToPowerFunction( ...
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


% function factor = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
%     fieldStrengths,relaxationRates,exponentValue,startValueForFactor)
% 
% input(1,:) = fieldStrengths;
% input(2,:) = exponentValue;
% 
% powerFunction = @(factor,input)(factor ...
%     *input(1,:).^(-input(2,:)));
% 
% factor = lsqcurvefit(powerFunction,startValueForFactor,input ...
%     ,relaxationRates);
% 
% 
% end














