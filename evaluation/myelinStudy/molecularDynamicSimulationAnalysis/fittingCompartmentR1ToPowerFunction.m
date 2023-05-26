clear all; close all; fclose('all');


%% set dependencies
addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

%% set up
saving = 1;
highFieldDefinition = 1.4; %T
verySmallFieldDefinition = 0.4; %T

%% directories
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));

%% relaxation rates 
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
r1RatesFileName  ...
    = "SM_MW_SP_histCompartmentAndCrossR1";
r1Data =  load(resultsDir + r1RatesFolder + r1RatesFileName);

%% fitting field strength behavior of comparments to power function
% Fitting the MD-Sim. results to a power function.

[factorSM,exponentSM] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Eff_SM,[5.98 0.45]);

[factorMW,exponentMW] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Hist_MW,[1.53 0.35]);

r1Data.smFitParams = [factorSM exponentSM];
r1Data.mwFitParams = [factorMW exponentMW];

fprintf("<strong>My fitting results: a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSM,exponentSM);
fprintf("  MW: a = %.4f, b = %.4f\n",factorMW,exponentMW);

%% Fitting R1_SP in GM, WM and their average at different field strengths
% As first approach all datapoints determined by the fast-exchange model
% are used to fit the power-function. The low field values (<= 0.55T) are
% calculated without the lipid pool for GM and WM. It seems reasonable that
% for GM this low field application works best due to the low lipid
% content.

availableProteinFieldStrengths = r1Data.litFieldStrengths;
r1_ProteinWM = r1Data.r1Prot_WM;
r1_ProteinGM = r1Data.r1Prot_GM;
r1_AvgProtein = (r1_ProteinWM + r1_ProteinGM)/2;

[factorProteinWM,exponentProteinWM] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    availableProteinFieldStrengths,r1_ProteinWM ...
    ,[factorSM exponentSM]);

[factorProteinGM,exponentProteinGM] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    availableProteinFieldStrengths,r1_ProteinGM ...
    ,[factorSM exponentSM]);

[avgFactorProtein,avgExponentProtein] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    availableProteinFieldStrengths,r1_AvgProtein ...
    ,[factorSM exponentSM]);

%% Fit R1_SP only at high field strengths
% This approach uses the high field fitting of GM, WM and their average
% data to determine the power function. The goal of this investigation is
% to find out whether the high field fitting can predict the same values
% like the low field calculations are giving.

% high field fits
fieldStrengthsForHighWmFitting = availableProteinFieldStrengths( ...
    availableProteinFieldStrengths > highFieldDefinition);
r1_ProteinWmForHighFieldFitting = r1_ProteinWM( ...
    availableProteinFieldStrengths > highFieldDefinition);
[highFieldWmFactor,highFieldWmExponent] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForHighWmFitting,r1_ProteinWmForHighFieldFitting ...
    ,[factorSM,exponentMW]);

fieldStrengthsForHighGmFitting = availableProteinFieldStrengths( ...
    availableProteinFieldStrengths > highFieldDefinition);
r1_ProteinGmForHighFieldFitting = r1_ProteinGM( ...
    availableProteinFieldStrengths > highFieldDefinition);
[highFieldGmFactor,highFieldGmExponent] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForHighGmFitting,r1_ProteinGmForHighFieldFitting ...
    ,[factorSM,exponentSM]);

fieldStrengthsForHighAvgFitting = availableProteinFieldStrengths( ...
    availableProteinFieldStrengths > highFieldDefinition);
r1_ProteinAvgForHighFieldFitting = r1_AvgProtein( ...
    availableProteinFieldStrengths > highFieldDefinition);
[highFieldAvgFactor,highFieldAvgExponent] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForHighAvgFitting,r1_ProteinAvgForHighFieldFitting ...
    ,[factorSM,factorMW]);

% low field fits ignoring very low field strengths
fieldStrengthsForLowGmFitting = availableProteinFieldStrengths( ...
    availableProteinFieldStrengths > verySmallFieldDefinition);
r1_ProteinGmForLowFitting = r1_ProteinGM(availableProteinFieldStrengths ...
    > verySmallFieldDefinition);
[lowFieldGmFactor,lowFieldGmExponent] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForLowGmFitting,r1_ProteinGmForLowFitting ...
    ,[factorSM,factorMW]);

fig1 = initializeFigure();
initializeSubplot(fig1,2,2,3:4);
% WM
plot(availableProteinFieldStrengths,r1_ProteinWM,'*');
plot(r1Data.fieldStrengths ...
    ,factorProteinWM*r1Data.fieldStrengths.^(-exponentProteinWM));
plot(r1Data.fieldStrengths ...
    ,highFieldWmFactor*r1Data.fieldStrengths.^(-highFieldWmExponent),'--');
% GM
plot(availableProteinFieldStrengths,r1_ProteinGM,'*');
plot(r1Data.fieldStrengths ...
    ,factorProteinGM*r1Data.fieldStrengths.^(-exponentProteinGM));
plot(r1Data.fieldStrengths ...
    ,highFieldGmFactor*r1Data.fieldStrengths.^(-highFieldGmExponent),'--');
% avgerage
plot(availableProteinFieldStrengths,r1_AvgProtein,'*');
plot(r1Data.fieldStrengths ...
    ,avgFactorProtein*r1Data.fieldStrengths.^(-avgExponentProtein));
plot(r1Data.fieldStrengths ...
    ,highFieldAvgFactor*r1Data.fieldStrengths.^(-highFieldAvgExponent) ...
    ,'--');
axis([0 inf 0 25]);
legend("WM","WM Fit","WM Fit high field","GM","GM Fit" ...
    ,"GM Fit high field","Average","Average fit" ...
    ,"Average fit high field");
xlabel("Field strength [T]");
ylabel("R$_{1,SP}$");

% only GM
initializeSubplot(fig1,2,2,2);
plot(availableProteinFieldStrengths,r1_ProteinGM,'*');
plot(r1Data.fieldStrengths ...
    ,factorProteinGM*r1Data.fieldStrengths.^(-exponentProteinGM));
plot(r1Data.fieldStrengths ...
    ,highFieldGmFactor*r1Data.fieldStrengths.^(-highFieldGmExponent),'--');
plot(r1Data.fieldStrengths ...
    ,lowFieldGmFactor*r1Data.fieldStrengths.^(-lowFieldGmExponent),':');
title("Only GM")
legend("All data points","Fitted" ...
    ,sprintf("Fitted $>$ %.2f T",highFieldDefinition) ...
    ,sprintf("Fitted $>$ %.2f T",verySmallFieldDefinition));

% parameters
initializeSubplot(fig1,2,2,1);
writeIntoPlotWindow(fig1,[2,2,1],sprintf( ...
    "Fitted power functions:\n High fields $>$ %.2f T\n " ...
    + "Low field $>$ %.2f: \n\n" ...
    + "R$_{1,SP,WM}$ = %.3f*B$_{0}^{-%.3f}$ \n" ...
    + "R$_{1,SP,WM}^{high}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP,GM}$ = %.3f*B$_{0}^{-%.3f}$ \n" ...
    + "R$_{1,SP,GM}^{high}$ = %.3f*B$_{0}^{-%.3f}$ \n" ...
    + "R$_{1,SP,GM}^{low}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP,AVG}$ = %.3f*B$_{0}^{-%.3f}$ \n" ...
    + "R$_{1,SP,AVG}^{high}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    ,highFieldDefinition,verySmallFieldDefinition ...
    ,factorProteinWM,exponentProteinWM ...
    ,highFieldWmFactor,highFieldWmExponent ...
    ,factorProteinGM,exponentProteinGM ...
    ,highFieldGmFactor,highFieldGmExponent ...
    ,lowFieldGmFactor,lowFieldGmExponent ...
    ,avgFactorProtein,avgExponentProtein ...
    ,highFieldAvgFactor,highFieldAvgExponent));

r1Data.highFieldDefinition = highFieldDefinition;

% saving WM
r1Data.r1_SP_WM = factorProteinWM * r1Data.fieldStrengths ...
    .^(-exponentProteinWM);
r1Data.spFitParamsWM_AllFieldStrengths = [factorProteinWM exponentProteinWM];
r1Data.r1_SP_WmHighField = highFieldWmFactor * r1Data.fieldStrengths ...
    .^(-highFieldWmExponent);
r1Data.spFitParamsWM_HighFields = [highFieldWmFactor highFieldWmExponent];
% saving GM
r1Data.r1_SP_GM = factorProteinGM * r1Data.fieldStrengths ...
    .^(-exponentProteinGM);
r1Data.spFitParamsGM_AllFieldStrengths = [factorProteinGM exponentProteinGM];
r1Data.r1_SP_GmHighField = highFieldGmFactor * r1Data.fieldStrengths ...
    .^(-highFieldGmExponent);
r1Data.spFitParamsGM_HighFields = [highFieldGmFactor highFieldGmExponent];
% saving Avg
r1Data.r1_SP_Avg = avgFactorProtein * r1Data.fieldStrengths ...
    .^(-avgExponentProtein);
r1Data.spFitParamsAvg_AllFieldStrengths = [avgFactorProtein avgExponentProtein];
r1Data.r1_SP_AvgHighField = highFieldAvgFactor * r1Data.fieldStrengths ...
    .^(-highFieldAvgExponent);
r1Data.spFitPararmsAvg_highFields = [highFieldAvgFactor highFieldAvgExponent];

if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"","solidProteinR1FittedToPowerFunction",false); 
    save(resultsDir + r1RatesFolder + datestr(now,"yyyymmdd") ...
        + "_SM_WM_SP_histCompartmentAndCrossR1_fitted",'-struct','r1Data');
end

%% Validation of my fitting results compared to Schyboll2019

highFieldStrengthIndices = ... 
    r1Data.fieldStrengths > 1.4999 & r1Data.fieldStrengths < 1.5001 ...
    | r1Data.fieldStrengths > 2.999 & r1Data.fieldStrengths < 3.0001 ...
    | r1Data.fieldStrengths > 6.999 & r1Data.fieldStrengths < 7.0001;

[highFieldFactorSM,highFieldExponentSM] =  ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths(highFieldStrengthIndices ) ...
    ,r1Data.r1Eff_SM(highFieldStrengthIndices),[factorSM,exponentSM]);

[highFieldFactorMW,highFieldExponentMW] = ...
    fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths(highFieldStrengthIndices) ...
    ,r1Data.r1Hist_MW(highFieldStrengthIndices),[factorMW,exponentMW]);

fprintf("<strong>Reproducing Schyboll2019 with high field strengths :" ...
    + "a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",highFieldFactorSM,highFieldExponentSM);
fprintf("  MW: a = %.4f, b = %.4f\n",highFieldFactorMW,highFieldExponentMW);

factorSchybollSM = 5.98;
factorSchybollMW = 1.53;
exponentSchybollSM = 0.45;
exponentSchybollMW = 0.31;

lowFieldStrengthIndices = ...
    r1Data.fieldStrengths > 0.34999 & r1Data.fieldStrengths < 0.35001 ...
    | r1Data.fieldStrengths > 0.54999 & r1Data.fieldStrengths < 0.550001;

lowFieldsAppendedToSchyboll = [r1Data.fieldStrengths( ...
    lowFieldStrengthIndices) ...
    r1Data.fieldStrengths(highFieldStrengthIndices)];
schybollR1AppendWithLowFieldR1 = [r1Data.r1Eff_SM(lowFieldStrengthIndices) ...
    factorSchybollSM * [1.5 3 7].^(-exponentSchybollSM)];
[factorSM_SchybollExpanded,exponentSM_SchybollExpanded] ...
    = fitFieldStrengthBehaviorToPowerFunction(lowFieldsAppendedToSchyboll ...
    ,schybollR1AppendWithLowFieldR1,[factorSM,exponentSM]);
    
fprintf("<strong>Schyboll2019's findings: a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSchybollSM,exponentSchybollSM);
fprintf("  MW: a = %.4f, b = %.4f\n",factorSchybollMW,exponentSchybollMW);

fprintf("<strong>Schyboll2019's finding appended with low field:" ...
    + "a*B_0^(-b) </strong> \n");
fprintf("  SM: a = %.4f, b = %.4f\n",factorSM_SchybollExpanded ...
    ,exponentSM_SchybollExpanded);

fig2 = initializeFigure();
initializeSubplot(fig2,2,1,1);
title("Reproduced behavior");
plot(lowFieldsAppendedToSchyboll,schybollR1AppendWithLowFieldR1,'*');
plot(r1Data.fieldStrengths ...
    ,highFieldFactorSM*r1Data.fieldStrengths.^(-highFieldExponentSM));
plot(r1Data.fieldStrengths ...
    ,factorSchybollSM*r1Data.fieldStrengths.^(-exponentSchybollSM));
plot(r1Data.fieldStrengths ...
    ,factorSM_SchybollExpanded ...
    *r1Data.fieldStrengths.^(-exponentSM_SchybollExpanded));
legend("Schyboll + my data","Fitted high field","Schyboll2019 fit" ...
    ,"Fitted Schyboll expanded");
xlabel("Field strength [T]");
ylabel("R$_{1,SM}$");

initializeSubplot(fig2,2,1,2);
writeIntoPlotWindow(fig1,[2,1,2],sprintf( ...
    "Fitted power functions:\n\n" ...
    + "R$_{1,SM}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SM}^{high}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP}^{Schyboll}$ = %.3f*B$_{0}^{-%.3f}$ \n\n" ...
    + "R$_{1,SP}^{Schyboll,expanded}$ = %.3f*B$_{0}^{-%.3f}$" ...
    ,factorSM,exponentSM,highFieldFactorSM,highFieldExponentSM ...
    ,factorSchybollSM,exponentSchybollSM,factorSM_SchybollExpanded ...
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

title("Solid myelin (SM)");
legend("MD-Sim.","Fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

% myelin water (MW) fit
initializeSubplot(fig,2,2,2);
plot(r1Data.fieldStrengths,r1Data.r1Hist_MW,'*');
plot(r1Data.fieldStrengths,factorMW*r1Data.fieldStrengths.^(-exponentMW));

title("Myelin water (MW)");
legend("MD-Sim.","Fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

% solid protein (SP) fit
% -> GM only because most valid low field assumption
initializeSubplot(fig,2,2,3);
plot(availableProteinFieldStrengths,r1_ProteinGM,'*');
plot(r1Data.fieldStrengths,factorProteinGM*r1Data.fieldStrengths ...
    .^(-exponentProteinGM));

title("Solid protein (SP)");
legend("Caclulated in GM","Fitted in GM");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

writeIntoPlotWindow(fig,[2,2,4],sprintf( ...
    "R$_{1,SM}$ = %.3f*B$_{0}^{-%.3f}$\n\n" ...
    + "R$_{1,MW}$ = %.3f*B$_{0}^{-%.3f}$\n\n" ...
    + "R$_{1,SP,GM}$ = %.3f*B$_{0}^{-%.3f}$\n\n" ...
    + "R$_{1,SP,WM}$ = %.3f*B$_{0}^{-%.3f}$\n\n" ...
    + "R$_{1,SP,avg}$ = %.3f*B$_{0}^{-%.3f}$\n\n" ...
    ,factorSM,exponentSM,factorMW,exponentMW,factorProteinGM ...
    ,exponentProteinGM,factorProteinWM,exponentProteinWM ...
    ,avgFactorProtein,avgExponentProtein));

if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"solidMyelin_MyelinWater_SolidProtein" ...
        ,"R1FittedToPowerFunction",false); 
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














