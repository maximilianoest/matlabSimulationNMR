clc; clear all; close all; fclose('all');

%% set dependencies
addpath(genpath(sprintf('..%s..%s..%slibrary',filesep,filesep,filesep)));
addpath(genpath(sprintf('..%s..%s..%stxtFiles' ...
    ,createFilesepStringArray(3))));
constants = readConstantsFile('constants.txt');

%% set up
saving = 0;

%% directories
resultsDir = sprintf("..%s..%s..%sRESULTS%",createFilesepStringArray(3));

%% relaxation rates 
r1RatesFolder = sprintf("%swholeMyelin_relaxationRates%s" ...
    ,createFilesepStringArray(2));
r1RatesFileName  ...
    = "20230414_histCompartmentAndCrossR1WithProteinR1";
r1Data =  load(resultsDir + r1RatesFolder + r1RatesFileName);



%% fitting field strength behavior of comparments to power function

[factorSM,exponentSM] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Eff_SM,[5.98 0.45]);

[factorMW,exponentMW] = fitFieldStrengthBehaviorToPowerFunction( ...
    r1Data.fieldStrengths,r1Data.r1Hist_MW,[1.53 0.35]);



%% transfer behavior 
% transfer field strength behavior of solid myelin to field strength
% behavior of solide protein for high field strengths because for low field
% strengths a fast exchange model cannot be assumed andymore

smallestFieldStrengthForFitting = 1.4; %T

protFieldStrengths = r1Data.litFieldStrengths;
r1_ProteinWM = r1Data.r1Prot_WM;
r1_ProteinGM = r1Data.r1Prot_GM;

r1_AvgProtein = (r1_ProteinWM + r1_ProteinGM)/2;

fieldStrengthsForFitting = protFieldStrengths(protFieldStrengths ...
    > smallestFieldStrengthForFitting);
r1_ProteinWMForFitting = r1_ProteinWM(protFieldStrengths ...
    > smallestFieldStrengthForFitting);
r1_ProteinGMForFitting = r1_ProteinGM(protFieldStrengths ...
    > smallestFieldStrengthForFitting);

% Fitting calulcated R1s of Protein pool to power function and using only
% higher fields

[factorProtWM,exponentProtWM] = fitFieldStrengthBehaviorToPowerFunction(...
    fieldStrengthsForFitting,r1_ProteinWMForFitting,[factorSM exponentSM]);

[factorProtGM,exponentProtGM] = fitFieldStrengthBehaviorToPowerFunction(...
    fieldStrengthsForFitting,r1_ProteinGMForFitting,[factorSM exponentSM]);

% -> Fit works, but not so well. There is a large difference between GM and
% WM. The values for the fitted parameters work not too good. They are more
% like the values from Schyboll  2019 and less like the ones from solid
% myelin. The reason for that is that the low field is not in the fitting
% data which leads to results that are not so well.

% Can fitting the factor and not the exponent be a helpfull approach?
% Assume that change in exponent in fit from caluclated data is due to
% missing low field strengths (compare to Schyboll2019)

factorProtWMBasedOnSM = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForFitting,r1_ProteinWMForFitting,exponentSM,factorSM);

factorProtGMBasedOnSM = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForFitting,r1_ProteinGMForFitting,exponentSM,factorSM);

% The curves from SM and SP look very similiar now. Maybe try to fit the
% a curve to the datapoints where the exponent can be fitted with limited
% range.

exponentBound = exponentSM * [0.8 1.2];

[factorProtWM_boundedExp,exponentProtWM_boundedExp]  ...
    = fitFieldStrengthBehaviorToPowerFunction(...
    fieldStrengthsForFitting,r1_ProteinWMForFitting ...
    ,[factorSM exponentSM],exponentBound);

[factorProtGM_boundedExp,exponentProtGM_boundedExp] ...
    = fitFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengthsForFitting,r1_ProteinWMForFitting ...
    ,[factorSM exponentSM],exponentBound);

% The problem with this approach is that the exponential fitting parameter
% is running into the exponentBound. Therefore this approach is not so
% suitable either.


%% plotting
fig = initializeFigure();
initializeSubplot(fig,2,2,1);
plot(r1Data.fieldStrengths,r1Data.r1Eff_SM,'*');
plot(r1Data.fieldStrengths,factorSM*r1Data.fieldStrengths.^(-exponentSM));

title(sprintf("R$_{1}$: R$_{1,SM}$ = %.3f * B$_{0}^{-%.3f}$" ...
    ,factorSM,exponentSM));
legend("Caclulated SM","Fitted SM");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

initializeSubplot(fig,2,2,3);
plot(r1Data.fieldStrengths,r1Data.r1Hist_MW,'*');
plot(r1Data.fieldStrengths,factorMW*r1Data.fieldStrengths.^(-exponentMW));
title(sprintf("Solid myelin R$_{1}$: R$_{1,MW}$ = %.3f * B$_{0}^{-%.3f}$" ...
    ,factorMW,exponentMW));
legend("Caclulated","Fitted");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");


initializeSubplot(fig,2,2,2);
plot(fieldStrengthsForFitting,r1_ProteinWMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtWM*r1Data.fieldStrengths ...
    .^(-exponentProtWM));
plot(fieldStrengthsForFitting,r1_ProteinGMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtGM*r1Data.fieldStrengths ...
    .^(-exponentProtGM));
title(sprintf("R$_{1,SP,WM(GM)}$ = %.3f (%.3f) * B$_{0}^{-%.3f (-%.3f)}$" ...
    ,factorProtWM,factorProtGM,exponentProtWM,exponentProtGM));
legend("Caclulated SP on WM","Fitted SP on WM","Caclulated SP on GM" ...
    ,"Fitted SP on GM");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

initializeSubplot(fig,2,2,4);
plot(fieldStrengthsForFitting,r1_ProteinWMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtWMBasedOnSM*r1Data.fieldStrengths ...
    .^(-exponentSM));
plot(fieldStrengthsForFitting,r1_ProteinGMForFitting,'*');
plot(r1Data.fieldStrengths,factorProtGMBasedOnSM*r1Data.fieldStrengths ...
    .^(-exponentSM));

plot(r1Data.fieldStrengths,factorProtWM_boundedExp*r1Data.fieldStrengths ...
    .^(-exponentProtWM_boundedExp));
plot(r1Data.fieldStrengths,factorProtGM_boundedExp*r1Data.fieldStrengths ...
    .^(-exponentProtGM_boundedExp));
title(sprintf("R$_{1,SP,WM(GM)}$ = %.3f (%.3f) * B$_{0}^{-%.3f (-%.3f)}$" ...
    ,factorProtWMBasedOnSM,factorProtGMBasedOnSM,exponentSM,exponentSM));
legend("Caclulated SP on WM","Fitted SP on WM","Caclulated SP on GM" ...
    ,"Fitted SP on GM","bounded WM","bounded GM");
ylabel("R$_1$ [Hz]");
xlabel("Field strength [T]");

if saving
    saveFigureTo(resultsDir + r1RatesFolder,datestr(now,"yyyymmdd") ...
        ,"solidMyelinAndMyelinWater","R1FittedToPowerFunction",false); 
end

%% functions

function [factor,exponent] = fitFieldStrengthBehaviorToPowerFunction( ...
   fieldStrengths,relaxationRates,startValues,varargin)

powerFunction = @(parameters,fieldStrengths)(parameters(1) ...
        *fieldStrengths.^(-parameters(2)));

if isempty(varargin)

    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,relaxationRates);
    
elseif length(varargin{:}) == 2
    
    lowerLimitFactor = 0;
    upperLimitFactor = inf;
    
    varargin = sort(varargin{:});
    lowerLimitExponent = varargin(1);
    upperLimitExponent = varargin(2);
    
    lowerLimits = [lowerLimitFactor lowerLimitExponent];
    upperLimits = [upperLimitFactor upperLimitExponent];
    
    parameters = lsqcurvefit(powerFunction,startValues,fieldStrengths ...
        ,relaxationRates,lowerLimits,upperLimits);
    
else
    error("Not implemented yet.");
    
end

factor = parameters(1);
    exponent = parameters(2);

end


function factor = fitFactorOfFieldStrengthBehaviorToPowerFunction( ...
    fieldStrengths,relaxationRates,exponentValue,startValueForFactor)

input(1,:) = fieldStrengths;
input(2,:) = exponentValue;

powerFunction = @(factor,input)(factor ...
    *input(1,:).^(-input(2,:)));

factor = lsqcurvefit(powerFunction,startValueForFactor,input ...
    ,relaxationRates);


end














