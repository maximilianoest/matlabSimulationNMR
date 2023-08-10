clc; clear all; close all;

fieldStrengthInLit = [1.5 3 7];
r1SLInLit = [4.95+0.112 3.72+0.090 2.43+0.065];
r1LWInLit = [1.175+0.165 0.978+0.132 0.718+0.095];
fprintf("Schyboll2019:\n");
[litFact,litExp,~] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengthInLit,r1SLInLit,[1 1]);

myR1Data = load("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\wholeMyelin_relaxationRates\SM_MW_SP_histCompartmentAndCrossR1.mat");
fieldStrengthsToFit = [0.35 fieldStrengthInLit];
fieldStrengthsULF = myR1Data.fieldStrengths; 

r1SL = [];
r1LW = [];
for fieldStrength = fieldStrengthsToFit
r1SL(end+1) = myR1Data.r1Eff_SM( ...
    myR1Data.fieldStrengths > fieldStrength - 0.0001 ...
    & myR1Data.fieldStrengths < fieldStrength + 0.0001); %#ok<SAGROW>
r1LW(end+1) = myR1Data.r1Hist_MW( ...
    myR1Data.fieldStrengths > fieldStrength - 0.0001 ...
    & myR1Data.fieldStrengths < fieldStrength + 0.0001); %#ok<SAGROW>
end
fprintf("\nMy data:\n");
[simFact,simExp,~] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengthsToFit,r1SL,[1 1]);

fprintf("\nMy data high field: \n");
[simFactHighField,simExpHighField,~] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengthsToFit(2:end),r1SL(2:end),[1 1]);

fprintf("\nSchyboll2019 with my low field data:\n");
[schybollWithLowFieldFact,schybollWithLowFieldExp,~] ...
    = fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengthsToFit,[r1SL(1) r1SLInLit],[1 1]);

r1SL_ULF = [];
for fieldStrength = fieldStrengthsULF
r1SL_ULF(end+1) = myR1Data.r1Eff_SM( ...
    myR1Data.fieldStrengths > fieldStrength - 0.0001 ...
    & myR1Data.fieldStrengths < fieldStrength + 0.0001); %#ok<SAGROW>
end
fprintf("\nMy ULF data:\n");
[simFactULF,simExpULF,covariance_ULF] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunctionWithoutWeights( ...
    fieldStrengthsULF,r1SL_ULF,[1 1]);

fprintf("\nWM lit data\n");
[WMLitFact,WMLitExp,~] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    myR1Data.wmFieldStrengthArray,myR1Data.r1WMAvg,[1 1] ...
    ,myR1Data.r1WMSTD);

fprintf("GM lit data\n");
[GMLitFact,GMLitExp,~] = ...
    fitFieldStrengthBehaviorOfR1ToPowerFunction( ...
    myR1Data.gmFieldStrengthArray,myR1Data.r1GMAVg,[1 1] ... 
    ,myR1Data.r1GMSTD);


initializeFigure();
plot(myR1Data.fieldStrengths,litFact * myR1Data.fieldStrengths.^(-litExp));
plot(myR1Data.fieldStrengths,simFact * myR1Data.fieldStrengths.^(-simExp));
plot(myR1Data.fieldStrengths,simFactHighField  ...
    * myR1Data.fieldStrengths.^(-simFactHighField));
plot(myR1Data.fieldStrengths,schybollWithLowFieldFact ...
    * myR1Data.fieldStrengths.^(-schybollWithLowFieldExp));
plot(myR1Data.fieldStrengths,simFactULF * myR1Data.fieldStrengths.^( ...
    -simExpULF));

legend("Schyboll2019","my data","my data high field" ...
    ,"my data + Schyboll2019","ULF");
axis([0 inf 0 100])





%% functions
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

fprintf("R1(B_0) = %.4f B_0 ^(-%.4f) \n",factor,exponent);

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
fprintf("R1(B_0) = %.4f B_0 ^(-%.4f) \n",factor,exponent);

end