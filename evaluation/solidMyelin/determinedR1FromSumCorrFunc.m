clc; clear all; close all; fclose('all');

addpath(genpath("../../library/"));
constants = readConstantsFile("../../txtFiles/constants.txt");
results = load("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\solidMyelin_20230109_DetermineCorrelationFunctions\20230127_Results_MYELINsolidMyelin_20221222_MYELIN_TIP4_Bilayer_50water_solidMyelin_H_whole_dt4ps_simTime1000ns.mat");

nearestNeighbourCase = 5500;
nNIndex = results.nearestNeighbourCases == nearestNeighbourCase;
if isempty(nNIndex)
    error("Not determined nearestNeighbourCase");
end

zeroPaddingLength = 100000;
corrFuncFirstOrder =  results.sumCorrFuncFirstOrder(nNIndex,:) ...
    /results.atomCounter;
corrFuncFirstOrder(end+1:end+zeroPaddingLength) ...
    = zeros(1,zeroPaddingLength);
corrFuncSecondOrder = results.sumCorrFuncSecondOrder(nNIndex,:) ...
    /results.atomCounter;
corrFuncSecondOrder(end+1:end+zeroPaddingLength) ...
    = zeros(1,zeroPaddingLength);

avgRegion = round(length(corrFuncFirstOrder)*[0.99 1]);


fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;

fig = initializeFigure();
ax = initializeSubplot(fig,2,2,1);
timeAxis = [0:length(corrFuncFirstOrder)-1]*results.deltaTInS;
plot(timeAxis,real(corrFuncFirstOrder));
plot(timeAxis,imag(corrFuncFirstOrder));
plot(timeAxis,real(corrFuncSecondOrder));
plot(timeAxis,imag(corrFuncSecondOrder));
legend("1st Order real","1st Order imag","2nd Order real" ...
    ,"2nd Order imag");
title(sprintf("NN: %i, zero pad: %i" ...
    ,nearestNeighbourCase(nNIndex),zeroPaddingLength));

[specDens1, specDens2] ...
    = calculateSpectralDensities(corrFuncFirstOrder ...
    ,corrFuncSecondOrder,omega0,results.deltaTInS,[0 1]);
initializeSubplot(fig,2,2,2);
plot(timeAxis,real(specDens1));
plot(timeAxis,real(specDens2));
legend("1st Order", "2nd Order");

r1 = calculateR1WithSpectralDensity( ...
    specDens1(avgRegion(1):avgRegion(2)) ...
    ,specDens2(avgRegion(1):avgRegion(2)),results.dipolDipolConstant);
title(sprintf("Real(specDens) zero pad R$_1$: %.4f",r1));


timeAxis = timeAxis(1:end-zeroPaddingLength); 
avgRegion = avgRegion - zeroPaddingLength;
[specDens1WithoutZeroPad,specDens2WithoutZeroPad] ...
    = calculateSpectralDensities( ...
    corrFuncFirstOrder(1:end-zeroPaddingLength) ...
    ,corrFuncSecondOrder(1:end-zeroPaddingLength) ...
    ,omega0,results.deltaTInS,[0 1]);

initializeSubplot(fig,2,2,3);
plot(timeAxis,real(specDens1WithoutZeroPad));
plot(timeAxis,real(specDens2WithoutZeroPad));
legend("1st Order", "2nd Order");

r1WithoutZeroPad = calculateR1WithSpectralDensity( ...
    specDens1WithoutZeroPad(avgRegion(1):avgRegion(2)) ...
    ,specDens2WithoutZeroPad(avgRegion(1):avgRegion(2)) ...
    ,results.dipolDipolConstant);
title(sprintf("Real(specDens) R$_1$: %.4f" ...
    ,r1WithoutZeroPad));

timeStepSkips = [1 2 3 4 5 6];
for timeStep = timeStepSkips
    [specDens1WithoutZeroPad,specDens2WithoutZeroPad] ...
        = calculateSpectralDensities( ...
        corrFuncFirstOrder(1:timeStep:end-zeroPaddingLength) ...
        ,corrFuncSecondOrder(1:timeStep:end-zeroPaddingLength) ...
        ,omega0,results.deltaTInS,[0 1]);
    avgRegion = round(avgRegion./timeStep);
    r1AtTimeStep(timeStep) = calculateR1WithSpectralDensity( ...
        specDens1WithoutZeroPad(avgRegion(1):avgRegion(2)) ...
        ,specDens2WithoutZeroPad(avgRegion(1):avgRegion(2)) ...
        ,results.dipolDipolConstant);
end

initializeSubplot(fig,2,2,4);
plot(timeStepSkips,r1AtTimeStep,"*-");
xlabel(sprintf("Multiples of dT : %.4d sec",results.deltaTInS))
ylabel("R$_1$");
lgd = legend();
lgd.Visible = 'off';


