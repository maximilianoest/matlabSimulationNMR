clc; clear all; close all; fclose('all');

addpath(genpath("../../library/"));
constants = readConstantsFile("../../txtFiles/constants.txt");
results = load("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\solidMyelin_20230109_DetermineCorrelationFunctions\20230127_Results_MYELINsolidMyelin_20221222_MYELIN_TIP4_Bilayer_50water_solidMyelin_H_whole_dt4ps_simTime1000ns.mat");

corrFuncFirstOrder =  results.sumCorrFuncFirstOrder(2,:)/results.atomCounter;
% corrFuncFirstOrder(end+1:end+100000) = zeros(1,100000);
corrFuncSecondOrder = results.sumCorrFuncSecondOrder(2,:)/results.atomCounter;
% corrFuncSecondOrder(end+1:end+100000) = zeros(1,100000);

fieldStrength = 3; % Tesla
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrength;

fig = initializeFigure();
ax = initializeSubplot(fig,2,2,1);
timeAxis = [0:length(corrFuncFirstOrder)-1]*results.deltaTInS;
plot(timeAxis,abs(corrFuncFirstOrder));
plot(timeAxis,abs(corrFuncSecondOrder));
legend("1st Order", "2nd Order");

[specDens1, specDens2] ...
    = calculateSpectralDensities(corrFuncFirstOrder ...
    ,corrFuncSecondOrder,omega0,results.deltaTInS,[0 1]);
ax = initializeSubplot(fig,2,2,2);
plot(timeAxis,real(specDens1));
plot(timeAxis,real(specDens2));
legend("1st Order", "2nd Order");

r1 = calculateR1WithSpectralDensity(specDens1,specDens2 ...
    ,results.dipolDipolConstant);


[averageSpectralDensityFirstOrder,averageSpectralDensitySecondOrder] ...
    = calculateSpectralDensities(results.sumCorrFuncFirstOrder(1,:)/results.atomCounter ...
    ,results.sumCorrFuncSecondOrder(1,:)/results.atomCounter,omega0,results.deltaTInS ...
    ,[0 1]);
ax = initializeSubplot(fig,2,2,3);
plot(timeAxis,real(averageSpectralDensityFirstOrder));
plot(timeAxis,real(averageSpectralDensitySecondOrder));
legend("1st Order", "2nd Order");

r1Estimation = calculateR1WithSpectralDensity( ...
    averageSpectralDensityFirstOrder ...
    ,averageSpectralDensitySecondOrder,results.dipolDipolConstant);





