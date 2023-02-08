clc; clear all; close all; fclose('all');

addpath(genpath("../../library/"));
constants = readConstantsFile("../../txtFiles/constants.txt");
results = load("C:\Users\maxoe\Google Drive\Promotion\Simulation\RESULTS\myelinWater_20230206_DetermineCorrelationFunctions\20230208_Results_MYELINmyelinWater_20230202_MYELIN_TIP4_Monolayer_50water_myelinWater_H_whole_dt02ps_simTime20ns.mat");

nearestNeighbourCase = 15000;
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

avgPercent = [0.99 1];
avgRegion = round(length(corrFuncFirstOrder)*avgPercent);


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
title("Real(specDens)");


avgRegionWithoutPadding = int32((length(corrFuncFirstOrder)-zeroPaddingLength) ...
    *avgPercent);
[specDens1WithoutZeroPad,specDens2WithoutZeroPad] ...
    = calculateSpectralDensities( ...
    corrFuncFirstOrder(1:end-zeroPaddingLength) ...
    ,corrFuncSecondOrder(1:end-zeroPaddingLength) ...
    ,omega0,results.deltaTInS,[0 1]);
r1WithoutZeroPad = calculateR1WithSpectralDensity( ...
    specDens1WithoutZeroPad(avgRegionWithoutPadding(1):avgRegionWithoutPadding(2)) ...
    ,specDens2WithoutZeroPad(avgRegionWithoutPadding(1):avgRegionWithoutPadding(2)) ...
    ,results.dipolDipolConstant);




timeStepSkips = [1 2 3 4 5 6];
for timeStep = timeStepSkips
    corrFuncFirstOrderTimeSkips ...
        = corrFuncFirstOrder(1:timeStep:end-zeroPaddingLength);
    corrFuncSecondOrderTimeSkips ...
        = corrFuncSecondOrder(1:timeStep:end-zeroPaddingLength);
    
    [specDens1JustSkipped,specDens2JustSkipped] ...
        = calculateSpectralDensities( ...
        corrFuncFirstOrderTimeSkips,corrFuncSecondOrderTimeSkips...
        ,omega0,results.deltaTInS*timeStep,avgPercent);
    r1AtTimeStep(timeStep) = calculateR1WithSpectralDensity( ...
        specDens1JustSkipped,specDens2JustSkipped ...
        ,results.dipolDipolConstant);
    
    interpolCorrFuncFirstOrder = complex(zeros(1 ...
        ,length(corrFuncFirstOrderTimeSkips)*timeStep));
    interpolCorrFuncFirstOrder(1:timeStep:end) = corrFuncFirstOrderTimeSkips;
    interpolCorrFuncSecondOrder = zeros(1 ...
        ,length(corrFuncSecondOrderTimeSkips)*timeStep);
    interpolCorrFuncSecondOrder(1:timeStep:end) = corrFuncSecondOrderTimeSkips;
    
    for index = 1:timeStep:length(interpolCorrFuncFirstOrder)-timeStep
        for skip = 1:timeStep-1
            diffFirstOrder = interpolCorrFuncFirstOrder(index+timeStep) ...
                - interpolCorrFuncFirstOrder(index);
            interpolCorrFuncFirstOrder(index+skip) = ...
                interpolCorrFuncFirstOrder(index) ...
                + diffFirstOrder*skip/timeStep;
            
            diffSecondOrder = interpolCorrFuncSecondOrder(index+timeStep) ...
                - interpolCorrFuncSecondOrder(index);
            interpolCorrFuncSecondOrder(index+skip) = ...
                interpolCorrFuncSecondOrder(index) ...
                + diffSecondOrder*skip/timeStep;
            
        end
    end
    interpolCorrFuncFirstOrder = interpolCorrFuncFirstOrder( ...
        1:length(corrFuncFirstOrder)-zeroPaddingLength);
    interpolCorrFuncSecondOrder = interpolCorrFuncSecondOrder( ...
        1:length(corrFuncSecondOrder)-zeroPaddingLength);
    
    [interpolSpecDens1,interpolSpecDens2] = calculateSpectralDensities( ...
        interpolCorrFuncFirstOrder,interpolCorrFuncSecondOrder,omega0 ...
        ,results.deltaTInS,avgPercent);
    
    interpolatedR1(timeStep) = calculateR1WithSpectralDensity( ...
        interpolSpecDens1,interpolSpecDens2,results.dipolDipolConstant);
    
end

initializeSubplot(fig,2,2,3);
plot(timeStepSkips,r1AtTimeStep,"*-");
plot(timeStepSkips,interpolatedR1,"*-");
xlabel(sprintf("Multiples of dT : %.4d sec",results.deltaTInS))
ylabel("R$_1$");
lgd = legend("just skipped","linearly interpolated");
lgd.Location = "northwest";
title("Skip time steps");


calculatedAtomsAxis = 1:results.atomCounter;
r1Estimation = results.r1Estimation(1:results.atomCounter);
initializeSubplot(fig,2,2,4);
plot(calculatedAtomsAxis,r1Estimation);
% axis([0 inf 0.75 1.25] * mean(r1Estimation(1:end)))
lgd = legend();
lgd.Visible = 'Off';
ylabel("R$_1$ [Hz]");
xlabel("Number of evaluateda atoms");
title(sprintf("R$_1$ pad/no pad: %.4f / %.4f" ...
    ,r1,r1WithoutZeroPad));



