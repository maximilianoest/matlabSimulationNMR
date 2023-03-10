%% initialize
clc; clear all; close all;

addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile("constants.txt");
saving = 1;

fieldStrength = 3; % Tesla
gyromagneticRatio = constants.gyromagneticRatioOfHydrogenAtom;
omega0 = fieldStrength * gyromagneticRatio;

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);

%% set up dependencies
resultsDirectory = "C:\Users\maxoe\Google Drive\Promotion\Simulation" ...
    + "\RESULTS\angleDependency_solidLipid_bigSimulation\";

filesInDirectory = [ ...
    "20221112_Results_DOPSlipid_20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns" ...
    "20221104_Results_PLPClipid_20220401_PLPC_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns" ...
    "20221108_Results_PSMlipid_20220509_PSM_TIP4_Bilayer_50water_lipid_H_whole_dt4ps_simTime1000ns"];
filesInDirectory = filesInDirectory + ".mat";

%% load data and determine R1

datasetNr = 1;
partOfSpecDensToPlot = [0.8 1];
partOfSpecDensToCalculateR1 = 0.99;
allLipidNames = {};

for datasetNr = 1:length(filesInDirectory)
    fig = initializeFigure();
    binWidth = 0.1;
    results = load(resultsDirectory + filesInDirectory(datasetNr));
    lipidName = results.whichLipid;
    deltaTInS = results.deltaTInS;
    atomCounter = results.atomCounter;
    simDuration = results.simulationDurationInS;
    
    % correlation function
    corrFuncFirstOrder = results.sumCorrFuncFirstOrder / atomCounter;
    corrFuncSecondOrder = results.sumCorrFuncSecondOrder / atomCounter;
    timeAxis = 0:deltaTInS:(length(corrFuncFirstOrder) - 1)*deltaTInS;
    initializeSubplot(fig,2,2,1);
    plot(timeAxis,real(corrFuncFirstOrder));
    plot(timeAxis,imag(corrFuncFirstOrder));
    plot(timeAxis,real(corrFuncSecondOrder));
    plot(timeAxis,imag(corrFuncSecondOrder));
    legend("real 1st order","imag 1st order" ...
        ,"real 2nd order","imag 2nd order");
    title(lipidName + "  corrFunc");
    xlabel("lag time [s]");
    
    % spectral density
    initializeSubplot(fig,2,2,2);
    shortenedTimeAxis = partOfSpecDensToPlot(1)*length(corrFuncFirstOrder)*deltaTInS ...
        :deltaTInS:(length(corrFuncFirstOrder) - 1)*deltaTInS*partOfSpecDensToPlot(2);
    [specDensFirstOrder, specDensSecondOrder] = ...
        calculateSpectralDensities(corrFuncFirstOrder,corrFuncSecondOrder ...
        ,omega0,deltaTInS,[0 1]);
    plot(shortenedTimeAxis,real(specDensFirstOrder(round(partOfSpecDensToPlot(1)*end) + 1 ...
        :round(partOfSpecDensToPlot(2)*end))));
    plot(shortenedTimeAxis,real(specDensSecondOrder(round(partOfSpecDensToPlot(1)*end) + 1 ...
        :round(partOfSpecDensToPlot(2)*end))));
    legend("1st order", "2nd order","location","east");
    title("real(Spectral density)");
    xlabel("Upper integration limit [s]");
    
    % determine R1
    % 1. with real and imag part of corr funcs
    overallPoolR1(datasetNr) = calculateR1WithSpectralDensity( ...
        specDensFirstOrder(round(partOfSpecDensToCalculateR1*end):end) ...
        ,specDensSecondOrder(round(partOfSpecDensToCalculateR1*end):end) ...
        ,dipolDipolConstant); %#ok<SAGROW>
    
    % 2. only with real part of corr funcs
    [realSpecDensFirstOrder, realSpecDensSecondOrder] = ...
        calculateSpectralDensities(real(corrFuncFirstOrder) ...
        ,real(corrFuncSecondOrder),omega0,deltaTInS ...
        ,[partOfSpecDensToCalculateR1 1]);
    r1OnlyRealCorrFunc = calculateR1WithSpectralDensity( ...
        realSpecDensFirstOrder,realSpecDensSecondOrder,dipolDipolConstant);
    
    % 3. location dependent R1
    numberOfAtomsInLocation = sum(results.xPosMatrix ~= 0,2);
    
    for locationNr = 1:length(numberOfAtomsInLocation)
        numberOfAtomsInThisLocation = numberOfAtomsInLocation(locationNr);
        location = mean(results.locationsForLocationDep( ...
            locationNr:locationNr+1)) * constants.nanoMeter;
        locDepCorrFuncFirstOrder = results.sumCorrFuncFirstOrderLocDep( ...
            locationNr,:) / numberOfAtomsInThisLocation;
        locDepCorrFuncSecondOrder = results.sumCorrFuncSecondOrderLocDep( ...
            locationNr,:) / numberOfAtomsInThisLocation;
        
        [locDepSpecDensFirstOrder, locDepSpecDensSecondOrder] = ...
            calculateSpectralDensities(locDepCorrFuncFirstOrder,locDepCorrFuncSecondOrder,omega0,deltaTInS,[partOfSpecDensToCalculateR1 1]);
        
        locDepR1(1,locationNr) = location; %#ok<SAGROW>
        locDepR1(2,locationNr) = calculateR1WithSpectralDensity( ...
            locDepSpecDensFirstOrder,locDepSpecDensSecondOrder ...
            ,dipolDipolConstant); %#ok<SAGROW>
    end
    locationDepR1ForAllLipids(datasetNr,:,:) = locDepR1;
    allLipidNames{end+1} = lipidName;
    
    % plotting location dependency of R1 in lipids
    initializeSubplot(fig,2,2,3);
    plot(locDepR1(1,:),locDepR1(2,:));
    title("Location dependent R1")
    ylabel("R1 [Hz]");
    xlabel("Location [m]");
    lgd = legend();
    lgd.Visible = "off";
    
    % plotting results into fourth subplot slot
    fourthSubplotText = sprintf("Lipid: %s\n overall R1: %.4f \n R1 based" ...
        + " only real corr func: %.4f \n" ...
        ,lipidName,overallPoolR1(datasetNr),r1OnlyRealCorrFunc);
    
    ax = initializeSubplot(fig,2,2,4);
    text(0.05,0.5,fourthSubplotText,'interpreter','latex');
    set(ax,'visible','off');
    lgd = legend();
    lgd.Visible = 'Off';
    
    if saving
        saveFigureTo(resultsDirectory,"SolidLipidR1",lipidName ...
            ,results.matlabSimulationDate);
    end
    
end

%%
initializeFigure();
title("Location dependent $R_1$");
for datasetNr = 1:size(locationDepR1ForAllLipids,1)
   plot(squeeze(locationDepR1ForAllLipids(datasetNr,1,:)) ...
       ,squeeze(locationDepR1ForAllLipids(datasetNr,2,:)),"*-");
end
xlabel("Location $[m]$");
ylabel("Relaxation Rate $R_1$ $[Hz]$");
legend(allLipidNames);

if saving
    saveFigureTo(resultsDirectory,"LocationDependentR1","allLipids"...
            ,"");
end

initializeFigure();
title("Location dependent $R_1$");
for datasetNr = 1:size(locationDepR1ForAllLipids,1)
   plot(squeeze(locationDepR1ForAllLipids(datasetNr,1,2:end)) ...
       ,squeeze(locationDepR1ForAllLipids(datasetNr,2,2:end)),"*-");
end
xlabel("Location $[m]$");
ylabel("Relaxation Rate $R_1$ $[Hz]$");
legend(allLipidNames);

if saving
    saveFigureTo(resultsDirectory,"LocationDependentR1","allLipids"...
            ,"withoutFirstDatapoint");
    lipidNames = string(allLipidNames);
    save(resultsDirectory + "lipidRelaxationRates.mat",'lipidNames' ...
        ,'overallPoolR1');
end

%% 
% results = load(resultsDirectory + filesInDirectory(datasetNr));
% lipidName = results.whichLipid;
% deltaTInS = results.deltaTInS;
% atomCounter = results.atomCounter;
% simDuration = results.simulationDurationInS;
% 
% % correlation function
% corrFuncFirstOrder = results.sumCorrFuncFirstOrder / atomCounter;
% corrFuncSecondOrder = results.sumCorrFuncSecondOrder / atomCounter;
% timeAxis = 0:deltaTInS:(length(corrFuncFirstOrder) - 1)*deltaTInS;
% initializeSubplot(fig,2,2,1);
% plot(timeAxis,real(corrFuncFirstOrder));
% plot(timeAxis,imag(corrFuncFirstOrder));
% plot(timeAxis,real(corrFuncSecondOrder));
% plot(timeAxis,imag(corrFuncSecondOrder));
% legend("real 1st order","imag 1st order" ...
%     ,"real 2nd order","imag 2nd order");
% title(lipidName + "  corrFunc");
% xlabel("lag time [s]");
% 
% % spectral density
% initializeSubplot(fig,2,2,2);
% shortenedTimeAxis = partOfSpecDensToPlot(1)*length(corrFuncFirstOrder)*deltaTInS ...
%     :deltaTInS:(length(corrFuncFirstOrder) - 1)*deltaTInS*partOfSpecDensToPlot(2);
% [specDensFirstOrder, specDensSecondOrder] = ...
%     calculateSpectralDensities(corrFuncFirstOrder,corrFuncSecondOrder ...
%     ,omega0,deltaTInS,[0 1]);
% plot(shortenedTimeAxis,real(specDensFirstOrder(round(partOfSpecDensToPlot(1)*end) + 1 ...
%     :round(partOfSpecDensToPlot(2)*end))));
% plot(shortenedTimeAxis,real(specDensSecondOrder(round(partOfSpecDensToPlot(1)*end) + 1 ...
%     :round(partOfSpecDensToPlot(2)*end))));
% legend("1st order", "2nd order","location","east");
% title("real(Spectral density)");
% xlabel("Upper integration limit [s]");
% 
% % determine R1
% % 1. with real and imag part of corr funcs
% overallPoolR1 = calculateR1WithSpectralDensity( ...
%     specDensFirstOrder(round(partOfSpecDensToCalculateR1*end):end) ...
%     ,specDensSecondOrder(round(partOfSpecDensToCalculateR1*end):end) ...
%     ,dipolDipolConstant);
% 
% % 2. only with real part of corr funcs
% [realSpecDensFirstOrder, realSpecDensSecondOrder] = ...
%     calculateSpectralDensities(real(corrFuncFirstOrder) ...
%     ,real(corrFuncSecondOrder),omega0,deltaTInS ...
%     ,[partOfSpecDensToCalculateR1 1]);
% r1OnlyRealCorrFunc = calculateR1WithSpectralDensity( ...
%     realSpecDensFirstOrder,realSpecDensSecondOrder,dipolDipolConstant);
% 
% % 3. location dependent R1
% numberOfAtomsInLocation = sum(results.xPosMatrix ~= 0,2);
% 
% for locationNr = 1:length(numberOfAtomsInLocation)
%     numberOfAtomsInThisLocation = numberOfAtomsInLocation(locationNr);
%     location = mean(results.locationsForLocationDep( ...
%         locationNr:locationNr+1)) * constants.nanoMeter;
%     locDepCorrFuncFirstOrder = results.sumCorrFuncFirstOrderLocDep( ...
%         locationNr,:) / numberOfAtomsInThisLocation;
%     locDepCorrFuncSecondOrder = results.sumCorrFuncSecondOrderLocDep( ...
%         locationNr,:) / numberOfAtomsInThisLocation;
%     
%     [locDepSpecDensFirstOrder, locDepSpecDensSecondOrder] = ...
%         calculateSpectralDensities(locDepCorrFuncFirstOrder,locDepCorrFuncSecondOrder,omega0,deltaTInS,[partOfSpecDensToCalculateR1 1]);
%     
%     locDepR1(1,locationNr) = location; %#ok<SAGROW>
%     locDepR1(2,locationNr) = calculateR1WithSpectralDensity( ...
%         locDepSpecDensFirstOrder,locDepSpecDensSecondOrder ...
%         ,dipolDipolConstant); %#ok<SAGROW>
% end
% 
% % plotting location dependency of R1 in lipids
% initializeSubplot(fig,2,2,3);
% plot(locDepR1(1,:),locDepR1(2,:));
% title("Location dependent R1")
% ylabel("R1 [Hz]");
% xlabel("Location [m]");
% lgd = legend();
% lgd.Visible = "off";
% 
% % plotting results into fourth subplot slot
% fourthSubplotText = sprintf("Lipid: %s\n overall R1: %.4f \n R1 based" ...
%     + " only real corr func: %.4f \n" ...
%     ,lipidName,overallPoolR1,r1OnlyRealCorrFunc);
% 
% ax = initializeSubplot(fig,2,2,4);
% text(0.05,0.5,fourthSubplotText,'interpreter','latex');
% set(ax,'visible','off');
% lgd = legend();
% lgd.Visible = 'Off';






