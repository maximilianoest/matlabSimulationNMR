clc; clear all; close all;
addpath(genpath(sprintf('..%s..%slibrary',filesep,filesep)));
addpath(genpath(sprintf('..%s..%stxtFiles',createFilesepStringArray(2))));
constants = readConstantsFile('constants.txt');

resultsDir = sprintf('..%s..%sRESULTS%sangleDependency_lipidWater_bigSimulation' ...
    ,createFilesepStringArray(3));
if ~exist(resultsDir,'dir')
    error('The results directory does not exist.');
end
addpath(genpath(resultsDir));

savingDirectory = sprintf('%s%sPlots%s',resultsDir ...
    ,createFilesepStringArray(2));

dipoleDipoleConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
fieldStrengthInT = 3; % main magnetic field strength
omega0 = constants.gyromagneticRatioOfHydrogenAtom * fieldStrengthInT;
saving = 0;

nNCase = 'nearestNeighbours8000';
allMatFilesInResultsDirectory = dir(sprintf( ...
    '%s%s*.mat',resultsDir,filesep));
averagingRegion = [0.98 1.0];
% [left bottom width height]
figPosAndSize = [50 50 1000 700];
widthAndHeight = [0.43 0.40];
subFigPosAndSize = [ ...
    0.07 0.55 widthAndHeight ; ...
    0.55 0.55 widthAndHeight ; ...
    0.07 0.06 widthAndHeight ; ...
    0.55 0.06 widthAndHeight];

configVariables = who;

%% get relevant results from files
sumCorrFuncNames = ["sumCorrFuncZerothOrder", "sumCorrFuncFirstOrder" ...
    ,"sumCorrFuncSecondOrder"]; 
corrFuncNames = strrep(sumCorrFuncNames,"sumC","c");
relevantResults = struct();

for lipidNr = 1:length(allMatFilesInResultsDirectory)
    data = load(allMatFilesInResultsDirectory(lipidNr).name);
    relevantResults(lipidNr).lipid = data.whichLipid;
    relevantResults(lipidNr).atomCounter = data.atomCounter;
    relevantResults(lipidNr).simDur = data.simulationDurationInS;
    relevantResults(lipidNr).dT = data.deltaTInS;
    relevantResults(lipidNr).theta = rad2deg(data.fibreAnglesTheta);
    relevantResults(lipidNr).phi = rad2deg(data.fibreAnglesPhi);
    relevantResults(lipidNr).matlabSimulationDate = ...
        data.matlabSimulationDate;
    relevantResults(lipidNr).meanPositions = data.meanPositions;
    for name = sumCorrFuncNames
        corrFuncName = strrep(name,"sumC","c");
        relevantResults(lipidNr).(corrFuncName) = data.(name).(nNCase) ...
            /relevantResults(lipidNr).atomCounter;
    end
end


%% analyse results with plots

for lipidNr = 1:length(relevantResults)
   lipidData = relevantResults(lipidNr);
   fig = initializeFigure();
   corrFuncs = [];
   % correlation function plot
   initializeSubplot(fig,2,1,1);
   lgd = legend();
   legendEntries = {};
   timeAxis = 0 : lipidData.dT : lipidData.simDur;
   for name = corrFuncNames
      switch name
          case "corrFuncFirstOrder"
              corrFuncs(1,:) = lipidData.(name);
          case "corrFuncSecondOrder"
              corrFuncs(2,:) = lipidData.(name);
          otherwise
              disp("Correlation function zeroth Order.")
      end
      plot(timeAxis,abs(lipidData.(name))); 
      legendEntries{end+1} = name; %#ok<SAGROW>
   end
   lgd.String = legendEntries;
   axis([-inf inf 0 6e3]);
   title(sprintf("Lipid: %s", lipidData.lipid));
   ylabel("abs(Corr.Func)");
   xlabel("lag time [s]");
   
   % spectral density plot
   [specDensFirstOrder,specDensSecondOrder] = ...
       calculateSpectralDensities(corrFuncs(1,:),corrFuncs(2,:) ...
       ,omega0,lipidData.dT,[0 1]);
   initializeSubplot(fig,2,1,2);
   plot(timeAxis,real(specDensFirstOrder));
   plot(timeAxis,real(specDensSecondOrder));
   legend("first order", "second order");
   xlabel("Upper integration limit [s]")
   ylabel("real(spec. Dens.)");
   relevantResults(lipidNr).R1 = calculateR1WithSpectralDensity( ...
       mean(specDensFirstOrder(round(0.9*end):end)) ...
       ,mean(specDensSecondOrder(round(0.9*end):end)) ...
       ,dipoleDipoleConstant);
   title(sprintf("R1: %.4f, $\\theta$: %.2f, $\\phi$: %.2f" ...
       ,relevantResults(lipidNr).R1,lipidData.theta ...
       ,lipidData.phi));
   if saving
       saveFigureTo(savingDirectory,lipidData.lipid ...
           ,lipidData.matlabSimulationDate ...
           ,sprintf("lipidWater_R1analysis_nrAtoms%i" ...
           ,lipidData.atomCounter));
   end
   
end
