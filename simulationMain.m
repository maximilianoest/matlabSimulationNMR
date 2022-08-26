clc
clearvars -except trajectoryX trajectoryY trajectoryZ ...
    simulationConfiguration

%% Initialize system

addpath(genpath('library'));
addpath(genpath('txtFiles'));
addpath(genpath('scripts'));
configuration = readConfigurationFile('configMain.txt');
[dataDirectory,resultsFileSavingPath,logFilePath] = ...
    setUpDependenciesBasedOn(configuration);

logMessage('System initialized. Dependendencies were added',logFilePath)

%% Analyse data set

analyseFileNameAndCreateVariablesInBaseWorkspace(configuration ...
    ,logFilePath);
dataFilePath = sprintf('%s%s.mat',dataDirectory,fileName);
checkIfFileExists(dataFilePath,logFilePath);

%% Load Data

if ~configuration.dataLoaded
    logMessage('Start loading data.',dataFilePath);
    [trajectoryX,trajectoryY,trajectoryZ,simulationConfiguration] = ...
        loadTrajectoriesAndSimConfig(dataFilePath);
    logMessage('Loading data finished',dataFilePath)
else
    logMessage('Data was already loaded in run before.',dataFilePath)
end

%% Define constants
logMessage('Defining constants.',logFilePath,true);
constants = readConstantsFile('constants.txt');

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *configuration.mainMagneticField;

%% Define simulation parameters
% logMessage('Defining simulation parameters.',logFilePath,true);
% [numberOfHs,timeSteps] = size(trajectoryX);
% 
% logMessage(sprintf(['Found %d hydrogen atoms at %d time ' ...
%     'steps of %.3d s'],numberOfHs,timeSteps,timeBetweenTimeStepsInS) ...
%     ,logFilePath,false);
% 
% lags = round(configuration.fractionForLags*timeSteps);
% logMessage(sprintf(['The lag is set to %d time steps, which ' ...
%     'is equivalent to %.2f %%. This configuration only shortens the '...
%     'correlation functions and NOT the simulation time.'],lags ...
%     ,(configuration.fractionForLags)*100),logFilePath,false);
% 
% nearestNeighbours = configuration.nearestNeighboursCase;
% if nearestNeighbours >= numberOfHs
%     logMessage(['The number of nearest neighbours is higher than '
%         'the number of possible atoms. PLEASE CHECK YOUR CONFIG ' ...
%         'FILE!'],path2LogFile);
%     error(['The number of nearest neighbours is higher than the ' ...
%         'number of possible atoms. Please check your config file!']);
% end
% 
% logMessage(sprintf(['Analysing %.f nearst neighbours of ' ...
%     'overall %.f hydrogen atoms'],nearestNeighbours,numberOfHs) ...
%     ,path2LogFile,false);
% 
% atomsToCalculate = configuration.atomsToCalculate;
% startDateOfSimulation = datestr(now,'yyyymmdd');
% logMessage(sprintf('Calculate %i atoms',atomsToCalculate) ...
%     ,path2LogFile,false);
% 
% orientationAngles = deg2rad(getValuesFromStringEnumeration( ...
%     configuration.fibreOrientations,';','numeric'));
% fibreOrientationsCount = size(orientationAngles,2);
% logMessage(['Found the orientations' sprintf(' %.f',rad2deg( ...
%     orientationAngles))],path2LogFile,false);
% 
% positionAngles = deg2rad(getValuesFromStringEnumeration( ...
%     configuration.myelinPositions,';','numeric'));
% positionsInMyelinCount = size(positionAngles,2);
% logMessage(['Found the positions' sprintf(' %.f',rad2deg( ...
%     positionAngles))],path2LogFile,false);









