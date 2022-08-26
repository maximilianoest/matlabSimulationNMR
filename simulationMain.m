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

logMessage('Defining simulation parameters.',logFilePath,true);
[numberOfHs,timeSteps] = size(trajectoryX);
logMessage(sprintf(['    Found %d hydrogen atoms at %d time ' ...
    'steps of %.3d s'],numberOfHs,timeSteps,timeBetweenTimeStepsInS) ...
    ,logFilePath,false);

atomsToCalculate = configuration.atomsToCalculate;
startDateOfSimulation = datestr(now,'yyyymmdd');
logMessage(sprintf('    Calculate %i atoms',atomsToCalculate) ...
    ,logFilePath,false);

%% Starting external script
orientationDependency_lipidWater_longSimulationTime








