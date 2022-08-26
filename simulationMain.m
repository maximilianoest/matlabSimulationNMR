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
dataFilePath = sprintf('%s%s%s',dataDirectory,filesep,fileName);
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


% logMessage(sprintf('     The system has %i hydrogen atoms.' ...
%     ,size(trajectoryX,1)),path2LogFile,false);
% 
% %% Define constants
% logMessage('Defining constants.',path2LogFile,false);
% constants = readConstantsFile(path2ConstantsFile);
% 
% dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
%     *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
%     /(constants.nanoMeter^6);
% omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
%     *configuration.mainMagneticField;